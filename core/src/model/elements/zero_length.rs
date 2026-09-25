use nalgebra::{SMatrix, SVector};

use super::super::{
    FiberSection, FiberSection3, Material, Node, Node3, Node3Id, NodeId, ELEMENT_DOF, NDF, SPATIAL_NDF,
};
use super::truss::{SpatialElementMatrix, SpatialElementVector};

/// Couples one "shear" DOF to a *different* DOF's own material — the
/// normal direction's current force — via a bilinear-kinematic-hardening
/// Coulomb return map: sticks (stiffness `k0`) until `|V|` reaches
/// `mu * yield-eligible normal force`, then slides with post-slip tangent
/// `b*k0` (see the `b` field doc comment for why `b` is required, not
/// optional).
///
/// This is deliberately *not* a `Material` variant: `Material::evaluate` is
/// a pure function of one scalar strain (every call site — `Truss`,
/// `ZeroLength`, `Fiber`, the `Parallel`/`Series`/`MinMax` composites —
/// depends on that signature), and friction's yield force needs a second
/// DOF's current force, which a single `Material` has no way to see. So the
/// coupling lives here, in the element, not the material catalog.
///
/// Sign convention: this codebase is compression-negative (see
/// `Material::Ent`/`Gap` — `strain <= 0.0` is the compression/engaged
/// branch). Friction is mobilized only when the normal direction is in
/// compression: yield force is `mu * (-normal_force).max(0.0)`.
///
/// The return map: `trial = k0*(shear_relative - slip)` is the elastic
/// predictor (same shape as `ElasticPP`'s `trial_stress = e*(strain - ep)`,
/// `elastic_like.rs`). Inside the yield surface, `force = trial`. Past it,
/// this is standard 1D linear-kinematic-hardening plasticity with
/// (fixed, for one call) yield force `yield_force` and hardening modulus
/// chosen so the post-yield tangent is exactly `b*k0`:
/// `force = b*trial + (1-b)*yield_force*sign(trial)`, which by construction
/// satisfies `k0*(shear_relative - next_slip) == force` — so, exactly like
/// `ElasticPP`'s `next_ep = strain - stress/e`, `next_slip =
/// shear_relative - force/k0`.
#[derive(Debug, Clone, Copy)]
pub struct Friction {
    pub normal_dof: usize,
    pub shear_dof: usize,
    pub mu: f64,
    pub k0: f64,
    /// Post-slip hardening ratio (fraction of `k0`) — required at
    /// construction, not defaulted, so callers make an active choice. A
    /// small nonzero value (recommended starting range `1e-3`-`1e-4`) keeps
    /// the post-slip tangent away from exactly zero: while sliding, this
    /// DOF is commonly the *only* stiffness path on its shear direction, so
    /// an exactly-zero tangent risks a structurally singular row/column in
    /// the global system (`AnalysisError::SingularSystem`) — the same
    /// reasoning behind `Material::MinMax`'s nonzero failed-state tangent
    /// floor, except here `b` is part of the actual force law (rounding the
    /// stick/slip corner, reducing Newton bouncing at the transition), not
    /// a tangent-only fudge. Set to `0.0` for an exact textbook Coulomb
    /// comparison and accept the singularity risk that comes back with it.
    pub b: f64,
    slip: f64,
}

impl Friction {
    pub fn new(normal_dof: usize, shear_dof: usize, mu: f64, k0: f64, b: f64) -> Self {
        Friction { normal_dof, shear_dof, mu, k0, b, slip: 0.0 }
    }

    /// `(shear_force, d(force)/d(shear_relative), d(force)/d(normal_force),
    /// next_slip)` — same `(value, tangent, next_committed_state)` shape as
    /// `Material::trial_stress_tangent`/committed-state pair, plus the extra
    /// coupling derivative friction needs so the caller can chain-rule it
    /// through the normal material's own tangent. Pure — safe to call every
    /// Newton iteration, same contract as `Material::trial_stress_tangent`.
    fn evaluate(&self, normal_force: f64, shear_relative: f64) -> (f64, f64, f64, f64) {
        let yield_force = self.mu * (-normal_force).max(0.0);
        let trial = self.k0 * (shear_relative - self.slip);
        if trial.abs() <= yield_force {
            (trial, self.k0, 0.0, self.slip)
        } else {
            let force = (self.b * trial.abs() + (1.0 - self.b) * yield_force).copysign(trial);
            // d(yield_force)/d(normal_force): -mu in the engaged
            // (compression) branch, 0 past it — matches the `max(.., 0.0)`
            // clamp above, including at the boundary (subgradient
            // convention, same as `Gap`/`Ent`'s `<=`).
            let d_yield_d_normal = if normal_force < 0.0 { -self.mu } else { 0.0 };
            let d_force_d_normal = (1.0 - self.b) * trial.signum() * d_yield_d_normal;
            let next_slip = shear_relative - force / self.k0;
            (force, self.b * self.k0, d_force_d_normal, next_slip)
        }
    }
}

/// `Friction`'s spatial counterpart: one normal DOF, two shear DOFs, each
/// run independently through `Friction`'s return map against the same
/// normal force (a square, not circular, interaction surface between the
/// two shear directions — the physically correct circular friction cone is
/// a materially harder shared-slip-vector return-map problem, deferred
/// until needed).
#[derive(Debug, Clone, Copy)]
pub struct Friction3 {
    pub normal_dof: usize,
    pub shear_dofs: [usize; 2],
    pub mu: f64,
    pub k0: f64,
    /// See `Friction::b`'s doc comment.
    pub b: f64,
    slip: [f64; 2],
}

impl Friction3 {
    pub fn new(normal_dof: usize, shear_dofs: [usize; 2], mu: f64, k0: f64, b: f64) -> Self {
        Friction3 { normal_dof, shear_dofs, mu, k0, b, slip: [0.0, 0.0] }
    }

    /// Runs `Friction::evaluate`'s return map independently for shear axis
    /// `i` (`0` or `1`) against the shared normal force.
    fn evaluate(&self, axis: usize, normal_force: f64, shear_relative: f64) -> (f64, f64, f64, f64) {
        let leg = Friction {
            normal_dof: self.normal_dof,
            shear_dof: self.shear_dofs[axis],
            mu: self.mu,
            k0: self.k0,
            b: self.b,
            slip: self.slip[axis],
        };
        leg.evaluate(normal_force, shear_relative)
    }
}

/// A 2-node, zero-length connector: no geometry or integration, just direct
/// per-DOF material evaluation (§3.1) — each direction with a material
/// assigned independently relates that DOF's relative displacement between
/// the two nodes to a force along that same (global) direction. No
/// orientation vectors (unlike OpenSees' general `ZeroLength`, which can
/// evaluate materials along arbitrary local axes) — directions are the
/// global DOF axes (including rotation, since M3's DOF bump), which is all
/// the current scope needs.
///
/// `friction` adds one cross-DOF coupling term on top of the independent
/// `materials` (see `Friction`'s doc comment for why this lives here rather
/// than as a `Material` variant): the shear DOF's own slot in `materials`
/// must stay `None` — `friction` owns that direction's resistance entirely.
#[derive(Debug, Clone)]
pub struct ZeroLength {
    pub node_i: NodeId,
    pub node_j: NodeId,
    materials: [Option<Material>; NDF],
    friction: Option<Friction>,
}

impl ZeroLength {
    pub fn new(node_i: NodeId, node_j: NodeId) -> Self {
        ZeroLength {
            node_i,
            node_j,
            materials: std::array::from_fn(|_| None),
            friction: None,
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    /// `friction.normal_dof` must already have a material set via
    /// `with_material` (debug-asserted, since evaluating friction reads the
    /// normal direction's current force from it); `friction.shear_dof`
    /// must *not*.
    pub fn with_friction(mut self, friction: Friction) -> Self {
        debug_assert!(
            self.materials[friction.shear_dof].is_none(),
            "friction owns its shear DOF's resistance; don't also assign it a material"
        );
        debug_assert!(
            self.materials[friction.normal_dof].is_some(),
            "friction reads its normal DOF's force from that DOF's own material"
        );
        self.friction = Some(friction);
        self
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let mut k = SMatrix::<f64, ELEMENT_DOF, ELEMENT_DOF>::zeros();
        let mut resistance = SVector::<f64, ELEMENT_DOF>::zeros();

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            // Local DOF order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; this
            // direction only couples node_i's and node_j's copy of the same
            // dof.
            let mut b = SVector::<f64, ELEMENT_DOF>::zeros();
            b[dof] = -1.0;
            b[NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        if let Some(fr) = &self.friction {
            let normal_material = self.materials[fr.normal_dof].as_ref().expect("friction needs a normal-direction material");
            let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
            let (normal_force, normal_tangent) = normal_material.trial_stress_tangent(normal_rel);
            let shear_rel = node_j.displacement[fr.shear_dof] - node_i.displacement[fr.shear_dof];
            let (force, dv_dshear, dv_dnormal, _) = fr.evaluate(normal_force, shear_rel);

            let mut b_v = SVector::<f64, ELEMENT_DOF>::zeros();
            b_v[fr.shear_dof] = -1.0;
            b_v[NDF + fr.shear_dof] = 1.0;
            let mut b_n = SVector::<f64, ELEMENT_DOF>::zeros();
            b_n[fr.normal_dof] = -1.0;
            b_n[NDF + fr.normal_dof] = 1.0;

            // Only a `b_v * b_nᵀ` block, not its transpose — `k` is not
            // symmetric in general. Fine here: the solver uses a general
            // sparse LU, not a symmetric-only factorization.
            k += dv_dshear * (b_v * b_v.transpose()) + (dv_dnormal * normal_tangent) * (b_v * b_n.transpose());
            resistance += force * b_v;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }

        if let Some(fr) = &mut self.friction {
            // The normal DOF's own slot in `materials` already advanced its
            // state in the loop above; read (not commit) it here via
            // `trial_stress_tangent` — every leaf material in this catalog
            // reproduces the same force whether read just before or just
            // after its own `commit()`, so this ordering is safe, but
            // committing it a second time from here would not be.
            let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
            let normal_force = self.materials[fr.normal_dof].as_ref().unwrap().trial_stress_tangent(normal_rel).0;
            let shear_rel = node_j.displacement[fr.shear_dof] - node_i.displacement[fr.shear_dof];
            let (.., next_slip) = fr.evaluate(normal_force, shear_rel);
            fr.slip = next_slip;
        }
    }

    /// No orientation vectors, no separate local frame — see `Truss::
    /// local_force`'s doc comment for why this is just the resistance
    /// vector, recomputed fresh from each direction's current committed
    /// material state.
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// `ZeroLength`'s spatial counterpart: independent per-DOF materials along
/// the six global directions `[ux, uy, uz, rx, ry, rz]` — no orientation
/// vectors, same simplification `ZeroLength` already makes relative to
/// OpenSees' general `ZeroLength` (which can evaluate materials along
/// arbitrary local axes via a user-supplied orientation). Arbitrarily
/// oriented local springs are a later, explicit transform feature (see
/// `docs/spatial-architecture.md`'s "Elements and transforms" section).
/// `friction` — see `ZeroLength`'s doc comment — adds two independent
/// (Phase 2, see `Friction3`) shear-DOF couplings on top of the
/// independent `materials`.
#[derive(Debug, Clone)]
pub struct ZeroLength3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    materials: [Option<Material>; SPATIAL_NDF],
    friction: Option<Friction3>,
}

impl ZeroLength3 {
    pub fn new(node_i: Node3Id, node_j: Node3Id) -> Self {
        ZeroLength3 {
            node_i,
            node_j,
            materials: std::array::from_fn(|_| None),
            friction: None,
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    /// See `ZeroLength::with_friction`'s doc comment — same preconditions,
    /// checked for both entries of `friction.shear_dofs`.
    pub fn with_friction(mut self, friction: Friction3) -> Self {
        for shear_dof in friction.shear_dofs {
            debug_assert!(
                self.materials[shear_dof].is_none(),
                "friction owns its shear DOFs' resistance; don't also assign them a material"
            );
        }
        debug_assert!(
            self.materials[friction.normal_dof].is_some(),
            "friction reads its normal DOF's force from that DOF's own material"
        );
        self.friction = Some(friction);
        self
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let mut k = SpatialElementMatrix::zeros();
        let mut resistance = SpatialElementVector::zeros();

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut b = SpatialElementVector::zeros();
            b[dof] = -1.0;
            b[SPATIAL_NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        if let Some(fr) = &self.friction {
            let normal_material = self.materials[fr.normal_dof].as_ref().expect("friction needs a normal-direction material");
            let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
            let (normal_force, normal_tangent) = normal_material.trial_stress_tangent(normal_rel);

            let mut b_n = SpatialElementVector::zeros();
            b_n[fr.normal_dof] = -1.0;
            b_n[SPATIAL_NDF + fr.normal_dof] = 1.0;

            for (axis, &shear_dof) in fr.shear_dofs.iter().enumerate() {
                let shear_rel = node_j.displacement[shear_dof] - node_i.displacement[shear_dof];
                let (force, dv_dshear, dv_dnormal, _) = fr.evaluate(axis, normal_force, shear_rel);

                let mut b_v = SpatialElementVector::zeros();
                b_v[shear_dof] = -1.0;
                b_v[SPATIAL_NDF + shear_dof] = 1.0;

                k += dv_dshear * (b_v * b_v.transpose()) + (dv_dnormal * normal_tangent) * (b_v * b_n.transpose());
                resistance += force * b_v;
            }
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }

        if let Some(fr) = &mut self.friction {
            // See `ZeroLength::commit`'s doc comment for why reading (not
            // committing) the normal material here is correct.
            let normal_rel = node_j.displacement[fr.normal_dof] - node_i.displacement[fr.normal_dof];
            let normal_force = self.materials[fr.normal_dof].as_ref().unwrap().trial_stress_tangent(normal_rel).0;

            let mut next_slip = fr.slip;
            for (axis, &shear_dof) in fr.shear_dofs.iter().enumerate() {
                let shear_rel = node_j.displacement[shear_dof] - node_i.displacement[shear_dof];
                let (.., slip) = fr.evaluate(axis, normal_force, shear_rel);
                next_slip[axis] = slip;
            }
            fr.slip = next_slip;
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// A 2-node, zero-length connector driven by a coupled `FiberSection`
/// instead of `ZeroLength`'s independent per-DOF materials — the axial
/// (`ux`) and flexural (`rz`) response come out of the same fiber
/// discretization, so axial force and bending moment interact through
/// shared fiber strain, the way a concentrated-plasticity "fiber hinge"
/// needs (e.g. at a beam-column end) and `ZeroLength` structurally cannot
/// produce. Modeled on OpenSees/Xara's `ZeroLengthSection`
/// (`SRC/element/Point/ZeroLengthSection.cpp`): the section's basic
/// deformation is the direct nodal DOF difference (no shape functions or
/// integration — there's no length to integrate over), and `K`/`P` are
/// `Bᵀ*k_section*B`/`Bᵀ*q` with `B` selecting `[ux, rz]`, the same pattern
/// `DispBeamColumn::form_tangent_and_resistance` uses per integration point
/// (here with a single "point" and unit scale).
///
/// `uy` — the one DOF `FiberSection` has no resultant for — can still carry
/// an independent uniaxial `Material` spring via `with_material`, exactly
/// like `ZeroLength`, additive with (and uncoupled from) the section's
/// axial-moment response.
#[derive(Debug, Clone)]
pub struct ZeroLengthSection {
    pub node_i: NodeId,
    pub node_j: NodeId,
    section: FiberSection,
    materials: [Option<Material>; NDF],
}

impl ZeroLengthSection {
    pub fn new(node_i: NodeId, node_j: NodeId, section: FiberSection) -> Self {
        ZeroLengthSection {
            node_i,
            node_j,
            section,
            materials: std::array::from_fn(|_| None),
        }
    }

    /// Independent spring for a DOF the section doesn't drive — intended
    /// for `dof = 1` (`uy`, shear), since the section already owns `ux`
    /// (dof 0) and `rz` (dof 2).
    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa = node_j.displacement[2] - node_i.displacement[2];
        let (n, m, k_section) = self.section.trial(eps0, kappa);

        let b_eps0 = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[-1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        let b_kappa = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[0.0, 0.0, -1.0, 0.0, 0.0, 1.0]);

        let mut k = k_section[0][0] * (b_eps0 * b_eps0.transpose())
            + k_section[0][1] * (b_eps0 * b_kappa.transpose())
            + k_section[1][0] * (b_kappa * b_eps0.transpose())
            + k_section[1][1] * (b_kappa * b_kappa.transpose());
        let mut resistance = b_eps0 * n + b_kappa * m;

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut b = SVector::<f64, ELEMENT_DOF>::zeros();
            b[dof] = -1.0;
            b[NDF + dof] = 1.0;

            k += tangent_modulus * (b * b.transpose());
            resistance += force * b;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa = node_j.displacement[2] - node_i.displacement[2];
        self.section.commit(eps0, kappa);

        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

/// `ZeroLengthSection`'s spatial counterpart: coupled axial-biaxial-moment
/// response (`ux`, `ry`, `rz`) from a `FiberSection3`, plus optional
/// independent springs (via `with_material`) for `uy`, `uz` (shear) and
/// `rx` (torsion) — the three DOFs the section has no resultant for.
/// `FiberSection3` deliberately excludes torsion from its fiber loop (see
/// its doc comment), and there's no meaningful `G*J/length` term at zero
/// length the way `DispBeamColumn3` has, so an explicit `rx` material is
/// the only way to give this element torsional stiffness at all.
#[derive(Debug, Clone)]
pub struct ZeroLengthSection3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    section: FiberSection3,
    materials: [Option<Material>; SPATIAL_NDF],
}

impl ZeroLengthSection3 {
    pub fn new(node_i: Node3Id, node_j: Node3Id, section: FiberSection3) -> Self {
        ZeroLengthSection3 {
            node_i,
            node_j,
            section,
            materials: std::array::from_fn(|_| None),
        }
    }

    /// Independent spring for a DOF the section doesn't drive — `uy`/`uz`
    /// (shear, dofs 1/2) or `rx` (torsion, dof 3) — since the section
    /// already owns `ux` (dof 0), `ry` (dof 4) and `rz` (dof 5).
    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa_z = node_j.displacement[5] - node_i.displacement[5];
        let kappa_y = node_j.displacement[4] - node_i.displacement[4];
        let (n, mz, my, k_section) = self.section.trial(eps0, kappa_z, kappa_y);

        let mut b_eps0 = SpatialElementVector::zeros();
        b_eps0[0] = -1.0;
        b_eps0[SPATIAL_NDF] = 1.0;

        let mut b_kappa_z = SpatialElementVector::zeros();
        b_kappa_z[5] = -1.0;
        b_kappa_z[SPATIAL_NDF + 5] = 1.0;

        let mut b_kappa_y = SpatialElementVector::zeros();
        b_kappa_y[4] = -1.0;
        b_kappa_y[SPATIAL_NDF + 4] = 1.0;

        let mut k = SpatialElementMatrix::zeros();
        let mut resistance = b_eps0 * n + b_kappa_z * mz + b_kappa_y * my;

        let b = [&b_eps0, &b_kappa_z, &b_kappa_y];
        for (bi, row_i) in b.iter().enumerate() {
            for (bj, row_j) in b.iter().enumerate() {
                k += k_section[bi][bj] * (*row_i * row_j.transpose());
            }
        }

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.trial_stress_tangent(relative);

            let mut bm = SpatialElementVector::zeros();
            bm[dof] = -1.0;
            bm[SPATIAL_NDF + dof] = 1.0;

            k += tangent_modulus * (bm * bm.transpose());
            resistance += force * bm;
        }

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        let eps0 = node_j.displacement[0] - node_i.displacement[0];
        let kappa_z = node_j.displacement[5] - node_i.displacement[5];
        let kappa_y = node_j.displacement[4] - node_i.displacement[4];
        self.section.commit(eps0, kappa_z, kappa_y);

        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }

    /// See `ZeroLength::local_force`'s doc comment.
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        self.form_tangent_and_resistance(node_i, node_j).1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Fiber, Fiber3, Material, Node};

    /// Proves the `Material` enum's dispatch generalizes to `ZeroLength`
    /// across every M2 variant and regime (elastic, plastic plateau, open
    /// gap, closed gap, tension, compression) by driving the element
    /// directly off manually-set node displacements — independent of
    /// `Analysis`/`Algorithm`, which for M2 is still linear-only (no
    /// Newton iteration until M4) and so can't itself resolve equilibrium
    /// across a material's nonlinear regimes in one step.
    #[test]
    fn zero_length_dispatches_across_material_regimes() {
        let node_i = Node::new([0.0, 0.0]);

        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -0.02; // compression, past both EPP yield and Gap closure

        // `resistance = force * b` with `b[dof_i] = -1`, so node_i's
        // component of the resistance vector is `-force`.
        let epp = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::elastic_pp(100.0, 0.01));
        let (_, r) = epp.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * 0.01, "should be clamped to yield force");

        let gap = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::Gap { e: 100.0, gap: 0.01 });
        let (_, r) = gap.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * (0.02 - 0.01), "gap engaged past closure");

        let ent = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::Ent { e: 100.0 });
        let (_, r) = ent.form_tangent_and_resistance(&node_i, &node_j);
        assert_eq!(r[0], 100.0 * 0.02, "ENT carries full compression");

        let mut node_tension = Node::new([0.0, 0.0]);
        node_tension.displacement[0] = 0.02; // tension
        let (k, r) = ent.form_tangent_and_resistance(&node_i, &node_tension);
        assert_eq!(r[0], 0.0, "ENT carries no tension");
        assert_eq!(k[(0, 0)], 0.0, "no tangent stiffness while open");
    }

    /// `ZeroLength3`'s spatial counterpart, exercising a rotational DOF
    /// (index 3, `SpatialDof::Rx`) in addition to a translational one — the
    /// planar `ZeroLength` test above only ever touches dof 0.
    #[test]
    fn zero_length3_dispatches_independently_per_global_direction() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.01; // Ux
        node_j.displacement[3] = -0.02; // Rx

        let zl = ZeroLength3::new(Node3Id::default(), Node3Id::default())
            .with_material(0, Material::Elastic { e: 100.0 })
            .with_material(3, Material::elastic_pp(50.0, 0.01));

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);
        // `resistance = force * b` with `b[dof_i] = -1`, so node_i's
        // component of the resistance vector is `-force` (see the planar
        // `ZeroLength` test's doc comment).
        assert_eq!(r[0], -(100.0 * 0.01), "Ux: linear elastic");
        assert_eq!(r[SPATIAL_NDF], 100.0 * 0.01);
        assert_eq!(r[3], 50.0 * 0.01, "Rx: clamped to yield force (sign follows relative disp)");
        assert_eq!(r[SPATIAL_NDF + 3], -(50.0 * 0.01));
        assert_eq!(k[(0, 0)], 100.0, "Ux tangent unaffected by Rx material");
        assert_eq!(k[(3, 3)], 0.0, "Rx tangent is zero past yield");
        for dof in [1, 2, 4, 5] {
            assert_eq!(r[dof], 0.0, "unassigned directions carry no force");
        }
    }

    /// `ZeroLengthSection`: two symmetric elastic fibers reproduce the
    /// closed-form `EA`/`EI` (already validated directly against
    /// `FiberSection::trial` in `fiber_section.rs`'s
    /// `two_symmetric_elastic_fibers_reproduce_ea_ei_exactly`), now routed
    /// through the element's `B`-matrix mapping onto `[ux, rz]`.
    #[test]
    fn zero_length_section_reproduces_symmetric_ea_ei() {
        let (e, area, iz): (f64, f64, f64) = (30000.0, 2.0, 1000.0);
        let h = (iz / area).sqrt();
        let section = FiberSection::new(vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ]);

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001; // ux -> eps0
        node_j.displacement[2] = 0.0002; // rz -> kappa

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);
        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA mismatch: {}", k[(0, 0)]);
        assert!((k[(2, 2)] - e * iz).abs() < 1e-6, "EI mismatch: {}", k[(2, 2)]);
        assert!(k[(0, 2)].abs() < 1e-9, "no axial-moment coupling for a symmetric section");

        let n = e * area * 0.001;
        let m = e * iz * 0.0002;
        assert!((r[3] - n).abs() < 1e-6, "node_j axial force");
        assert!((r[0] + n).abs() < 1e-6, "node_i axial force");
        assert!((r[5] - m).abs() < 1e-6, "node_j moment");
        assert!((r[2] + m).abs() < 1e-6, "node_i moment");
        assert_eq!(r[1], 0.0, "no material assigned to uy");
    }

    /// An off-centroid single fiber couples axial and moment response
    /// (`k_section[0][1] != 0`) — something `ZeroLength`'s independent
    /// per-DOF materials can never produce, since each of its directions is
    /// evaluated against its own material in isolation.
    #[test]
    fn zero_length_section_couples_axial_and_moment_for_an_asymmetric_section() {
        let (e, area, y): (f64, f64, f64) = (1000.0, 2.0, 3.0);
        let section = FiberSection::new(vec![Fiber::new(y, area, Material::Elastic { e })]);

        let (_n, _m, k_section) = section.trial(0.001, 0.0002);
        assert!(k_section[0][1].abs() > 1e-9, "test setup: section must actually couple");

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[2] = 0.0002;

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);
        let (k, _r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        // b_eps0[0] = -1, b_kappa[2] = -1, so k(0,2) = k_section[0][1].
        assert!(
            (k[(0, 2)] - k_section[0][1]).abs() < 1e-9,
            "axial-moment coupling should carry through to the element tangent"
        );
    }

    /// `ZeroLengthSection3`: four corner elastic fibers reproduce `EA`,
    /// `EIy`, `EIz` (see `fiber_section.rs`'s
    /// `four_corner_elastic_fibers_reproduce_ea_eiy_eiz_with_zero_cross_coupling`)
    /// routed through the element's `B`-matrix mapping onto `[ux, ry, rz]`;
    /// `uy`/`uz`/`rx` — DOFs the section has no resultant for — stay zero.
    #[test]
    fn zero_length_section3_reproduces_biaxial_ea_eiy_eiz() {
        let (e, area, iy, iz): (f64, f64, f64, f64) = (30000.0, 4.0, 500.0, 2000.0);
        let hz = (iy / area).sqrt();
        let hy = (iz / area).sqrt();
        let a4 = area / 4.0;
        let section = FiberSection3::new(vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]);

        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.001; // ux -> eps0
        node_j.displacement[5] = 0.0002; // rz -> kappa_z
        node_j.displacement[4] = 0.0003; // ry -> kappa_y

        let zl = ZeroLengthSection3::new(Node3Id::default(), Node3Id::default(), section);
        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA mismatch");
        assert!((k[(5, 5)] - e * iz).abs() < 1e-6, "EIz mismatch");
        assert!((k[(4, 4)] - e * iy).abs() < 1e-6, "EIy mismatch");

        for dof in [1usize, 2, 3] {
            for local in [dof, SPATIAL_NDF + dof] {
                assert_eq!(r[local], 0.0, "unmodeled DOF carries no force");
                for col in 0..2 * SPATIAL_NDF {
                    assert_eq!(k[(local, col)], 0.0, "unmodeled DOF has no stiffness");
                }
            }
        }
    }

    /// A section commit propagates into the element's subsequent trial
    /// response, the same commit/trial split `Material` and `ZeroLength`
    /// already rely on — proves `ZeroLengthSection::commit` actually drives
    /// `FiberSection::commit` (not just recomputing `eps0`/`kappa` and
    /// discarding them).
    #[test]
    fn zero_length_section_commit_propagates_fiber_history() {
        let mat = || Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);

        // Pure axial fiber (y = 0) isolates the material's own history from
        // any bending coupling.
        let section = FiberSection::new(vec![Fiber::new(0.0, 1.0, mat())]);
        let mut zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section);

        // Independent "by hand" copy, committed in lockstep, to compute the
        // expected post-yield residual stress from scratch (same pattern as
        // `fiber_section.rs`'s RC test).
        let mut hand = mat();

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.01; // well past yield (epsy = 60/29000 ~= 0.00207)

        zl.commit(&node_i, &node_j);
        hand = hand.commit(0.01);

        // Fully unload back to zero relative displacement — with committed
        // plastic strain, elastic unload from 0.01 leaves residual stress.
        node_j.displacement[0] = 0.0;
        let (_k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        let (expected_n, _tangent) = hand.trial_stress_tangent(0.0);
        assert!(expected_n.abs() > 1e-6, "test setup: material should have residual stress after yield");
        assert!(
            (r[3] - expected_n).abs() < 1e-9,
            "committed plastic history should carry into subsequent trial"
        );
    }

    /// A `uy` spring alongside the section is additive and uncoupled: it
    /// contributes independently to `k`/`r` on top of (not instead of) the
    /// section's coupled axial-moment response.
    #[test]
    fn zero_length_section_combines_with_an_independent_shear_spring() {
        let (e, area, iz): (f64, f64, f64) = (30000.0, 2.0, 1000.0);
        let h = (iz / area).sqrt();
        let section = FiberSection::new(vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ]);
        let shear_e = 500.0;

        let zl = ZeroLengthSection::new(NodeId::default(), NodeId::default(), section)
            .with_material(1, Material::Elastic { e: shear_e });

        let node_i = Node::new([0.0, 0.0]);
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[1] = 0.002; // uy -> shear spring
        node_j.displacement[2] = 0.0002;

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(1, 1)] - shear_e).abs() < 1e-9, "shear spring stiffness");
        assert_eq!(k[(1, 0)], 0.0, "shear spring uncoupled from axial");
        assert_eq!(k[(1, 2)], 0.0, "shear spring uncoupled from moment");
        assert!((r[4] - shear_e * 0.002).abs() < 1e-9, "node_j shear force");

        // Section response unchanged by the added spring.
        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA unaffected by shear spring");
        assert!((k[(2, 2)] - e * iz).abs() < 1e-6, "EI unaffected by shear spring");
    }

    /// `ZeroLengthSection3` analogue: an `rx` torsion spring alongside the
    /// section's biaxial axial-moment response — the only way this element
    /// gets torsional stiffness at all (see the type's doc comment).
    #[test]
    fn zero_length_section3_combines_with_an_independent_torsion_spring() {
        let (e, area, iy, iz): (f64, f64, f64, f64) = (30000.0, 4.0, 500.0, 2000.0);
        let hz = (iy / area).sqrt();
        let hy = (iz / area).sqrt();
        let a4 = area / 4.0;
        let section = FiberSection3::new(vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]);
        let torsion_e = 700.0;

        let zl = ZeroLengthSection3::new(Node3Id::default(), Node3Id::default(), section)
            .with_material(3, Material::Elastic { e: torsion_e });

        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = 0.001;
        node_j.displacement[3] = 0.004; // rx -> torsion spring
        node_j.displacement[5] = 0.0002;
        node_j.displacement[4] = 0.0003;

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((k[(3, 3)] - torsion_e).abs() < 1e-9, "torsion spring stiffness");
        assert!((r[9] - torsion_e * 0.004).abs() < 1e-9, "node_j torque");
        assert!((k[(0, 0)] - e * area).abs() < 1e-6, "EA unaffected by torsion spring");
        assert!((k[(5, 5)] - e * iz).abs() < 1e-6, "EIz unaffected by torsion spring");
        assert!((k[(4, 4)] - e * iy).abs() < 1e-6, "EIy unaffected by torsion spring");
    }

    /// Common friction setup for the 2D tests below: dof 0 is the normal
    /// direction (an `Elastic` spring, `e = 1000`), dof 1 is the friction
    /// shear direction (`mu = 0.3`, `k0 = 500`, `b = 0.01`).
    fn friction_fixture() -> (Node, ZeroLength) {
        let node_i = Node::new([0.0, 0.0]);
        let zl = ZeroLength::new(NodeId::default(), NodeId::default())
            .with_material(0, Material::Elastic { e: 1000.0 })
            .with_friction(Friction::new(0, 1, 0.3, 500.0, 0.01));
        (node_i, zl)
    }

    /// Below the yield surface, friction is a plain elastic spring:
    /// `V == k0 * shear_rel`, tangent `== k0`.
    #[test]
    fn friction_sticks_below_yield_surface() {
        let (node_i, zl) = friction_fixture();
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -1.0; // compression -> normal_force = -1000, yield_force = 300
        node_j.displacement[1] = 0.1; // trial = 500*0.1 = 50, well under 300

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);
        assert!((r[NDF + 1] - 50.0).abs() < 1e-9, "V == k0*shear_rel while sticking");
        assert!((k[(1, 1)] - 500.0).abs() < 1e-9, "tangent == k0 while sticking");
    }

    /// Past the yield surface, `V` follows the bilinear-kinematic-hardening
    /// return map (`Friction::evaluate`'s doc comment): `V = b*trial +
    /// (1-b)*yield_force*sign(trial)`, tangent `== b*k0`.
    #[test]
    fn friction_slides_past_yield_surface() {
        let (node_i, zl) = friction_fixture();
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -1.0; // normal_force = -1000, yield_force = 300
        node_j.displacement[1] = 5.0; // trial = 500*5 = 2500, well past 300

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        let (mu, k0, b) = (0.3, 500.0, 0.01);
        let yield_force = mu * 1000.0;
        let trial = k0 * 5.0;
        let expected_v = b * trial + (1.0 - b) * yield_force;

        assert!((r[NDF + 1] - expected_v).abs() < 1e-9, "V matches the hardening-branch return map");
        assert!((k[(1, 1)] - b * k0).abs() < 1e-9, "tangent == b*k0 while sliding");
    }

    /// The same shear displacement against two different normal-direction
    /// states: shear resistance depends on axial load level, not just a
    /// generic plasticity check — a lighter normal force slides at the same
    /// shear displacement a heavier one still sticks at.
    #[test]
    fn friction_yield_surface_scales_with_normal_force() {
        let (node_i, zl) = friction_fixture();
        let shear_rel = 1.0; // trial = k0*1.0 = 500

        let mut node_j_light = Node::new([0.0, 0.0]);
        node_j_light.displacement[0] = -0.5; // normal_force = -500, yield_force = 150 < 500: slides
        node_j_light.displacement[1] = shear_rel;

        let mut node_j_heavy = Node::new([0.0, 0.0]);
        node_j_heavy.displacement[0] = -2.0; // normal_force = -2000, yield_force = 600 > 500: sticks
        node_j_heavy.displacement[1] = shear_rel;

        let (_, r_light) = zl.form_tangent_and_resistance(&node_i, &node_j_light);
        let (_, r_heavy) = zl.form_tangent_and_resistance(&node_i, &node_j_heavy);

        let (mu, k0, b) = (0.3, 500.0, 0.01);
        let expected_v_light = b * (k0 * shear_rel) + (1.0 - b) * (mu * 500.0);
        assert!((r_light[NDF + 1] - expected_v_light).abs() < 1e-9, "lighter normal force: sliding");
        assert!((r_heavy[NDF + 1] - k0 * shear_rel).abs() < 1e-9, "heavier normal force: still sticking, V == trial");
    }

    /// Mirrors `elastic_pp_remembers_permanent_set_after_commit`
    /// (`elastic_like.rs`): slide past the yield surface, commit, partially
    /// unload — the response must be elastic from the new committed `slip`,
    /// not from zero.
    #[test]
    fn friction_remembers_permanent_slip_after_commit() {
        let (node_i, mut zl) = friction_fixture();
        let (mu, k0, b) = (0.3, 500.0, 0.01);

        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -1.0; // normal_force = -1000, yield_force = 300
        node_j.displacement[1] = 5.0; // well past yield -> slides, accrues slip
        zl.commit(&node_i, &node_j);

        let yield_force = mu * 1000.0;
        let trial = k0 * 5.0;
        let force = b * trial + (1.0 - b) * yield_force;
        let committed_slip = 5.0 - force / k0;

        // Small elastic excursion from the newly committed slip.
        node_j.displacement[1] = committed_slip + 0.05;
        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        assert!((r[NDF + 1] - k0 * 0.05).abs() < 1e-9, "elastic unload from committed slip, not from zero");
        assert!((k[(1, 1)] - k0).abs() < 1e-9, "tangent back to k0 (elastic) after unload");
    }

    /// Out of compression (tension, per this codebase's sign convention),
    /// `yield_force` clamps to `0`: any nonzero shear displacement slides
    /// immediately, with only the post-slip hardening term surviving.
    #[test]
    fn friction_drops_capacity_when_normal_dof_is_in_tension() {
        let (node_i, zl) = friction_fixture();
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = 1.0; // tension -> normal_force = +1000 -> yield_force = 0
        node_j.displacement[1] = 0.001;

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        let (k0, b) = (500.0, 0.01);
        let expected_v = b * (k0 * 0.001); // yield_force == 0, so only the hardening term survives
        assert!((r[NDF + 1] - expected_v).abs() < 1e-9, "near-zero force: no normal-direction capacity");
        assert!((k[(1, 1)] - b * k0).abs() < 1e-9, "post-slip tangent even at zero capacity");
    }

    /// The most important test here, and the easiest place to get a sign or
    /// chain-rule factor wrong: perturb the normal DOF's relative
    /// displacement by a small `h` and confirm the analytical
    /// `k[shear_dof, normal_dof]` block matches `dV/d(normal_rel)` by
    /// finite difference while sliding — validates `Friction::evaluate`'s
    /// `d_force_d_normal` (and the chain rule through the normal material's
    /// own tangent in `form_tangent_and_resistance`) directly, rather than
    /// trusting the derivation by inspection.
    #[test]
    fn friction_off_diagonal_tangent_matches_finite_difference() {
        let (node_i, zl) = friction_fixture();
        let mut node_j = Node::new([0.0, 0.0]);
        node_j.displacement[0] = -1.0; // normal_force = -1000
        node_j.displacement[1] = 5.0; // well past yield -> sliding, d(force)/d(normal_force) != 0

        let (k, r0) = zl.form_tangent_and_resistance(&node_i, &node_j);
        let analytical = k[(1, 0)];

        let h = 1e-6;
        let mut node_j_plus = node_j.clone();
        node_j_plus.displacement[0] += h;
        let (_, r_plus) = zl.form_tangent_and_resistance(&node_i, &node_j_plus);

        let finite_diff = (r_plus[NDF + 1] - r0[NDF + 1]) / h;

        assert!(
            (analytical - finite_diff).abs() < 1e-6,
            "analytical dV/d(normal_rel) = {analytical}, finite difference = {finite_diff}"
        );
    }

    /// Phase 2: `ZeroLength3` with one normal DOF and two independent shear
    /// DOFs — each slides/sticks independently against the same normal
    /// force with no cross-talk between the two shear axes (the absence of
    /// coupling that distinguishes Phase 2 from the deferred, circular-
    /// interaction-surface Phase 3).
    #[test]
    fn friction3_shear_axes_slide_independently_against_shared_normal_force() {
        let node_i = Node3::new([0.0, 0.0, 0.0]);
        let mut node_j = Node3::new([0.0, 0.0, 0.0]);
        node_j.displacement[0] = -1.0; // normal DOF (Ux): normal_force = -1000
        node_j.displacement[1] = 5.0; // shear axis 0 (Uy): well past yield -> slides
        node_j.displacement[2] = 0.1; // shear axis 1 (Uz): well under yield -> sticks

        let (mu, k0, b) = (0.3, 500.0, 0.01);
        let zl = ZeroLength3::new(Node3Id::default(), Node3Id::default())
            .with_material(0, Material::Elastic { e: 1000.0 })
            .with_friction(Friction3::new(0, [1, 2], mu, k0, b));

        let (k, r) = zl.form_tangent_and_resistance(&node_i, &node_j);

        let yield_force = mu * 1000.0;
        let expected_v_uy = b * (k0 * 5.0) + (1.0 - b) * yield_force;
        let expected_v_uz = k0 * 0.1;

        assert!((r[SPATIAL_NDF + 1] - expected_v_uy).abs() < 1e-9, "Uy slides independently");
        assert!((r[SPATIAL_NDF + 2] - expected_v_uz).abs() < 1e-9, "Uz sticks independently");
        assert!((k[(1, 1)] - b * k0).abs() < 1e-9, "Uy tangent == b*k0 (sliding)");
        assert!((k[(2, 2)] - k0).abs() < 1e-9, "Uz tangent == k0 (sticking)");
        assert_eq!(k[(1, 2)], 0.0, "no cross-talk between the two shear axes");
        assert_eq!(k[(2, 1)], 0.0, "no cross-talk between the two shear axes");
    }
}
