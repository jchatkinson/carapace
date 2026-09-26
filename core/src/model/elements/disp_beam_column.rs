use nalgebra::{SMatrix, SVector};

use super::super::transform::Corotational2d;
use super::super::{
    BeamIntegration, Fiber, Fiber3, FiberSection, FiberSection3, GeomTransf, GeomTransf3, Node,
    Node3, Node3Id, NodeId,
};
use super::truss::{SpatialElementMatrix, SpatialElementVector};

/// A 2-node, displacement-based, fiber-discretized 2D beam-column (§3.1):
/// nodal displacements directly give the strain/curvature field along the
/// length (cubic Hermite transverse + linear axial shape functions, same
/// as `ElasticBeamColumn`'s), each `BeamIntegration` point's `FiberSection`
/// converts that to stress-resultants, and one length-integral gives the
/// element's resisting force/tangent — a single, direct evaluation, unlike
/// `ForceBeamColumn` (M8), which needs its own internal equilibrium
/// iteration per element.
///
/// Prismatic-member assumption: every integration point starts from the
/// same fiber layout (`fibers`, replicated `points` times) but evolves
/// independently thereafter, since each point's material history depends
/// on that point's own strain path. A tapered/non-prismatic member (a
/// different section per point) is a natural extension, not built until
/// needed.
///
/// Small-displacement by default, with an opt-in `GeomTransf::Corotational`
/// formulation that integrates fiber response from the three objective
/// basic deformations. `PDelta` remains unsupported for this element (a
/// fiber section's force state does not reduce to the scalar axial force
/// used by `ElasticBeamColumn`'s closed-form correction).
/// No element loads (`ElasticBeamColumn`'s uniform transverse load) either
/// — a fiber element's consistent load vector needs the same per-section
/// integration machinery core to this element, not a closed-form formula,
/// so it's deferred until `pysees` actually needs distributed loads on a
/// fiber-section member.
#[derive(Debug, Clone)]
pub struct DispBeamColumn {
    pub node_i: NodeId,
    pub node_j: NodeId,
    integration: BeamIntegration,
    sections: Vec<FiberSection>,
    transform: GeomTransf,
    /// Mass per unit volume, applied uniformly over the section's total
    /// fiber area. Zero (the default) means massless.
    pub density: f64,
}

impl DispBeamColumn {
    pub fn new(node_i: NodeId, node_j: NodeId, fibers: Vec<Fiber>, integration: BeamIntegration) -> Self {
        let n_points = integration.points().len();
        let sections = (0..n_points).map(|_| FiberSection::new(fibers.clone())).collect();
        DispBeamColumn {
            node_i,
            node_j,
            integration,
            sections,
            transform: GeomTransf::Linear,
            density: 0.0,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Use finite-rotation, chord-following kinematics for this member.
    pub fn with_corotational(mut self) -> Self {
        self.transform = GeomTransf::Corotational;
        self
    }

    fn geometry(&self, node_i: &Node, node_j: &Node) -> (f64, f64, f64) {
        let dx = node_j.coords[0] - node_i.coords[0];
        let dy = node_j.coords[1] - node_i.coords[1];
        let length = (dx * dx + dy * dy).sqrt();
        (length, dx / length, dy / length)
    }

    fn transformation(&self, cx: f64, cy: f64) -> SMatrix<f64, 6, 6> {
        #[rustfmt::skip]
        let block = SMatrix::<f64, 3, 3>::new(
             cx,  cy, 0.0,
            -cy,  cx, 0.0,
            0.0, 0.0, 1.0,
        );
        let mut t = SMatrix::<f64, 6, 6>::zeros();
        t.fixed_view_mut::<3, 3>(0, 0).copy_from(&block);
        t.fixed_view_mut::<3, 3>(3, 3).copy_from(&block);
        t
    }

    /// Strain-displacement vectors `(b_eps0, b_kappa)` at `xi` in `[0,1]`
    /// along the length, relating local nodal displacements
    /// `[u1,v1,th1,u2,v2,th2]` to `[eps0, kappa]` — `b_eps0` from the
    /// linear axial shape functions (constant), `b_kappa` the second
    /// derivative of the same cubic Hermite transverse shape functions
    /// `ElasticBeamColumn` uses, derived directly (not copied from a
    /// reference — this is standard Euler-Bernoulli/Hermite beam theory,
    /// verified against `ElasticBeamColumn`'s closed-form stiffness in
    /// `core/tests/m7_disp_beam_column.rs`: for constant `EA`/`EI` a rule
    /// of 2 or more points must reproduce it exactly, since `b_kappa` is
    /// linear in `xi` and its square is then only degree 2).
    fn strain_displacement(xi: f64, length: f64) -> (SVector<f64, 6>, SVector<f64, 6>) {
        let l = length;
        let b_eps0 = SVector::<f64, 6>::from_column_slice(&[-1.0 / l, 0.0, 0.0, 1.0 / l, 0.0, 0.0]);
        #[rustfmt::skip]
        let b_kappa = SVector::<f64, 6>::from_column_slice(&[
            0.0,
            (-6.0 + 12.0 * xi) / (l * l),
            (-4.0 + 6.0 * xi) / l,
            0.0,
            (6.0 - 12.0 * xi) / (l * l),
            (-2.0 + 6.0 * xi) / l,
        ]);
        (b_eps0, b_kappa)
    }

    fn local_displacement(&self, node_i: &Node, node_j: &Node) -> (f64, SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let d_global = SVector::<f64, 6>::from_row_slice(&[
            node_i.displacement[0],
            node_i.displacement[1],
            node_i.displacement[2],
            node_j.displacement[0],
            node_j.displacement[1],
            node_j.displacement[2],
        ]);
        (length, t, t * d_global)
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        let corotational = (self.transform == GeomTransf::Corotational)
            .then(|| Corotational2d::new(node_i, node_j));
        let (length, t, d_local) = if let Some(state) = &corotational {
            (
                state.initial_length(),
                SMatrix::<f64, 6, 6>::identity(),
                Corotational2d::local_basic_modes() * state.basic_deformation(),
            )
        } else {
            self.local_displacement(node_i, node_j)
        };

        let mut k_local = SMatrix::<f64, 6, 6>::zeros();
        let mut r_local = SVector::<f64, 6>::zeros();

        for ((xi, w), section) in self.integration.points().iter().zip(&self.sections) {
            let (b_eps0, b_kappa) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa = b_kappa.dot(&d_local);
            let (n, m, k_section) = section.trial(eps0, kappa);

            let scale = w * length;
            r_local += scale * (b_eps0 * n + b_kappa * m);
            k_local += scale
                * (k_section[0][0] * (b_eps0 * b_eps0.transpose())
                    + k_section[0][1] * (b_eps0 * b_kappa.transpose())
                    + k_section[1][0] * (b_kappa * b_eps0.transpose())
                    + k_section[1][1] * (b_kappa * b_kappa.transpose()));
        }

        if let Some(state) = &corotational {
            let modes = Corotational2d::local_basic_modes();
            let basic_tangent = modes.transpose() * k_local * modes;
            let basic_resistance = modes.transpose() * r_local;
            (
                state.global_tangent(&basic_tangent, &basic_resistance),
                state.global_resistance(&basic_resistance),
            )
        } else {
            (t.transpose() * k_local * t, t.transpose() * r_local)
        }
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let (length, d_local) = if self.transform == GeomTransf::Corotational {
            let state = Corotational2d::new(node_i, node_j);
            (
                state.initial_length(),
                Corotational2d::local_basic_modes() * state.basic_deformation(),
            )
        } else {
            let (length, _t, d_local) = self.local_displacement(node_i, node_j);
            (length, d_local)
        };
        let points = self.integration.points();
        for ((xi, _w), section) in points.iter().zip(&mut self.sections) {
            let (b_eps0, b_kappa) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa = b_kappa.dot(&d_local);
            section.commit(eps0, kappa);
        }
    }

    /// Local nodal force — the same `r_local`/`basic_resistance` value
    /// `form_tangent_and_resistance` already computes and discards,
    /// recomputed fresh here (cheaper than that method, in fact: no
    /// tangent accumulation needed). Always exactly consistent with the
    /// current committed section state — no cache to keep in sync, since
    /// this element has no Newton loop between a target and the committed
    /// answer the way `ForceBeamColumn` does.
    pub(super) fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        let corotational = (self.transform == GeomTransf::Corotational).then(|| Corotational2d::new(node_i, node_j));
        let (length, d_local) = if let Some(state) = &corotational {
            (state.initial_length(), Corotational2d::local_basic_modes() * state.basic_deformation())
        } else {
            let (length, _t, d_local) = self.local_displacement(node_i, node_j);
            (length, d_local)
        };

        let mut r_local = SVector::<f64, 6>::zeros();
        for ((xi, w), section) in self.integration.points().iter().zip(&self.sections) {
            let (b_eps0, b_kappa) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa = b_kappa.dot(&d_local);
            let (n, m, _k_section) = section.trial(eps0, kappa);
            r_local += (w * length) * (b_eps0 * n + b_kappa * m);
        }

        // Non-corotational: `r_local` is already the local nodal force
        // (`d_local` came from the fixed-orientation `t`, same as
        // `form_tangent_and_resistance`'s pre-transform value) — no further
        // rotation needed, unlike that method's `t.transpose() * r_local`
        // step, which converts *to* global.
        match &corotational {
            Some(state) => state.local_resistance(&(Corotational2d::local_basic_modes().transpose() * r_local)),
            None => r_local,
        }
    }

    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        // Prismatic assumption (see the type doc comment): every section
        // has the same fiber layout, so any one of them gives the total
        // area.
        let total_area = self.sections[0].total_area();
        let half = self.density * total_area * length / 2.0;
        SVector::<f64, 6>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }

    /// Every integration point's per-fiber `(strain, stress)`, freshly
    /// recomputed from `node_i`/`node_j`'s current committed displacement —
    /// same "recompute, don't cache" reasoning as `local_force` (this
    /// element has no Newton loop between a target and the committed
    /// answer, unlike `ForceBeamColumn`). Outer `Vec` is one entry per
    /// integration point (`BeamIntegration::points()` order), inner `Vec`
    /// one entry per fiber (the order originally passed to `new`).
    pub fn fiber_responses(&self, node_i: &Node, node_j: &Node) -> Vec<Vec<(f64, f64)>> {
        let corotational = (self.transform == GeomTransf::Corotational).then(|| Corotational2d::new(node_i, node_j));
        let (length, d_local) = if let Some(state) = &corotational {
            (state.initial_length(), Corotational2d::local_basic_modes() * state.basic_deformation())
        } else {
            let (length, _t, d_local) = self.local_displacement(node_i, node_j);
            (length, d_local)
        };

        self.integration
            .points()
            .iter()
            .zip(&self.sections)
            .map(|((xi, _w), section)| {
                let (b_eps0, b_kappa) = Self::strain_displacement(*xi, length);
                let eps0 = b_eps0.dot(&d_local);
                let kappa = b_kappa.dot(&d_local);
                section.fiber_responses(eps0, kappa)
            })
            .collect()
    }
}

/// `DispBeamColumn`'s spatial (biaxial) counterpart: a 2-node,
/// displacement-based, fiber-discretized 3D beam-column. Bending in both
/// planes comes directly from `FiberSection3` the same way `DispBeamColumn`
/// gets its single bending plane from `FiberSection`; torsion is *not*
/// fiber-derived (see `FiberSection3`'s doc comment) and is instead a
/// simple decoupled elastic `G*J/L` term, folded directly into the local
/// tangent/resistance alongside the fiber-integrated axial/biaxial-bending
/// terms — the same decoupling `ElasticBeamColumn3` uses for torsion.
///
/// `vec_xz` plus `Node3`'s coordinates give the local axes via the same
/// `GeomTransf3::local_axes`/`rotation_matrix` machinery `ElasticBeamColumn3`
/// uses, but only ever in its `Linear3` sense — no `PDelta3` option is
/// exposed, matching planar `DispBeamColumn`'s "`GeomTransf::Linear` only"
/// scope (see its doc comment: a fiber section's state isn't a single
/// scalar axial force, so it doesn't plug into the closed-form geometric-
/// stiffness formula directly, and building the right generalization is
/// deferred, not silently wrong).
///
/// Local DOF order per node `[u, v, w, rx, ry, rz]`, matching
/// `ElasticBeamColumn3`/`SpatialDof` exactly.
#[derive(Debug, Clone)]
pub struct DispBeamColumn3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    pub g: f64,
    pub j: f64,
    vec_xz: [f64; 3],
    integration: BeamIntegration,
    sections: Vec<FiberSection3>,
    /// Mass per unit volume, same convention as `DispBeamColumn::density`.
    pub density: f64,
}

impl DispBeamColumn3 {
    pub fn new(
        node_i: Node3Id,
        node_j: Node3Id,
        g: f64,
        j: f64,
        vec_xz: [f64; 3],
        fibers: Vec<Fiber3>,
        integration: BeamIntegration,
    ) -> Self {
        let n_points = integration.points().len();
        let sections = (0..n_points).map(|_| FiberSection3::new(fibers.clone())).collect();
        DispBeamColumn3 {
            node_i,
            node_j,
            g,
            j,
            vec_xz,
            integration,
            sections,
            density: 0.0,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    /// Strain-displacement vectors `(b_eps0, b_kappa_z, b_kappa_y)` at `xi`
    /// in `[0,1]`, relating local nodal displacements
    /// `[u1,v1,w1,rx1,ry1,rz1,u2,v2,w2,rx2,ry2,rz2]` to `[eps0, kappa_z,
    /// kappa_y]`. `b_eps0`/`b_kappa_z` are `DispBeamColumn::strain_
    /// displacement`'s `b_eps0`/`b_kappa` re-indexed onto `[u,v,rz]`'s
    /// spatial slots (`0,1,5` / `6,7,11`); `b_kappa_y` is the same Hermite
    /// curvature operator re-derived for `[w, ry]`'s slots (`2,4` / `8,10`)
    /// with its rotation-column signs flipped, from substituting
    /// `theta_equiv = -ry` (the standard-Hermite formula assumes
    /// `theta = dv/dx`, but this element's convention is `ry = -dw/dx`, see
    /// `ElasticBeamColumn3::local_elastic_stiffness`'s doc comment) into the
    /// unmodified `b_kappa_z` formula — the same substitution
    /// `ElasticBeamColumn3::geometric_stiffness`'s `y_block` uses relative
    /// to `z_block`. Verified, not just plausible-by-analogy, in
    /// `core/tests/m17_disp_beam_column3.rs` against `ElasticBeamColumn3`'s
    /// exact closed-form biaxial-bending stiffness.
    fn strain_displacement(xi: f64, length: f64) -> (SpatialElementVector, SpatialElementVector, SpatialElementVector) {
        let l = length;
        let mut b_eps0 = SpatialElementVector::zeros();
        b_eps0[0] = -1.0 / l;
        b_eps0[6] = 1.0 / l;

        let mut b_kappa_z = SpatialElementVector::zeros();
        b_kappa_z[1] = (-6.0 + 12.0 * xi) / (l * l);
        b_kappa_z[5] = (-4.0 + 6.0 * xi) / l;
        b_kappa_z[7] = (6.0 - 12.0 * xi) / (l * l);
        b_kappa_z[11] = (-2.0 + 6.0 * xi) / l;

        let mut b_kappa_y = SpatialElementVector::zeros();
        b_kappa_y[2] = (-6.0 + 12.0 * xi) / (l * l);
        b_kappa_y[4] = (4.0 - 6.0 * xi) / l;
        b_kappa_y[8] = (6.0 - 12.0 * xi) / (l * l);
        b_kappa_y[10] = (2.0 - 6.0 * xi) / l;

        (b_eps0, b_kappa_z, b_kappa_y)
    }

    fn local_displacement(&self, node_i: &Node3, node_j: &Node3) -> (f64, SpatialElementMatrix, SpatialElementVector) {
        let transform = GeomTransf3::linear(self.vec_xz);
        let (length, r) = transform.local_axes(node_i, node_j);
        let t = GeomTransf3::rotation_matrix(&r);
        let d_global = SpatialElementVector::from_iterator(node_i.displacement.iter().chain(node_j.displacement.iter()).copied());
        (length, t, t * d_global)
    }

    /// Decoupled elastic torsion contribution `[rx1, rx2]` — see the type
    /// doc comment for why torsion isn't fiber-derived.
    fn torsion_stiffness(&self, length: f64) -> SpatialElementMatrix {
        let gj_l = self.g * self.j / length;
        let mut k = SpatialElementMatrix::zeros();
        k[(3, 3)] = gj_l;
        k[(9, 9)] = gj_l;
        k[(3, 9)] = -gj_l;
        k[(9, 3)] = -gj_l;
        k
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let (length, t, d_local) = self.local_displacement(node_i, node_j);

        let mut k_local = self.torsion_stiffness(length);
        let mut r_local = k_local * d_local;

        for ((xi, w), section) in self.integration.points().iter().zip(&self.sections) {
            let (b_eps0, b_kappa_z, b_kappa_y) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa_z = b_kappa_z.dot(&d_local);
            let kappa_y = b_kappa_y.dot(&d_local);
            let (n, mz, my, k_section) = section.trial(eps0, kappa_z, kappa_y);

            let scale = w * length;
            r_local += scale * (b_eps0 * n + b_kappa_z * mz + b_kappa_y * my);

            let b = [&b_eps0, &b_kappa_z, &b_kappa_y];
            for (bi, row_i) in b.iter().enumerate() {
                for (bj, row_j) in b.iter().enumerate() {
                    k_local += scale * k_section[bi][bj] * (*row_i * row_j.transpose());
                }
            }
        }

        (t.transpose() * k_local * t, t.transpose() * r_local)
    }

    /// Local nodal force (including the decoupled torsion term) — the same
    /// pre-transform `r_local` value above, recomputed fresh (cheap: no
    /// tangent accumulation needed, and no Newton loop to avoid re-running).
    pub(super) fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        let mut r_local = self.torsion_stiffness(length) * d_local;
        for ((xi, w), section) in self.integration.points().iter().zip(&self.sections) {
            let (b_eps0, b_kappa_z, b_kappa_y) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa_z = b_kappa_z.dot(&d_local);
            let kappa_y = b_kappa_y.dot(&d_local);
            let (n, mz, my, _k_section) = section.trial(eps0, kappa_z, kappa_y);
            r_local += (w * length) * (b_eps0 * n + b_kappa_z * mz + b_kappa_y * my);
        }
        r_local
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        let points = self.integration.points();
        for ((xi, _w), section) in points.iter().zip(&mut self.sections) {
            let (b_eps0, b_kappa_z, b_kappa_y) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa_z = b_kappa_z.dot(&d_local);
            let kappa_y = b_kappa_y.dot(&d_local);
            section.commit(eps0, kappa_z, kappa_y);
        }
    }

    pub(super) fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SpatialElementVector {
        let transform = GeomTransf3::linear(self.vec_xz);
        let (length, _r) = transform.local_axes(node_i, node_j);
        let total_area = self.sections[0].total_area();
        let half = self.density * total_area * length / 2.0;
        let mut mass = SpatialElementVector::zeros();
        for dof in [0, 1, 2, 6, 7, 8] {
            mass[dof] = half;
        }
        mass
    }

    /// See `DispBeamColumn::fiber_responses`'s doc comment — same
    /// contract, biaxial strain field.
    pub fn fiber_responses(&self, node_i: &Node3, node_j: &Node3) -> Vec<Vec<(f64, f64)>> {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        self.integration
            .points()
            .iter()
            .zip(&self.sections)
            .map(|((xi, _w), section)| {
                let (b_eps0, b_kappa_z, b_kappa_y) = Self::strain_displacement(*xi, length);
                let eps0 = b_eps0.dot(&d_local);
                let kappa_z = b_kappa_z.dot(&d_local);
                let kappa_y = b_kappa_y.dot(&d_local);
                section.fiber_responses(eps0, kappa_z, kappa_y)
            })
            .collect()
    }
}
