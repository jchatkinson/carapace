use nalgebra::{Matrix3, SMatrix, SVector, Vector3};

use super::super::{BeamIntegration, Fiber, Fiber3, FiberSection, FiberSection3, GeomTransf3, Node, Node3, Node3Id, NodeId};
use super::truss::{SpatialElementMatrix, SpatialElementVector};

/// A 2-node, force-based (flexibility-method) 2D beam-column (§3.1, M8):
/// unlike every other element in the catalog, its section forces are
/// determined directly from the element's basic (reduced, rigid-body-mode-
/// free) end forces via statics-exact force-interpolation functions —
/// always in equilibrium, regardless of material state — while *section*
/// deformations (and hence the element's tangent) come from inverting each
/// section's nonlinear force-deformation response. That inversion is why
/// this element needs its own internal Newton iteration
/// (`state_determination`) every time `form_tangent_and_resistance` is
/// called, instead of a single direct evaluation like `DispBeamColumn`.
///
/// **Basic system.** 3 dof, rigid-body-mode-free: `v = [v1, v2, v3]` where
/// `v1` is axial elongation and `v2`/`v3` are the end rotations *relative to
/// the chord* (`θ_i - ψ`, `θ_j - ψ`, `ψ = (v2_local - v1_local)/L` the
/// chord's rigid-rotation) — the same objective quantity `GeomTransf`
/// implicitly builds into `ElasticBeamColumn`'s closed-form stiffness, made
/// explicit here because the flexibility method needs an isolated,
/// statically-determinate reduced system to define force-interpolation
/// functions on. Basic forces `q = [q1, q2, q3]` are q1 = axial force, and
/// `q2`/`q3` are *not* simply the local nodal moments — the force
/// interpolation `b(x)` below (and the `a_matrix`/`a_matrix_transpose` map
/// between basic and local systems) were derived to be consistent with
/// `DispBeamColumn`'s already-verified curvature-displacement convention
/// (`FiberSection`'s `M = -sum(stress*area*y)`, `ElasticBeamColumn`'s
/// closed-form `4EI/L`/`2EI/L` pattern) and checked against
/// `ElasticBeamColumn`'s exact elastic stiffness in
/// `core/tests/m8_force_beam_column.rs` — that test is the actual source of
/// truth for these signs, not a textbook convention taken on faith.
///
/// **State determination algorithm** — ported directly from
/// `xara/SRC/element/Frame/Other/Force/ForceBeamColumn2d.cpp`'s
/// `update()` (the Neuenhofer-Filippou/Spacone-Ciampi-Filippou algorithm
/// that file's header cites), not re-derived by hand; see
/// [`state_determination`](ForceBeamColumn::state_determination) and
/// [`try_state_determination`](ForceBeamColumn::try_state_determination)'s
/// own doc comments for the mechanics (a two-stage per-section correction
/// each outer iteration, an energy-based convergence check, and a
/// bisection fallback when a single jump from `v_commit` to the target
/// `v` doesn't converge directly) — done fresh on every call from
/// `q_commit`/`e_commit`/`v_commit` (see those fields' doc comments for
/// why there's no mutable trial cache). At convergence, the basic tangent
/// stiffness is `F^-1` (`F` the element flexibility `integral(b^T * f_s *
/// b) dx`) — the standard result for force-based elements, from
/// differentiating the equilibrium `r(q, v) = 0` implicit relation.
///
/// `GeomTransf::Linear` only, no element loads, `BeamIntegration::Lobatto`
/// the usual choice (endpoints included — plastic hinges concentrate at
/// member ends) — same M7-deferred scope as `DispBeamColumn` for both;
/// nothing here needs them yet.
#[derive(Debug, Clone)]
pub struct ForceBeamColumn {
    pub node_i: NodeId,
    pub node_j: NodeId,
    integration: BeamIntegration,
    sections: Vec<FiberSection>,
    /// Mass per unit volume, applied uniformly over the section's total
    /// fiber area, same convention as `DispBeamColumn`. Zero (the default)
    /// means massless.
    pub density: f64,
    /// Last-committed basic force `[q1, q2, q3]` — the state-determination
    /// Newton iteration's starting guess for `q`, refreshed by `commit`.
    q_commit: SVector<f64, 3>,
    /// Last-committed per-section deformation `(eps0, kappa)` — the
    /// starting guess for each section's `e_s`, refreshed by `commit`.
    /// Needed because `FiberSection` has no "given force, find deformation"
    /// operation of its own (only the reverse); the element itself must
    /// track the deformation iterate that state determination corrects.
    e_commit: Vec<(f64, f64)>,
    /// Last-committed basic deformation — the endpoint `q_commit`/
    /// `e_commit` are consistent with, needed so `state_determination` can
    /// subdivide a too-large jump to `v` into smaller intermediate targets
    /// (see its doc comment) rather than only its own start and end point.
    v_commit: SVector<f64, NBD>,
    max_iters: usize,
    tolerance: f64,
}

/// `2` — number of section force/deformation components (`N`, `M`; §3.1's
/// 2D scope, no shear/torsion).
const NSD: usize = 2;
/// `3` — number of basic (reduced) element dof.
const NBD: usize = 3;

/// `(basic force, per-section (eps0, kappa), basic tangent stiffness)` —
/// [`ForceBeamColumn::state_determination`] and
/// [`ForceBeamColumn::try_state_determination`]'s common result shape.
type StateDeterminationResult = (SVector<f64, NBD>, Vec<(f64, f64)>, SMatrix<f64, NBD, NBD>);

impl ForceBeamColumn {
    pub fn new(
        node_i: NodeId,
        node_j: NodeId,
        fibers: Vec<Fiber>,
        integration: BeamIntegration,
    ) -> Self {
        let n_points = integration.points().len();
        let sections = (0..n_points)
            .map(|_| FiberSection::new(fibers.clone()))
            .collect();
        ForceBeamColumn {
            node_i,
            node_j,
            integration,
            sections,
            density: 0.0,
            q_commit: SVector::<f64, NBD>::zeros(),
            e_commit: vec![(0.0, 0.0); n_points],
            v_commit: SVector::<f64, NBD>::zeros(),
            max_iters: 50,
            tolerance: 1e-6,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
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

    fn local_displacement(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (f64, SMatrix<f64, 6, 6>, SVector<f64, 6>) {
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

    /// Basic (rigid-body-mode-free) deformation `v = A * d_local` — see the
    /// type doc comment for `A`'s definition (chord-relative end rotations).
    fn basic_deformation_matrix(length: f64) -> SMatrix<f64, NBD, 6> {
        let l_inv = 1.0 / length;
        #[rustfmt::skip]
        let a = SMatrix::<f64, NBD, 6>::from_row_slice(&[
            -1.0,    0.0, 0.0, 1.0,    0.0, 0.0,
             0.0, l_inv, 1.0, 0.0, -l_inv, 0.0,
             0.0, l_inv, 0.0, 0.0, -l_inv, 1.0,
        ]);
        a
    }

    /// Force interpolation `s(x) = b(xi) * q`, `s = [N, M]`, `xi` in
    /// `[0,1]` along the length — see the type doc comment for the sign
    /// convention (`M(0) = -q2`, `M(1) = q3`, matched to `DispBeamColumn`'s
    /// curvature convention and verified by
    /// `core/tests/m8_force_beam_column.rs`).
    fn b_matrix(xi: f64) -> SMatrix<f64, NSD, NBD> {
        #[rustfmt::skip]
        let b = SMatrix::<f64, NSD, NBD>::from_row_slice(&[
            1.0,      0.0,  0.0,
            0.0, xi - 1.0,   xi,
        ]);
        b
    }

    fn invert_2x2(m: [[f64; 2]; 2]) -> SMatrix<f64, NSD, NSD> {
        let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        SMatrix::<f64, NSD, NSD>::new(m[1][1] / det, -m[0][1] / det, -m[1][0] / det, m[0][0] / det)
    }

    /// Outer driver: try [`try_state_determination`]'s single-shot Newton
    /// solve straight from `v_commit` to `v`; if it doesn't converge,
    /// subdivide the gap into progressively more, smaller intermediate
    /// targets (each sub-step starting from the previous one's *converged*
    /// state) — the source's own `numSubdivide` fallback (`update()`'s
    /// outer `while` loop), simplified to bisection alone (the source also
    /// tries an initial-tangent iteration variant before subdividing; not
    /// needed here since bisection alone is enough to keep every
    /// sub-step's jump small). A too-large jump in one shot can walk a
    /// section briefly through an unreachable, exactly-singular
    /// fully-plastic state on the way to a perfectly reasonable final
    /// answer — subdividing keeps each step's target close enough to its
    /// start that this doesn't happen.
    fn state_determination(&self, v: SVector<f64, NBD>, length: f64) -> StateDeterminationResult {
        let mut divisions = 1;
        loop {
            let mut q = self.q_commit;
            let mut e = self.e_commit.clone();
            let mut k_basic = SMatrix::<f64, NBD, NBD>::identity();
            let mut all_converged = true;

            for step in 1..=divisions {
                let v_step = self.v_commit + (v - self.v_commit) * (step as f64 / divisions as f64);
                let ((q2, e2, k2), converged) = self.try_state_determination(v_step, length, q, e);
                q = q2;
                e = e2;
                k_basic = k2;
                if !converged {
                    all_converged = false;
                    break;
                }
            }

            if all_converged || divisions >= 256 {
                return (q, e, k_basic);
            }
            divisions *= 4;
        }
    }

    /// Single-shot Newton solve from `(q0, e0)` to the target `v` —
    /// ported directly from
    /// `xara/SRC/element/Frame/Other/Force/ForceBeamColumn2d.cpp`'s
    /// `update()` (the Neuenhofer-Filippou/Spacone-Ciampi-Filippou
    /// algorithm that file's header cites), not re-derived by hand.
    /// Simplified from the source: no distributed element loads (§3.1
    /// scope), no initial-tangent iteration variant (the source's
    /// `Newton`/`InitialFirst`/`Initial` algorithm ladder — an additional
    /// robustness measure on top of subdivision, not needed for
    /// correctness, addable later if a model actually demands it).
    ///
    /// Each outer iteration corrects every section's deformation `e_s`
    /// *twice* against the target section force `b(xi_s) * q`, both using
    /// `FiberSection::trial` (never mutating): once to advance `e_s`
    /// itself (`vsSubdivide[i] += dvs` in the source), and again — after
    /// re-evaluating the section at that advanced `e_s` — to get the
    /// "leftover" unbalance folded into the compatibility integral (`vr`)
    /// and flexibility (`f`) without being written back into `e_s`. This
    /// distinction is the actual fix for a bug an earlier, hand-derived
    /// version of this function had: correcting `e_s` only *once* per
    /// outer iteration and checking convergence via the bare compatibility
    /// residual `v - integral(b^T * e)` is a trap — that residual is
    /// driven to *exactly* zero after a single outer iteration by the
    /// linear algebra of the update alone, regardless of whether any
    /// section is actually in equilibrium with `q`, so the loop
    /// "converges" immediately every time and silently freezes `q`/`e` at
    /// a materially wrong, too-stiff answer once any section is nonlinear
    /// (caught by `core/tests/m8_force_beam_column.rs`'s past-yield test).
    /// The source's two-stage correction plus its energy-based convergence
    /// check (`dv . dSe`, not the raw residual norm) avoid that trap;
    /// match its structure, not a simplification of it.
    ///
    /// Returns whether it actually converged within `self.max_iters`;
    /// [`state_determination`] is what decides what to do if not.
    fn try_state_determination(
        &self,
        v: SVector<f64, NBD>,
        length: f64,
        q0: SVector<f64, NBD>,
        e0: Vec<(f64, f64)>,
    ) -> (StateDeterminationResult, bool) {
        let points = self.integration.points();
        let mut q = q0;
        let mut e = e0;
        let mut k_basic = SMatrix::<f64, NBD, NBD>::identity();
        let mut converged = false;

        for _ in 0..self.max_iters {
            let mut flexibility = SMatrix::<f64, NBD, NBD>::zeros();
            let mut vr = SVector::<f64, NBD>::zeros();

            for (i, (xi, w)) in points.iter().enumerate() {
                let b = Self::b_matrix(*xi);
                let target = b * q;
                let scale = w * length;

                // First correction: advance this section's own stored
                // deformation estimate toward the target, using its
                // tangent from *before* this correction (matches the
                // source's first `dvs.addMatrixVector(..fsSubdivide[i]..)`
                // before `setTrialSectionDeformation`).
                let (n0, m0, k_sec0) = self.sections[i].trial(e[i].0, e[i].1);
                let f_sec0 = Self::invert_2x2(k_sec0);
                let d0 = f_sec0 * (target - SVector::<f64, NSD>::new(n0, m0));
                e[i].0 += d0[0];
                e[i].1 += d0[1];

                // Re-evaluate at the advanced deformation (the source's
                // `setTrialSectionDeformation` + `getStressResultant`/
                // `getSectionFlexibility`) — this is the section's
                // genuine nonlinear response, generally still not
                // exactly matching `target`.
                let (n1, m1, k_sec1) = self.sections[i].trial(e[i].0, e[i].1);
                let f_sec1 = Self::invert_2x2(k_sec1);

                // Second correction, using the *fresh* tangent against
                // the *fresh* response — folded into `vr`/`flexibility`
                // only, never written back into `e[i]` (matches the
                // source computing this into a fresh local `dvs`, added
                // to `vsSubdivide[i]` only for the `vr` integration).
                let d1 = f_sec1 * (target - SVector::<f64, NSD>::new(n1, m1));
                let e_star = SVector::<f64, NSD>::new(e[i].0 + d1[0], e[i].1 + d1[1]);

                flexibility += scale * (b.transpose() * f_sec1 * b);
                vr += scale * (b.transpose() * e_star);
            }

            k_basic = match flexibility.try_inverse() {
                Some(k) if k.iter().all(|x| x.is_finite()) => k,
                _ => break,
            };
            let dv = v - vr;
            let dq = k_basic * dv;
            if !dq.iter().all(|x| x.is_finite()) {
                break;
            }
            q += dq;

            if dv.dot(&dq).abs() < self.tolerance {
                converged = true;
                break;
            }
        }

        ((q, e, k_basic), converged)
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        let (length, t, d_local) = self.local_displacement(node_i, node_j);
        let a = Self::basic_deformation_matrix(length);
        let v = a * d_local;

        let (q, _e, k_basic) = self.state_determination(v, length);

        let k_local = a.transpose() * k_basic * a;
        let r_local = a.transpose() * q;

        (t.transpose() * k_local * t, t.transpose() * r_local)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        let a = Self::basic_deformation_matrix(length);
        let v = a * d_local;

        let (q, e, _k_basic) = self.state_determination(v, length);

        for (i, (eps0, kappa)) in e.iter().enumerate() {
            self.sections[i].commit(*eps0, *kappa);
        }
        self.q_commit = q;
        self.e_commit = e;
        self.v_commit = v;
    }

    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let total_area = self.sections[0].total_area();
        let half = self.density * total_area * length / 2.0;
        SVector::<f64, 6>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}

/// `3` — number of section force/deformation components for
/// `ForceBeamColumn3` (`N`, `Mz`, `My` — biaxial bending, no torsion; see
/// `FiberSection3`'s doc comment for why).
const NSD3: usize = 3;
/// `5` — number of basic (reduced) dof for `ForceBeamColumn3`: axial
/// elongation plus two chord-relative end rotations in *each* bending
/// plane. Torsion is excluded from the basic system entirely (unlike
/// `NSD3`'s `N`/`Mz`/`My`, there's no `q6`/`v6` pair) — it's a decoupled
/// elastic term folded onto the local tangent/resistance directly, the
/// same treatment `DispBeamColumn3` uses.
const NBD3: usize = 5;

/// `(basic force, per-section (eps0, kappa_z, kappa_y), basic tangent
/// stiffness)` — `ForceBeamColumn3`'s counterpart to `StateDeterminationResult`.
type StateDeterminationResult3 = (SVector<f64, NBD3>, Vec<(f64, f64, f64)>, SMatrix<f64, NBD3, NBD3>);

/// `ForceBeamColumn`'s spatial (biaxial) counterpart: a 2-node, force-based
/// 3D beam-column. Basic system, force interpolation, and state-
/// determination algorithm are the direct biaxial generalization of
/// `ForceBeamColumn`'s (see that type's doc comment for the algorithm
/// itself, ported unchanged in structure — only `NSD`/`NBD` and the
/// matrices below grow); torsion is excluded from the basic system and
/// handled as a decoupled elastic `G*J/L` term, same as `DispBeamColumn3`.
///
/// **Basic system**, `v = [v1..v5]`: `v1` axial elongation (same as
/// `ForceBeamColumn`'s `v1`); `v2`/`v3` the z-bending chord-relative end
/// rotations, numerically identical in form to `ForceBeamColumn`'s
/// `v2`/`v3` since local z-bending uses the same `[v, rz]` convention;
/// `v4`/`v5` the y-bending chord-relative end rotations, built by
/// substituting `theta_equiv = -ry` into the same construction (matching
/// `DispBeamColumn3::strain_displacement`'s `b_kappa_y` substitution) —
/// this flips the sign of `v4`/`v5`'s *rotation* coefficient relative to
/// `v2`/`v3`'s (their translation coefficients are unchanged, since the
/// chord slope `psi` is purely geometric and doesn't depend on the
/// rotation sign convention). Basic forces `q = [q1..q5]` are then
/// automatically work-conjugate to `v` by construction (`A`/`A^T` share the
/// same matrix — the standard force-method contragredience), with `q1=N`,
/// `(q2,q3)` conjugate to `(Mz`'s end values`)`, `(q4,q5)` conjugate to
/// `My`'s. This whole construction — not just guessed by analogy — is
/// verified against `ElasticBeamColumn3`'s exact biaxial-bending elastic
/// stiffness in `core/tests/m17_force_beam_column3.rs`, the same role
/// `core/tests/m8_force_beam_column.rs` plays for the planar element.
#[derive(Debug, Clone)]
pub struct ForceBeamColumn3 {
    pub node_i: Node3Id,
    pub node_j: Node3Id,
    pub g: f64,
    pub j: f64,
    vec_xz: [f64; 3],
    integration: BeamIntegration,
    sections: Vec<FiberSection3>,
    /// Mass per unit volume, same convention as `ForceBeamColumn::density`.
    pub density: f64,
    q_commit: SVector<f64, NBD3>,
    e_commit: Vec<(f64, f64, f64)>,
    v_commit: SVector<f64, NBD3>,
    max_iters: usize,
    tolerance: f64,
}

impl ForceBeamColumn3 {
    #[allow(clippy::too_many_arguments)]
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
        ForceBeamColumn3 {
            node_i,
            node_j,
            g,
            j,
            vec_xz,
            integration,
            sections,
            density: 0.0,
            q_commit: SVector::<f64, NBD3>::zeros(),
            e_commit: vec![(0.0, 0.0, 0.0); n_points],
            v_commit: SVector::<f64, NBD3>::zeros(),
            max_iters: 50,
            tolerance: 1e-6,
        }
    }

    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }

    fn local_displacement(&self, node_i: &Node3, node_j: &Node3) -> (f64, SpatialElementMatrix, SpatialElementVector) {
        let transform = GeomTransf3::linear(self.vec_xz);
        let (length, r) = transform.local_axes(node_i, node_j);
        let t = GeomTransf3::rotation_matrix(&r);
        let d_global = SpatialElementVector::from_iterator(node_i.displacement.iter().chain(node_j.displacement.iter()).copied());
        (length, t, t * d_global)
    }

    fn torsion_stiffness(&self, length: f64) -> SpatialElementMatrix {
        let gj_l = self.g * self.j / length;
        let mut k = SpatialElementMatrix::zeros();
        k[(3, 3)] = gj_l;
        k[(9, 9)] = gj_l;
        k[(3, 9)] = -gj_l;
        k[(9, 3)] = -gj_l;
        k
    }

    /// Basic deformation `v = A * d_local` — see the type doc comment for
    /// `A`'s rows. `d_local` is `[u,v,w,rx,ry,rz]` per node; rows here skip
    /// the torsion columns (`3`, `9`) entirely, since torsion isn't part of
    /// the basic system.
    fn basic_deformation_matrix(length: f64) -> SMatrix<f64, NBD3, 12> {
        let l_inv = 1.0 / length;
        let mut a = SMatrix::<f64, NBD3, 12>::zeros();
        // v1: axial.
        a[(0, 0)] = -1.0;
        a[(0, 6)] = 1.0;
        // v2, v3: z-bending ([v, rz] at indices 1,5,7,11).
        a[(1, 1)] = l_inv;
        a[(1, 5)] = 1.0;
        a[(1, 7)] = -l_inv;
        a[(2, 1)] = l_inv;
        a[(2, 7)] = -l_inv;
        a[(2, 11)] = 1.0;
        // v4, v5: y-bending ([w, ry] at indices 2,4,8,10), rotation
        // coefficient sign-flipped relative to z-bending (see doc comment).
        a[(3, 2)] = l_inv;
        a[(3, 4)] = -1.0;
        a[(3, 8)] = -l_inv;
        a[(4, 2)] = l_inv;
        a[(4, 8)] = -l_inv;
        a[(4, 10)] = -1.0;
        a
    }

    /// Force interpolation `s(x) = b(xi) * q`, `s = [N, Mz, My]` — `Mz`'s
    /// row is `ForceBeamColumn::b_matrix`'s unchanged (same `[v, rz]`
    /// convention); `My`'s row has the identical `[xi-1, xi]` form against
    /// `q4`/`q5`, with no further sign change needed — the rotation-sign
    /// flip is already fully absorbed into `v4`/`v5`'s definition, and
    /// work-conjugacy (`q` conjugate to `v`, `My` conjugate to `kappa_y`,
    /// both by the same construction) carries it through consistently.
    fn b_matrix(xi: f64) -> SMatrix<f64, NSD3, NBD3> {
        let mut b = SMatrix::<f64, NSD3, NBD3>::zeros();
        b[(0, 0)] = 1.0;
        b[(1, 1)] = xi - 1.0;
        b[(1, 2)] = xi;
        b[(2, 3)] = xi - 1.0;
        b[(2, 4)] = xi;
        b
    }

    /// See `ForceBeamColumn::state_determination`'s doc comment — identical
    /// bisection-subdivision structure, `NSD3`/`NBD3` in place of `NSD`/`NBD`.
    fn state_determination(&self, v: SVector<f64, NBD3>, length: f64) -> StateDeterminationResult3 {
        let mut divisions = 1;
        loop {
            let mut q = self.q_commit;
            let mut e = self.e_commit.clone();
            let mut k_basic = SMatrix::<f64, NBD3, NBD3>::identity();
            let mut all_converged = true;

            for step in 1..=divisions {
                let v_step = self.v_commit + (v - self.v_commit) * (step as f64 / divisions as f64);
                let ((q2, e2, k2), converged) = self.try_state_determination(v_step, length, q, e);
                q = q2;
                e = e2;
                k_basic = k2;
                if !converged {
                    all_converged = false;
                    break;
                }
            }

            if all_converged || divisions >= 256 {
                return (q, e, k_basic);
            }
            divisions *= 4;
        }
    }

    /// See `ForceBeamColumn::try_state_determination`'s doc comment —
    /// identical two-stage-per-section-correction, energy-based-convergence
    /// algorithm, `NSD3`/`NBD3` in place of `NSD`/`NBD` and a 3x3 section
    /// flexibility (`nalgebra`'s `try_inverse`) in place of the hand-rolled
    /// `invert_2x2`.
    fn try_state_determination(
        &self,
        v: SVector<f64, NBD3>,
        length: f64,
        q0: SVector<f64, NBD3>,
        e0: Vec<(f64, f64, f64)>,
    ) -> (StateDeterminationResult3, bool) {
        let points = self.integration.points();
        let mut q = q0;
        let mut e = e0;
        let mut k_basic = SMatrix::<f64, NBD3, NBD3>::identity();
        let mut converged = false;

        for _ in 0..self.max_iters {
            let mut flexibility = SMatrix::<f64, NBD3, NBD3>::zeros();
            let mut vr = SVector::<f64, NBD3>::zeros();

            for (i, (xi, w)) in points.iter().enumerate() {
                let b = Self::b_matrix(*xi);
                let target = b * q;
                let scale = w * length;

                let (n0, mz0, my0, k_sec0) = self.sections[i].trial(e[i].0, e[i].1, e[i].2);
                let f_sec0 = Matrix3::from_row_slice(&[
                    k_sec0[0][0], k_sec0[0][1], k_sec0[0][2], k_sec0[1][0], k_sec0[1][1], k_sec0[1][2], k_sec0[2][0], k_sec0[2][1],
                    k_sec0[2][2],
                ])
                .try_inverse()
                .expect("fiber section tangent should be nonsingular for a real cross-section");
                let d0 = f_sec0 * (target - Vector3::new(n0, mz0, my0));
                e[i].0 += d0[0];
                e[i].1 += d0[1];
                e[i].2 += d0[2];

                let (n1, mz1, my1, k_sec1) = self.sections[i].trial(e[i].0, e[i].1, e[i].2);
                let f_sec1 = Matrix3::from_row_slice(&[
                    k_sec1[0][0], k_sec1[0][1], k_sec1[0][2], k_sec1[1][0], k_sec1[1][1], k_sec1[1][2], k_sec1[2][0], k_sec1[2][1],
                    k_sec1[2][2],
                ])
                .try_inverse()
                .expect("fiber section tangent should be nonsingular for a real cross-section");

                let d1 = f_sec1 * (target - Vector3::new(n1, mz1, my1));
                let e_star = Vector3::new(e[i].0 + d1[0], e[i].1 + d1[1], e[i].2 + d1[2]);

                flexibility += scale * (b.transpose() * f_sec1 * b);
                vr += scale * (b.transpose() * e_star);
            }

            k_basic = match flexibility.try_inverse() {
                Some(k) if k.iter().all(|x| x.is_finite()) => k,
                _ => break,
            };
            let dv = v - vr;
            let dq = k_basic * dv;
            if !dq.iter().all(|x| x.is_finite()) {
                break;
            }
            q += dq;

            if dv.dot(&dq).abs() < self.tolerance {
                converged = true;
                break;
            }
        }

        ((q, e, k_basic), converged)
    }

    pub(super) fn form_tangent_and_resistance(&self, node_i: &Node3, node_j: &Node3) -> (SpatialElementMatrix, SpatialElementVector) {
        let (length, t, d_local) = self.local_displacement(node_i, node_j);
        let a = Self::basic_deformation_matrix(length);
        let v = a * d_local;

        let (q, _e, k_basic) = self.state_determination(v, length);

        let k_torsion = self.torsion_stiffness(length);
        let k_local = a.transpose() * k_basic * a + k_torsion;
        let r_local = a.transpose() * q + k_torsion * d_local;

        (t.transpose() * k_local * t, t.transpose() * r_local)
    }

    pub(super) fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        let a = Self::basic_deformation_matrix(length);
        let v = a * d_local;

        let (q, e, _k_basic) = self.state_determination(v, length);

        for (i, (eps0, kappa_z, kappa_y)) in e.iter().enumerate() {
            self.sections[i].commit(*eps0, *kappa_z, *kappa_y);
        }
        self.q_commit = q;
        self.e_commit = e;
        self.v_commit = v;
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
}
