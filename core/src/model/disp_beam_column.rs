use nalgebra::{SMatrix, SVector};

use super::{BeamIntegration, Fiber, FiberSection, Node, NodeId};

/// A 2-node, displacement-based, fiber-discretized 2D beam-column (§4.1):
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
/// `GeomTransf::Linear` only — small-displacement, no `PDelta` geometric-
/// stiffness correction (unlike `ElasticBeamColumn`, whose axial force is a
/// single scalar; a fiber section's isn't nearly as direct to plug into a
/// geometric-stiffness formula, so this is deferred, not silently wrong).
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
            density: 0.0,
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
        let (length, t, d_local) = self.local_displacement(node_i, node_j);

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

        (t.transpose() * k_local * t, t.transpose() * r_local)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let (length, _t, d_local) = self.local_displacement(node_i, node_j);
        let points = self.integration.points();
        for ((xi, _w), section) in points.iter().zip(&mut self.sections) {
            let (b_eps0, b_kappa) = Self::strain_displacement(*xi, length);
            let eps0 = b_eps0.dot(&d_local);
            let kappa = b_kappa.dot(&d_local);
            section.commit(eps0, kappa);
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
}
