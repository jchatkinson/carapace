use nalgebra::{SMatrix, SVector};

use super::super::{GeomTransf, Node, NodeId};

/// A 2-node, prismatic, linear-elastic 2D beam-column (Euler-Bernoulli, no
/// shear deformation): closed-form stiffness, no iteration (§3.1) — unlike
/// `Truss`/`ZeroLength`, there's no `Material` dispatch here, since the
/// element's response is fully determined by `e`/`a`/`iz` with no nonlinear
/// stress-strain law to evaluate.
#[derive(Debug, Clone)]
pub struct ElasticBeamColumn {
    pub node_i: NodeId,
    pub node_j: NodeId,
    pub e: f64,
    pub a: f64,
    pub iz: f64,
    pub transform: GeomTransf,
    /// Uniform transverse load (force/length) in the element's local
    /// +y direction, converted to equivalent nodal loads (§3.4). Zero means
    /// no element load. A single field, not a `Vec<ElementLoad>` — M3's
    /// scope is exactly this one load case; generalize only when a second
    /// element-load type is actually needed.
    pub w_transverse: f64,
    /// Mass per unit volume. Zero (the default) means massless — existing
    /// models are unaffected unless they opt in via `with_density`.
    pub density: f64,
}

impl ElasticBeamColumn {
    pub fn new(node_i: NodeId, node_j: NodeId, e: f64, a: f64, iz: f64, transform: GeomTransf) -> Self {
        ElasticBeamColumn {
            node_i,
            node_j,
            e,
            a,
            iz,
            transform,
            w_transverse: 0.0,
            density: 0.0,
        }
    }

    pub fn with_uniform_load(mut self, w_transverse: f64) -> Self {
        self.w_transverse = w_transverse;
        self
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

    /// Local-DOF-order `[u1, v1, θ1, u2, v2, θ2]` elastic stiffness —
    /// standard Euler-Bernoulli prismatic beam (Przemieniecki / Cook §2).
    fn local_elastic_stiffness(&self, length: f64) -> SMatrix<f64, 6, 6> {
        let (e, a, i) = (self.e, self.a, self.iz);
        let l = length;
        let ea_l = e * a / l;
        let ei = e * i;

        #[rustfmt::skip]
        let k = SMatrix::<f64, 6, 6>::from_row_slice(&[
             ea_l,             0.0,             0.0,  -ea_l,             0.0,             0.0,
              0.0,  12.0*ei/l.powi(3),  6.0*ei/l.powi(2),   0.0, -12.0*ei/l.powi(3),  6.0*ei/l.powi(2),
              0.0,   6.0*ei/l.powi(2),      4.0*ei/l,       0.0,  -6.0*ei/l.powi(2),      2.0*ei/l,
            -ea_l,             0.0,             0.0,   ea_l,             0.0,             0.0,
              0.0, -12.0*ei/l.powi(3), -6.0*ei/l.powi(2),   0.0,  12.0*ei/l.powi(3), -6.0*ei/l.powi(2),
              0.0,   6.0*ei/l.powi(2),      2.0*ei/l,       0.0,  -6.0*ei/l.powi(2),      4.0*ei/l,
        ]);
        k
    }

    /// First-order (linearized) geometric-stiffness correction from the
    /// current axial force `p` (tension positive), added to the elastic
    /// bending block for `GeomTransf::PDelta` — the standard consistent
    /// geometric stiffness matrix (e.g. Przemieniecki §2), *not* a full
    /// corotational update. Only the tangent used to solve for the next
    /// displacement increment is corrected this way; the resisting-force
    /// recovery (`form_tangent_and_resistance`'s second return value)
    /// stays purely elastic. Consistent path-following P-Delta (where the
    /// resisting force itself must reflect the same correction) needs
    /// Newton iteration to track correctly and is deferred alongside M4.
    fn geometric_stiffness(&self, p: f64, length: f64) -> SMatrix<f64, 6, 6> {
        let l = length;
        let mut kg = SMatrix::<f64, 6, 6>::zeros();
        #[rustfmt::skip]
        let block = [
            [ 6.0/5.0,       l/10.0,      -6.0/5.0,       l/10.0],
            [   l/10.0,  2.0*l*l/15.0,      -l/10.0,   -l*l/30.0],
            [-6.0/5.0,      -l/10.0,       6.0/5.0,      -l/10.0],
            [   l/10.0,    -l*l/30.0,      -l/10.0,  2.0*l*l/15.0],
        ];
        let indices = [1, 2, 4, 5];
        for (bi, &gi) in indices.iter().enumerate() {
            for (bj, &gj) in indices.iter().enumerate() {
                kg[(gi, gj)] = (p / l) * block[bi][bj];
            }
        }
        kg
    }

    /// Block-diagonal rotation from global to local DOFs: `d_local = T * d_global`.
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

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let k_local = self.local_elastic_stiffness(length);

        let d_global = SVector::<f64, 6>::from_row_slice(&[
            node_i.displacement[0],
            node_i.displacement[1],
            node_i.displacement[2],
            node_j.displacement[0],
            node_j.displacement[1],
            node_j.displacement[2],
        ]);
        let d_local = t * d_global;

        let resistance_local = k_local * d_local;
        let resistance = t.transpose() * resistance_local;

        let k_total_local = match self.transform {
            GeomTransf::Linear => k_local,
            GeomTransf::PDelta => {
                let axial_force = resistance_local[3]; // tension-positive axial force at node j
                k_local + self.geometric_stiffness(axial_force, length)
            }
        };
        let k = t.transpose() * k_total_local * t;

        (k, resistance)
    }

    /// Equivalent nodal load (global coordinates) from the element's
    /// uniform transverse load, via consistent (virtual-work) Hermite
    /// cubic shape-function integration — see implementation-plan §3.4.
    /// Zero when `w_transverse` is zero.
    pub(super) fn form_load_vector(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        if self.w_transverse == 0.0 {
            return SVector::<f64, 6>::zeros();
        }
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let t = self.transformation(cx, cy);
        let w = self.w_transverse;
        let l = length;
        let local = SVector::<f64, 6>::from_row_slice(&[
            0.0,
            w * l / 2.0,
            w * l * l / 12.0,
            0.0,
            w * l / 2.0,
            -w * l * l / 12.0,
        ]);
        t.transpose() * local
    }

    /// Lumped mass: half the element's total mass (`density * a * length`)
    /// at each node's translational DOFs, zero rotational contribution —
    /// the simplest standard lumped-mass model (§3.4: "lumped, to start");
    /// a consistent (non-diagonal) mass matrix or a nonzero rotational
    /// lumped inertia is a further refinement, not built until needed.
    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, 6> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let half = self.density * self.a * length / 2.0;
        SVector::<f64, 6>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}
