use nalgebra::{SMatrix, SVector};

use super::super::{Material, Node, NodeId, ELEMENT_DOF};

/// A 2-node axial truss: fixed-size element-local linear algebra (§2.4), no
/// heap allocation in the hot path. Only ever populates the translational
/// DOF entries of the (now 3-DOF-per-node) local/global matrices — the
/// rotation entries stay zero, so a node connected only to `Truss`/
/// `ZeroLength` elements needs its rotation DOF fixed by the model, or the
/// global system is singular in that DOF.
#[derive(Debug, Clone)]
pub struct Truss {
    pub node_i: NodeId,
    pub node_j: NodeId,
    pub area: f64,
    pub material: Material,
    /// Mass per unit volume. Zero (the default) means massless — existing
    /// models are unaffected unless they opt in via `with_density`.
    pub density: f64,
}

impl Truss {
    pub fn new(node_i: NodeId, node_j: NodeId, area: f64, material: Material) -> Self {
        Truss {
            node_i,
            node_j,
            area,
            material,
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

    /// Axial elongation from current nodal displacements, projected onto
    /// the (undeformed) bar axis — small-displacement theory. Shared by
    /// `form_tangent_and_resistance` (trial) and `commit`.
    fn strain(&self, node_i: &Node, node_j: &Node) -> f64 {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        ((node_j.displacement[0] - node_i.displacement[0]) * cx
            + (node_j.displacement[1] - node_i.displacement[1]) * cy)
            / length
    }

    pub(super) fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);
        let strain = self.strain(node_i, node_j);
        let (stress, tangent_modulus) = self.material.trial_stress_tangent(strain);

        // Local DOF order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; rotation
        // entries (2, 5) stay zero — a truss carries no moment.
        let b = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[-cx, -cy, 0.0, cx, cy, 0.0]);

        let k = (tangent_modulus * self.area / length) * (b * b.transpose());
        let resistance = (stress * self.area) * b;

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        let strain = self.strain(node_i, node_j);
        self.material = self.material.commit(strain);
    }

    /// Lumped mass: half the element's total mass (`density * area * length`)
    /// at each node, split equally between that node's translational DOFs —
    /// a point mass has no directional preference. Zero rotational
    /// contribution (index 2, 5) — a truss carries no moment either.
    pub(super) fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let half = self.density * self.area * length / 2.0;
        SVector::<f64, ELEMENT_DOF>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}
