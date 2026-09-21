use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

use super::{Material, Node, NodeId};

new_key_type! {
    /// Generational index into `Domain`'s element store (§3.2).
    pub struct ElementId;
}

/// Element catalog. Closed enum, `match`-based dispatch, no `Box<dyn Trait>`
/// (§3.1). `Truss` is the only variant for M1; `ZeroLength` /
/// `ElasticBeamColumn` / ... land at later milestones per the plan's §4.1
/// table.
#[derive(Debug, Clone)]
pub enum Element {
    Truss(Truss),
}

impl Element {
    pub fn nodes(&self) -> [NodeId; 2] {
        match self {
            Element::Truss(t) => [t.node_i, t.node_j],
        }
    }

    /// Form this element's contribution to the global tangent stiffness and
    /// internal resisting force, given its two nodes' current state.
    /// Returned in element-local DOF order `[ux_i, uy_i, ux_j, uy_j]`; the
    /// caller (`Domain::form_tangent_and_residual`) scatters these into the
    /// global system using each node's equation numbers.
    pub fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 4, 4>, SVector<f64, 4>) {
        match self {
            Element::Truss(t) => t.form_tangent_and_resistance(node_i, node_j),
        }
    }
}

/// A 2-node axial truss: fixed-size element-local linear algebra (§3.4), no
/// heap allocation in the hot path.
#[derive(Debug, Clone)]
pub struct Truss {
    pub node_i: NodeId,
    pub node_j: NodeId,
    pub area: f64,
    pub material: Material,
}

impl Truss {
    pub fn new(node_i: NodeId, node_j: NodeId, area: f64, material: Material) -> Self {
        Truss {
            node_i,
            node_j,
            area,
            material,
        }
    }

    fn geometry(&self, node_i: &Node, node_j: &Node) -> (f64, f64, f64) {
        let dx = node_j.coords[0] - node_i.coords[0];
        let dy = node_j.coords[1] - node_i.coords[1];
        let length = (dx * dx + dy * dy).sqrt();
        (length, dx / length, dy / length)
    }

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 4, 4>, SVector<f64, 4>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);

        // Axial elongation from current nodal displacements, projected onto
        // the (undeformed) bar axis — small-displacement theory.
        let strain = ((node_j.displacement[0] - node_i.displacement[0]) * cx
            + (node_j.displacement[1] - node_i.displacement[1]) * cy)
            / length;
        let (stress, tangent_modulus) = self.material.stress_tangent(strain);

        // Direction cosine vector, local DOF order [ux_i, uy_i, ux_j, uy_j].
        let b = SVector::<f64, 4>::from_column_slice(&[-cx, -cy, cx, cy]);

        let k = (tangent_modulus * self.area / length) * (b * b.transpose());
        let resistance = (stress * self.area) * b;

        (k, resistance)
    }
}
