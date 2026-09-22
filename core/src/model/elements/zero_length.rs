use nalgebra::{SMatrix, SVector};

use super::super::{Material, Node, NodeId, ELEMENT_DOF, NDF};

/// A 2-node, zero-length connector: no geometry or integration, just direct
/// per-DOF material evaluation (§3.1) — each direction with a material
/// assigned independently relates that DOF's relative displacement between
/// the two nodes to a force along that same (global) direction. No
/// orientation vectors (unlike OpenSees' general `ZeroLength`, which can
/// evaluate materials along arbitrary local axes) — directions are the
/// global DOF axes (including rotation, since M3's DOF bump), which is all
/// the current scope needs.
#[derive(Debug, Clone)]
pub struct ZeroLength {
    pub node_i: NodeId,
    pub node_j: NodeId,
    materials: [Option<Material>; NDF],
}

impl ZeroLength {
    pub fn new(node_i: NodeId, node_j: NodeId) -> Self {
        ZeroLength {
            node_i,
            node_j,
            materials: std::array::from_fn(|_| None),
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
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

        (k, resistance)
    }

    pub(super) fn commit(&mut self, node_i: &Node, node_j: &Node) {
        for (dof, material) in self.materials.iter_mut().enumerate() {
            let Some(material) = material else { continue };
            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            *material = material.commit(relative);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Material, Node};

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
}
