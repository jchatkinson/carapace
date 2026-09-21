use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

use super::{ElasticBeamColumn, Material, Node, NodeId, ELEMENT_DOF, NDF};

new_key_type! {
    /// Generational index into `Domain`'s element store (§3.2).
    pub struct ElementId;
}

/// Element catalog. Closed enum, `match`-based dispatch, no `Box<dyn Trait>`
/// (§3.1). `DispBeamColumn` / `ForceBeamColumn` / ... land at later
/// milestones per the plan's §4.1 table.
#[derive(Debug, Clone)]
pub enum Element {
    Truss(Truss),
    ZeroLength(ZeroLength),
    ElasticBeamColumn(ElasticBeamColumn),
}

impl Element {
    pub fn nodes(&self) -> [NodeId; 2] {
        match self {
            Element::Truss(t) => [t.node_i, t.node_j],
            Element::ZeroLength(z) => [z.node_i, z.node_j],
            Element::ElasticBeamColumn(b) => [b.node_i, b.node_j],
        }
    }

    /// Form this element's contribution to the global tangent stiffness and
    /// internal resisting force, given its two nodes' current state.
    /// Returned in element-local DOF order
    /// `[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]`; the caller
    /// (`Domain::form_tangent_and_residual`) scatters these into the global
    /// system using each node's equation numbers.
    pub fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        match self {
            Element::Truss(t) => t.form_tangent_and_resistance(node_i, node_j),
            Element::ZeroLength(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element::ElasticBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j),
        }
    }

    /// This element's equivalent nodal load vector (global coordinates,
    /// same DOF order as above) from any element load applied to it (§4.4)
    /// — e.g. a beam-column's distributed transverse load. Zero for
    /// elements with no element-load support (`Truss`, `ZeroLength`).
    pub fn form_load_vector(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(_) | Element::ZeroLength(_) => SVector::<f64, ELEMENT_DOF>::zeros(),
            Element::ElasticBeamColumn(b) => b.form_load_vector(node_i, node_j),
        }
    }

    /// This element's lumped-mass contribution (diagonal only — see §4.4:
    /// "lumped, to start") in the same local DOF order, geometry-dependent
    /// (`length`) for `Truss`/`ElasticBeamColumn` so it needs both nodes.
    /// Zero for `ZeroLength` (a spring/connector, not a mass-bearing
    /// member) and for `Truss`/`ElasticBeamColumn` with the default
    /// `density = 0.0`.
    pub fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(t) => t.form_mass(node_i, node_j),
            Element::ZeroLength(_) => SVector::<f64, ELEMENT_DOF>::zeros(),
            Element::ElasticBeamColumn(b) => b.form_mass(node_i, node_j),
        }
    }
}

/// A 2-node axial truss: fixed-size element-local linear algebra (§3.4), no
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

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let (length, cx, cy) = self.geometry(node_i, node_j);

        // Axial elongation from current nodal displacements, projected onto
        // the (undeformed) bar axis — small-displacement theory.
        let strain = ((node_j.displacement[0] - node_i.displacement[0]) * cx
            + (node_j.displacement[1] - node_i.displacement[1]) * cy)
            / length;
        let (stress, tangent_modulus) = self.material.stress_tangent(strain);

        // Local DOF order [ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]; rotation
        // entries (2, 5) stay zero — a truss carries no moment.
        let b = SVector::<f64, ELEMENT_DOF>::from_column_slice(&[-cx, -cy, 0.0, cx, cy, 0.0]);

        let k = (tangent_modulus * self.area / length) * (b * b.transpose());
        let resistance = (stress * self.area) * b;

        (k, resistance)
    }

    /// Lumped mass: half the element's total mass (`density * area * length`)
    /// at each node, split equally between that node's translational DOFs —
    /// a point mass has no directional preference. Zero rotational
    /// contribution (index 2, 5) — a truss carries no moment either.
    fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        let (length, _cx, _cy) = self.geometry(node_i, node_j);
        let half = self.density * self.area * length / 2.0;
        SVector::<f64, ELEMENT_DOF>::from_column_slice(&[half, half, 0.0, half, half, 0.0])
    }
}

/// A 2-node, zero-length connector: no geometry or integration, just direct
/// per-DOF material evaluation (§4.1) — each direction with a material
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
            materials: [None; NDF],
        }
    }

    pub fn with_material(mut self, dof: usize, material: Material) -> Self {
        self.materials[dof] = Some(material);
        self
    }

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        let mut k = SMatrix::<f64, ELEMENT_DOF, ELEMENT_DOF>::zeros();
        let mut resistance = SVector::<f64, ELEMENT_DOF>::zeros();

        for (dof, material) in self.materials.iter().enumerate() {
            let Some(material) = material else { continue };

            let relative = node_j.displacement[dof] - node_i.displacement[dof];
            let (force, tangent_modulus) = material.stress_tangent(relative);

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Node;

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
            .with_material(0, Material::ElasticPP { e: 100.0, eyp: 0.01 });
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
