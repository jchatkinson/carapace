use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

use super::{Node, NodeId, ELEMENT_DOF};

mod disp_beam_column;
mod elastic_beam_column;
mod force_beam_column;
mod truss;
mod zero_length;

pub use disp_beam_column::DispBeamColumn;
pub use elastic_beam_column::ElasticBeamColumn;
pub use force_beam_column::ForceBeamColumn;
pub use truss::Truss;
pub use zero_length::ZeroLength;

new_key_type! {
    /// Generational index into `Domain`'s element store (§2.2).
    pub struct ElementId;
}

/// Element catalog. Closed enum, `match`-based dispatch, no `Box<dyn Trait>`
/// (§2.1).
#[derive(Debug, Clone)]
pub enum Element {
    Truss(Truss),
    ZeroLength(ZeroLength),
    ElasticBeamColumn(ElasticBeamColumn),
    DispBeamColumn(DispBeamColumn),
    ForceBeamColumn(ForceBeamColumn),
}

impl Element {
    pub fn nodes(&self) -> [NodeId; 2] {
        match self {
            Element::Truss(t) => [t.node_i, t.node_j],
            Element::ZeroLength(z) => [z.node_i, z.node_j],
            Element::ElasticBeamColumn(b) => [b.node_i, b.node_j],
            Element::DispBeamColumn(b) => [b.node_i, b.node_j],
            Element::ForceBeamColumn(b) => [b.node_i, b.node_j],
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
            Element::DispBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element::ForceBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j),
        }
    }

    /// This element's equivalent nodal load vector (global coordinates,
    /// same DOF order as above) from any element load applied to it (§3.4)
    /// — e.g. a beam-column's distributed transverse load. Zero for
    /// elements with no element-load support (`Truss`, `ZeroLength`,
    /// `DispBeamColumn` — see its doc comment for why).
    pub fn form_load_vector(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(_) | Element::ZeroLength(_) | Element::DispBeamColumn(_) | Element::ForceBeamColumn(_) => {
                SVector::<f64, ELEMENT_DOF>::zeros()
            }
            Element::ElasticBeamColumn(b) => b.form_load_vector(node_i, node_j),
        }
    }

    /// This element's lumped-mass contribution (diagonal only — see §3.4:
    /// "lumped, to start") in the same local DOF order, geometry-dependent
    /// (`length`) so it needs both nodes. Zero for `ZeroLength` (a spring/
    /// connector, not a mass-bearing member) and for any element with the
    /// default `density = 0.0`.
    pub fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(t) => t.form_mass(node_i, node_j),
            Element::ZeroLength(_) => SVector::<f64, ELEMENT_DOF>::zeros(),
            Element::ElasticBeamColumn(b) => b.form_mass(node_i, node_j),
            Element::DispBeamColumn(b) => b.form_mass(node_i, node_j),
            Element::ForceBeamColumn(b) => b.form_mass(node_i, node_j),
        }
    }

    /// Commit this element's material(s) at the given (final, converged)
    /// node state — see `Material`'s doc comment for why this is the only
    /// place a `Material` ever mutates. A no-op for `ElasticBeamColumn`
    /// (no `Material` — its response is closed-form, §3.1) and for any
    /// `ZeroLength` direction with no material assigned.
    pub fn commit(&mut self, node_i: &Node, node_j: &Node) {
        match self {
            Element::Truss(t) => t.commit(node_i, node_j),
            Element::ZeroLength(z) => z.commit(node_i, node_j),
            Element::ElasticBeamColumn(_) => {}
            Element::DispBeamColumn(b) => b.commit(node_i, node_j),
            Element::ForceBeamColumn(b) => b.commit(node_i, node_j),
        }
    }
}
