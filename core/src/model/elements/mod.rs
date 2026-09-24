use std::convert::Infallible;

use nalgebra::{SMatrix, SVector};
use slotmap::{new_key_type, Key};

use super::{ElementLoad, Node, Node3, Node3Id, NodeId, ELEMENT_DOF, NDF, PLANAR_NDIM, SPATIAL_ELEMENT_DOF, SPATIAL_NDF, SPATIAL_NDIM};

/// What `Domain`'s generic assembly/state plumbing needs from an element
/// catalog — `Element`/`Element3`'s shared shape, not their (genuinely
/// different) physics. A trait rather than a shared enum: `Element` and
/// `Element3` hold different variant sets (`Truss` vs `Truss3`, etc.) and a
/// single enum spanning both would either duplicate every variant anyway or
/// require `Box<dyn Trait>` — heap allocation and vtable dispatch in the
/// hot assembly loop, which the spatial-architecture plan rules out for a
/// browser/wasm target. `NDIM`/`NDOF`/`ELEMENT_DOF` are separate const
/// generic parameters (not computed from each other) because stable Rust
/// can't evaluate `2 * NDOF` inside a generic item (no `generic_const_exprs`
/// yet) — the same reason top-level `ELEMENT_DOF`/`SPATIAL_ELEMENT_DOF`
/// constants exist instead of inline expressions.
pub trait ElementOps<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId: Copy> {
    /// An element-load kind this catalog supports (`ElementLoad` for
    /// `Element`). `Element3` uses `Infallible` — no spatial element
    /// carries a load yet, so a `HashMap<_, Infallible>` is always empty
    /// and `add_element_load` is uncallable (no `Infallible` value exists)
    /// until a real spatial load type replaces it.
    type Load;

    /// This catalog's own `Domain` element-store key (`ElementId`/
    /// `Element3Id`) — an associated type rather than a further generic
    /// parameter of `Domain` itself, so it's always determined by `E` and
    /// never a separate, independently-unconstrained type variable: a
    /// `Domain::add_element(element)` call pins `E` concretely from real
    /// usage (the `Element`/`Element3` value passed in), but nothing in an
    /// ordinary test ever *uses* the returned element ID for anything, so a
    /// bare `EId` parameter with nothing forcing it — even with a default —
    /// left `Domain::new()` unable to infer it.
    type Id: Key;

    fn nodes(&self) -> [NId; 2];

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node<NDIM, NDOF>,
        node_j: &Node<NDIM, NDOF>,
    ) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>);

    fn form_load_vector(&self, node_i: &Node<NDIM, NDOF>, node_j: &Node<NDIM, NDOF>, load: Option<&Self::Load>) -> SVector<f64, ELEMENT_DOF>;

    fn form_mass(&self, node_i: &Node<NDIM, NDOF>, node_j: &Node<NDIM, NDOF>) -> SVector<f64, ELEMENT_DOF>;

    fn commit(&mut self, node_i: &Node<NDIM, NDOF>, node_j: &Node<NDIM, NDOF>);
}

mod disp_beam_column;
mod elastic_beam_column;
mod force_beam_column;
mod truss;
mod zero_length;

pub use disp_beam_column::{DispBeamColumn, DispBeamColumn3};
pub use elastic_beam_column::{ElasticBeamColumn, ElasticBeamColumn3};
pub use force_beam_column::{ForceBeamColumn, ForceBeamColumn3};
pub use truss::{SpatialElementMatrix, SpatialElementVector, Truss, Truss3};
pub use zero_length::{ZeroLength, ZeroLength3};

new_key_type! {
    /// Generational index into `Domain`'s element store (§2.2).
    pub struct ElementId;

    /// Generational index into `Domain3`'s element store. Distinct from
    /// `ElementId` for the same reason `Node3Id` is distinct from `NodeId`
    /// (see its doc comment).
    pub struct Element3Id;
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
    /// same DOF order as above) from `load` — the `ElementLoad` (if any)
    /// that whichever `LoadPattern` is currently being assembled has on
    /// this element (§3.4) — e.g. a beam-column's distributed transverse
    /// load. Zero when `load` is `None`, and for elements with no
    /// element-load support at all (`Truss`, `ZeroLength`, `DispBeamColumn`
    /// — see its doc comment for why).
    pub fn form_load_vector(&self, node_i: &Node, node_j: &Node, load: Option<&ElementLoad>) -> SVector<f64, ELEMENT_DOF> {
        match (self, load) {
            (Element::ElasticBeamColumn(b), Some(ElementLoad::UniformTransverse(w))) => b.form_load_vector(node_i, node_j, *w),
            _ => SVector::<f64, ELEMENT_DOF>::zeros(),
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

impl ElementOps<PLANAR_NDIM, NDF, ELEMENT_DOF, NodeId> for Element {
    type Load = ElementLoad;
    type Id = ElementId;

    fn nodes(&self) -> [NodeId; 2] {
        Element::nodes(self)
    }

    fn form_tangent_and_resistance(&self, node_i: &Node, node_j: &Node) -> (SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>, SVector<f64, ELEMENT_DOF>) {
        Element::form_tangent_and_resistance(self, node_i, node_j)
    }

    fn form_load_vector(&self, node_i: &Node, node_j: &Node, load: Option<&ElementLoad>) -> SVector<f64, ELEMENT_DOF> {
        Element::form_load_vector(self, node_i, node_j, load)
    }

    fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        Element::form_mass(self, node_i, node_j)
    }

    fn commit(&mut self, node_i: &Node, node_j: &Node) {
        Element::commit(self, node_i, node_j)
    }
}

/// Spatial element catalog — `Domain3`'s counterpart to `Element`. Closed
/// enum, `match`-based dispatch, same reasoning as `Element` (§2.1).
#[derive(Debug, Clone)]
pub enum Element3 {
    Truss3(Truss3),
    ZeroLength3(ZeroLength3),
    ElasticBeamColumn3(ElasticBeamColumn3),
    DispBeamColumn3(DispBeamColumn3),
    ForceBeamColumn3(ForceBeamColumn3),
}

impl Element3 {
    pub fn nodes(&self) -> [Node3Id; 2] {
        match self {
            Element3::Truss3(t) => [t.node_i, t.node_j],
            Element3::ZeroLength3(z) => [z.node_i, z.node_j],
            Element3::ElasticBeamColumn3(b) => [b.node_i, b.node_j],
            Element3::DispBeamColumn3(b) => [b.node_i, b.node_j],
            Element3::ForceBeamColumn3(b) => [b.node_i, b.node_j],
        }
    }

    pub fn form_tangent_and_resistance(
        &self,
        node_i: &Node3,
        node_j: &Node3,
    ) -> (SMatrix<f64, SPATIAL_ELEMENT_DOF, SPATIAL_ELEMENT_DOF>, SVector<f64, SPATIAL_ELEMENT_DOF>) {
        match self {
            Element3::Truss3(t) => t.form_tangent_and_resistance(node_i, node_j),
            Element3::ZeroLength3(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element3::ElasticBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element3::DispBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j),
        }
    }

    /// This element's lumped-mass contribution — see `Element::form_mass`.
    pub fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match self {
            Element3::Truss3(t) => t.form_mass(node_i, node_j),
            Element3::ZeroLength3(_) => SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros(),
            Element3::ElasticBeamColumn3(b) => b.form_mass(node_i, node_j),
            Element3::DispBeamColumn3(b) => b.form_mass(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.form_mass(node_i, node_j),
        }
    }

    pub fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        match self {
            Element3::Truss3(t) => t.commit(node_i, node_j),
            Element3::ZeroLength3(z) => z.commit(node_i, node_j),
            Element3::ElasticBeamColumn3(_) => {}
            Element3::DispBeamColumn3(b) => b.commit(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.commit(node_i, node_j),
        }
    }
}

impl ElementOps<SPATIAL_NDIM, SPATIAL_NDF, SPATIAL_ELEMENT_DOF, Node3Id> for Element3 {
    /// No spatial element carries a load yet — see this trait's doc comment.
    type Load = Infallible;
    type Id = Element3Id;

    fn nodes(&self) -> [Node3Id; 2] {
        Element3::nodes(self)
    }

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node3,
        node_j: &Node3,
    ) -> (SMatrix<f64, SPATIAL_ELEMENT_DOF, SPATIAL_ELEMENT_DOF>, SVector<f64, SPATIAL_ELEMENT_DOF>) {
        Element3::form_tangent_and_resistance(self, node_i, node_j)
    }

    /// Always zero: `Option<&Infallible>` can only ever be `None`.
    fn form_load_vector(&self, _node_i: &Node3, _node_j: &Node3, _load: Option<&Infallible>) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros()
    }

    fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        Element3::form_mass(self, node_i, node_j)
    }

    fn commit(&mut self, node_i: &Node3, node_j: &Node3) {
        Element3::commit(self, node_i, node_j)
    }
}
