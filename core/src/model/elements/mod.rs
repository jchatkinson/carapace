use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

use super::{
    ElementLoad, ElementLoad3, Node, Node3, Node3Id, NodeId, ELEMENT_DOF, NDF, PLANAR_NDIM,
    SPATIAL_ELEMENT_DOF, SPATIAL_NDF, SPATIAL_NDIM,
};

mod ops;
pub use ops::{
    two_node_dofs, DofMask, DofRef, ElementForce, ElementOps, NodeList, NodeView, TangentSink,
    VectorSink, MAX_ELEMENT_NODES,
};

mod disp_beam_column;
mod elastic_beam_column;
mod force_beam_column;
mod truss;
mod uniform_load;
mod zero_length;

pub use disp_beam_column::{DispBeamColumn, DispBeamColumn3};
pub use elastic_beam_column::{ElasticBeamColumn, ElasticBeamColumn3};
pub use force_beam_column::{ForceBeamColumn, ForceBeamColumn3};
pub use truss::{SpatialElementMatrix, SpatialElementVector, Truss, Truss3};
pub use zero_length::{
    Friction, Friction3, Orientation, OrientationError, ZeroLength, ZeroLength3, ZeroLengthSection,
    ZeroLengthSection3,
};

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
    ZeroLengthSection(ZeroLengthSection),
    ElasticBeamColumn(ElasticBeamColumn),
    DispBeamColumn(DispBeamColumn),
    ForceBeamColumn(ForceBeamColumn),
}

impl Element {
    pub fn nodes(&self) -> [NodeId; 2] {
        match self {
            Element::Truss(t) => [t.node_i, t.node_j],
            Element::ZeroLength(z) => [z.node_i, z.node_j],
            Element::ZeroLengthSection(z) => [z.node_i, z.node_j],
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
    fn form_tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
        load: Option<&ElementLoad>,
    ) -> (
        SMatrix<f64, ELEMENT_DOF, ELEMENT_DOF>,
        SVector<f64, ELEMENT_DOF>,
    ) {
        match self {
            Element::Truss(t) => t.form_tangent_and_resistance(node_i, node_j),
            Element::ZeroLength(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element::ZeroLengthSection(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element::ElasticBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element::DispBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element::ForceBeamColumn(b) => b.form_tangent_and_resistance(node_i, node_j, load),
        }
    }

    /// This element's equivalent nodal load vector (global coordinates,
    /// same DOF order as above) from `load` — the `ElementLoad` (if any)
    /// that whichever `LoadPattern` is currently being assembled has on
    /// this element (§3.4) — e.g. a beam-column's distributed transverse
    /// load. Zero when `load` is `None`, and for elements with no
    /// element-load support at all (`Truss`, `ZeroLength`).
    fn form_load_vector(
        &self,
        node_i: &Node,
        node_j: &Node,
        load: Option<&ElementLoad>,
    ) -> SVector<f64, ELEMENT_DOF> {
        match (self, load) {
            (Element::ElasticBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy)
            }
            (Element::DispBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy)
            }
            (Element::ForceBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy)
            }
            _ => SVector::<f64, ELEMENT_DOF>::zeros(),
        }
    }

    /// See `ElementOps::form_local_load_vector`'s doc comment.
    fn form_local_load_vector(
        &self,
        node_i: &Node,
        node_j: &Node,
        load: Option<&ElementLoad>,
    ) -> SVector<f64, ELEMENT_DOF> {
        match (self, load) {
            (Element::ElasticBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy)
            }
            (Element::DispBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy)
            }
            (Element::ForceBeamColumn(b), Some(ElementLoad::Uniform { wx, wy })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy)
            }
            _ => SVector::<f64, ELEMENT_DOF>::zeros(),
        }
    }

    /// This element's lumped-mass contribution (diagonal only — see §3.4:
    /// "lumped, to start") in the same local DOF order, geometry-dependent
    /// (`length`) so it needs both nodes. Zero for `ZeroLength` (a spring/
    /// connector, not a mass-bearing member) and for any element with the
    /// default `density = 0.0`.
    fn form_mass(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(t) => t.form_mass(node_i, node_j),
            Element::ZeroLength(_) => SVector::<f64, ELEMENT_DOF>::zeros(),
            Element::ZeroLengthSection(_) => SVector::<f64, ELEMENT_DOF>::zeros(),
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
    fn commit(&mut self, node_i: &Node, node_j: &Node, load: Option<&ElementLoad>) {
        match self {
            Element::Truss(t) => t.commit(node_i, node_j),
            Element::ZeroLength(z) => z.commit(node_i, node_j),
            Element::ZeroLengthSection(z) => z.commit(node_i, node_j),
            Element::ElasticBeamColumn(_) => {}
            Element::DispBeamColumn(b) => b.commit(node_i, node_j),
            Element::ForceBeamColumn(b) => b.commit(node_i, node_j, load),
        }
    }

    /// See `ElementOps::local_force`'s doc comment. `ForceBeamColumn`'s
    /// local force is cached at `commit` time (its own `local_force()`
    /// takes no node arguments — reading the cache, not recomputing);
    /// every other variant recomputes fresh from `node_i`/`node_j`.
    fn local_force(&self, node_i: &Node, node_j: &Node) -> SVector<f64, ELEMENT_DOF> {
        match self {
            Element::Truss(t) => t.local_force(node_i, node_j),
            Element::ZeroLength(z) => z.local_force(node_i, node_j),
            Element::ZeroLengthSection(z) => z.local_force(node_i, node_j),
            Element::ElasticBeamColumn(b) => b.local_force(node_i, node_j),
            Element::DispBeamColumn(b) => b.local_force(node_i, node_j),
            Element::ForceBeamColumn(b) => b.local_force(),
        }
    }

    /// See `ElementOps::fiber_responses`'s doc comment.
    fn fiber_responses(&self, node_i: &Node, node_j: &Node) -> Option<Vec<Vec<(f64, f64)>>> {
        match self {
            Element::DispBeamColumn(b) => Some(b.fiber_responses(node_i, node_j)),
            Element::ForceBeamColumn(b) => Some(b.fiber_responses()),
            Element::Truss(_)
            | Element::ZeroLength(_)
            | Element::ZeroLengthSection(_)
            | Element::ElasticBeamColumn(_) => None,
        }
    }
}

impl ElementOps<PLANAR_NDIM, NDF, NodeId> for Element {
    type Load = ElementLoad;
    type Id = ElementId;

    fn nodes(&self) -> NodeList<NodeId> {
        Element::nodes(self).into_iter().collect()
    }

    fn dof_mask(&self) -> DofMask {
        DofMask::all(NDF)
    }

    fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let [i, j] = Element::nodes(self);
        let (k, r) = Element::form_tangent_and_resistance(self, nodes.get(i), nodes.get(j), load);
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &k, &r);
    }

    fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let [i, j] = Element::nodes(self);
        let v = Element::form_load_vector(self, nodes.get(i), nodes.get(j), load);
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &v);
    }

    fn assemble_mass<S: VectorSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        sink: &mut S,
    ) {
        let [i, j] = Element::nodes(self);
        let v = Element::form_mass(self, nodes.get(i), nodes.get(j));
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &v);
    }

    fn commit(
        &mut self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
    ) {
        let [i, j] = Element::nodes(self);
        Element::commit(self, nodes.get(i), nodes.get(j), load)
    }

    fn local_force_width(&self) -> usize {
        ELEMENT_DOF
    }

    fn local_force(&self, nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>) -> ElementForce {
        let [i, j] = Element::nodes(self);
        Element::local_force(self, nodes.get(i), nodes.get(j)).into()
    }

    fn local_load_force(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
    ) -> ElementForce {
        let [i, j] = Element::nodes(self);
        Element::form_local_load_vector(self, nodes.get(i), nodes.get(j), load).into()
    }

    fn fiber_responses(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
    ) -> Option<Vec<Vec<(f64, f64)>>> {
        let [i, j] = Element::nodes(self);
        Element::fiber_responses(self, nodes.get(i), nodes.get(j))
    }
}
/// Spatial element catalog — `Domain3`'s counterpart to `Element`. Closed
/// enum, `match`-based dispatch, same reasoning as `Element` (§2.1).
#[derive(Debug, Clone)]
pub enum Element3 {
    Truss3(Truss3),
    ZeroLength3(ZeroLength3),
    ZeroLengthSection3(ZeroLengthSection3),
    ElasticBeamColumn3(ElasticBeamColumn3),
    DispBeamColumn3(DispBeamColumn3),
    ForceBeamColumn3(ForceBeamColumn3),
}

impl Element3 {
    pub fn nodes(&self) -> [Node3Id; 2] {
        match self {
            Element3::Truss3(t) => [t.node_i, t.node_j],
            Element3::ZeroLength3(z) => [z.node_i, z.node_j],
            Element3::ZeroLengthSection3(z) => [z.node_i, z.node_j],
            Element3::ElasticBeamColumn3(b) => [b.node_i, b.node_j],
            Element3::DispBeamColumn3(b) => [b.node_i, b.node_j],
            Element3::ForceBeamColumn3(b) => [b.node_i, b.node_j],
        }
    }

    fn form_tangent_and_resistance(
        &self,
        node_i: &Node3,
        node_j: &Node3,
        load: Option<&ElementLoad3>,
    ) -> (
        SMatrix<f64, SPATIAL_ELEMENT_DOF, SPATIAL_ELEMENT_DOF>,
        SVector<f64, SPATIAL_ELEMENT_DOF>,
    ) {
        match self {
            Element3::Truss3(t) => t.form_tangent_and_resistance(node_i, node_j),
            Element3::ZeroLength3(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element3::ZeroLengthSection3(z) => z.form_tangent_and_resistance(node_i, node_j),
            Element3::ElasticBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element3::DispBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.form_tangent_and_resistance(node_i, node_j, load),
        }
    }

    /// This element's equivalent nodal load vector from `load` — the
    /// spatial counterpart of `Element::form_load_vector`. Zero when
    /// `load` is `None`, and for elements with no element-load support
    /// (`Truss3`, `ZeroLength3`, `DispBeamColumn3`, `ForceBeamColumn3` —
    /// matching their planar counterparts, see `ElementOps::Load`'s doc
    /// comment).
    fn form_load_vector(
        &self,
        node_i: &Node3,
        node_j: &Node3,
        load: Option<&ElementLoad3>,
    ) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match (self, load) {
            (Element3::ElasticBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            (Element3::DispBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            (Element3::ForceBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            _ => SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros(),
        }
    }

    /// See `ElementOps::form_local_load_vector`'s doc comment.
    fn form_local_load_vector(
        &self,
        node_i: &Node3,
        node_j: &Node3,
        load: Option<&ElementLoad3>,
    ) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match (self, load) {
            (Element3::ElasticBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            (Element3::DispBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            (Element3::ForceBeamColumn3(b), Some(ElementLoad3::Uniform { wx, wy, wz })) => {
                b.form_local_load_vector(node_i, node_j, *wx, *wy, *wz)
            }
            _ => SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros(),
        }
    }

    /// This element's lumped-mass contribution — see `Element::form_mass`.
    fn form_mass(&self, node_i: &Node3, node_j: &Node3) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match self {
            Element3::Truss3(t) => t.form_mass(node_i, node_j),
            Element3::ZeroLength3(_) => SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros(),
            Element3::ZeroLengthSection3(_) => SVector::<f64, SPATIAL_ELEMENT_DOF>::zeros(),
            Element3::ElasticBeamColumn3(b) => b.form_mass(node_i, node_j),
            Element3::DispBeamColumn3(b) => b.form_mass(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.form_mass(node_i, node_j),
        }
    }

    fn commit(&mut self, node_i: &Node3, node_j: &Node3, load: Option<&ElementLoad3>) {
        match self {
            Element3::Truss3(t) => t.commit(node_i, node_j),
            Element3::ZeroLength3(z) => z.commit(node_i, node_j),
            Element3::ZeroLengthSection3(z) => z.commit(node_i, node_j),
            Element3::ElasticBeamColumn3(_) => {}
            Element3::DispBeamColumn3(b) => b.commit(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.commit(node_i, node_j, load),
        }
    }

    /// See `Element::local_force`'s doc comment.
    fn local_force(&self, node_i: &Node3, node_j: &Node3) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match self {
            Element3::Truss3(t) => t.local_force(node_i, node_j),
            Element3::ZeroLength3(z) => z.local_force(node_i, node_j),
            Element3::ZeroLengthSection3(z) => z.local_force(node_i, node_j),
            Element3::ElasticBeamColumn3(b) => b.local_force(node_i, node_j),
            Element3::DispBeamColumn3(b) => b.local_force(node_i, node_j),
            Element3::ForceBeamColumn3(b) => b.local_force(),
        }
    }

    /// See `ElementOps::fiber_responses`'s doc comment.
    fn fiber_responses(&self, node_i: &Node3, node_j: &Node3) -> Option<Vec<Vec<(f64, f64)>>> {
        match self {
            Element3::DispBeamColumn3(b) => Some(b.fiber_responses(node_i, node_j)),
            Element3::ForceBeamColumn3(b) => Some(b.fiber_responses()),
            Element3::Truss3(_)
            | Element3::ZeroLength3(_)
            | Element3::ZeroLengthSection3(_)
            | Element3::ElasticBeamColumn3(_) => None,
        }
    }
}

impl ElementOps<SPATIAL_NDIM, SPATIAL_NDF, Node3Id> for Element3 {
    type Load = ElementLoad3;
    type Id = Element3Id;

    fn nodes(&self) -> NodeList<Node3Id> {
        Element3::nodes(self).into_iter().collect()
    }

    fn dof_mask(&self) -> DofMask {
        DofMask::all(SPATIAL_NDF)
    }

    fn assemble_tangent<S: TangentSink<Node3Id>>(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
        sink: &mut S,
    ) {
        let [i, j] = Element3::nodes(self);
        let (k, r) = Element3::form_tangent_and_resistance(self, nodes.get(i), nodes.get(j), load);
        sink.add(
            &two_node_dofs::<_, SPATIAL_ELEMENT_DOF>(i, j, SPATIAL_NDF),
            &k,
            &r,
        );
    }

    fn assemble_load<S: VectorSink<Node3Id>>(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
        sink: &mut S,
    ) {
        let [i, j] = Element3::nodes(self);
        let v = Element3::form_load_vector(self, nodes.get(i), nodes.get(j), load);
        sink.add(
            &two_node_dofs::<_, SPATIAL_ELEMENT_DOF>(i, j, SPATIAL_NDF),
            &v,
        );
    }

    fn assemble_mass<S: VectorSink<Node3Id>>(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        sink: &mut S,
    ) {
        let [i, j] = Element3::nodes(self);
        let v = Element3::form_mass(self, nodes.get(i), nodes.get(j));
        sink.add(
            &two_node_dofs::<_, SPATIAL_ELEMENT_DOF>(i, j, SPATIAL_NDF),
            &v,
        );
    }

    fn commit(
        &mut self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
    ) {
        let [i, j] = Element3::nodes(self);
        Element3::commit(self, nodes.get(i), nodes.get(j), load)
    }

    fn local_force_width(&self) -> usize {
        SPATIAL_ELEMENT_DOF
    }

    fn local_force(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> ElementForce {
        let [i, j] = Element3::nodes(self);
        Element3::local_force(self, nodes.get(i), nodes.get(j)).into()
    }

    fn local_load_force(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
    ) -> ElementForce {
        let [i, j] = Element3::nodes(self);
        Element3::form_local_load_vector(self, nodes.get(i), nodes.get(j), load).into()
    }

    fn fiber_responses(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> Option<Vec<Vec<(f64, f64)>>> {
        let [i, j] = Element3::nodes(self);
        Element3::fiber_responses(self, nodes.get(i), nodes.get(j))
    }
}
