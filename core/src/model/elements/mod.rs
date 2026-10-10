use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

use super::{
    ElementLoad, ElementLoad3, Node, Node3, Node3Id, NodeId, ELEMENT_DOF, NDF, PLANAR_NDIM,
    SPATIAL_ELEMENT_DOF, SPATIAL_NDF, SPATIAL_NDIM,
};

mod ops;
pub use ops::{
    two_node_dofs, DofMask, DofRef, ElementForce, ElementOps, GaussResponse, NodeList, NodeView,
    ShellResponse, TangentSink, VectorSink, MAX_ELEMENT_NODES,
};

mod disp_beam_column;
mod elastic_beam_column;
mod force_beam_column;
mod plane_common;
mod quad4;
mod shell3;
mod shell4;
mod tri3;
mod truss;
mod uniform_load;
mod zero_length;

pub use disp_beam_column::{DispBeamColumn, DispBeamColumn3};
pub use elastic_beam_column::{ElasticBeamColumn, ElasticBeamColumn3};
pub use force_beam_column::{ForceBeamColumn, ForceBeamColumn3};
pub use quad4::{Quad4, Quad4Formulation};
pub use shell3::Shell3;
pub use shell4::Shell4;
pub use tri3::Tri3;
pub use truss::{SpatialElementMatrix, SpatialElementVector, Truss, Truss3};
pub use zero_length::{
    Friction, Friction3, Orientation, OrientationError, ZeroLength, ZeroLength3, ZeroLengthSection,
    ZeroLengthSection3,
};

new_key_type! {
    /// Generational index into `Domain`'s element store.
    pub struct ElementId;

    /// Generational index into `Domain3`'s element store. Distinct from
    /// `ElementId` for the same reason `Node3Id` is distinct from `NodeId`
    /// (see its doc comment).
    pub struct Element3Id;
}

/// Element catalog. Closed enum, `match`-based dispatch, no `Box<dyn Trait>`
///.
#[derive(Debug, Clone)]
pub enum Element {
    Truss(Truss),
    ZeroLength(ZeroLength),
    ZeroLengthSection(ZeroLengthSection),
    ElasticBeamColumn(ElasticBeamColumn),
    DispBeamColumn(DispBeamColumn),
    ForceBeamColumn(ForceBeamColumn),
    Tri3(Tri3),
    Quad4(Quad4),
}

const CONTINUUM: &str = "continuum elements assemble through ElementOps, not the two-node path";

impl Element {
    /// The nodes this element connects, in local order.
    pub fn nodes(&self) -> NodeList<NodeId> {
        match self {
            Element::Tri3(t) => t.nodes.into_iter().collect(),
            Element::Quad4(q) => q.nodes.into_iter().collect(),
            _ => self.pair().into_iter().collect(),
        }
    }

    /// The two end nodes of a two-node element.
    fn pair(&self) -> [NodeId; 2] {
        match self {
            Element::Tri3(_) | Element::Quad4(_) => unreachable!("{CONTINUUM}"),
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
            Element::Tri3(_) | Element::Quad4(_) => unreachable!("{CONTINUUM}"),
        }
    }

    /// This element's equivalent nodal load vector (global coordinates,
    /// same DOF order as above) from `load` — the `ElementLoad` (if any)
    /// that whichever `LoadPattern` is currently being assembled has on
    /// this element — e.g. a beam-column's distributed transverse
    /// load. Zero when `load` is `None`, and for elements with no
    /// element-load support at all (`Truss`, `ZeroLength`).
    fn form_load_vector(
        &self,
        node_i: &Node,
        node_j: &Node,
        load: Option<&ElementLoad>,
    ) -> SVector<f64, ELEMENT_DOF> {
        match (self, load) {
            (Element::ElasticBeamColumn(b), Some(load)) => {
                b.form_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
            }
            (Element::DispBeamColumn(b), Some(load)) => {
                b.form_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
            }
            (Element::ForceBeamColumn(b), Some(load)) => {
                b.form_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
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
            (Element::ElasticBeamColumn(b), Some(load)) => {
                b.form_local_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
            }
            (Element::DispBeamColumn(b), Some(load)) => {
                b.form_local_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
            }
            (Element::ForceBeamColumn(b), Some(load)) => {
                b.form_local_load_vector(node_i, node_j, load.beam_uniform[0], load.beam_uniform[1])
            }
            _ => SVector::<f64, ELEMENT_DOF>::zeros(),
        }
    }

    /// This element's lumped-mass contribution (diagonal only, lumped) in the same local DOF
    /// order, geometry-dependent (`length`) so it needs both nodes. Zero for `ZeroLength` (a spring/
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
            Element::Tri3(_) | Element::Quad4(_) => unreachable!("{CONTINUUM}"),
        }
    }

    /// Commit this element's material(s) at the given (final, converged)
    /// node state — see `Material`'s doc comment for why this is the only
    /// place a `Material` ever mutates. A no-op for `ElasticBeamColumn`
    /// (no `Material` — its response is closed-form) and for any
    /// `ZeroLength` direction with no material assigned.
    fn commit(&mut self, node_i: &Node, node_j: &Node, load: Option<&ElementLoad>) {
        match self {
            Element::Truss(t) => t.commit(node_i, node_j),
            Element::ZeroLength(z) => z.commit(node_i, node_j),
            Element::ZeroLengthSection(z) => z.commit(node_i, node_j),
            Element::ElasticBeamColumn(_) => {}
            Element::DispBeamColumn(b) => b.commit(node_i, node_j),
            Element::ForceBeamColumn(b) => b.commit(node_i, node_j, load),
            Element::Tri3(_) | Element::Quad4(_) => unreachable!("{CONTINUUM}"),
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
            Element::Tri3(_) | Element::Quad4(_) => unreachable!("{CONTINUUM}"),
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
            | Element::ElasticBeamColumn(_)
            | Element::Tri3(_)
            | Element::Quad4(_) => None,
        }
    }
}

impl ElementOps<PLANAR_NDIM, NDF, NodeId> for Element {
    type Load = ElementLoad;
    type Id = ElementId;

    fn nodes(&self) -> NodeList<NodeId> {
        Element::nodes(self)
    }

    fn dof_mask(&self) -> DofMask {
        match self {
            // A truss carries no moment: translations only.
            Element::Truss(_) => DofMask::none().with(0).with(1),
            Element::ZeroLength(z) => z.dof_mask(),
            Element::ZeroLengthSection(z) => z.dof_mask(),
            Element::ElasticBeamColumn(_)
            | Element::DispBeamColumn(_)
            | Element::ForceBeamColumn(_) => DofMask::all(NDF),
            Element::Tri3(t) => t.dof_mask(),
            Element::Quad4(q) => q.dof_mask(),
        }
    }

    fn validate(&self, nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>) -> Result<(), &'static str> {
        match self {
            Element::Tri3(t) => t.validate(nodes),
            Element::Quad4(q) => q.validate(nodes),
            _ => Ok(()),
        }
    }

    fn accepts_load(&self, load: &ElementLoad) -> bool {
        match self {
            Element::ElasticBeamColumn(_)
            | Element::DispBeamColumn(_)
            | Element::ForceBeamColumn(_) => !load.has_body() && !load.has_edge(),
            Element::Tri3(_) => !load.has_beam() && !load.has_edge_from(3),
            Element::Quad4(_) => !load.has_beam(),
            Element::Truss(_) | Element::ZeroLength(_) | Element::ZeroLengthSection(_) => {
                !load.has_beam() && !load.has_body() && !load.has_edge()
            }
        }
    }

    fn prepare(&mut self, nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>) {
        match self {
            Element::Tri3(t) => t.prepare(nodes),
            Element::Quad4(q) => q.prepare(nodes),
            _ => {}
        }
    }

    fn gauss_point_count(&self) -> usize {
        match self {
            Element::Tri3(_) => 1,
            Element::Quad4(_) => 4,
            _ => 0,
        }
    }

    fn gauss_responses(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
    ) -> Option<Vec<GaussResponse>> {
        match self {
            Element::Tri3(t) => Some(t.gauss_responses(nodes)),
            Element::Quad4(q) => Some(q.gauss_responses(nodes)),
            _ => None,
        }
    }

    fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        match self {
            Element::Tri3(t) => return t.assemble_tangent(nodes, sink),
            Element::Quad4(q) => return q.assemble_tangent(nodes, sink),
            _ => {}
        }
        let [i, j] = Element::pair(self);
        let (k, r) = Element::form_tangent_and_resistance(self, nodes.get(i), nodes.get(j), load);
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &k, &r);
    }

    fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        match (self, load) {
            (Element::Tri3(t), load) => {
                return load.map_or((), |load| t.assemble_load(nodes, load, sink))
            }
            (Element::Quad4(q), load) => {
                return load.map_or((), |load| q.assemble_load(nodes, load, sink))
            }
            _ => {}
        }
        let [i, j] = Element::pair(self);
        let v = Element::form_load_vector(self, nodes.get(i), nodes.get(j), load);
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &v);
    }

    fn assemble_mass<S: VectorSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        sink: &mut S,
    ) {
        match self {
            Element::Tri3(t) => return t.assemble_mass(nodes, sink),
            Element::Quad4(q) => return q.assemble_mass(nodes, sink),
            _ => {}
        }
        let [i, j] = Element::pair(self);
        let v = Element::form_mass(self, nodes.get(i), nodes.get(j));
        sink.add(&two_node_dofs::<_, ELEMENT_DOF>(i, j, NDF), &v);
    }

    fn commit(
        &mut self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
    ) {
        match self {
            Element::Tri3(t) => return t.commit(nodes),
            Element::Quad4(q) => return q.commit(nodes),
            _ => {}
        }
        let [i, j] = Element::pair(self);
        Element::commit(self, nodes.get(i), nodes.get(j), load)
    }

    fn local_force_width(&self) -> usize {
        match self {
            Element::Tri3(_) => 6,
            Element::Quad4(_) => 8,
            _ => ELEMENT_DOF,
        }
    }

    fn local_force(&self, nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>) -> ElementForce {
        match self {
            Element::Tri3(t) => return t.local_force(nodes),
            Element::Quad4(q) => return q.local_force(nodes),
            _ => {}
        }
        let [i, j] = Element::pair(self);
        Element::local_force(self, nodes.get(i), nodes.get(j)).into()
    }

    fn local_load_force(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
        load: Option<&ElementLoad>,
    ) -> ElementForce {
        match self {
            Element::Tri3(t) => {
                return load
                    .map_or(SVector::<f64, 6>::zeros(), |l| t.load_vector(nodes, l))
                    .into()
            }
            Element::Quad4(q) => {
                return load
                    .map_or(SVector::<f64, 8>::zeros(), |l| q.load_vector(nodes, l))
                    .into()
            }
            _ => {}
        }
        let [i, j] = Element::pair(self);
        Element::form_local_load_vector(self, nodes.get(i), nodes.get(j), load).into()
    }

    fn fiber_responses(
        &self,
        nodes: &NodeView<'_, PLANAR_NDIM, NDF, NodeId>,
    ) -> Option<Vec<Vec<(f64, f64)>>> {
        if let Element::Tri3(_) | Element::Quad4(_) = self {
            return None;
        }
        let [i, j] = Element::pair(self);
        Element::fiber_responses(self, nodes.get(i), nodes.get(j))
    }
}
/// 3D element catalog — `Domain3`'s counterpart to `Element`. Closed
/// enum, `match`-based dispatch, same reasoning as `Element`.
#[derive(Debug, Clone)]
pub enum Element3 {
    Truss3(Truss3),
    ZeroLength3(ZeroLength3),
    ZeroLengthSection3(ZeroLengthSection3),
    ElasticBeamColumn3(ElasticBeamColumn3),
    DispBeamColumn3(DispBeamColumn3),
    ForceBeamColumn3(ForceBeamColumn3),
    Shell3(Shell3),
    Shell4(Shell4),
}

const SHELL: &str = "shell elements assemble through ElementOps, not the two-node path";

impl Element3 {
    /// The nodes this element connects, in local order.
    pub fn nodes(&self) -> NodeList<Node3Id> {
        match self {
            Element3::Shell3(s) => s.nodes.into_iter().collect(),
            Element3::Shell4(s) => s.nodes.into_iter().collect(),
            _ => self.pair().into_iter().collect(),
        }
    }

    /// The two end nodes of a two-node element.
    fn pair(&self) -> [Node3Id; 2] {
        match self {
            Element3::Shell3(_) | Element3::Shell4(_) => unreachable!("{SHELL}"),
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
            Element3::Shell3(_) | Element3::Shell4(_) => unreachable!("{SHELL}"),
        }
    }

    /// This element's equivalent nodal load vector from `load` — the
    /// 3D counterpart of `Element::form_load_vector`. Zero when
    /// `load` is `None`, and for elements with no element-load support
    /// (`Truss3`, `ZeroLength3`).
    fn form_load_vector(
        &self,
        node_i: &Node3,
        node_j: &Node3,
        load: Option<&ElementLoad3>,
    ) -> SVector<f64, SPATIAL_ELEMENT_DOF> {
        match (self, load) {
            (Element3::ElasticBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_load_vector(node_i, node_j, wx, wy, wz)
            }
            (Element3::DispBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_load_vector(node_i, node_j, wx, wy, wz)
            }
            (Element3::ForceBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_load_vector(node_i, node_j, wx, wy, wz)
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
            (Element3::ElasticBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_local_load_vector(node_i, node_j, wx, wy, wz)
            }
            (Element3::DispBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_local_load_vector(node_i, node_j, wx, wy, wz)
            }
            (Element3::ForceBeamColumn3(b), Some(l)) => {
                let [wx, wy, wz] = l.uniform;
                b.form_local_load_vector(node_i, node_j, wx, wy, wz)
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
            Element3::Shell3(_) | Element3::Shell4(_) => unreachable!("{SHELL}"),
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
            Element3::Shell3(_) | Element3::Shell4(_) => unreachable!("{SHELL}"),
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
            Element3::Shell3(_) | Element3::Shell4(_) => unreachable!("{SHELL}"),
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
            | Element3::ElasticBeamColumn3(_)
            | Element3::Shell3(_)
            | Element3::Shell4(_) => None,
        }
    }
}

impl ElementOps<SPATIAL_NDIM, SPATIAL_NDF, Node3Id> for Element3 {
    type Load = ElementLoad3;
    type Id = Element3Id;

    fn nodes(&self) -> NodeList<Node3Id> {
        Element3::nodes(self)
    }

    fn dof_mask(&self) -> DofMask {
        match self {
            Element3::Truss3(_) => DofMask::none().with(0).with(1).with(2),
            Element3::ZeroLength3(z) => z.dof_mask(),
            Element3::ZeroLengthSection3(z) => z.dof_mask(),
            Element3::ElasticBeamColumn3(_)
            | Element3::DispBeamColumn3(_)
            | Element3::ForceBeamColumn3(_) => DofMask::all(SPATIAL_NDF),
            Element3::Shell3(s) => s.dof_mask(),
            Element3::Shell4(s) => s.dof_mask(),
        }
    }

    fn validate(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> Result<(), &'static str> {
        match self {
            Element3::Shell3(s) => s.validate(nodes),
            Element3::Shell4(s) => s.validate(nodes),
            _ => Ok(()),
        }
    }

    fn accepts_load(&self, load: &ElementLoad3) -> bool {
        match self {
            Element3::Shell3(_) | Element3::Shell4(_) => !load.has_uniform(),
            _ => !load.has_body() && !load.has_pressure(),
        }
    }

    fn prepare(&mut self, nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>) {
        match self {
            Element3::Shell3(s) => s.prepare(nodes),
            Element3::Shell4(s) => s.prepare(nodes),
            _ => {}
        }
    }

    fn gauss_point_count(&self) -> usize {
        match self {
            Element3::Shell3(_) | Element3::Shell4(_) => 4,
            _ => 0,
        }
    }

    fn gauss_component_count(&self) -> usize {
        match self {
            Element3::Shell3(_) | Element3::Shell4(_) => 8,
            _ => 3,
        }
    }

    fn shell_responses(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> Option<Vec<ShellResponse>> {
        match self {
            Element3::Shell3(s) => Some(s.shell_responses(nodes)),
            Element3::Shell4(s) => Some(s.shell_responses(nodes)),
            _ => None,
        }
    }

    fn assemble_tangent<S: TangentSink<Node3Id>>(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
        sink: &mut S,
    ) {
        match self {
            Element3::Shell3(s) => return s.assemble_tangent(nodes, sink),
            Element3::Shell4(s) => return s.assemble_tangent(nodes, sink),
            _ => {}
        }
        let [i, j] = Element3::pair(self);
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
        match self {
            Element3::Shell3(s) => {
                return load.map_or((), |load| s.assemble_load(nodes, load, sink))
            }
            Element3::Shell4(s) => {
                return load.map_or((), |load| s.assemble_load(nodes, load, sink))
            }
            _ => {}
        }
        let [i, j] = Element3::pair(self);
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
        match self {
            Element3::Shell3(s) => return s.assemble_mass(nodes, sink),
            Element3::Shell4(s) => return s.assemble_mass(nodes, sink),
            _ => {}
        }
        let [i, j] = Element3::pair(self);
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
        match self {
            Element3::Shell3(s) => return s.commit(nodes),
            Element3::Shell4(s) => return s.commit(nodes),
            _ => {}
        }
        let [i, j] = Element3::pair(self);
        Element3::commit(self, nodes.get(i), nodes.get(j), load)
    }

    fn local_force_width(&self) -> usize {
        match self {
            Element3::Shell3(_) => 18,
            Element3::Shell4(_) => 24,
            _ => SPATIAL_ELEMENT_DOF,
        }
    }

    fn local_force(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> ElementForce {
        match self {
            Element3::Shell3(s) => return s.local_force(nodes),
            Element3::Shell4(s) => return s.local_force(nodes),
            _ => {}
        }
        let [i, j] = Element3::pair(self);
        Element3::local_force(self, nodes.get(i), nodes.get(j)).into()
    }

    fn local_load_force(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
        load: Option<&ElementLoad3>,
    ) -> ElementForce {
        match self {
            Element3::Shell3(s) => {
                return load
                    .map_or(SVector::<f64, 18>::zeros(), |l| s.load_vector(nodes, l))
                    .into()
            }
            Element3::Shell4(s) => {
                return load
                    .map_or(SVector::<f64, 24>::zeros(), |l| s.load_vector(nodes, l))
                    .into()
            }
            _ => {}
        }
        let [i, j] = Element3::pair(self);
        Element3::form_local_load_vector(self, nodes.get(i), nodes.get(j), load).into()
    }

    fn fiber_responses(
        &self,
        nodes: &NodeView<'_, SPATIAL_NDIM, SPATIAL_NDF, Node3Id>,
    ) -> Option<Vec<Vec<(f64, f64)>>> {
        if let Element3::Shell3(_) | Element3::Shell4(_) = self {
            return None;
        }
        let [i, j] = Element3::pair(self);
        Element3::fiber_responses(self, nodes.get(i), nodes.get(j))
    }
}

/// `ElementOps::dof_mask` must cover every node DOF slot on which an element
/// produces stiffness or resistance: a slot outside the mask that nothing
/// else uses is not an equation, so the contribution would be silently
/// dropped. This drives every element kind at a generic nonzero state and
/// compares the slots it actually stiffens with its declared mask. Any new
/// element kind needs an entry here.
#[cfg(test)]
mod mask_conformance {
    use std::collections::BTreeSet;

    use slotmap::{Key, SlotMap};

    use super::*;
    use crate::model::{
        BeamIntegration, ElasticBeamColumn, ElasticBeamColumn3, Fiber, Fiber3, FiberSection,
        FiberSection3, ForceBeamColumn, ForceBeamColumn3, Friction, Friction3, GeomTransf,
        GeomTransf3, Material, Orientation, Truss, Truss3, ZeroLength, ZeroLength3,
        ZeroLengthSection, ZeroLengthSection3,
    };

    /// Records which `(node, slot)` an element's tangent or resistance is
    /// nonzero on, relative to the element's own largest entry.
    struct Touched<NId> {
        slots: BTreeSet<u8>,
        _marker: std::marker::PhantomData<NId>,
    }

    impl<NId: Key> TangentSink<NId> for Touched<NId> {
        fn add<const N: usize>(
            &mut self,
            dofs: &[DofRef<NId>; N],
            k: &SMatrix<f64, N, N>,
            r: &SVector<f64, N>,
        ) {
            let k_scale = (0..N).map(|i| k[(i, i)].abs()).fold(0.0, f64::max);
            let r_scale = r.iter().map(|v| v.abs()).fold(0.0, f64::max);
            for (a, &(_, slot)) in dofs.iter().enumerate() {
                if k[(a, a)].abs() > 1e-10 * k_scale || r[a].abs() > 1e-10 * r_scale {
                    self.slots.insert(slot);
                }
            }
        }
    }

    fn touched<const NDIM: usize, const NDOF: usize, NId: Key, E: ElementOps<NDIM, NDOF, NId>>(
        element: &E,
        nodes: &SlotMap<NId, Node<NDIM, NDOF>>,
    ) -> BTreeSet<u8> {
        let mut sink = Touched::<NId> {
            slots: BTreeSet::new(),
            _marker: std::marker::PhantomData,
        };
        element.assemble_tangent(&NodeView::new(nodes), None, &mut sink);
        sink.slots
    }

    fn mask_slots<const NDOF: usize>(mask: DofMask) -> BTreeSet<u8> {
        (0..NDOF)
            .filter(|&s| mask.contains(s))
            .map(|s| s as u8)
            .collect()
    }

    /// `exact`: the mask must equal the touched slots (no needless
    /// activation); otherwise it only has to cover them.
    fn check<const NDIM: usize, const NDOF: usize, NId: Key, E: ElementOps<NDIM, NDOF, NId>>(
        name: &str,
        element: &E,
        nodes: &SlotMap<NId, Node<NDIM, NDOF>>,
        exact: bool,
    ) {
        let touched = touched::<NDIM, NDOF, NId, E>(element, nodes);
        let mask = mask_slots::<NDOF>(element.dof_mask());
        assert!(
            touched.is_subset(&mask),
            "{name}: stiffens slots {touched:?} but declares only {mask:?}"
        );
        if exact {
            assert_eq!(
                touched, mask,
                "{name}: declared mask is wider than what it stiffens"
            );
        }
    }

    /// Two 2D nodes with a generic nonzero state on every DOF, skew to the axes
    /// (so a truss uses both translations).
    fn nodes_2d() -> (SlotMap<NodeId, Node>, NodeId, NodeId) {
        let mut nodes = SlotMap::with_key();
        let a = nodes.insert(Node::new([0.0, 0.0]));
        let mut j = Node::new([300.0, 200.0]);
        j.displacement = [0.4, -0.7, 0.003];
        let b = nodes.insert(j);
        (nodes, a, b)
    }

    fn nodes_3d() -> (SlotMap<Node3Id, Node3>, Node3Id, Node3Id) {
        let mut nodes = SlotMap::with_key();
        let a = nodes.insert(Node3::new([0.0, 0.0, 0.0]));
        let mut j = Node3::new([300.0, 200.0, 150.0]);
        j.displacement = [0.4, -0.7, 0.3, 0.002, -0.003, 0.004];
        let b = nodes.insert(j);
        (nodes, a, b)
    }

    fn steel() -> Material {
        Material::Elastic { e: 200_000.0 }
    }

    /// A skew four-node patch with a generic nonzero displacement, for the continuum elements.
    fn nodes_plane() -> (SlotMap<NodeId, Node>, [NodeId; 4]) {
        let mut nodes = SlotMap::with_key();
        let coords = [[0.0, 0.0], [2.0, 0.2], [2.3, 1.7], [-0.1, 1.4]];
        let ids = coords.map(|c| {
            let mut node = Node::new(c);
            node.displacement = [0.01 * (c[1] + 1.0), -0.02 * (c[0] + 0.5), 0.003];
            nodes.insert(node)
        });
        (nodes, ids)
    }

    #[test]
    fn continuum_elements_declare_exactly_the_translations_they_stiffen() {
        let (nodes, [a, b, c, d]) = nodes_plane();
        let material = || crate::model::PlaneMaterial::plane_stress(30e3, 0.2).unwrap();
        check(
            "Tri3",
            &Element::Tri3(Tri3::new(a, b, c, 0.5, material())),
            &nodes,
            true,
        );
        check(
            "Quad4",
            &Element::Quad4(Quad4::new([a, b, c, d], 0.5, material())),
            &nodes,
            true,
        );
        check(
            "Quad4 enhanced",
            &Element::Quad4(
                Quad4::new([a, b, c, d], 0.5, material())
                    .with_formulation(crate::model::Quad4Formulation::Enhanced),
            ),
            &nodes,
            true,
        );
    }

    #[test]
    fn shell_declares_all_six_slots_it_stiffens() {
        let mut nodes = SlotMap::with_key();
        let coords = [
            [0.0, 0.0, 0.0],
            [2.0, 0.2, 0.3],
            [2.3, 1.7, 0.9],
            [-0.1, 1.4, 0.2],
        ];
        let ids = coords.map(|c| {
            let mut node = Node3::new(c);
            node.displacement = [0.01, -0.02, 0.03, 0.002, -0.003, 0.004];
            nodes.insert(node)
        });
        let section =
            crate::model::ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 0.0).unwrap();
        check(
            "Shell4",
            &Element3::Shell4(Shell4::new(ids, section.clone())),
            &nodes,
            true,
        );
        check(
            "Shell3",
            &Element3::Shell3(Shell3::new([ids[0], ids[1], ids[2]], section)),
            &nodes,
            true,
        );
    }

    #[test]
    fn planar_elements_declare_the_slots_they_stiffen() {
        let (nodes, a, b) = nodes_2d();
        let fibers = || {
            vec![
                Fiber::new(10.0, 100.0, steel()),
                Fiber::new(-10.0, 100.0, steel()),
            ]
        };
        let rot = 30.0_f64.to_radians();
        let oriented =
            || Orientation::new([rot.cos(), rot.sin(), 0.0], [-rot.sin(), rot.cos(), 0.0]).unwrap();

        check(
            "Truss",
            &Element::Truss(Truss::new(a, b, 10.0, steel())),
            &nodes,
            true,
        );
        check(
            "ZeroLength ux only",
            &Element::ZeroLength(ZeroLength::new(a, b).with_material(0, steel())),
            &nodes,
            true,
        );
        check(
            "ZeroLength rz only",
            &Element::ZeroLength(ZeroLength::new(a, b).with_material(2, steel())),
            &nodes,
            true,
        );
        check(
            "ZeroLength oriented",
            &Element::ZeroLength(
                ZeroLength::new(a, b)
                    .with_material(0, steel())
                    .with_material(1, steel())
                    .with_orientation(oriented())
                    .unwrap(),
            ),
            &nodes,
            true,
        );
        check(
            "ZeroLength friction",
            &Element::ZeroLength(
                ZeroLength::new(a, b)
                    .with_material(0, steel())
                    .with_friction(Friction::new(0, 1, 0.3, 500.0, 0.01)),
            ),
            &nodes,
            true,
        );
        check(
            "ZeroLengthSection",
            &Element::ZeroLengthSection(ZeroLengthSection::new(a, b, FiberSection::new(fibers()))),
            &nodes,
            true,
        );
        check(
            "ZeroLengthSection + shear spring",
            &Element::ZeroLengthSection(
                ZeroLengthSection::new(a, b, FiberSection::new(fibers())).with_material(1, steel()),
            ),
            &nodes,
            true,
        );
        check(
            "ElasticBeamColumn",
            &Element::ElasticBeamColumn(ElasticBeamColumn::new(
                a,
                b,
                200_000.0,
                1e4,
                1e8,
                GeomTransf::PDelta,
            )),
            &nodes,
            true,
        );
        check(
            "DispBeamColumn",
            &Element::DispBeamColumn(DispBeamColumn::new(
                a,
                b,
                fibers(),
                BeamIntegration::Lobatto { points: 3 },
            )),
            &nodes,
            true,
        );
        check(
            "ForceBeamColumn",
            &Element::ForceBeamColumn(ForceBeamColumn::new(
                a,
                b,
                fibers(),
                BeamIntegration::Lobatto { points: 3 },
            )),
            &nodes,
            true,
        );
    }

    #[test]
    fn spatial_elements_declare_the_slots_they_stiffen() {
        let (nodes, a, b) = nodes_3d();
        let fibers = || {
            vec![
                Fiber3::new(10.0, 10.0, 50.0, steel()),
                Fiber3::new(10.0, -10.0, 50.0, steel()),
                Fiber3::new(-10.0, 10.0, 50.0, steel()),
                Fiber3::new(-10.0, -10.0, 50.0, steel()),
            ]
        };
        let tilt = Orientation::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();

        check(
            "Truss3",
            &Element3::Truss3(Truss3::new(a, b, 10.0, steel())),
            &nodes,
            true,
        );
        check(
            "ZeroLength3 translations",
            &Element3::ZeroLength3(
                ZeroLength3::new(a, b)
                    .with_material(0, steel())
                    .with_material(1, steel())
                    .with_material(2, steel()),
            ),
            &nodes,
            true,
        );
        check(
            "ZeroLength3 torsion only",
            &Element3::ZeroLength3(ZeroLength3::new(a, b).with_material(3, steel())),
            &nodes,
            true,
        );
        check(
            "ZeroLength3 oriented",
            &Element3::ZeroLength3(
                ZeroLength3::new(a, b)
                    .with_material(0, steel())
                    .with_orientation(tilt.clone()),
            ),
            &nodes,
            false,
        );
        check(
            "ZeroLength3 friction",
            &Element3::ZeroLength3(
                ZeroLength3::new(a, b)
                    .with_material(0, steel())
                    .with_friction(Friction3::new(0, [1, 2], 0.3, 500.0, 0.01)),
            ),
            &nodes,
            true,
        );
        check(
            "ZeroLengthSection3",
            &Element3::ZeroLengthSection3(ZeroLengthSection3::new(
                a,
                b,
                FiberSection3::new(fibers()),
            )),
            &nodes,
            true,
        );
        check(
            "ZeroLengthSection3 + torsion spring",
            &Element3::ZeroLengthSection3(
                ZeroLengthSection3::new(a, b, FiberSection3::new(fibers()))
                    .with_material(3, steel()),
            ),
            &nodes,
            true,
        );
        check(
            "ElasticBeamColumn3",
            &Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
                a,
                b,
                200_000.0,
                80_000.0,
                1e4,
                1e7,
                5e7,
                1e8,
                GeomTransf3::Linear3 {
                    vec_xz: [0.0, 0.0, 1.0],
                },
            )),
            &nodes,
            true,
        );
        check(
            "DispBeamColumn3",
            &Element3::DispBeamColumn3(DispBeamColumn3::new(
                a,
                b,
                80_000.0,
                1e7,
                [0.0, 0.0, 1.0],
                fibers(),
                BeamIntegration::Lobatto { points: 3 },
            )),
            &nodes,
            true,
        );
        check(
            "ForceBeamColumn3",
            &Element3::ForceBeamColumn3(ForceBeamColumn3::new(
                a,
                b,
                80_000.0,
                1e7,
                [0.0, 0.0, 1.0],
                fibers(),
                BeamIntegration::Lobatto { points: 3 },
            )),
            &nodes,
            true,
        );
    }
}
