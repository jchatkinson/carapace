//! The interface between `Domain` and an element catalog.
//!
//! An element reports which nodes it connects and which node DOF slots it
//! stiffens, and *pushes* its own contributions (tangent and resistance,
//! equivalent loads, lumped mass) through a sink that is generic over the
//! element's matrix size. `Domain` therefore never assumes a node count or
//! an element DOF count: a 6x6 truss and an 8x8 quadrilateral can live in
//! one catalog with no padding, and every element keeps fixed-size
//! (`SMatrix`) math and a closed-enum dispatch (no `dyn`, no heap in the
//! assembly loop).

use arrayvec::ArrayVec;
use nalgebra::{SMatrix, SVector};
use slotmap::{Key, SlotMap};

use crate::model::Node;

/// The most nodes any element in this crate connects: two-node frame
/// elements, three- and four-node planar elements (2D membranes, 3D shell
/// triangles and quadrilaterals). 3D solid elements are not planned.
pub const MAX_ELEMENT_NODES: usize = 4;

/// The nodes an element connects, in the element's own local order.
pub type NodeList<NId> = ArrayVec<NId, MAX_ELEMENT_NODES>;

/// One element DOF: `(node, slot)`, where `slot` indexes the node's
/// profile DOFs (`ux, uy, rz` in 2D; `ux..rz` in 3D).
pub type DofRef<NId> = (NId, u8);

/// Which node DOF slots an element touches, one bit per slot. One mask per
/// element, applied to every node the element connects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DofMask(u8);

impl DofMask {
    /// Every one of a node's first `ndof` slots.
    pub const fn all(ndof: usize) -> Self {
        assert!(ndof <= 8);
        DofMask(((1u16 << ndof) - 1) as u8)
    }

    /// No slots.
    pub const fn none() -> Self {
        DofMask(0)
    }

    /// This mask plus `slot`.
    pub const fn with(self, slot: usize) -> Self {
        assert!(slot < 8);
        DofMask(self.0 | (1 << slot))
    }

    /// Every slot in either mask.
    pub const fn union(self, other: DofMask) -> Self {
        DofMask(self.0 | other.0)
    }

    pub const fn contains(self, slot: usize) -> bool {
        slot < 8 && (self.0 >> slot) & 1 == 1
    }
}

/// Read access to the domain's nodes, handed to every element method that
/// needs node coordinates or state. A thin wrapper so the element interface
/// does not expose how `Domain` stores nodes.
pub struct NodeView<'a, const NDIM: usize, const NDOF: usize, NId: Key> {
    nodes: &'a SlotMap<NId, Node<NDIM, NDOF>>,
}

impl<'a, const NDIM: usize, const NDOF: usize, NId: Key> NodeView<'a, NDIM, NDOF, NId> {
    pub(crate) fn new(nodes: &'a SlotMap<NId, Node<NDIM, NDOF>>) -> Self {
        NodeView { nodes }
    }

    pub fn get(&self, id: NId) -> &'a Node<NDIM, NDOF> {
        &self.nodes[id]
    }
}

/// Receives an element's tangent stiffness `k` and internal resisting force
/// `r`, indexed by the element's DOFs `dofs` (row `a` of `k` and entry `a`
/// of `r` belong to `dofs[a]`).
pub trait TangentSink<NId> {
    fn add<const N: usize>(
        &mut self,
        dofs: &[DofRef<NId>; N],
        k: &SMatrix<f64, N, N>,
        r: &SVector<f64, N>,
    );
}

/// Receives a vector indexed by an element's DOFs: an equivalent nodal load
/// or a lumped mass.
pub trait VectorSink<NId> {
    fn add<const N: usize>(&mut self, dofs: &[DofRef<NId>; N], v: &SVector<f64, N>);
}

/// DOFs of a two-node element laid out as `[all of node i's slots, all of
/// node j's slots]`; `N` must equal `2 * ndof`.
pub fn two_node_dofs<NId: Copy, const N: usize>(i: NId, j: NId, ndof: usize) -> [DofRef<NId>; N] {
    debug_assert_eq!(N, 2 * ndof);
    std::array::from_fn(|a| {
        if a < ndof {
            (i, a as u8)
        } else {
            (j, (a - ndof) as u8)
        }
    })
}

/// An element's local nodal force (or equivalent-load) vector, for results
/// recording. Its length is `ElementOps::local_force_width`. Heap-backed:
/// this is only ever used on recorder paths, never in assembly.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementForce(Vec<f64>);

impl ElementForce {
    pub fn from_svector<const N: usize>(v: &SVector<f64, N>) -> Self {
        ElementForce(v.iter().copied().collect())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.0.iter()
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.0
    }

    /// `self -= factor * other`, component by component.
    pub fn subtract_scaled(&mut self, factor: f64, other: &ElementForce) {
        assert_eq!(self.0.len(), other.0.len());
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            *a -= factor * b;
        }
    }
}

impl<const N: usize> From<SVector<f64, N>> for ElementForce {
    fn from(v: SVector<f64, N>) -> Self {
        ElementForce::from_svector(&v)
    }
}

impl std::ops::Index<usize> for ElementForce {
    type Output = f64;
    fn index(&self, index: usize) -> &f64 {
        &self.0[index]
    }
}

/// What `Domain`'s generic assembly/state plumbing needs from an element
/// catalog. A trait rather than a shared enum because `Element` (2D) and
/// `Element3` (3D) hold different variant sets, and a single enum would
/// either duplicate every variant or need `Box<dyn Trait>` in the assembly
/// loop. `NDIM`/`NDOF` are separate const generics (not computed from each
/// other) because stable Rust cannot evaluate `2 * NDOF` in a generic item.
pub trait ElementOps<const NDIM: usize, const NDOF: usize, NId: Key> {
    /// An element-load kind this catalog supports. `Add` so several loads on
    /// one element in one pattern accumulate.
    type Load: Copy
        + crate::model::ElementLoadComponents
        + std::ops::Add<Output = Self::Load>
        + std::ops::Mul<f64, Output = Self::Load>;

    /// This catalog's own `Domain` element-store key. An associated type
    /// rather than a further generic parameter of `Domain`, so it is always
    /// determined by `E` and never an independently unconstrained type
    /// variable (which broke inference at every `Domain::new()`).
    type Id: Key;

    /// The nodes this element connects, in local order.
    fn nodes(&self) -> NodeList<NId>;

    /// The node DOF slots this element stiffens (the same for each node).
    /// Must cover every slot on which `assemble_tangent` can produce a
    /// nonzero stiffness or resistance: a slot outside the mask that no
    /// other element, constraint or boundary condition uses is not an
    /// equation at all, so such a contribution would be silently dropped
    /// (debug builds assert against it).
    fn dof_mask(&self) -> DofMask;

    /// Rejects an element whose geometry or parameters make it unusable
    /// (coincident nodes, inverted orientation, ...); called by
    /// `Domain::validate`. The default accepts everything.
    fn validate(&self, _nodes: &NodeView<'_, NDIM, NDOF, NId>) -> Result<(), &'static str> {
        Ok(())
    }

    /// Tangent and internal resisting force at the current nodal state.
    /// `load` is the effective element load at the pseudo-time being
    /// assembled (every pattern's load on this element, scaled by its
    /// factor). Only elements whose internal state depends on the load
    /// itself (`ForceBeamColumn`) read it; for every other element loads
    /// enter only through `assemble_load`.
    fn assemble_tangent<S: TangentSink<NId>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, NId>,
        load: Option<&Self::Load>,
        sink: &mut S,
    );

    /// The equivalent nodal load vector of `load` (zero if `load` is
    /// `None` or of a kind this element does not carry).
    fn assemble_load<S: VectorSink<NId>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, NId>,
        load: Option<&Self::Load>,
        sink: &mut S,
    );

    /// Lumped (diagonal) mass, per element DOF.
    fn assemble_mass<S: VectorSink<NId>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, NId>,
        sink: &mut S,
    );

    /// Commit this element's state at the converged nodal state. `load` is
    /// the effective element load at the committed pseudo-time.
    fn commit(&mut self, nodes: &NodeView<'_, NDIM, NDOF, NId>, load: Option<&Self::Load>);

    /// Length of `local_force`/`local_load_force`; the bound a recorder's
    /// `component` index is checked against.
    fn local_force_width(&self) -> usize;

    /// This element's local nodal force at its current committed state,
    /// in the element's own axis frame (the original fixed orientation for
    /// `Linear`/`PDelta`, the current chord for `Corotational`). Never
    /// mutates state.
    fn local_force(&self, nodes: &NodeView<'_, NDIM, NDOF, NId>) -> ElementForce;

    /// `assemble_load`'s equivalent nodal load in the same local frame as
    /// `local_force`; what `Domain::element_end_force` subtracts so a
    /// loaded member's end forces include the load's fixed-end effect.
    fn local_load_force(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, NId>,
        load: Option<&Self::Load>,
    ) -> ElementForce;

    /// Every integration point's per-fiber `(strain, stress)`; `None` for
    /// every element that is not fiber-discretized.
    fn fiber_responses(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, NId>,
    ) -> Option<Vec<Vec<(f64, f64)>>>;
}
