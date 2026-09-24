use slotmap::new_key_type;

/// Planar coordinate components `[x, y]`.
pub const PLANAR_NDIM: usize = 2;

/// Spatial coordinate components `[x, y, z]`.
pub const SPATIAL_NDIM: usize = 3;

/// Spatial frame degrees of freedom `[ux, uy, uz, rx, ry, rz]`.
pub const SPATIAL_NDF: usize = 6;

/// Two-node spatial element DOFs: six at node i followed by six at node j.
pub const SPATIAL_ELEMENT_DOF: usize = 2 * SPATIAL_NDF;

/// Named spatial DOF indices, in the order used by every future spatial
/// element/result/worker record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum SpatialDof {
    Ux = 0,
    Uy = 1,
    Uz = 2,
    Rx = 3,
    Ry = 4,
    Rz = 5,
}

/// A global coordinate axis (`0..SPATIAL_NDIM`) — used where a spatial API
/// needs to name one of `[x, y, z]` itself rather than a DOF (e.g.
/// `Domain3::rigid_diaphragm_about`'s diaphragm-normal parameter), where a
/// bare `usize` would silently also accept a rotational `SpatialDof` index.
/// `as usize` gives the matching translational `SpatialDof`/coordinate
/// index directly (`Axis3::Y as usize == SpatialDof::Uy as usize == 1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum Axis3 {
    X = 0,
    Y = 1,
    Z = 2,
}

new_key_type! {
    /// Generational index into the planar `Domain` node store. Stays valid
    /// across removals of other nodes; using a stale key against a different
    /// domain is a logic error the type cannot catch, but a use-after-free
    /// within one domain's lifetime is.
    pub struct NodeId;

    /// Generational index for a future spatial `Domain3` node store. It is
    /// distinct from `NodeId`, so planar elements cannot connect to spatial
    /// nodes at compile time.
    pub struct Node3Id;
}

/// A node in a concrete kinematic profile. Coordinates have `NDIM`
/// components; state, boundary conditions, and mass have `NDOF` components.
/// `Node` defaults to the existing planar `(2, 3)` profile, so its public API
/// remains source-compatible; `Node3` is the spatial `(3, 6)` specialization.
///
/// `velocity`/`acceleration` are only meaningful for/written by
/// `TransientAnalysis` — static `Analysis` never touches them. Applied load
/// is not here; it belongs to `Domain` load patterns because one node can
/// carry independently scaled loads across phases.
#[derive(Debug, Clone)]
pub struct Node<const NDIM: usize = PLANAR_NDIM, const NDOF: usize = 3> {
    pub coords: [f64; NDIM],
    pub fixed: [bool; NDOF],
    pub displacement: [f64; NDOF],
    pub velocity: [f64; NDOF],
    pub acceleration: [f64; NDOF],
    /// Lumped nodal mass per DOF — a user-assigned point mass, additive with
    /// any element-consistent lumped mass.
    pub mass: [f64; NDOF],
    /// Equation number for each free DOF, assigned by the owning domain.
    pub(crate) equation: [Option<usize>; NDOF],
}

/// Explicit planar node specialization, useful where a profile must be named
/// rather than inferred from `Node`'s default generic parameters.
pub type Node2 = Node<PLANAR_NDIM, 3>;

/// Explicit spatial node specialization. It shares the node implementation
/// with planar `Node`; only fixed array dimensions differ.
pub type Node3 = Node<SPATIAL_NDIM, SPATIAL_NDF>;

impl<const NDIM: usize, const NDOF: usize> Node<NDIM, NDOF> {
    pub fn new(coords: [f64; NDIM]) -> Self {
        Node {
            coords,
            fixed: [false; NDOF],
            displacement: [0.0; NDOF],
            velocity: [0.0; NDOF],
            acceleration: [0.0; NDOF],
            mass: [0.0; NDOF],
            equation: [None; NDOF],
        }
    }

    pub fn fix(mut self, dof: usize) -> Self {
        self.fixed[dof] = true;
        self
    }

    pub fn with_mass(mut self, dof: usize, value: f64) -> Self {
        self.mass[dof] = value;
        self
    }

    /// Sets an initial displacement without solving for it.
    pub fn with_initial_displacement(mut self, dof: usize, value: f64) -> Self {
        self.displacement[dof] = value;
        self
    }

    pub fn with_initial_velocity(mut self, dof: usize, value: f64) -> Self {
        self.velocity[dof] = value;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planar_and_spatial_aliases_share_one_node_implementation() {
        let planar: Node2 = Node::new([1.0, 2.0]).fix(2);
        assert_eq!(planar.coords, [1.0, 2.0]);
        assert!(planar.fixed[2]);
        assert_eq!(planar.displacement.len(), 3);

        let spatial = Node3::new([1.0, 2.0, 3.0])
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Uz as usize, 4.0);
        assert_eq!(spatial.coords, [1.0, 2.0, 3.0]);
        assert!(spatial.fixed[SpatialDof::Rz as usize]);
        assert_eq!(spatial.mass[SpatialDof::Uz as usize], 4.0);
        assert_eq!(spatial.displacement.len(), 6);
    }
}
