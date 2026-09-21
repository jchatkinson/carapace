use slotmap::new_key_type;

use super::NDF;

new_key_type! {
    /// Generational index into `Domain`'s node store. Stays valid across
    /// removals of *other* nodes; using a stale key against a different
    /// domain is a logic error the type can't catch, but a use-after-free
    /// within one domain's lifetime is (§3.2 of the implementation plan).
    pub struct NodeId;
}

/// A node: fixed 2D coordinates, `NDF` degrees of freedom, and the mutable
/// state (displacement, velocity, acceleration, applied load, boundary
/// conditions) an analysis reads and writes each step. `velocity`/
/// `acceleration` are only meaningful for/written by `TransientAnalysis`
/// (M6) — static `Analysis` never touches them.
#[derive(Debug, Clone)]
pub struct Node {
    pub coords: [f64; 2],
    pub fixed: [bool; NDF],
    pub load: [f64; NDF],
    pub displacement: [f64; NDF],
    pub velocity: [f64; NDF],
    pub acceleration: [f64; NDF],
    /// Lumped nodal mass per DOF — a user-assigned point mass, additive
    /// with any element-consistent lumped mass (`Element::form_mass`, M6).
    /// Needed by modal analysis (M5) and time-history analysis (M6).
    pub mass: [f64; NDF],
    /// Equation number for each free DOF, or `None` if fixed. Assigned by
    /// `Domain::number_dofs` during `AnalysisBuilder::build`.
    pub(crate) equation: [Option<usize>; NDF],
}

impl Node {
    pub fn new(coords: [f64; 2]) -> Self {
        Node {
            coords,
            fixed: [false; NDF],
            load: [0.0; NDF],
            displacement: [0.0; NDF],
            velocity: [0.0; NDF],
            acceleration: [0.0; NDF],
            mass: [0.0; NDF],
            equation: [None; NDF],
        }
    }

    pub fn fix(mut self, dof: usize) -> Self {
        self.fixed[dof] = true;
        self
    }

    pub fn with_load(mut self, dof: usize, value: f64) -> Self {
        self.load[dof] = value;
        self
    }

    pub fn with_mass(mut self, dof: usize, value: f64) -> Self {
        self.mass[dof] = value;
        self
    }

    /// Sets an initial condition for `TransientAnalysis` (M6) — bypasses
    /// the usual equation-driven displacement path (static `Analysis` only
    /// ever *adds* increments via `Domain::apply_displacement_increment`),
    /// since an initial condition is a given, not something solved for.
    pub fn with_initial_displacement(mut self, dof: usize, value: f64) -> Self {
        self.displacement[dof] = value;
        self
    }

    pub fn with_initial_velocity(mut self, dof: usize, value: f64) -> Self {
        self.velocity[dof] = value;
        self
    }
}
