/// Constraint-handling strategy. Closed enum (§2.1) — small, homogeneous,
/// user-selected per analysis.
#[derive(Debug, Clone, Copy)]
pub enum ConstraintHandler {
    /// Single-point constraints only (fixed DOFs), applied directly via
    /// `Domain`'s DOF numbering (fixed DOFs get no equation number, so they
    /// never enter the free-DOF system at all). No Lagrange multipliers, no
    /// penalty terms, no multi-point constraints — those would need a
    /// different variant if/when needed.
    Plain,
}
