/// Constraint-handling strategy. Closed enum (§2.1) — small, homogeneous,
/// user-selected per analysis.
#[derive(Debug, Clone, Copy)]
pub enum ConstraintHandler {
    /// Single-point constraints only (fixed DOFs), applied directly via
    /// `Domain`'s DOF numbering (fixed DOFs get no equation number, so they
    /// never enter the free-DOF system at all). No Lagrange multipliers, no
    /// penalty terms. Using this with a `Domain` that has multi-point
    /// constraints (`Domain::equal_dof`/`rigid_diaphragm`) is a model
    /// construction error — see `AnalysisBuilder<Ready>::build`; use
    /// `Transformation` instead.
    Plain,
    /// Single-point *and* multi-point constraints (`Domain::equal_dof`,
    /// `Domain::rigid_diaphragm`). Multi-point ties are resolved by DOF
    /// aliasing at `Domain::number_dofs` time: a constrained DOF gets no
    /// equation number of its own and instead reuses its retained node's
    /// equation number for that DOF, so tied DOFs are literally the same
    /// unknown in the assembled system — no coefficient/transformation
    /// matrix, no extra unknowns (Lagrange), no penalty stiffness. This
    /// only works because Carapace's `equal_dof`/`rigid_diaphragm` are
    /// identity ties (coefficient 1, same DOF index on both sides); a
    /// general affine multi-point constraint (arbitrary coefficients,
    /// cross-DOF terms) would need a real transformation matrix instead —
    /// out of scope until a concrete need arises (see
    /// `docs/implementation-plan.md`).
    Transformation,
}
