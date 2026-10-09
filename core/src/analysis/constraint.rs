/// Constraint-handling strategy. Closed enum — small, homogeneous,
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
    /// Single-point *and* multi-point constraints: general linear ties
    /// `u_slave = sum(coeff * u_master)` with masters on any nodes
    /// (`Domain::add_constraint`, and the helpers built on it:
    /// `equal_dof`, `rigid_diaphragm`, `rigid_link`). They are resolved at
    /// `Domain::number_dofs` time by substitution, with no extra unknowns
    /// (no Lagrange multipliers) and no penalty: a slave that reduces to
    /// exactly one master with coefficient 1 is literally the same unknown
    /// (it shares the master's equation number); any other slave is a
    /// stored linear combination of unknowns that assembly and state
    /// scatter apply (`Domain::value_at` reads one back). Chains (a slave
    /// that is another constraint's master) resolve to any depth. The
    /// constraints are homogeneous, so a fixed master must stay at zero.
    Transformation,
}
