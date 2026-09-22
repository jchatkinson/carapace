/// Real error enum, not sentinel integers (§2.8). Carries the context an
/// OpenSees-style `-1`..`-5` return code would otherwise leave to
/// out-of-band documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalysisError {
    FailedToConverge { step: usize },
    SingularSystem,
    /// An `Integrator::DisplacementControl`'s controlled DOF is fixed (not
    /// part of the free-DOF system), so there's no equation to control.
    InvalidConstraint,
    /// `modal_analysis` was asked for more modes than the model has free
    /// DOFs (or zero modes) — there's no such Krylov subspace to build.
    InvalidModeCount { requested: usize, free_dofs: usize },
}
