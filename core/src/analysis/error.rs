/// Real error enum, not sentinel integers (§3.8). Carries the context an
/// OpenSees-style `-1`..`-5` return code would otherwise leave to
/// out-of-band documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalysisError {
    FailedToConverge { step: usize },
    SingularSystem,
    /// An `Integrator::DisplacementControl`'s controlled DOF is fixed (not
    /// part of the free-DOF system), so there's no equation to control.
    InvalidConstraint,
}
