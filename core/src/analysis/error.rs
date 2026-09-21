/// Real error enum, not sentinel integers (§3.8). Carries the context an
/// OpenSees-style `-1`..`-5` return code would otherwise leave to
/// out-of-band documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalysisError {
    FailedToConverge { step: usize },
    SingularSystem,
}
