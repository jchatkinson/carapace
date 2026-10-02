/// Real error enum, not sentinel integers (§2.8). Carries the context an
/// OpenSees-style `-1`..`-5` return code would otherwise leave to
/// out-of-band documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalysisError {
    FailedToConverge {
        step: usize,
    },
    SingularSystem,
    /// An `Integrator::DisplacementControl`'s controlled DOF is fixed (not
    /// part of the free-DOF system), so there's no equation to control.
    InvalidConstraint,
    /// `modal_analysis` was asked for more modes than the model has free
    /// DOFs (or zero modes) — there's no such Krylov subspace to build.
    InvalidModeCount {
        requested: usize,
        free_dofs: usize,
    },
    /// A configuration value is out of range, nonfinite, or not supported
    /// in combination with the rest of the analysis. `field` names the
    /// option in its wire-format spelling (e.g. `"integrator.minRadius"`).
    InvalidOption {
        field: &'static str,
    },
    /// Arc-length continuation found an unfrozen `LoadSeries::Path` pattern
    /// (its slope isn't valid across a load reversal or breakpoint).
    UnsupportedLoadSeries,
    /// Arc-length continuation's combined reference-load sensitivity is
    /// zero on every free equation: there is no load to continue.
    ZeroLoadSensitivity,
    /// Arc-length continuation must start from equilibrium. `measure` is the
    /// force-test measure at entry (`||r||` for `NormUnbalance`, the largest
    /// scaled component for `Combined`).
    InitialStateNotInEquilibrium {
        measure: f64,
    },
    /// The first arc-length tangent is singular and no `seed` direction was
    /// configured to orient it.
    MissingSeedDirection,
    /// Every arc-length attempt for this step failed and the radius can't
    /// be reduced further (or the retry budget is spent). `last_failure`
    /// says why the final attempt failed.
    CutbacksExhausted {
        step: usize,
        attempts: usize,
        radius: f64,
        last_failure: ArcFailure,
    },
    /// An arc-length stop criterion was already met; the phase is over and
    /// further steps are refused without changing any state.
    ContinuationComplete,
}

/// Why one arc-length attempt was rejected. Attempts are retried with a
/// smaller radius; this only surfaces through
/// `AnalysisError::CutbacksExhausted`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcFailure {
    /// The iteration budget ran out before every criterion passed.
    NotConverged,
    /// Several iterations made no meaningful progress on the merit function.
    Stagnated,
    /// Coupled backtracking found no step length that decreases the merit
    /// function.
    NoDescent,
    /// A trial state, residual or correction was not finite.
    NonFinite,
    /// The bordered system was singular.
    SingularSystem,
    /// The bordered solve's scaled backward error stayed above target after
    /// refinement.
    InaccurateSolve,
    /// The converged increment pointed backward along the path.
    BranchOrientation,
    /// The iterate returned to the committed point, where the constraint
    /// row vanishes.
    DegenerateConstraint,
}
