/// Which tangent stiffness an `Algorithm::Newton`/`KrylovNewton` iteration
/// solves against. Closed enum — mirrors Xara's `TangentFlagType`
/// (`CURRENT_TANGENT` / `INITIAL_TANGENT` / `PREDICTOR_TANGENT`), which
/// every Newton-family algorithm class there takes as a constructor
/// argument rather than hard-coding (`docs/algorithms.md` §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangentStrategy {
    /// Re-form and re-factor the tangent at the *current* trial state,
    /// every iteration. Full Newton-Raphson — today's only behavior before
    /// this type existed.
    Current,
    /// Form the tangent once, at the state the step *started* from, and
    /// reuse it (still factored once) for every iteration of that step.
    /// Xara's `ModifiedNewton` with `CURRENT_TANGENT` passed at
    /// construction (misleadingly named from the caller's perspective —
    /// it means "current as of step start," not "current every
    /// iteration").
    ReuseAtStepStart,
    /// Form the tangent once, from the model's very first (undeformed)
    /// state, ever — reused for every iteration of every step. Xara's
    /// `ModifiedNewton` with `INITIAL_TANGENT`. Cheapest, worst
    /// convergence rate (linear, not quadratic); pairs naturally with
    /// `Algorithm::KrylovNewton`'s acceleration to claw the convergence
    /// rate back without paying for re-factorization.
    Initial,
}
