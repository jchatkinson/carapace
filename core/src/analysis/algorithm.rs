/// Iteration strategy for resolving each step's equilibrium. Closed enum
/// (§2.1).
#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    /// One tangent formation, one solve, no re-iteration and no convergence
    /// check — correct (and exact, to solver tolerance) whenever the
    /// tangent doesn't change with the solution, i.e. linear response.
    /// Matches OpenSees' `Linear` algorithm semantics.
    Linear,
    /// Full Newton-Raphson: re-forms the tangent and residual every
    /// iteration, at a load factor held fixed for the whole step (only the
    /// `Integrator`'s predictor changes it), until `ConvergenceTest`
    /// passes or its `max_iter` is exceeded (`AnalysisError::FailedToConverge`).
    NewtonRaphson,
}
