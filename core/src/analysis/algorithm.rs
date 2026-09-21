/// Iteration strategy for resolving each step's equilibrium. Closed enum
/// (§3.1). `NewtonRaphson` lands at M4.
#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    /// One tangent formation, one solve, no re-iteration and no convergence
    /// check — correct (and exact, to solver tolerance) whenever the
    /// tangent doesn't change with the solution, i.e. linear response.
    /// Matches OpenSees' `Linear` algorithm semantics.
    Linear,
}
