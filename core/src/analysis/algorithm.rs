use super::{LineSearch, TangentStrategy};

/// Iteration strategy for resolving each step's equilibrium. Closed enum
/// (§2.1). `docs/algorithms.md` designs `Newton`/`KrylovNewton` as three
/// composable axes (tangent strategy, line search, Krylov acceleration)
/// rather than a combinatorial pile of variants — see that doc for the
/// full design and Xara reference map.
#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    /// One tangent formation, one solve, no re-iteration and no convergence
    /// check — correct (and exact, to solver tolerance) whenever the
    /// tangent doesn't change with the solution, i.e. linear response.
    /// Matches OpenSees' `Linear` algorithm semantics.
    Linear,
    /// Newton-family iteration: re-forms the residual every iteration and
    /// solves against whichever tangent `tangent` says to form/factor
    /// (`TangentStrategy`), until `ConvergenceTest` passes or its
    /// `max_iter` is exceeded (`AnalysisError::FailedToConverge`).
    /// `tangent: TangentStrategy::Current, line_search: None` is today's
    /// old bare `NewtonRaphson` variant, exactly.
    Newton {
        /// Which tangent gets formed/factored at each iteration.
        tangent: TangentStrategy,
        /// `None`: accept each raw Newton correction outright. `Some(_)`:
        /// rescale a correction that doesn't already satisfy the
        /// convergence test at `eta = 1` (`docs/algorithms.md` §4).
        line_search: Option<LineSearch>,
    },
    /// Modified Newton (solving against a *stale* tangent — `tangent`,
    /// almost always `Initial`/`ReuseAtStepStart` in practice) with a
    /// Krylov-subspace correction (`docs/algorithms.md` §5) applied to
    /// each raw solve, extrapolating from the history of how the last
    /// `max_dimension` corrections actually moved the residual to recover
    /// most of full Newton's convergence rate without paying for a fresh
    /// factorization every iteration.
    KrylovNewton { tangent: TangentStrategy, max_dimension: usize },
}
