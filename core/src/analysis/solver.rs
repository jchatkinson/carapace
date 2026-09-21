use nalgebra::{DMatrix, DVector};

use super::AnalysisError;

/// The one linear solver Carapace ships (§3.6) — a concrete struct, not a
/// trait or enum, since there is exactly one implementation to commit to.
///
/// Resolves implementation-plan §6 open decision #1: `nalgebra` over
/// `faer`, using its dense LU factorization for now. At M1's DOF counts
/// (single-digit) dense vs. sparse is not the bottleneck either way — the
/// real sparse-vs-dense and nalgebra-vs-faer comparison is deferred to a
/// milestone where problem sizes actually stress it (M3+, once multi-
/// element frame models exist). Revisit then; until this struct's
/// implementation changes, it's the one place that needs to.
pub struct SparseSolver;

impl SparseSolver {
    pub fn new() -> Self {
        SparseSolver
    }

    pub fn solve(&self, k: &DMatrix<f64>, r: &DVector<f64>) -> Result<DVector<f64>, AnalysisError> {
        k.clone().lu().solve(r).ok_or(AnalysisError::SingularSystem)
    }
}

impl Default for SparseSolver {
    fn default() -> Self {
        Self::new()
    }
}
