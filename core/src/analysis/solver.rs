use faer::prelude::*;
use faer::sparse::linalg::solvers::Lu;
use nalgebra::DVector;

use crate::model::SparseMatrix;

use super::AnalysisError;

/// The one linear solver Carapace ships (§2.6) — a concrete struct, not a
/// trait or enum, since there is exactly one implementation to commit to.
///
/// Resolves implementation-plan §5 open decision #1 for real: `faer`'s
/// sparse LU (with built-in COLAMD/AMD fill-reducing ordering), not a dense
/// `nalgebra` solve. The dense placeholder from M1 was fine while every
/// test model had single-digit DOF counts, but a dense O(n^3) LU doesn't
/// scale to real models — see the implementation plan's §5 decision #1 note
/// for the wasm32 build/solve verification this choice is based on.
pub struct SparseSolver;

impl SparseSolver {
    pub fn new() -> Self {
        SparseSolver
    }

    /// Factor `k`, without solving — the half of a solve that
    /// `TangentStrategy::ReuseAtStepStart`/`Initial` (`docs/algorithms.md`
    /// §3.1) let a caller skip repeating across iterations. `faer`'s `Lu`
    /// is a fully owned value (no lifetime tied back to `k`), so the result
    /// outlives the matrix it was formed from — see `SparseFactorization`.
    pub fn factor(&self, k: &SparseMatrix) -> Result<SparseFactorization, AnalysisError> {
        Ok(SparseFactorization(
            k.sp_lu().map_err(|_| AnalysisError::SingularSystem)?,
        ))
    }

    /// Factor-then-solve in one call — for `Algorithm::Linear` and anywhere
    /// else a single ad-hoc solve is enough. Equivalent to `self.factor(k)?
    /// .solve(r)`, kept as its own method so a one-off caller doesn't have
    /// to name `SparseFactorization` just to discard it immediately.
    pub fn solve(&self, k: &SparseMatrix, r: &DVector<f64>) -> Result<DVector<f64>, AnalysisError> {
        Ok(self.factor(k)?.solve(r))
    }
}

impl Default for SparseSolver {
    fn default() -> Self {
        Self::new()
    }
}

/// An LU factorization held independently of the matrix it was formed
/// from — cached across Newton iterations by `TangentStrategy::
/// ReuseAtStepStart`/`Initial` so re-solving against the same tangent
/// doesn't re-factor it (`docs/algorithms.md` §3.1's real prerequisite for
/// those strategies to actually save anything).
pub struct SparseFactorization(Lu<usize, f64>);

impl SparseFactorization {
    pub fn solve(&self, r: &DVector<f64>) -> DVector<f64> {
        let n = r.len();
        let rhs = faer::col::Col::from_fn(n, |i| r[i]);
        let x = self.0.solve(&rhs);
        DVector::from_fn(n, |i, _| x[i])
    }
}
