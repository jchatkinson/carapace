use std::sync::Mutex;

use faer::prelude::*;
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
use faer::sparse::{SymbolicSparseColMat, SymbolicSparseColMatRef};
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
pub struct SparseSolver {
    // Keep factor/solve callable through &self and preserve Send + Sync.
    // The lock is released before numeric factorization and solving.
    symbolic_cache: Mutex<Option<CachedSymbolic>>,
    #[cfg(test)]
    symbolic_analyses: std::sync::atomic::AtomicUsize,
}

struct CachedSymbolic {
    pattern: SymbolicSparseColMat<usize>,
    factorization: SymbolicLu<usize>,
}

impl CachedSymbolic {
    fn matches(&self, pattern: SymbolicSparseColMatRef<'_, usize>) -> bool {
        self.pattern.nrows() == pattern.nrows()
            && self.pattern.ncols() == pattern.ncols()
            && self.pattern.col_ptr() == pattern.col_ptr()
            && self.pattern.col_nnz() == pattern.col_nnz()
            && self.pattern.row_idx() == pattern.row_idx()
    }
}

impl SparseSolver {
    pub fn new() -> Self {
        Self {
            symbolic_cache: Mutex::new(None),
            #[cfg(test)]
            symbolic_analyses: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// Factor `k`, without solving — the half of a solve that
    /// `TangentStrategy::ReuseAtStepStart`/`Initial` (`docs/algorithms.md`
    /// §3.1) let a caller skip repeating across iterations. `faer`'s `Lu`
    /// is a fully owned value (no lifetime tied back to `k`), so the result
    /// outlives the matrix it was formed from — see `SparseFactorization`.
    /// Symbolic analysis is reused while the full CSC pattern matches;
    /// values (including zeros) never participate in cache invalidation.
    /// Numeric factorization is always recomputed here, independently of
    /// callers' numeric tangent-reuse caches.
    pub fn factor(&self, k: &SparseMatrix) -> Result<SparseFactorization, AnalysisError> {
        let symbolic = {
            let mut cache = self
                .symbolic_cache
                .lock()
                .expect("symbolic LU cache poisoned");
            let pattern = k.symbolic();
            if !cache.as_ref().is_some_and(|cached| cached.matches(pattern)) {
                let factorization =
                    SymbolicLu::try_new(pattern).map_err(|_| AnalysisError::SingularSystem)?;
                let pattern = pattern
                    .to_owned()
                    .map_err(|_| AnalysisError::SingularSystem)?;
                *cache = Some(CachedSymbolic {
                    pattern,
                    factorization,
                });
                #[cfg(test)]
                self.symbolic_analyses
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            cache
                .as_ref()
                .expect("symbolic cache populated above")
                .factorization
                .clone()
        };
        Ok(SparseFactorization(
            Lu::try_new_with_symbolic(symbolic, k.as_ref())
                .map_err(|_| AnalysisError::SingularSystem)?,
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

#[cfg(test)]
mod tests {
    use super::*;
    use faer::sparse::Triplet;
    use std::sync::atomic::Ordering;

    fn matrix(n: usize, entries: &[(usize, usize, f64)]) -> SparseMatrix {
        let triplets: Vec<_> = entries
            .iter()
            .map(|&(i, j, v)| Triplet::new(i, j, v))
            .collect();
        SparseMatrix::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    fn assert_solution(factor: &SparseFactorization, rhs: &[f64], expected: &[f64]) {
        let solution = factor.solve(&DVector::from_column_slice(rhs));
        assert!((solution - DVector::from_column_slice(expected)).norm() < 1e-12);
    }

    #[test]
    fn changing_values_and_pivots_reuses_symbolic_but_not_numeric_factors() {
        let solver = SparseSolver::new();
        let first = matrix(2, &[(0, 0, 4.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 5.0)]);
        let factor = solver.factor(&first).unwrap();
        assert_solution(&factor, &[8.0, 13.0], &[1.0, 2.0]);

        // A structural entry becomes zero and numeric LU must change pivots.
        let second = matrix(2, &[(0, 0, 0.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)]);
        let updated = solver.factor(&second).unwrap();
        assert_solution(&updated, &[4.0, 11.0], &[1.0, 2.0]);
        assert_eq!(solver.symbolic_analyses.load(Ordering::Relaxed), 1);
        // Previously returned numeric factors remain independently owned.
        assert_solution(&factor, &[8.0, 13.0], &[1.0, 2.0]);
    }

    #[test]
    fn changed_pattern_with_same_size_and_nnz_and_changed_dimensions_rebuild_cache() {
        let solver = SparseSolver::new();
        let upper = matrix(2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 1, 3.0)]);
        assert_solution(&solver.factor(&upper).unwrap(), &[6.0, 6.0], &[1.0, 2.0]);
        let lower = matrix(2, &[(0, 0, 4.0), (1, 0, 1.0), (1, 1, 3.0)]);
        assert_solution(&solver.factor(&lower).unwrap(), &[4.0, 7.0], &[1.0, 2.0]);
        assert_eq!(solver.symbolic_analyses.load(Ordering::Relaxed), 2);

        let larger = matrix(3, &[(0, 0, 2.0), (1, 1, 3.0), (2, 2, 4.0)]);
        assert_solution(
            &solver.factor(&larger).unwrap(),
            &[2.0, 6.0, 12.0],
            &[1.0, 2.0, 3.0],
        );
        assert_eq!(solver.symbolic_analyses.load(Ordering::Relaxed), 3);
        assert_solution(&solver.factor(&lower).unwrap(), &[4.0, 7.0], &[1.0, 2.0]);
        assert_eq!(solver.symbolic_analyses.load(Ordering::Relaxed), 4);
    }

    #[test]
    fn singular_pattern_does_not_poison_cache_or_invalidate_owned_factors() {
        let solver = SparseSolver::new();
        let nonsingular = matrix(2, &[(0, 0, 2.0), (1, 1, 3.0)]);
        let factor = solver.factor(&nonsingular).unwrap();
        let singular = matrix(2, &[(1, 1, 3.0)]);
        assert!(matches!(
            solver.factor(&singular),
            Err(AnalysisError::SingularSystem)
        ));
        assert_solution(&factor, &[2.0, 6.0], &[1.0, 2.0]);
        assert_solution(
            &solver.factor(&nonsingular).unwrap(),
            &[2.0, 6.0],
            &[1.0, 2.0],
        );
    }
}
