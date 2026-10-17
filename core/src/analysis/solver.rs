use std::sync::Mutex;

use faer::prelude::*;
use faer::sparse::linalg::solvers::{Llt, Lu, SymbolicLlt, SymbolicLu};
use faer::sparse::{SymbolicSparseColMat, SymbolicSparseColMatRef};
use faer::Side;
use nalgebra::DVector;

use crate::model::SparseMatrix;

use super::AnalysisError;

/// The one linear solver Carapace ships — a concrete struct, not a
/// trait or enum, since there is exactly one implementation to commit to.
///
/// `faer`'s sparse direct solvers (with built-in COLAMD/AMD fill-reducing
/// ordering), not a dense `nalgebra` solve. A dense O(n^3) LU doesn't scale
/// to real models. They build and solve correctly under wasm32.
///
/// A symmetric matrix is first tried with sparse Cholesky (`LLT`), roughly
/// an order of magnitude faster than general LU on the stiffness matrices
/// of stable structures. The matrix must be symmetric (Cholesky reads only
/// its lower triangle) and positive definite; if either fails, e.g. past a
/// limit point or for an unsymmetric tangent, this falls back to LU, which
/// handles any nonsingular matrix.
pub struct SparseSolver {
    // Keep factor/solve callable through &self and preserve Send + Sync.
    // The lock is released before numeric factorization and solving.
    symbolic_cache: Mutex<Option<CachedSymbolic>>,
    cholesky_cache: Mutex<Option<CachedCholesky>>,
    #[cfg(test)]
    symbolic_analyses: std::sync::atomic::AtomicUsize,
}

struct CachedSymbolic {
    pattern: SymbolicSparseColMat<usize>,
    factorization: SymbolicLu<usize>,
}

struct CachedCholesky {
    pattern: SymbolicSparseColMat<usize>,
    factorization: SymbolicLlt<usize>,
}

/// Whether two CSC patterns are identical (same shape, column pointers, counts and row indices).
fn same_pattern(
    a: SymbolicSparseColMatRef<'_, usize>,
    b: SymbolicSparseColMatRef<'_, usize>,
) -> bool {
    a.nrows() == b.nrows()
        && a.ncols() == b.ncols()
        && a.col_ptr() == b.col_ptr()
        && a.col_nnz() == b.col_nnz()
        && a.row_idx() == b.row_idx()
}

/// Whether `k` is square and symmetric to rounding: every off-diagonal entry matches its transpose to a
/// relative `1e-10` of the largest entry magnitude, a missing counterpart counting as zero. Assembly
/// sums the two triangles in different orders, so exact equality would reject sound stiffness matrices.
///
/// Linear in the number of entries: columns are visited in order, so the lookup of row `j` in column `i`
/// only ever moves forward (`cursor[i]`).
fn is_symmetric(k: &SparseMatrix) -> bool {
    let (symbolic, values) = (k.symbolic(), k.val());
    let n = symbolic.ncols();
    if symbolic.nrows() != n {
        return false;
    }
    let rows = symbolic.row_idx();
    let tolerance = 1e-10 * values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut cursor: Vec<usize> = (0..n).map(|col| symbolic.col_range(col).start).collect();
    for col in 0..n {
        for at in symbolic.col_range(col) {
            let row = rows[at];
            if row == col {
                continue;
            }
            let end = symbolic.col_range(row).end;
            let next = &mut cursor[row];
            while *next < end && rows[*next] < col {
                *next += 1;
            }
            let mirrored = if *next < end && rows[*next] == col {
                values[*next]
            } else {
                0.0
            };
            if (values[at] - mirrored).abs() > tolerance {
                return false;
            }
        }
    }
    true
}

impl CachedSymbolic {
    fn matches(&self, pattern: SymbolicSparseColMatRef<'_, usize>) -> bool {
        same_pattern(self.pattern.as_ref(), pattern)
    }
}

impl CachedCholesky {
    fn matches(&self, pattern: SymbolicSparseColMatRef<'_, usize>) -> bool {
        same_pattern(self.pattern.as_ref(), pattern)
    }
}

impl SparseSolver {
    pub fn new() -> Self {
        Self {
            symbolic_cache: Mutex::new(None),
            cholesky_cache: Mutex::new(None),
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
        if is_symmetric(k) {
            if let Some(llt) = self.try_cholesky(k) {
                return Ok(SparseFactorization(Factors::Cholesky(llt)));
            }
        }
        self.factor_lu(k)
    }

    /// Cholesky factors of `k`, or `None` when it is not positive definite (or its ordering fails), in which
    /// case the caller falls back to LU. The symbolic analysis is cached like the LU one.
    fn try_cholesky(&self, k: &SparseMatrix) -> Option<Llt<usize, f64>> {
        let symbolic = {
            let mut cache = self
                .cholesky_cache
                .lock()
                .expect("symbolic LLT cache poisoned");
            let pattern = k.symbolic();
            if !cache.as_ref().is_some_and(|cached| cached.matches(pattern)) {
                let factorization = SymbolicLlt::try_new(pattern, Side::Lower).ok()?;
                let pattern = pattern.to_owned().ok()?;
                *cache = Some(CachedCholesky {
                    pattern,
                    factorization,
                });
            }
            cache
                .as_ref()
                .expect("symbolic cache populated above")
                .factorization
                .clone()
        };
        Llt::try_new_with_symbolic(symbolic, k.as_ref(), Side::Lower).ok()
    }

    fn factor_lu(&self, k: &SparseMatrix) -> Result<SparseFactorization, AnalysisError> {
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
        Ok(SparseFactorization(Factors::Lu(Box::new(
            Lu::try_new_with_symbolic(symbolic, k.as_ref())
                .map_err(|_| AnalysisError::SingularSystem)?,
        ))))
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
pub struct SparseFactorization(Factors);

enum Factors {
    Cholesky(Llt<usize, f64>),
    Lu(Box<Lu<usize, f64>>),
}

impl SparseFactorization {
    pub fn solve(&self, r: &DVector<f64>) -> DVector<f64> {
        let n = r.len();
        let rhs = faer::col::Col::from_fn(n, |i| r[i]);
        let x = match &self.0 {
            Factors::Cholesky(llt) => llt.solve(&rhs),
            Factors::Lu(lu) => lu.solve(&rhs),
        };
        DVector::from_fn(n, |i, _| x[i])
    }

    #[cfg(test)]
    fn is_cholesky(&self) -> bool {
        matches!(self.0, Factors::Cholesky(_))
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

        // Unsymmetric on purpose: a symmetric positive-definite matrix would take the Cholesky path and
        // not exercise this LU cache.
        let larger = matrix(3, &[(0, 0, 2.0), (0, 1, 1.0), (1, 1, 3.0), (2, 2, 4.0)]);
        assert_solution(
            &solver.factor(&larger).unwrap(),
            &[4.0, 6.0, 12.0],
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

    /// A symmetric positive-definite matrix is factored by Cholesky and solves correctly, including for a
    /// second matrix with the same pattern but new values (the cached symbolic analysis is reused).
    #[test]
    fn symmetric_positive_definite_uses_cholesky() {
        let solver = SparseSolver::new();
        let k = matrix(
            3,
            &[
                (0, 0, 4.0),
                (0, 1, 1.0),
                (1, 0, 1.0),
                (1, 1, 3.0),
                (1, 2, 1.0),
                (2, 1, 1.0),
                (2, 2, 2.0),
            ],
        );
        let factor = solver.factor(&k).unwrap();
        assert!(factor.is_cholesky());
        // K * [1, 2, 3] = [6, 10, 8]
        assert_solution(&factor, &[6.0, 10.0, 8.0], &[1.0, 2.0, 3.0]);

        let stiffer = matrix(
            3,
            &[
                (0, 0, 8.0),
                (0, 1, 1.0),
                (1, 0, 1.0),
                (1, 1, 6.0),
                (1, 2, 1.0),
                (2, 1, 1.0),
                (2, 2, 4.0),
            ],
        );
        let again = solver.factor(&stiffer).unwrap();
        assert!(again.is_cholesky());
        assert_solution(&again, &[10.0, 16.0, 14.0], &[1.0, 2.0, 3.0]);
    }

    /// A symmetric matrix that is not positive definite (a softened tangent) makes Cholesky fail, so the
    /// solve falls back to LU and still returns the right answer.
    #[test]
    fn symmetric_indefinite_falls_back_to_lu() {
        let k = matrix(2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0)]);
        let factor = SparseSolver::new().factor(&k).unwrap();
        assert!(!factor.is_cholesky());
        // K * [1, 2] = [5, 4]
        assert_solution(&factor, &[5.0, 4.0], &[1.0, 2.0]);
    }

    /// An unsymmetric matrix, or one whose pattern lacks the transpose entry, is never given to Cholesky
    /// (which reads only the lower triangle and would solve a different system).
    #[test]
    fn unsymmetric_matrices_use_lu() {
        let solver = SparseSolver::new();
        let values = matrix(2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 2.0), (1, 1, 3.0)]);
        let factor = solver.factor(&values).unwrap();
        assert!(!factor.is_cholesky());
        assert_solution(&factor, &[6.0, 8.0], &[1.0, 2.0]);

        let pattern = matrix(2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 1, 3.0)]);
        let factor = solver.factor(&pattern).unwrap();
        assert!(!factor.is_cholesky());
        assert_solution(&factor, &[6.0, 6.0], &[1.0, 2.0]);
    }

    /// Symmetry is judged to rounding: a 1e-14 mismatch (assembly summation order) still uses Cholesky, a
    /// 1e-3 mismatch does not.
    #[test]
    fn symmetry_tolerance_accepts_rounding_but_not_real_asymmetry() {
        let nearly = matrix(
            2,
            &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0 + 1e-14), (1, 1, 3.0)],
        );
        assert!(SparseSolver::new().factor(&nearly).unwrap().is_cholesky());
        let skewed = matrix(2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.001), (1, 1, 3.0)]);
        assert!(!SparseSolver::new().factor(&skewed).unwrap().is_cholesky());
    }

    /// The linear-time symmetry check agrees with the definition on a scattered pattern: mirrored entries in
    /// several columns are accepted, and a single entry without its mirror (or with a different one) is not.
    #[test]
    fn symmetry_check_handles_scattered_patterns() {
        let mut entries = vec![
            (0, 0, 10.0),
            (1, 1, 10.0),
            (2, 2, 10.0),
            (3, 3, 10.0),
            (4, 4, 10.0),
        ];
        for &(i, j, v) in &[
            (0, 3, 1.0),
            (1, 4, 2.0),
            (0, 4, 0.5),
            (2, 3, 1.5),
            (3, 4, 0.25),
            (1, 2, 0.75),
        ] {
            entries.push((i, j, v));
            entries.push((j, i, v));
        }
        assert!(is_symmetric(&matrix(5, &entries)));
        assert!(SparseSolver::new()
            .factor(&matrix(5, &entries))
            .unwrap()
            .is_cholesky());

        let mut lonely = entries.clone();
        lonely.push((2, 4, 1.0));
        assert!(!is_symmetric(&matrix(5, &lonely)));
        let mut lonely_lower = entries.clone();
        lonely_lower.push((4, 2, 1.0));
        assert!(!is_symmetric(&matrix(5, &lonely_lower)));
        let mut different = entries.clone();
        different.push((2, 4, 1.0));
        different.push((4, 2, 1.5));
        assert!(!is_symmetric(&matrix(5, &different)));
    }
}
