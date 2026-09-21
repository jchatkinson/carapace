use faer::prelude::*;
use nalgebra::DVector;

use crate::model::SparseMatrix;

use super::AnalysisError;

/// The one linear solver Carapace ships (§3.6) — a concrete struct, not a
/// trait or enum, since there is exactly one implementation to commit to.
///
/// Resolves implementation-plan §6 open decision #1 for real: `faer`'s
/// sparse LU (with built-in COLAMD/AMD fill-reducing ordering), not a dense
/// `nalgebra` solve. The dense placeholder from M1 was fine while every
/// test model had single-digit DOF counts, but a dense O(n^3) LU doesn't
/// scale to real models — see the implementation plan's §6 decision #1 note
/// for the wasm32 build/solve verification this choice is based on.
///
/// Re-factors from scratch on every call — the sparsity *pattern* is fixed
/// for a given `Domain` (only values change between Newton iterations), so
/// caching the symbolic factorization across solves within one `Analysis`
/// is a real future optimization, just not one this milestone needs yet.
pub struct SparseSolver;

impl SparseSolver {
    pub fn new() -> Self {
        SparseSolver
    }

    pub fn solve(&self, k: &SparseMatrix, r: &DVector<f64>) -> Result<DVector<f64>, AnalysisError> {
        let n = r.len();
        let rhs = faer::col::Col::from_fn(n, |i| r[i]);
        let lu = k.sp_lu().map_err(|_| AnalysisError::SingularSystem)?;
        let x = lu.solve(&rhs);
        Ok(DVector::from_fn(n, |i, _| x[i]))
    }
}

impl Default for SparseSolver {
    fn default() -> Self {
        Self::new()
    }
}
