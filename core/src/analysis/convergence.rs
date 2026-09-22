use nalgebra::DVector;

/// Convergence criterion for `Algorithm::NewtonRaphson`. Closed enum
/// (§2.1). Unused by `Algorithm::Linear` (which never iterates).
#[derive(Debug, Clone, Copy)]
pub enum ConvergenceTest {
    /// Converged once the unbalanced-force residual's norm drops below `tol`.
    NormUnbalance { tol: f64, max_iter: usize },
    /// Converged once the displacement correction's norm drops below `tol`.
    NormDispIncr { tol: f64, max_iter: usize },
    /// Converged once the (absolute) incremental work `|du . residual|`
    /// drops below `tol` — scale-invariant across mixed translational/
    /// rotational DOFs in a way `NormDispIncr` isn't.
    EnergyIncr { tol: f64, max_iter: usize },
}

impl ConvergenceTest {
    pub(crate) fn max_iter(&self) -> usize {
        match self {
            ConvergenceTest::NormUnbalance { max_iter, .. }
            | ConvergenceTest::NormDispIncr { max_iter, .. }
            | ConvergenceTest::EnergyIncr { max_iter, .. } => *max_iter,
        }
    }

    /// Whether the iteration that just computed `residual` (the unbalance
    /// solved against) and `du` (the resulting displacement correction,
    /// already applied) has converged.
    pub(crate) fn check(&self, residual: &DVector<f64>, du: &DVector<f64>) -> bool {
        match self {
            ConvergenceTest::NormUnbalance { tol, .. } => residual.norm() < *tol,
            ConvergenceTest::NormDispIncr { tol, .. } => du.norm() < *tol,
            ConvergenceTest::EnergyIncr { tol, .. } => du.dot(residual).abs() < *tol,
        }
    }
}
