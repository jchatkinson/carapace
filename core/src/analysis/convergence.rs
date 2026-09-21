/// Convergence criterion for iterative algorithms. Closed enum (§3.1).
/// Unused by `Algorithm::Linear` (which never iterates), but required by
/// `AnalysisBuilder` so the typestate wiring matches real usage once
/// `NewtonRaphson` lands at M4 and actually needs it. `NormDispIncr` /
/// `EnergyIncr` join `NormUnbalance` at M4.
#[derive(Debug, Clone, Copy)]
pub enum ConvergenceTest {
    NormUnbalance { tol: f64, max_iter: usize },
}
