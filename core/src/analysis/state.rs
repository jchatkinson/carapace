use crate::model::Domain;

use super::{Algorithm, AnalysisError, ConstraintHandler, ConvergenceTest, Integrator, SparseSolver};

/// A fully-wired analysis (only buildable via `AnalysisBuilder<Ready>::build`).
pub struct Analysis {
    pub(crate) domain: Domain,
    #[allow(dead_code)] // read once multi-point constraints exist
    pub(crate) constraint_handler: ConstraintHandler,
    pub(crate) integrator: Integrator,
    pub(crate) algorithm: Algorithm,
    #[allow(dead_code)] // read once Algorithm::NewtonRaphson (M4) needs it
    pub(crate) test: ConvergenceTest,
    pub(crate) solver: SparseSolver,
    pub(crate) step_count: usize,
    pub(crate) load_factor: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct StepResult {
    pub step: usize,
    pub load_factor: f64,
}

impl Analysis {
    pub fn domain(&self) -> &Domain {
        &self.domain
    }

    /// Advance one step: form the tangent and residual, solve, update state.
    /// See implementation-plan §5.4 for the sketch this follows.
    pub fn step(&mut self) -> Result<StepResult, AnalysisError> {
        self.step_count += 1;
        self.load_factor = self.integrator.next_load_factor(self.load_factor);

        let (k, residual) = self.domain.form_tangent_and_residual(self.load_factor);
        let du = self.solver.solve(&k, &residual)?;
        self.domain.apply_displacement_increment(&du);

        match self.algorithm {
            // A single tangent formation + solve is definitionally converged
            // for `Linear` — see algorithm.rs.
            Algorithm::Linear => {}
        }

        Ok(StepResult {
            step: self.step_count,
            load_factor: self.load_factor,
        })
    }
}
