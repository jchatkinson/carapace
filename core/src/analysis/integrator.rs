use crate::model::{Domain, NodeId};

use super::{AnalysisError, SparseSolver};

/// Static/pseudo-static integration strategy. Closed enum (§2.1). Both
/// variants only ever change how this *step's* load factor is chosen (the
/// "predictor") — the iteration that follows (§ `Algorithm`) proceeds
/// identically regardless of which one picked it.
#[derive(Debug, Clone, Copy)]
pub enum Integrator {
    LoadControl {
        increment: f64,
    },
    /// Chooses this step's load factor so that the controlled node/DOF's
    /// displacement changes by exactly `increment`, via the standard
    /// unit-load technique: solve the current tangent against the
    /// (unscaled) reference load pattern once, then scale by whatever
    /// factor makes the controlled DOF's response match `increment`.
    DisplacementControl {
        node: NodeId,
        dof: usize,
        increment: f64,
    },
}

impl Integrator {
    /// Compute this step's new load factor, given the domain's state as of
    /// the end of the *previous* (converged) step.
    pub(crate) fn predict(
        &self,
        domain: &Domain,
        solver: &SparseSolver,
        current_load_factor: f64,
    ) -> Result<f64, AnalysisError> {
        match self {
            Integrator::LoadControl { increment } => Ok(current_load_factor + increment),
            Integrator::DisplacementControl { node, dof, increment } => {
                let eq = domain.equation_of(*node, *dof).ok_or(AnalysisError::InvalidConstraint)?;
                let (k, _resistance) = domain.assemble_tangent_and_resistance();
                let reference_load = domain.assemble_reference_load();
                let unit_response = solver.solve(&k, &reference_load)?;
                let delta_lambda = increment / unit_response[eq];
                Ok(current_load_factor + delta_lambda)
            }
        }
    }
}
