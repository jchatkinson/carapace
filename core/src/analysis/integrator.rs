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
    /// Chooses this step's pseudo-time increment so that the controlled
    /// node/DOF's displacement changes by exactly `increment`, via the
    /// standard unit-load technique: solve the current tangent against the
    /// *sensitivity* of the total applied load to a unit pseudo-time
    /// increment (`Domain::assemble_reference_load_sensitivity` — every
    /// non-frozen `LoadPattern`'s contribution, weighted by its series'
    /// slope at the current pseudo-time; a `hold_pattern_constant`-frozen
    /// or `LoadSeries::Constant` pattern contributes nothing, since it
    /// doesn't respond to pseudo-time at all) once, then scale by whatever
    /// factor makes the controlled DOF's response match `increment`.
    DisplacementControl {
        node: NodeId,
        dof: usize,
        increment: f64,
    },
}

impl Integrator {
    /// Compute this step's new pseudo-time, given the domain's state as of
    /// the end of the *previous* (converged) step.
    pub(crate) fn predict(
        &self,
        domain: &Domain,
        solver: &SparseSolver,
        current_pseudo_time: f64,
    ) -> Result<f64, AnalysisError> {
        match self {
            Integrator::LoadControl { increment } => Ok(current_pseudo_time + increment),
            Integrator::DisplacementControl { node, dof, increment } => {
                let eq = domain.equation_of(*node, *dof).ok_or(AnalysisError::InvalidConstraint)?;
                let (k, _resistance) = domain.assemble_tangent_and_resistance();
                let sensitivity = domain.assemble_reference_load_sensitivity(current_pseudo_time);
                let unit_response = solver.solve(&k, &sensitivity)?;
                let delta_lambda = increment / unit_response[eq];
                Ok(current_pseudo_time + delta_lambda)
            }
        }
    }
}
