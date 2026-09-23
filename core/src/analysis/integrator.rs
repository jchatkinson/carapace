use nalgebra::DVector;

use crate::model::{Domain, NodeId, SparseMatrix};

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

    /// Corrector for the *second and later* Newton iterations of a step
    /// (the caller must skip this on the first iteration — see below).
    /// `predict`'s pseudo-time only accounts for the tangent at the *start*
    /// of the step — exact for `LoadControl` (whose pseudo-time never
    /// changes mid-step, so `du_bar` unmodified is already the right
    /// increment every iteration).
    ///
    /// For `DisplacementControl`, the first Newton iteration uses that same
    /// starting tangent (nothing has moved yet), so its residual-driven
    /// `du_bar` already delivers `predict`'s target displacement at the
    /// controlled DOF exactly — the caller applies it unmodified, with no
    /// call to `correct` (equivalent to OpenSees' `newStep`, which predicts
    /// *and* applies that first increment together). But the tangent used
    /// by `predict`/iteration 1 can go stale for the *rest* of the step
    /// (e.g. a step landing exactly on a material's yield breakpoint sees
    /// the elastic tangent at `predict` time, then the true, much softer,
    /// post-yield tangent from iteration 2 on) — left uncorrected, ordinary
    /// fixed-load-factor Newton iteration would then keep driving the
    /// controlled DOF *past* the target already reached, chasing force
    /// equilibrium at the wrong load factor instead. This is the standard
    /// Yang & Shieh consistent Displacement Control corrector: solve the
    /// current tangent against the reference load sensitivity again
    /// (`unit_response`, this iteration's version of `predict`'s unit-load
    /// probe), then pick `delta_lambda` so that adding `delta_lambda *
    /// unit_response` to `du_bar` cancels `du_bar`'s contribution at the
    /// controlled DOF — i.e. this iteration leaves the controlled DOF
    /// exactly where it already was after iteration 1, no matter how far
    /// the tangent has drifted since, while every other DOF still gets
    /// `du_bar`'s residual-driven correction.
    pub(crate) fn correct(
        &self,
        domain: &Domain,
        solver: &SparseSolver,
        k: &SparseMatrix,
        du_bar: DVector<f64>,
        pseudo_time: f64,
    ) -> Result<(f64, DVector<f64>), AnalysisError> {
        match self {
            Integrator::LoadControl { .. } => Ok((0.0, du_bar)),
            Integrator::DisplacementControl { node, dof, .. } => {
                let eq = domain.equation_of(*node, *dof).ok_or(AnalysisError::InvalidConstraint)?;
                let sensitivity = domain.assemble_reference_load_sensitivity(pseudo_time);
                let unit_response = solver.solve(k, &sensitivity)?;
                let delta_lambda = -du_bar[eq] / unit_response[eq];
                let du = du_bar + unit_response * delta_lambda;
                Ok((delta_lambda, du))
            }
        }
    }
}
