use nalgebra::DVector;
use slotmap::Key;

use crate::model::{Domain, ElementOps, NodeId, SparseMatrix};

use super::{AnalysisError, ArcLength, SparseFactorization, SparseSolver};

/// Static/pseudo-static integration strategy. Closed enum.
/// `LoadControl`/`DisplacementControl` only change how this *step's* load
/// factor is chosen (the "predictor") — the iteration that follows
/// (see `Algorithm`) proceeds identically regardless of which one picked it.
/// `ArcLength` is different in kind: it predicts displacement *and* load
/// factor and corrects both together through a bordered system, so
/// `Analysis::step` hands it to its own continuation driver instead of
/// `predict`/`correct` (see `arclength.rs` and `docs/arclength.md`).
///
/// Generic over `NId` (the controlled node's ID type — `NodeId` or
/// `Node3Id`), defaulted to the planar profile so bare `Integrator` keeps
/// working, the same trick `Domain` uses. `predict`/`correct` take the rest
/// of the profile (`NDIM`/`NDOF`/`ELEMENT_DOF`/`E`) as their own
/// generic parameters instead of `Integrator`'s, since an integrator's
/// identity only depends on which node/DOF it controls, not which `Domain`
/// it's later used with.
#[derive(Debug, Clone, PartialEq)]
pub enum Integrator<NId = NodeId> {
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
        node: NId,
        dof: usize,
        increment: f64,
    },
    /// Arc-length continuation along the equilibrium path, through load
    /// limits and snap-back — see `ArcLength`.
    ArcLength(ArcLength<NId>),
}

impl<NId: Copy> Integrator<NId> {
    /// Compute this step's new pseudo-time, given the domain's state as of
    /// the end of the *previous* (converged) step, plus the displacement
    /// predictor the caller must apply before iterating (`None` for
    /// `LoadControl`). `DisplacementControl`'s predictor is the unit-load
    /// response scaled by the pseudo-time increment, so the controlled DOF
    /// starts the step exactly at its target and every Newton iteration
    /// (the first included) measures a true correction. That is what lets
    /// a linear step converge in one iteration instead of spending its
    /// first on the predictor and a second on proving `du` is small.
    pub(crate) fn predict<const NDIM: usize, const NDOF: usize, E>(
        &self,
        domain: &Domain<NDIM, NDOF, NId, E>,
        solver: &SparseSolver,
        current_pseudo_time: f64,
    ) -> Result<(f64, Option<DVector<f64>>), AnalysisError>
    where
        NId: Key,
        E: ElementOps<NDIM, NDOF, NId>,
    {
        match self {
            Integrator::LoadControl { increment } => Ok((current_pseudo_time + increment, None)),
            Integrator::DisplacementControl {
                node,
                dof,
                increment,
            } => {
                let eq = domain
                    .equation_of(*node, *dof)
                    .ok_or(AnalysisError::InvalidConstraint)?;
                let (k, _resistance) = domain.assemble_tangent_and_resistance(current_pseudo_time);
                let sensitivity = domain.assemble_reference_load_sensitivity(current_pseudo_time);
                let unit_response = solver.solve(&k, &sensitivity)?;
                let delta_lambda = increment / unit_response[eq];
                Ok((
                    current_pseudo_time + delta_lambda,
                    Some(unit_response * delta_lambda),
                ))
            }
            Integrator::ArcLength(_) => {
                unreachable!("Analysis::step routes ArcLength to its continuation driver")
            }
        }
    }

    /// Corrector for every Newton iteration of a step, applied after
    /// `predict`'s displacement predictor.
    /// `predict`'s pseudo-time only accounts for the tangent at the *start*
    /// of the step — exact for `LoadControl` (whose pseudo-time never
    /// changes mid-step, so `du_bar` unmodified is already the right
    /// increment every iteration).
    ///
    /// For `DisplacementControl`, `predict` already moved the controlled DOF
    /// to its target (OpenSees' `newStep` does the same), but the tangent it
    /// used can go stale as the iteration proceeds (e.g. a step landing
    /// exactly on a material's yield breakpoint sees the elastic tangent at
    /// `predict` time, then the true, much softer, post-yield tangent) —
    /// left uncorrected, ordinary fixed-load-factor Newton iteration would
    /// then keep driving the controlled DOF *past* the target already
    /// reached, chasing force equilibrium at the wrong load factor instead.
    /// This is the standard Yang & Shieh consistent Displacement Control
    /// corrector: solve the current tangent against the reference load
    /// sensitivity again (`unit_response`, this iteration's version of
    /// `predict`'s unit-load probe), then pick `delta_lambda` so that adding
    /// `delta_lambda * unit_response` to `du_bar` cancels `du_bar`'s
    /// contribution at the controlled DOF — i.e. this iteration leaves the
    /// controlled DOF exactly where it already is, no matter how far the
    /// tangent has drifted, while every other DOF still gets `du_bar`'s
    /// residual-driven correction.
    /// `current_factorization`, when supplied, must factor `k`. An older
    /// Newton tangent must not be passed here: the load-sensitivity probe
    /// continues to use the current tangent even under tangent reuse.
    pub(crate) fn correct<const NDIM: usize, const NDOF: usize, E>(
        &self,
        domain: &Domain<NDIM, NDOF, NId, E>,
        solver: &SparseSolver,
        k: &SparseMatrix,
        current_factorization: Option<&SparseFactorization>,
        du_bar: DVector<f64>,
        pseudo_time: f64,
    ) -> Result<(f64, DVector<f64>), AnalysisError>
    where
        NId: Key,
        E: ElementOps<NDIM, NDOF, NId>,
    {
        match self {
            Integrator::LoadControl { .. } => Ok((0.0, du_bar)),
            Integrator::DisplacementControl { node, dof, .. } => {
                let eq = domain
                    .equation_of(*node, *dof)
                    .ok_or(AnalysisError::InvalidConstraint)?;
                let sensitivity = domain.assemble_reference_load_sensitivity(pseudo_time);
                let unit_response = match current_factorization {
                    Some(factorization) => factorization.solve(&sensitivity),
                    None => solver.solve(k, &sensitivity)?,
                };
                let delta_lambda = -du_bar[eq] / unit_response[eq];
                let du = du_bar + unit_response * delta_lambda;
                Ok((delta_lambda, du))
            }
            Integrator::ArcLength(_) => {
                unreachable!("Analysis::step routes ArcLength to its continuation driver")
            }
        }
    }
}
