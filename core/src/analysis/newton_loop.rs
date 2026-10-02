use nalgebra::DVector;
use slotmap::Key;

use crate::model::{Domain, ElementOps, SparseMatrix};

use super::convergence::{ForceTolerance, IterationContext};
use super::{
    Algorithm, AnalysisError, ConvergenceTest, KrylovAccelerator, SparseFactorization,
    SparseSolver, TangentStrategy,
};

/// How many Newton iterations ran, and how many times the tangent was
/// actually factored — `docs/algorithms.md` §8's verification hook for
/// `TangentStrategy`/`KrylovNewton` (proving reuse/acceleration actually
/// saves factorizations, not just "also converges").
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct IterationOutcome {
    pub(crate) iterations: usize,
    pub(crate) factorizations: usize,
}

/// The shared Newton-family driver (`docs/algorithms.md` §0): both static
/// `Analysis` and `TransientAnalysis` reduce to "solve `Op(u)*du = r(u)`,
/// apply `du`, check convergence, repeat" against some *effective
/// operator* `Op` — `TangentStrategy`/`LineSearch`/`KrylovAccelerator`
/// only ever look at `Op`/`r`/`du`, never at what physically produced them
/// (bare static tangent `K`, or transient's `c1*K + c2*C + c3*M`). Only
/// `Algorithm::Newton`/`KrylovNewton` — `Algorithm::Linear`'s one-shot,
/// no-iteration case stays a trivial inline arm in each caller, not routed
/// through here.
///
/// `scalar` is a caller-owned auxiliary value threaded through
/// `form_system`/`correct_du` and updated by `correct_du`'s returned delta
/// — `Analysis` uses it for the load factor (`Integrator::
/// DisplacementControl`'s corrector adjusts it every iteration, per its
/// own doc comment); `TransientAnalysis` (§6.1) has no analogous
/// solved-for scalar (`time` is fixed for the whole step by `dt`), so it
/// passes `0.0` and a `correct_du` that always returns `(du, 0.0)`
/// unchanged. This is what keeps `Integrator::DisplacementControl` working
/// without smuggling `Analysis`-specific concepts into this otherwise
/// fully shared loop — `form_system`/`correct_du` are the only two hooks
/// that see `scalar`, and the driver itself never interprets it.
///
/// Convergence is checked at the accepted state of each iteration — after
/// the correction (and any line-search rescaling) is applied and the
/// system reassembled — and that assembly doubles as the next iteration's
/// system. Force-based tests also check the predictor before any solve.
/// `force_tolerance` must be supplied for `ConvergenceTest::Combined`
/// (built once per step by the caller, see `ForceTolerance::for_test`).
///
/// `form_system` reforms `(Op, r)` at the domain's *current* trial state —
/// called every iteration regardless of `TangentStrategy` (assembly cost
/// isn't saved by tangent reuse, only factorization cost is; `docs/
/// algorithms.md` §3.1's "ship the factorization split alone first"
/// decision). `after_increment` runs after every applied displacement
/// increment (including line-search trial evaluations, so a transient
/// caller's velocity/acceleration stay consistent with whatever
/// displacement is actually on the domain at every point `form_system`
/// might read it back — not just at the end of each outer iteration).
///
/// `cached_factorization` is owned by the caller (a field on `Analysis`/
/// `TransientAnalysis`, persisting across `step()` calls) so
/// `TangentStrategy::Initial` can actually skip re-factoring across every
/// step of an analysis's lifetime, not just within one step.
/// `correct_du` receives a factorization only when it was formed from
/// this iteration's operator, allowing additional right-hand sides to
/// reuse it without accidentally using an older tangent.
#[allow(clippy::too_many_arguments)]
pub(crate) fn iterate_to_equilibrium<
    const NDIM: usize,
    const NDOF: usize,
    const ELEMENT_DOF: usize,
    NId,
    E,
>(
    step: usize,
    algorithm: &Algorithm,
    test: &ConvergenceTest,
    solver: &SparseSolver,
    cached_factorization: &mut Option<SparseFactorization>,
    domain: &mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
    scalar0: f64,
    mut form_system: impl FnMut(
        &Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        f64,
    ) -> (SparseMatrix, DVector<f64>),
    mut correct_du: impl FnMut(
        &Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        &SparseMatrix,
        Option<&SparseFactorization>,
        DVector<f64>,
        f64,
        usize,
    ) -> Result<(DVector<f64>, f64), AnalysisError>,
    mut after_increment: impl FnMut(&mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>),
    force_tolerance: Option<&ForceTolerance>,
) -> Result<(IterationOutcome, f64), AnalysisError>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    test.validate()?;
    let (tangent, line_search, krylov_max_dimension) = match algorithm {
        Algorithm::Newton {
            tangent,
            line_search,
        } => (*tangent, *line_search, None),
        Algorithm::KrylovNewton {
            tangent,
            max_dimension,
        } => (*tangent, None, Some(*max_dimension)),
        Algorithm::Linear => {
            unreachable!("Algorithm::Linear never routes through iterate_to_equilibrium")
        }
    };

    // `ReuseAtStepStart` means "once per this step" — whatever was cached
    // from a *previous* step is stale for this one, so force iteration 0
    // to factor fresh. `Initial` leaves the cache as-is: only factor if
    // genuinely never factored before (this step or any previous one).
    if tangent == TangentStrategy::ReuseAtStepStart {
        *cached_factorization = None;
    }

    let mut krylov = krylov_max_dimension.map(KrylovAccelerator::new);
    let mut outcome = IterationOutcome::default();
    let mut scalar = scalar0;
    let mut converged = false;

    // Every convergence check below sees the residual reassembled at the
    // state the iteration actually accepted (after its correction and any
    // line-search rescaling), and that same assembly is reused as the
    // next iteration's system. The predictor itself is checked first by
    // the force-based tests, so an already-balanced state costs no solve.
    after_increment(domain);
    let (mut k, mut residual) = form_system(domain, scalar);
    if test.accepts_predictor() {
        let zero = DVector::zeros(residual.len());
        converged = test.check(&IterationContext {
            residual: &residual,
            du: &zero,
            force_tolerance,
        });
    }

    let mut iteration = 0;
    while !converged && iteration < test.max_iter() {
        let force_refactor = tangent == TangentStrategy::Current
            || cached_factorization.is_none()
            || krylov
                .as_ref()
                .is_some_and(KrylovAccelerator::should_refactor);
        if force_refactor {
            *cached_factorization = Some(solver.factor(&k)?);
            outcome.factorizations += 1;
            if let Some(dim) = krylov_max_dimension {
                krylov = Some(KrylovAccelerator::new(dim));
            }
        }
        let factorization = cached_factorization
            .as_ref()
            .expect("just factored or already cached above");
        let du_raw = factorization.solve(&residual);

        let du_pre_correct = match &mut krylov {
            Some(accelerator) => accelerator.accelerate(du_raw, &residual),
            None => du_raw,
        };
        let (du, delta_scalar) = correct_du(
            domain,
            &k,
            force_refactor.then_some(factorization),
            du_pre_correct,
            scalar,
            iteration,
        )?;
        scalar += delta_scalar;
        let residual_before = residual;

        domain.apply_displacement_increment(&du);
        after_increment(domain);
        outcome.iterations += 1;
        iteration += 1;
        (k, residual) = form_system(domain, scalar);

        if test.check(&IterationContext {
            residual: &residual,
            du: &du,
            force_tolerance,
        }) {
            converged = true;
            break;
        }

        if let Some(ls) = line_search {
            let s0 = du.dot(&residual_before);
            let s1 = du.dot(&residual);
            let mut applied_eta = 1.0_f64;
            let mut last = None;
            let eta = ls.resolve(s0, s1, |eta| {
                let delta = (eta - applied_eta) * &du;
                domain.apply_displacement_increment(&delta);
                after_increment(domain);
                applied_eta = eta;
                let system = form_system(domain, scalar);
                let s = du.dot(&system.1);
                last = Some(system);
                s
            });
            if let Some(system) = last {
                (k, residual) = system;
            }
            // `resolve` can return a bracket midpoint or endpoint it never
            // evaluated last; move the domain there so the accepted state,
            // its residual and the next iteration's system all agree.
            if eta != applied_eta {
                let delta = (eta - applied_eta) * &du;
                domain.apply_displacement_increment(&delta);
                after_increment(domain);
                (k, residual) = form_system(domain, scalar);
            }
            let final_du = eta * &du;
            if test.check(&IterationContext {
                residual: &residual,
                du: &final_du,
                force_tolerance,
            }) {
                converged = true;
                break;
            }
        }
    }

    if !converged {
        return Err(AnalysisError::FailedToConverge { step });
    }
    Ok((outcome, scalar))
}
