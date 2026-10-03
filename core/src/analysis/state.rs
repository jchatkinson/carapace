use slotmap::Key;

use crate::model::{
    Domain, Element, Element3, ElementOps, Node3Id, NodeId, NDF, PLANAR_NDIM, SPATIAL_ELEMENT_DOF,
    SPATIAL_NDF, SPATIAL_NDIM,
};

use super::arclength::{integrator_change_invalidates, ArcState};
use super::bordered::BorderedSolver;
use super::convergence::ForceTolerance;
use super::{
    iterate_to_equilibrium, Algorithm, AnalysisError, ArcStepInfo, ConstraintHandler,
    ConvergenceTest, Integrator, SparseFactorization, SparseSolver,
};

/// A fully-wired analysis (only buildable via `AnalysisBuilder<Ready>::build`).
///
/// Generic over the same kinematic profile as `Domain` (see its doc
/// comment) — `Analysis`'s step/iteration logic doesn't care about element
/// physics either, only about `Domain`'s generic assembly interface, so one
/// implementation covers both profiles. Defaults to the planar profile, so
/// bare `Analysis` keeps working unchanged; `Analysis3` (spatial) is a type
/// alias below.
pub struct Analysis<
    const NDIM: usize = PLANAR_NDIM,
    const NDOF: usize = NDF,
    const ELEMENT_DOF: usize = { crate::model::ELEMENT_DOF },
    NId = NodeId,
    E = Element,
> where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    pub(crate) domain: Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
    /// Only actually consulted at `AnalysisBuilder<Ready>::build` time (to
    /// reject `Plain` against a domain with multi-point constraints) —
    /// kept here rather than dropped after `build` so `Analysis` still
    /// records which strategy it was built with.
    #[allow(dead_code)]
    pub(crate) constraint_handler: ConstraintHandler,
    pub(crate) integrator: Integrator<NId>,
    pub(crate) algorithm: Algorithm,
    pub(crate) test: ConvergenceTest,
    pub(crate) solver: SparseSolver,
    pub(crate) step_count: usize,
    pub(crate) load_factor: f64,
    /// Cached across `step()` calls (not reset each time) so
    /// `TangentStrategy::Initial` (`docs/algorithms.md` §3) can actually
    /// skip re-factoring for the whole analysis's lifetime, not just
    /// within one step — `iterate_to_equilibrium` is the only thing that
    /// reads/writes this.
    pub(crate) cached_factorization: Option<SparseFactorization>,
    /// Accepted arc-length continuation history (`Integrator::ArcLength`
    /// only): the metric, last accepted direction, next radius, chord
    /// length and stop state. `None` until the first arc step, and cleared
    /// whenever the path definition may have changed (`set_integrator`,
    /// `domain_mut`) so the next step re-validates from scratch.
    pub(crate) arc: Option<Box<ArcState>>,
    /// The bordered system's ordering cache — separate from `solver`'s
    /// `n`-DOF cache so neither evicts the other.
    pub(crate) bordered: BorderedSolver,
}

#[derive(Debug, Clone, Copy)]
pub struct StepResult {
    pub step: usize,
    pub load_factor: f64,
    /// `1` for `Algorithm::Linear` (never iterates). `docs/algorithms.md`
    /// §8's verification hook for `TangentStrategy`/line search.
    pub iterations: usize,
    /// How many times the tangent was actually factored this step — `1`
    /// for `Algorithm::Linear`; proves `TangentStrategy::ReuseAtStepStart`/
    /// `Initial` actually save factorizations, not just "also converges"
    /// (`docs/algorithms.md` §8).
    pub factorizations: usize,
    /// Continuation diagnostics, for `Integrator::ArcLength` steps only.
    /// For those, `iterations` counts corrector solves in the accepted
    /// attempt and `factorizations` counts bordered factorizations across
    /// every attempt, predictor included.
    pub arc: Option<ArcStepInfo>,
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>
    Analysis<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    pub fn domain(&self) -> &Domain<NDIM, NDOF, ELEMENT_DOF, NId, E> {
        &self.domain
    }

    /// Ends this phase and hands back the `Domain` (nodal state, committed
    /// material history, load-pattern freeze state, all intact) so a
    /// caller can compose the next phase — a fresh `AnalysisBuilder::build`
    /// (a different `Integrator`/`Algorithm`/`ConvergenceTest` allowed) or
    /// `TransientAnalysis::new` for a dynamic phase. Typically preceded by
    /// `domain_mut().hold_pattern_constant(pattern, load_factor)` on
    /// whichever pattern(s) should stop ramping (e.g. gravity) before the
    /// next phase starts — see implementation-plan's load-pattern/phase-
    /// composition milestone.
    pub fn into_domain(self) -> Domain<NDIM, NDOF, ELEMENT_DOF, NId, E> {
        self.domain
    }

    /// Mutable access to the domain mid-analysis — needed to call
    /// `Domain::hold_pattern_constant` right before `into_domain()` (the
    /// pattern must be frozen while `self.load_factor`, the pseudo-time to
    /// freeze at, is still known).
    ///
    /// Granting mutable access discards arc-length continuation history and
    /// cached factorizations: the caller may change loads or the model, so
    /// the next arc step re-validates equilibrium and load sensitivity.
    pub fn domain_mut(&mut self) -> &mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E> {
        self.arc = None;
        self.cached_factorization = None;
        &mut self.domain
    }

    /// Swap the integrator without rebuilding `Analysis` — no
    /// `Domain::number_dofs`, no new `SparseSolver`, `step_count`/
    /// `load_factor` carry over unchanged. This is what makes a cyclic/
    /// quasi-static protocol (many small `DisplacementControl` legs with
    /// varying sign/magnitude) cheap: build one `Analysis`, then loop
    /// `analysis.set_integrator(Integrator::DisplacementControl { increment: next, .. }); analysis.step()?;`
    /// over the protocol's prescribed excursions — the protocol itself
    /// (a `Vec<f64>` of signed increments) is caller-side data, not
    /// something `Analysis` needs to model. For a bigger phase transition
    /// (e.g. gravity → pushover, where `Algorithm`/`ConvergenceTest` also
    /// typically change and a pattern needs freezing), use
    /// `into_domain()` + a fresh `AnalysisBuilder` instead.
    ///
    /// Arc-length continuation history survives only a switch between two
    /// `ArcLength` configurations with the same path definition (scales,
    /// direction, seed, predictor) — e.g. new radius bounds or stop
    /// criteria. Anything else starts the next arc step afresh.
    pub fn set_integrator(&mut self, integrator: Integrator<NId>)
    where
        NId: PartialEq,
    {
        if integrator_change_invalidates(&self.integrator, &integrator) {
            self.arc = None;
        } else if let Some(state) = &mut self.arc {
            state.reopen();
        }
        self.integrator = integrator;
    }

    /// Advance one step: the integrator predicts this step's load factor,
    /// then the algorithm resolves equilibrium at that (fixed) load factor.
    /// See implementation-plan §4.4 for the sketch this follows.
    ///
    /// On failure, the domain (nodal displacement — mutated eagerly every
    /// Newton iteration, which is correct Newton behavior, but shouldn't
    /// be left half-applied if the step as a whole doesn't converge) is
    /// restored to its state before this call. Materials never need this:
    /// see `Material`'s doc comment for why they're never mutated until
    /// `Domain::commit`, which only runs on the success path below.
    pub fn step(&mut self) -> Result<StepResult, AnalysisError> {
        self.step_prescribed(&[])
    }

    /// `step`, after first imposing `(node, dof, value)` displacements on
    /// fixed DOFs (non-homogeneous single-point constraints, see
    /// `Domain::prescribe_displacement`): the free DOFs then equilibrate
    /// against the displaced supports at this step's load factor. Pair it
    /// with `Integrator::LoadControl { increment: 0.0 }` for a pure
    /// prescribed-displacement step that needs no reference load and no
    /// tangent at the imposed DOF — unlike `DisplacementControl`, it is
    /// well-defined on a zero-tangent branch. Read the force back with
    /// `Domain::reaction`.
    ///
    /// Targets are validated before anything changes: an unfixed DOF or
    /// nonfinite value is `InvalidConstraint` and leaves the analysis
    /// untouched. Once validated, the update is atomic with the step — on
    /// failure the imposed displacements are rolled back along with the
    /// rest of the domain, so the last committed state stays intact for a
    /// retry. Arc-length continuation has no meaningful prescribed-value
    /// step and rejects any targets.
    pub fn step_prescribed(
        &mut self,
        targets: &[(NId, usize, f64)],
    ) -> Result<StepResult, AnalysisError> {
        if !targets
            .iter()
            .all(|&(node, dof, value)| self.domain.can_prescribe(node, dof, value))
        {
            return Err(AnalysisError::InvalidConstraint);
        }
        if let Integrator::ArcLength(config) = &self.integrator {
            if !targets.is_empty() {
                return Err(AnalysisError::InvalidConstraint);
            }
            let config = config.clone();
            return self.arc_step(&config);
        }
        let snapshot = self.domain.clone();
        for &(node, dof, value) in targets {
            self.domain.prescribe_displacement(node, dof, value);
        }
        if !targets.is_empty() {
            // The tangent at the displaced state is not the cached one.
            self.cached_factorization = None;
        }
        let (snapshot_step_count, snapshot_load_factor) = (self.step_count, self.load_factor);

        let result = self.try_step();

        if result.is_err() {
            self.domain = snapshot;
            self.step_count = snapshot_step_count;
            self.load_factor = snapshot_load_factor;
        } else {
            self.domain.commit(self.load_factor);
        }

        result
    }

    fn try_step(&mut self) -> Result<StepResult, AnalysisError> {
        self.step_count += 1;
        // `Combined`'s per-equation reference (committed and predicted
        // external force, committed internal force), frozen for the step.
        let committed_forces = matches!(self.test, ConvergenceTest::Combined { .. }).then(|| {
            let (_k, internal) = self.domain.assemble_tangent_and_resistance(self.load_factor);
            (self.domain.assemble_reference_load(self.load_factor), internal)
        });
        self.load_factor = self
            .integrator
            .predict(&self.domain, &self.solver, self.load_factor)?;
        let force_tolerance = committed_forces.and_then(|(external, internal)| {
            let predicted = self.domain.assemble_reference_load(self.load_factor);
            ForceTolerance::for_test(
                &self.test,
                &self.domain.rotational_equations(),
                &[&external, &predicted, &internal],
            )
        });

        let (iterations, factorizations) = match self.algorithm {
            // A single tangent formation + solve, unconditionally accepted
            // — definitionally converged for `Linear` (see algorithm.rs).
            // Kept as its own trivial arm rather than routed through
            // `iterate_to_equilibrium` — there's no iteration to share the
            // driver's machinery with (`docs/algorithms.md` §5's driver
            // doc comment).
            Algorithm::Linear => {
                let (k, residual) = self.domain.form_tangent_and_residual(self.load_factor);
                let du = self.solver.solve(&k, &residual)?;
                self.domain.apply_displacement_increment(&du);
                (1, 1)
            }
            Algorithm::Newton { .. } | Algorithm::KrylovNewton { .. } => {
                let integrator = &self.integrator;
                let solver = &self.solver;
                let (outcome, load_factor) = iterate_to_equilibrium(
                    self.step_count,
                    &self.algorithm,
                    &self.test,
                    solver,
                    &mut self.cached_factorization,
                    &mut self.domain,
                    self.load_factor,
                    |domain, load_factor| domain.form_tangent_and_residual(load_factor),
                    |domain, k, current_factorization, du_bar, load_factor, iteration| {
                        // The first iteration uses the same tangent
                        // `predict` used, so `du_bar` already delivers
                        // `predict`'s target displacement at the
                        // controlled DOF exactly (see `Integrator::
                        // correct`'s doc comment) — only from the second
                        // iteration on does the controlled DOF need to be
                        // actively held there while other DOFs still get
                        // corrected.
                        if iteration == 0 {
                            Ok((du_bar, 0.0))
                        } else {
                            let (delta_lambda, du) = integrator.correct(
                                domain, solver, k, current_factorization, du_bar, load_factor,
                            )?;
                            Ok((du, delta_lambda))
                        }
                    },
                    |_domain| {},
                    force_tolerance.as_ref(),
                )?;
                self.load_factor = load_factor;
                (outcome.iterations, outcome.factorizations)
            }
        };

        Ok(StepResult {
            step: self.step_count,
            load_factor: self.load_factor,
            iterations,
            factorizations,
            arc: None,
        })
    }
}

/// `Analysis`'s spatial instantiation — see `Domain3`'s doc comment for why
/// this is a type alias rather than a hand-duplicated struct.
pub type Analysis3 = Analysis<SPATIAL_NDIM, SPATIAL_NDF, SPATIAL_ELEMENT_DOF, Node3Id, Element3>;
