use slotmap::Key;

use crate::model::{
    Domain, Element, Element3, ElementOps, Node3Id, NodeId, PLANAR_NDIM, SPATIAL_ELEMENT_DOF, SPATIAL_NDF,
    SPATIAL_NDIM, NDF,
};

use super::{Algorithm, AnalysisError, ConstraintHandler, ConvergenceTest, Integrator, SparseSolver};

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
}

#[derive(Debug, Clone, Copy)]
pub struct StepResult {
    pub step: usize,
    pub load_factor: f64,
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E> Analysis<NDIM, NDOF, ELEMENT_DOF, NId, E>
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
    pub fn domain_mut(&mut self) -> &mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E> {
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
    pub fn set_integrator(&mut self, integrator: Integrator<NId>) {
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
        let snapshot = self.domain.clone();
        let (snapshot_step_count, snapshot_load_factor) = (self.step_count, self.load_factor);

        let result = self.try_step();

        if result.is_err() {
            self.domain = snapshot;
            self.step_count = snapshot_step_count;
            self.load_factor = snapshot_load_factor;
        } else {
            self.domain.commit();
        }

        result
    }

    fn try_step(&mut self) -> Result<StepResult, AnalysisError> {
        self.step_count += 1;
        self.load_factor = self.integrator.predict(&self.domain, &self.solver, self.load_factor)?;

        match self.algorithm {
            // A single tangent formation + solve, unconditionally accepted
            // — definitionally converged for `Linear` (see algorithm.rs).
            Algorithm::Linear => {
                let (k, residual) = self.domain.form_tangent_and_residual(self.load_factor);
                let du = self.solver.solve(&k, &residual)?;
                self.domain.apply_displacement_increment(&du);
            }
            Algorithm::NewtonRaphson => {
                let mut converged = false;
                for iteration in 0..self.test.max_iter() {
                    let (k, residual) = self.domain.form_tangent_and_residual(self.load_factor);
                    let du_bar = self.solver.solve(&k, &residual)?;
                    // The first iteration uses the same tangent `predict`
                    // used, so `du_bar` already delivers `predict`'s target
                    // displacement at the controlled DOF exactly (see
                    // `Integrator::correct`'s doc comment) — only from the
                    // second iteration on does the controlled DOF need to
                    // be actively held there while other DOFs still get
                    // corrected.
                    let (delta_lambda, du) = if iteration == 0 {
                        (0.0, du_bar)
                    } else {
                        self.integrator.correct(&self.domain, &self.solver, &k, du_bar, self.load_factor)?
                    };
                    self.load_factor += delta_lambda;
                    self.domain.apply_displacement_increment(&du);
                    if self.test.check(&residual, &du) {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    return Err(AnalysisError::FailedToConverge { step: self.step_count });
                }
            }
        }

        Ok(StepResult {
            step: self.step_count,
            load_factor: self.load_factor,
        })
    }
}

/// `Analysis`'s spatial instantiation — see `Domain3`'s doc comment for why
/// this is a type alias rather than a hand-duplicated struct.
pub type Analysis3 = Analysis<SPATIAL_NDIM, SPATIAL_NDF, SPATIAL_ELEMENT_DOF, Node3Id, Element3>;
