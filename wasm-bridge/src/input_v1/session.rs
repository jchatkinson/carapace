//! `Session`/`StageRunner`/`advance` — pysees-handoff.md's "Decode and
//! session model": `Session` is chosen once at decode time from the
//! header's `space` discriminant and never branches on it again; stepping
//! is driven by the caller via a step budget rather than run to completion
//! inside one call, the mechanism cooperative cancellation is built on.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest,
    Integrator,
};
use carapace_core::model::{Domain, ElementId, ElementLoad, LoadPatternId, NodeId};
use serde::Serialize;

/// `AnalysisError`'s fields, restated so `advance`'s result doesn't need to
/// name `carapace_core`'s error type directly — kept in the same tagged-
/// variant style (implementation-plan.md §2.8).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum AnalysisErrorDetail {
    FailedToConverge { step: usize },
    SingularSystem,
    InvalidConstraint,
    InvalidModeCount { requested: usize, free_dofs: usize },
}

impl From<AnalysisError> for AnalysisErrorDetail {
    fn from(error: AnalysisError) -> Self {
        match error {
            AnalysisError::FailedToConverge { step } => {
                AnalysisErrorDetail::FailedToConverge { step }
            }
            AnalysisError::SingularSystem => AnalysisErrorDetail::SingularSystem,
            AnalysisError::InvalidConstraint => AnalysisErrorDetail::InvalidConstraint,
            AnalysisError::InvalidModeCount {
                requested,
                free_dofs,
            } => AnalysisErrorDetail::InvalidModeCount {
                requested,
                free_dofs,
            },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolvedRecorder {
    pub node: NodeId,
    pub dof: u8,
}

pub(super) struct CompiledStage {
    pub id: String,
    pub steps: u32,
    pub kind: CompiledStageKind,
    /// Registered on the `Domain` only when this stage starts — see
    /// `tables::NodalLoadTable::stage`.
    pub pending_nodal_loads: Vec<(LoadPatternId, NodeId, usize, f64)>,
    pub pending_element_loads: Vec<(LoadPatternId, ElementId, ElementLoad)>,
}

pub(super) enum CompiledStageKind {
    Static {
        integrator: Integrator<NodeId>,
        algorithm: Algorithm,
        convergence: ConvergenceTest,
        hold_patterns_after: Vec<LoadPatternId>,
    },
}

enum StageRunner {
    Static {
        analysis: Analysis,
        steps_remaining: u32,
    },
}

/// `advance`'s result — pysees-handoff.md's `{ done, stageComplete,
/// stepsTaken, progressSnapshot, recorderBatch? }`, minus the parts that
/// are the worker/JS boundary's job (`progressSnapshot`'s throttling,
/// `recorderBatch`'s `response_blocks` byte layout): `load_factor` here is
/// the raw signal those would be built from. `error`, once set, is sticky —
/// later stages are not attempted, matching "a stage's `AnalysisError`
/// stops the sequence".
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct StepOutcome {
    pub done: bool,
    pub stage_complete: bool,
    pub steps_taken: u32,
    pub load_factor: f64,
    pub error: Option<AnalysisErrorDetail>,
}

/// One opaque handle exposed to the caller, chosen once from `space` at
/// decode time. A `Spatial` arm is deliberately not modeled yet — decode
/// rejects `space == 3` before any `Session` is constructed, so there is
/// nothing this enum needs to branch on downstream today.
#[derive(Debug)]
pub enum Session {
    Planar(PlanarSession),
}

impl Session {
    pub fn advance(&mut self, step_budget: u32) -> StepOutcome {
        match self {
            Session::Planar(session) => session.advance(step_budget),
        }
    }

    pub fn recorder_samples(&self, recorder_index: usize) -> &[(f64, f64)] {
        match self {
            Session::Planar(session) => session.recorder_samples(recorder_index),
        }
    }

    /// The `AnalysisSequence` stage id `advance` is currently in (or, once
    /// `done`, the id of whichever stage stopped it) — for surfacing which
    /// stage an `error` belongs to, matching pysees-handoff.md's "the run
    /// is marked failed with the error's structured detail".
    pub fn current_stage_id(&self) -> Option<&str> {
        match self {
            Session::Planar(session) => session.current_stage_id(),
        }
    }
}

pub struct PlanarSession {
    domain: Option<Domain>,
    stages: Vec<CompiledStage>,
    current_stage: usize,
    runner: Option<StageRunner>,
    load_factor: f64,
    error: Option<AnalysisErrorDetail>,
    recorders: Vec<ResolvedRecorder>,
    recorder_history: Vec<Vec<(f64, f64)>>,
}

impl std::fmt::Debug for PlanarSession {
    // `Analysis` (inside `StageRunner`, held via `runner`/`domain`) doesn't
    // derive `Debug`, so this reports the fields useful for a test/error
    // message rather than the full solver state.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PlanarSession")
            .field("current_stage", &self.current_stage)
            .field("load_factor", &self.load_factor)
            .field("error", &self.error)
            .finish()
    }
}

impl PlanarSession {
    pub(super) fn new(
        domain: Domain,
        stages: Vec<CompiledStage>,
        recorders: Vec<ResolvedRecorder>,
    ) -> Self {
        let recorder_history = vec![Vec::new(); recorders.len()];
        Self {
            domain: Some(domain),
            stages,
            current_stage: 0,
            runner: None,
            load_factor: 0.0,
            error: None,
            recorders,
            recorder_history,
        }
    }

    pub fn recorder_samples(&self, recorder_index: usize) -> &[(f64, f64)] {
        &self.recorder_history[recorder_index]
    }

    pub fn current_stage_id(&self) -> Option<&str> {
        self.stages
            .get(self.current_stage)
            .map(|stage| stage.id.as_str())
    }

    /// Steps up to `step_budget` times, stopping early if the current stage
    /// completes or a step fails. Each call only ever executes up to
    /// `step_budget` steps before returning — the mechanism a worker's
    /// cooperative-cancellation loop is built on (checking a pending cancel
    /// between calls), even though there is no such loop calling this yet.
    pub fn advance(&mut self, step_budget: u32) -> StepOutcome {
        if self.error.is_none() && self.runner.is_none() && self.current_stage < self.stages.len() {
            self.start_current_stage();
        }
        if self.error.is_some() || self.current_stage >= self.stages.len() {
            return StepOutcome {
                done: true,
                stage_complete: false,
                steps_taken: 0,
                load_factor: self.load_factor,
                error: self.error,
            };
        }

        let mut steps_taken = 0;
        let mut stage_complete = false;
        while steps_taken < step_budget {
            let StageRunner::Static {
                analysis,
                steps_remaining,
            } = self.runner.as_mut().expect("stage started above");
            if *steps_remaining == 0 {
                stage_complete = true;
                break;
            }
            match analysis.step() {
                Ok(result) => {
                    self.load_factor = result.load_factor;
                    self.record_sample();
                    let StageRunner::Static {
                        steps_remaining, ..
                    } = self.runner.as_mut().expect("still running");
                    *steps_remaining -= 1;
                    steps_taken += 1;
                    if *steps_remaining == 0 {
                        stage_complete = true;
                        break;
                    }
                }
                Err(error) => {
                    self.error = Some(error.into());
                    break;
                }
            }
        }

        if stage_complete {
            self.finish_current_stage();
        }

        StepOutcome {
            done: self.error.is_some() || self.current_stage >= self.stages.len(),
            stage_complete,
            steps_taken,
            load_factor: self.load_factor,
            error: self.error,
        }
    }

    fn start_current_stage(&mut self) {
        let stage = &self.stages[self.current_stage];
        let mut domain = self.domain.take().expect("domain held between stages");
        for &(pattern, node, dof, value) in &stage.pending_nodal_loads {
            domain.add_nodal_load(pattern, node, dof, value);
        }
        for &(pattern, element, load) in &stage.pending_element_loads {
            domain.add_element_load(pattern, element, load);
        }
        match &stage.kind {
            // `ConstraintHandler::Plain` is hardcoded: M10's wire format has
            // no multi-point-constraint tables yet (`equal_dof`/
            // `rigid_diaphragm` are spatial-profile concerns, M20), so
            // decode never produces a domain `Transformation` would be
            // needed for.
            CompiledStageKind::Static {
                integrator,
                algorithm,
                convergence,
                ..
            } => {
                let analysis = AnalysisBuilder::new()
                    .constraint_handler(ConstraintHandler::Plain)
                    .integrator(*integrator)
                    .algorithm(*algorithm)
                    .test(*convergence)
                    .build(domain);
                self.runner = Some(StageRunner::Static {
                    analysis,
                    steps_remaining: stage.steps,
                });
            }
        }
        self.load_factor = 0.0;
    }

    /// `Analysis::into_domain`/fresh-`AnalysisBuilder` multi-phase pattern
    /// (pysees-handoff.md / core/tests/m8_force_beam_column.rs's two-phase
    /// tests): freeze this stage's held patterns at its final load factor,
    /// then hand the same `Domain` (committed material state and all) to
    /// the next stage's fresh `Analysis`.
    fn finish_current_stage(&mut self) {
        let CompiledStageKind::Static {
            hold_patterns_after,
            ..
        } = &self.stages[self.current_stage].kind;
        let hold_patterns_after = hold_patterns_after.clone();
        let load_factor = self.load_factor;

        let StageRunner::Static { analysis, .. } = self.runner.take().expect("stage running");
        let mut domain = analysis.into_domain();
        for pattern in hold_patterns_after {
            domain.hold_pattern_constant(pattern, load_factor);
        }
        self.domain = Some(domain);
        self.current_stage += 1;
    }

    fn record_sample(&mut self) {
        let Some(StageRunner::Static { analysis, .. }) = &self.runner else {
            return;
        };
        let load_factor = self.load_factor;
        for (recorder, history) in self.recorders.iter().zip(self.recorder_history.iter_mut()) {
            let value = analysis.domain().node(recorder.node).displacement[recorder.dof as usize];
            history.push((load_factor, value));
        }
    }
}
