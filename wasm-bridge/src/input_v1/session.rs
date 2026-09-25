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

/// A `RecorderSpec` resolved against a decoded `Domain` — real `NodeId`/
/// `ElementId` handles instead of wire-format table indices. Growing this
/// by one more response kind (see `RecorderSpec`'s doc comment) is one more
/// variant plus one more `record_sample` match arm, not a new field
/// anywhere on `PlanarSession`/`StepOutcome`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResolvedRecorder {
    NodeDisp { node: NodeId, dof: u8 },
    ElementForce { element: ElementId, component: u8 },
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
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct StepOutcome {
    pub done: bool,
    pub stage_complete: bool,
    pub steps_taken: u32,
    pub load_factor: f64,
    pub error: Option<AnalysisErrorDetail>,
    /// Only the samples *this* `advance()` call produced, one entry per recorder that recorded
    /// at least one sample this call (every recorder samples every step today, so in practice
    /// this is either empty — no step taken — or has one entry per recorder). This is the
    /// results-storage plan's (docs/results-storage-indexeddb.md) `recorderBatch`, restated in
    /// field names that match pysees's frozen `ResultBlock` protocol type
    /// (`src/app/types/resultsStorage.ts`); the caller assigns `blockIndex` and packs `data`,
    /// since neither concept exists on this side of the wasm boundary.
    pub recorder_batches: Vec<RecorderBatch>,
}

/// One recorder's new samples from a single `advance()` call. `first_sample` plus
/// `samples.len()` gives the sample range this batch covers, so a caller can retry a failed
/// persist without renumbering — the "each `advance()` batch deterministic and identifies its
/// first sample/count" contract from results-storage-indexeddb.md's "Failure and cancellation
/// contract". Scalar-only (`(pseudo_time, value)` pairs, one channel) until a vector recorder
/// needs more.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RecorderBatch {
    pub recorder_index: usize,
    pub stage_index: usize,
    pub first_sample: u32,
    pub samples: Vec<(f64, f64)>,
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
    // Per-advance only: cleared at the top of every `advance()` call, drained into that call's
    // `StepOutcome::recorder_batches` at the end. Nothing here survives across calls except the
    // running counts below — this is results-storage-indexeddb.md's "Bound memory" step:
    // `recorder_history`'s old whole-run `Vec<Vec<(f64, f64)>>` is gone.
    current_batch: Vec<Vec<(f64, f64)>>,
    // Total samples emitted so far per recorder, so the next batch's `first_sample` is correct
    // without retaining the samples themselves.
    recorder_sample_counts: Vec<u32>,
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
        let current_batch = vec![Vec::new(); recorders.len()];
        let recorder_sample_counts = vec![0; recorders.len()];
        Self {
            domain: Some(domain),
            stages,
            current_stage: 0,
            runner: None,
            load_factor: 0.0,
            error: None,
            recorders,
            current_batch,
            recorder_sample_counts,
        }
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
        for batch in &mut self.current_batch {
            batch.clear();
        }

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
                recorder_batches: Vec::new(),
            };
        }

        // Captured before stepping: `finish_current_stage` below advances `current_stage`, but
        // every sample this call records (the stepping loop never starts a new stage mid-call)
        // belongs to the stage active when the call began.
        let stage_index_for_batch = self.current_stage;

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

        let recorder_batches = self
            .current_batch
            .iter()
            .zip(self.recorder_sample_counts.iter_mut())
            .enumerate()
            .filter_map(|(recorder_index, (samples, sample_count))| {
                if samples.is_empty() {
                    return None;
                }
                let first_sample = *sample_count;
                *sample_count += samples.len() as u32;
                Some(RecorderBatch {
                    recorder_index,
                    stage_index: stage_index_for_batch,
                    first_sample,
                    samples: samples.clone(),
                })
            })
            .collect();

        StepOutcome {
            done: self.error.is_some() || self.current_stage >= self.stages.len(),
            stage_complete,
            steps_taken,
            load_factor: self.load_factor,
            error: self.error,
            recorder_batches,
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
        for (recorder, batch) in self.recorders.iter().zip(self.current_batch.iter_mut()) {
            let value = match *recorder {
                ResolvedRecorder::NodeDisp { node, dof } => analysis.domain().node(node).displacement[dof as usize],
                ResolvedRecorder::ElementForce { element, component } => {
                    analysis.domain().element_local_force(element)[component as usize]
                }
            };
            batch.push((load_factor, value));
        }
    }
}
