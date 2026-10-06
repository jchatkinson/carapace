//! `Session`/`StageRunner`/`advance` — pysees-handoff.md's "Decode and
//! session model": `Session` is chosen once at decode time from the
//! header's `space` discriminant and never branches on it again; stepping
//! is driven by the caller via a step budget rather than run to completion
//! inside one call, the mechanism cooperative cancellation is built on.
//!
//! `ModelSession<NDIM, NDOF, ELEMENT_DOF, NId, E>` is generic over the same
//! kinematic profile `core`'s own `Domain`/`Analysis` are generic over
//! (`AnalysisBuilder::build` already infers that whole profile from its
//! `domain` argument — see its doc comment) — the planar/spatial split only
//! actually differs in element physics and DTO shapes (`decode.rs`/
//! `decode3.rs`, `tables.rs`/`tables3.rs`), never in this stepping/
//! recording bookkeeping, so `PlanarSession`/`SpatialSession` are just two
//! instantiations of one implementation rather than two hand-duplicated
//! ones.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, AnalysisError, ArcFailure, ArcStepInfo,
    ConstraintHandler, ConvergenceTest, GroundMotion, Integrator, Mode, RayleighDamping,
    StopReason, TransientAnalysis,
};
use carapace_core::model::{
    Domain, Element, Element3, ElementOps, LoadPatternId, Node3Id, NodeId, ELEMENT_DOF, NDF,
    PLANAR_NDIM, SPATIAL_ELEMENT_DOF, SPATIAL_NDF, SPATIAL_NDIM,
};
use serde::Serialize;
use slotmap::Key;

use super::sequence::FiberResponseKind;

/// Shared by `record_sample`'s `Static`/`Transient` arms: resolves a
/// `ResolvedRecorder::Fiber` against whatever `element`'s current fiber
/// state is, or `None` if that element has no fiber section at all
/// (`Domain::element_fiber_responses`'s doc comment) or `point`/`fiber`
/// is out of range for it — the same graceful-skip handling `ModeShape`'s
/// out-of-range `mode` gets, since neither can be checked before the
/// analysis actually runs.
fn fiber_value<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>(
    domain: &Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
    element: E::Id,
    point: u32,
    fiber: u32,
    response: FiberResponseKind,
) -> Option<f64>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    let responses = domain.element_fiber_responses(element)?;
    let &(strain, stress) = responses.get(point as usize)?.get(fiber as usize)?;
    Some(match response {
        FiberResponseKind::Strain => strain,
        FiberResponseKind::Stress => stress,
    })
}

/// `AnalysisError`'s fields, restated so `advance`'s result doesn't need to
/// name `carapace_core`'s error type directly — kept in the same tagged-
/// variant style (implementation-plan.md §2.8).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum AnalysisErrorDetail {
    FailedToConverge {
        step: usize,
    },
    SingularSystem,
    InvalidConstraint,
    InvalidModeCount {
        requested: usize,
        free_dofs: usize,
    },
    /// A configuration value is out of range or unsupported in combination
    /// (`field` is its wire spelling). Decode rejects these up front; this
    /// covers anything only detectable once the model is known (e.g. a
    /// missing rotation scale for a model with rotational DOFs).
    InvalidOption {
        field: &'static str,
    },
    /// Arc length found an unfrozen `path` load series.
    UnsupportedLoadSeries,
    /// Arc length has no reference load to continue.
    ZeroLoadSensitivity,
    /// Arc length must start from equilibrium; `measure` is the force-test
    /// measure at entry.
    InitialStateNotInEquilibrium {
        measure: f64,
    },
    /// The first arc-length tangent is singular and no seed was given.
    MissingSeedDirection,
    /// Every arc-length attempt for `step` failed.
    CutbacksExhausted {
        step: usize,
        attempts: usize,
        radius: f64,
        last_failure: ArcFailureDetail,
    },
    /// An arc-length stop criterion already ended the stage.
    ContinuationComplete,
}

/// `core::ArcFailure`, restated for the wire.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum ArcFailureDetail {
    NotConverged,
    Stagnated,
    NoDescent,
    NonFinite,
    SingularSystem,
    InaccurateSolve,
    BranchOrientation,
    DegenerateConstraint,
}

impl From<ArcFailure> for ArcFailureDetail {
    fn from(failure: ArcFailure) -> Self {
        match failure {
            ArcFailure::NotConverged => Self::NotConverged,
            ArcFailure::Stagnated => Self::Stagnated,
            ArcFailure::NoDescent => Self::NoDescent,
            ArcFailure::NonFinite => Self::NonFinite,
            ArcFailure::SingularSystem => Self::SingularSystem,
            ArcFailure::InaccurateSolve => Self::InaccurateSolve,
            ArcFailure::BranchOrientation => Self::BranchOrientation,
            ArcFailure::DegenerateConstraint => Self::DegenerateConstraint,
        }
    }
}

/// `core::ArcStepInfo` for the last accepted arc-length step of an
/// `advance()` call (`StepOutcome::continuation`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ContinuationDetail {
    pub radius: f64,
    pub next_radius: f64,
    pub retries: usize,
    pub corrector_solves: usize,
    pub factorizations: usize,
    pub force_measure: f64,
    pub arc_error: f64,
    pub chord_length: f64,
    pub det_sign: Option<i8>,
    pub bifurcation_suspected: bool,
    pub displacement_scale: f64,
    pub rotation_scale: Option<f64>,
    pub stop: Option<ContinuationStopDetail>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ContinuationStopDetail {
    pub reason: StopReasonDetail,
    pub landed_exactly: bool,
    pub overshoot: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum StopReasonDetail {
    DisplacementTarget,
    LoadFactorTarget,
    LoadFactorZeroCrossing,
    ChordLength,
    StepCount,
}

impl From<ArcStepInfo> for ContinuationDetail {
    fn from(info: ArcStepInfo) -> Self {
        Self {
            radius: info.radius,
            next_radius: info.next_radius,
            retries: info.retries,
            corrector_solves: info.corrector_solves,
            factorizations: info.factorizations,
            force_measure: info.force_measure,
            arc_error: info.arc_error,
            chord_length: info.chord_length,
            det_sign: info.det_sign,
            bifurcation_suspected: info.bifurcation_suspected,
            displacement_scale: info.displacement_scale,
            rotation_scale: info.rotation_scale,
            stop: info.stop.map(|stop| ContinuationStopDetail {
                reason: match stop.reason {
                    StopReason::DisplacementTarget => StopReasonDetail::DisplacementTarget,
                    StopReason::LoadFactorTarget => StopReasonDetail::LoadFactorTarget,
                    StopReason::LoadFactorZeroCrossing => StopReasonDetail::LoadFactorZeroCrossing,
                    StopReason::ChordLength => StopReasonDetail::ChordLength,
                    StopReason::StepCount => StopReasonDetail::StepCount,
                },
                landed_exactly: stop.landed_exactly,
                overshoot: stop.overshoot,
            }),
        }
    }
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
            AnalysisError::InvalidOption { field } => AnalysisErrorDetail::InvalidOption { field },
            AnalysisError::UnsupportedLoadSeries => AnalysisErrorDetail::UnsupportedLoadSeries,
            AnalysisError::ZeroLoadSensitivity => AnalysisErrorDetail::ZeroLoadSensitivity,
            AnalysisError::InitialStateNotInEquilibrium { measure } => {
                AnalysisErrorDetail::InitialStateNotInEquilibrium { measure }
            }
            AnalysisError::MissingSeedDirection => AnalysisErrorDetail::MissingSeedDirection,
            AnalysisError::CutbacksExhausted {
                step,
                attempts,
                radius,
                last_failure,
            } => AnalysisErrorDetail::CutbacksExhausted {
                step,
                attempts,
                radius,
                last_failure: last_failure.into(),
            },
            AnalysisError::ContinuationComplete => AnalysisErrorDetail::ContinuationComplete,
        }
    }
}

/// A `RecorderSpec`/`RecorderSpec3` resolved against a decoded `Domain` —
/// real node/element ID handles instead of wire-format table indices.
/// Generic over the profile's ID types so one definition serves both
/// `PlanarSession` (`NodeId`/`ElementId`) and `SpatialSession`
/// (`Node3Id`/`Element3Id`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ResolvedRecorder<NId, EId> {
    NodeDisp {
        node: NId,
        dof: u8,
    },
    /// Only meaningful during a `Transient` stage (`RecorderSpec::NodeVel`'s
    /// doc comment) — a `Static`/`Modal` stage's `record_sample` skips it.
    NodeVel {
        node: NId,
        dof: u8,
    },
    NodeAccel {
        node: NId,
        dof: u8,
    },
    ElementForce {
        element: EId,
        component: u8,
    },
    /// See `RecorderSpec::ElementLoad`; skipped during `Modal`.
    ElementLoad {
        element: EId,
        component: u8,
    },
    /// Only meaningful during a `Modal` stage (`RecorderSpec::ModeShape`'s
    /// doc comment).
    ModeShape {
        mode: u32,
        node: NId,
        dof: u8,
    },
    /// Valid during `Static`/`Transient` (`RecorderSpec::Reaction`'s doc
    /// comment); skipped during `Modal`.
    Reaction {
        node: NId,
        dof: u8,
    },
    /// Valid during `Static`/`Transient` (`RecorderSpec::Fiber`'s doc
    /// comment); skipped during `Modal`.
    Fiber {
        element: EId,
        point: u32,
        fiber: u32,
        response: super::sequence::FiberResponseKind,
    },
}

pub(super) struct CompiledStage<NId, EId, Load> {
    pub id: String,
    pub steps: u32,
    pub kind: CompiledStageKind<NId>,
    /// Registered on the `Domain` only when this stage starts — see
    /// `tables::NodalLoadTable::stage`.
    pub pending_nodal_loads: Vec<(LoadPatternId, NId, usize, f64)>,
    pub pending_element_loads: Vec<(LoadPatternId, EId, Load)>,
}

pub(super) enum CompiledStageKind<NId> {
    Static {
        integrator: Integrator<NId>,
        algorithm: Algorithm,
        convergence: ConvergenceTest,
        hold_patterns_after: Vec<LoadPatternId>,
    },
    /// A single eigensolve, not an iterative loop — `CompiledStage::steps`
    /// is always `1` for this kind (`decode.rs`'s `compile_stages`).
    Modal { num_modes: usize },
    /// Swaps the domain for its as-built copy (`ModelSession::pristine`); `steps` is `0`.
    Reset,
    Transient {
        algorithm: Algorithm,
        convergence: ConvergenceTest,
        damping: RayleighDamping,
        dt: f64,
        ground_motions: Vec<GroundMotion>,
    },
}

enum StageRunner<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    Static {
        analysis: Analysis<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        steps_remaining: u32,
    },
    /// `modal_analysis` already ran (in `start_current_stage`) by the time
    /// this variant exists — there's no `core` type to iteratively step
    /// the way `Analysis`/`TransientAnalysis` are, so this just holds the
    /// domain (handed back to `finish_current_stage`, exactly like
    /// `Analysis::into_domain` does for the other two kinds) alongside the
    /// already-computed modes.
    Modal {
        domain: Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        modes: Vec<Mode>,
        steps_remaining: u32,
    },
    Transient {
        analysis: TransientAnalysis<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        steps_remaining: u32,
    },
    /// The restored domain, handed to the next stage by `finish_current_stage`.
    Reset {
        domain: Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        steps_remaining: u32,
    },
}

/// `advance`'s result — pysees-handoff.md's `{ done, stageComplete,
/// stepsTaken, progressSnapshot, recorderBatch? }`, minus the parts that
/// are the worker/JS boundary's job (`progressSnapshot`'s throttling,
/// `recorderBatch`'s `response_blocks` byte layout): `load_factor` here is
/// the raw signal those would be built from — for a `Static` stage the
/// integrator's load factor, for a `Transient` stage the elapsed time
/// (`TransientAnalysis::time`), and for a `Modal` stage always `0.0` (a
/// single eigensolve has no comparable incremental progress scalar; read
/// `Mode::frequency` from a `ModeShape` recorder's batch instead). `error`,
/// once set, is sticky — later stages are not attempted, matching "a
/// stage's `AnalysisError` stops the sequence".
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct StepOutcome {
    pub done: bool,
    pub stage_complete: bool,
    pub steps_taken: u32,
    pub load_factor: f64,
    pub error: Option<AnalysisErrorDetail>,
    /// Diagnostics of the last accepted arc-length step this call took
    /// (radius, retries, solves, chord length, determinant sign, stop
    /// criterion), or `None` when it took none. With arc length the load
    /// factor can decrease and repeat, so samples are ordered by their
    /// sample index, never by load factor.
    pub continuation: Option<ContinuationDetail>,
    /// Only the samples *this* `advance()` call produced, one entry per recorder that recorded
    /// at least one sample this call (every recorder samples every step today, so in practice
    /// this is either empty — no step taken — or has one entry per recorder). This is the
    /// results-storage plan's (docs/results-storage-indexeddb.md) `recorderBatch`, restated in
    /// field names that match pysees's frozen `ResultBlock` protocol type
    /// (`src/app/types/resultsStorage.ts`); the caller assigns `blockIndex` and packs `data`,
    /// since neither concept exists on this side of the wasm boundary.
    pub recorder_batches: Vec<RecorderBatch>,
}

/// One computed mode of a finished `Modal` stage. `shape` is node-major over every node in the
/// input's node-table order (`shape[node * ndf + dof]`, fixed DOFs `0.0`), M-normalized and
/// sign-fixed so its largest-magnitude entry is positive. `participation[d]` is the modal
/// participation factor `phi^T M r_d` for global translation direction `d`, and `mass_ratio[d]`
/// the effective modal mass `participation^2` as a fraction of `ModalStageResult::total_mass[d]`.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ModeResult {
    /// Natural circular frequency, rad/time.
    pub frequency: f64,
    pub shape: Vec<f64>,
    pub participation: Vec<f64>,
    pub mass_ratio: Vec<f64>,
}

/// Everything a finished `Modal` stage computed.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ModalStageResult {
    pub stage_index: usize,
    pub stage_id: String,
    /// DOFs per node (the stride of `ModeResult::shape`).
    pub ndf: usize,
    pub modes: Vec<ModeResult>,
    /// `r_d^T M r_d` per translation direction `d`: the mass the mass ratios are relative to.
    pub total_mass: Vec<f64>,
}

/// `Session::modal_results`'s payload: every `Modal` stage that has finished so far, in order.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ModalResultsReport {
    pub stages: Vec<ModalStageResult>,
}

/// One recorder's new samples from a single `advance()` call. `first_sample` plus
/// `samples.len()` gives the sample range this batch covers, so a caller can retry a failed
/// persist without renumbering — the "each `advance()` batch deterministic and identifies its
/// first sample/count" contract from results-storage-indexeddb.md's "Failure and cancellation
/// contract". Scalar-only (`(pseudo_time, value)` pairs, one channel) until a vector recorder
/// needs more.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct RecorderBatch {
    pub recorder_index: usize,
    pub stage_index: usize,
    pub first_sample: u32,
    pub samples: Vec<(f64, f64)>,
}

/// One opaque handle exposed to the caller, chosen once from `space` at
/// decode time.
#[derive(Debug)]
pub enum Session {
    Planar(PlanarSession),
    Spatial(SpatialSession),
}

impl Session {
    pub fn advance(&mut self, step_budget: u32) -> StepOutcome {
        match self {
            Session::Planar(session) => session.advance(step_budget),
            Session::Spatial(session) => session.advance(step_budget),
        }
    }

    /// Every `Modal` stage finished so far (see `ModalStageResult`).
    pub fn modal_results(&self) -> ModalResultsReport {
        match self {
            Session::Planar(session) => session.modal_results(),
            Session::Spatial(session) => session.modal_results(),
        }
    }

    /// The `AnalysisSequence` stage id `advance` is currently in (or, once
    /// `done`, the id of whichever stage stopped it) — for surfacing which
    /// stage an `error` belongs to, matching pysees-handoff.md's "the run
    /// is marked failed with the error's structured detail".
    pub fn current_stage_id(&self) -> Option<&str> {
        match self {
            Session::Planar(session) => session.current_stage_id(),
            Session::Spatial(session) => session.current_stage_id(),
        }
    }
}

pub type PlanarSession = ModelSession<PLANAR_NDIM, NDF, ELEMENT_DOF, NodeId, Element>;
pub type SpatialSession =
    ModelSession<SPATIAL_NDIM, SPATIAL_NDF, SPATIAL_ELEMENT_DOF, Node3Id, Element3>;

pub struct ModelSession<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    domain: Option<Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>>,
    // Node handles in input node-table order, for `ModeResult::shape`.
    node_ids: Vec<NId>,
    modal_results: Vec<ModalStageResult>,
    // The domain as built, kept only when a `Reset` stage exists.
    pristine: Option<Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>>,
    // Every `holdPatternsAfter` freeze applied so far, in order, so a `Reset` can re-apply them.
    holds: Vec<(LoadPatternId, f64)>,
    record_initial: bool,
    stages: Vec<CompiledStage<NId, E::Id, E::Load>>,
    current_stage: usize,
    runner: Option<StageRunner<NDIM, NDOF, ELEMENT_DOF, NId, E>>,
    load_factor: f64,
    error: Option<AnalysisErrorDetail>,
    recorders: Vec<ResolvedRecorder<NId, E::Id>>,
    // Per-advance only: the last arc-length step's diagnostics.
    continuation: Option<ContinuationDetail>,
    // Per-advance only: cleared at the top of every `advance()` call, drained into that call's
    // `StepOutcome::recorder_batches` at the end. Nothing here survives across calls except the
    // running counts below — this is results-storage-indexeddb.md's "Bound memory" step:
    // `recorder_history`'s old whole-run `Vec<Vec<(f64, f64)>>` is gone.
    current_batch: Vec<Vec<(f64, f64)>>,
    // Total samples emitted so far per recorder, so the next batch's `first_sample` is correct
    // without retaining the samples themselves.
    recorder_sample_counts: Vec<u32>,
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E> std::fmt::Debug
    for ModelSession<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    // `Analysis` (inside `StageRunner`, held via `runner`/`domain`) doesn't
    // derive `Debug`, so this reports the fields useful for a test/error
    // message rather than the full solver state.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ModelSession")
            .field("current_stage", &self.current_stage)
            .field("load_factor", &self.load_factor)
            .field("error", &self.error)
            .finish()
    }
}

impl<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>
    ModelSession<NDIM, NDOF, ELEMENT_DOF, NId, E>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId> + Clone,
    E::Load: Clone,
{
    pub(super) fn new(
        domain: Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        node_ids: Vec<NId>,
        stages: Vec<CompiledStage<NId, E::Id, E::Load>>,
        recorders: Vec<ResolvedRecorder<NId, E::Id>>,
        record_initial: bool,
    ) -> Self {
        let pristine = stages
            .iter()
            .any(|stage| matches!(stage.kind, CompiledStageKind::Reset))
            .then(|| domain.clone());
        let current_batch = vec![Vec::new(); recorders.len()];
        let recorder_sample_counts = vec![0; recorders.len()];
        Self {
            domain: Some(domain),
            node_ids,
            modal_results: Vec::new(),
            pristine,
            holds: Vec::new(),
            record_initial,
            stages,
            current_stage: 0,
            runner: None,
            load_factor: 0.0,
            error: None,
            recorders,
            continuation: None,
            current_batch,
            recorder_sample_counts,
        }
    }

    pub fn modal_results(&self) -> ModalResultsReport {
        ModalResultsReport {
            stages: self.modal_results.clone(),
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
        self.continuation = None;

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
                continuation: None,
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
            if self.steps_remaining() == 0 {
                stage_complete = true;
                break;
            }
            match self.step_once() {
                Ok(()) => {
                    steps_taken += 1;
                    if self.steps_remaining() == 0 {
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
            continuation: self.continuation,
            recorder_batches,
        }
    }

    fn start_current_stage(&mut self) {
        let stage = &self.stages[self.current_stage];
        let mut domain = self.domain.take().expect("domain held between stages");
        if matches!(stage.kind, CompiledStageKind::Reset) {
            // As-built domain plus every load an earlier stage had registered (loads are
            // registered when their stage starts, so the pristine copy has none), and the freezes
            // (`loadConst`) already applied: like OpenSees' `reset`, this reverts the state but not
            // which patterns are held constant, so a held pattern is applied in full from the
            // first step after the reset.
            domain = self.pristine.clone().expect("kept for Reset stages");
            for earlier in &self.stages[..self.current_stage] {
                for &(pattern, node, dof, value) in &earlier.pending_nodal_loads {
                    domain.add_nodal_load(pattern, node, dof, value);
                }
                for (pattern, element, load) in &earlier.pending_element_loads {
                    domain.add_element_load(*pattern, *element, load.clone());
                }
            }
            for &(pattern, factor) in &self.holds {
                domain.hold_pattern_constant(pattern, factor);
            }
            self.runner = Some(StageRunner::Reset {
                domain,
                steps_remaining: 0,
            });
            self.load_factor = 0.0;
            return;
        }
        for &(pattern, node, dof, value) in &stage.pending_nodal_loads {
            domain.add_nodal_load(pattern, node, dof, value);
        }
        for (pattern, element, load) in &stage.pending_element_loads {
            domain.add_element_load(*pattern, *element, load.clone());
        }
        // `Transformation` resolves `equal_dof`/`rigid_diaphragm` ties;
        // `Plain` is a strictly cheaper no-op for a domain with none (and
        // `AnalysisBuilder::build` would reject `Plain` outright against a
        // domain that does have them — see its own doc comment) — so this
        // just asks the domain which one it actually needs rather than the
        // caller ever choosing.
        let constraint_handler = if domain.has_mp_constraints() {
            ConstraintHandler::Transformation
        } else {
            ConstraintHandler::Plain
        };
        match &stage.kind {
            CompiledStageKind::Static {
                integrator,
                algorithm,
                convergence,
                ..
            } => {
                let analysis = AnalysisBuilder::new()
                    .constraint_handler(constraint_handler)
                    .integrator(integrator.clone())
                    .algorithm(*algorithm)
                    .test(*convergence)
                    .build(domain);
                self.runner = Some(StageRunner::Static {
                    analysis,
                    steps_remaining: stage.steps,
                });
            }
            // The eigensolve runs immediately, here, rather than lazily on
            // the first `step_once` call — there's no cheaper "start" step
            // to defer it to, unlike `Static`/`Transient`'s first Newton/
            // Newmark solve, which naturally happens on the first `step()`.
            CompiledStageKind::Modal { num_modes } => match modal_analysis(&mut domain, *num_modes)
            {
                Ok(modes) => {
                    self.runner = Some(StageRunner::Modal {
                        domain,
                        modes,
                        steps_remaining: stage.steps,
                    });
                }
                Err(error) => self.error = Some(error.into()),
            },
            CompiledStageKind::Reset => unreachable!("handled above"),
            CompiledStageKind::Transient {
                algorithm,
                convergence,
                damping,
                dt,
                ground_motions,
            } => match TransientAnalysis::new(domain, *damping, *dt) {
                Ok(mut analysis) => {
                    analysis = analysis.with_algorithm(*algorithm, *convergence);
                    for motion in ground_motions.iter().cloned() {
                        analysis = analysis.with_ground_motion(motion);
                    }
                    self.runner = Some(StageRunner::Transient {
                        analysis,
                        steps_remaining: stage.steps,
                    });
                }
                Err(error) => self.error = Some(error.into()),
            },
        }
        self.load_factor = 0.0;
        if self.record_initial
            && matches!(
                self.stages[self.current_stage].kind,
                CompiledStageKind::Static { .. } | CompiledStageKind::Transient { .. }
            )
        {
            self.record_sample();
        }
    }

    fn steps_remaining(&self) -> u32 {
        match self.runner.as_ref().expect("stage started above") {
            StageRunner::Static {
                steps_remaining, ..
            }
            | StageRunner::Modal {
                steps_remaining, ..
            }
            | StageRunner::Reset {
                steps_remaining, ..
            }
            | StageRunner::Transient {
                steps_remaining, ..
            } => *steps_remaining,
        }
    }

    /// Advances whichever kind of stage is currently running by one unit —
    /// a Newton step (`Static`), a Newmark step (`Transient`), or (for
    /// `Modal`, whose eigensolve already ran in `start_current_stage`) just
    /// consuming its one always-available "step" — then records a sample
    /// for every recorder this runner kind supports.
    fn step_once(&mut self) -> Result<(), AnalysisError> {
        match self.runner.as_mut().expect("stage started above") {
            StageRunner::Static {
                analysis,
                steps_remaining,
            } => {
                let result = analysis.step()?;
                self.load_factor = result.load_factor;
                *steps_remaining -= 1;
                if let Some(info) = result.arc {
                    self.continuation = Some(info.into());
                    // A met stop criterion ends the stage early.
                    if info.stop.is_some() {
                        *steps_remaining = 0;
                    }
                }
            }
            StageRunner::Transient {
                analysis,
                steps_remaining,
            } => {
                let result = analysis.step()?;
                self.load_factor = result.time;
                *steps_remaining -= 1;
            }
            StageRunner::Modal {
                steps_remaining, ..
            } => {
                *steps_remaining -= 1;
            }
            StageRunner::Reset { .. } => {}
        }
        self.record_sample();
        Ok(())
    }

    /// `Analysis::into_domain`/fresh-`AnalysisBuilder` multi-phase pattern
    /// (pysees-handoff.md / core/tests/m8_force_beam_column.rs's two-phase
    /// tests): freeze this stage's held patterns (a `Static`-only concept —
    /// `Modal`/`Transient` stages hold none) at its final load factor, then
    /// hand the same `Domain` (committed material state and all) to the
    /// next stage's fresh `Analysis`/`TransientAnalysis`/`modal_analysis`
    /// call.
    fn finish_current_stage(&mut self) {
        let hold_patterns_after = match &self.stages[self.current_stage].kind {
            CompiledStageKind::Static {
                hold_patterns_after,
                ..
            } => hold_patterns_after.clone(),
            CompiledStageKind::Modal { .. }
            | CompiledStageKind::Transient { .. }
            | CompiledStageKind::Reset => Vec::new(),
        };
        let load_factor = self.load_factor;

        let mut domain = match self.runner.take().expect("stage running") {
            StageRunner::Static { analysis, .. } => analysis.into_domain(),
            StageRunner::Modal { mut domain, modes, .. } => {
                let result = self.modal_stage_result(&mut domain, &modes);
                self.modal_results.push(result);
                domain
            }
            StageRunner::Transient { analysis, .. } => analysis.into_domain(),
            StageRunner::Reset { domain, .. } => domain,
        };
        for pattern in hold_patterns_after {
            domain.hold_pattern_constant(pattern, load_factor);
            self.holds.push((pattern, load_factor));
        }
        self.domain = Some(domain);
        self.current_stage += 1;
    }

    fn modal_stage_result(
        &self,
        domain: &mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
        modes: &[Mode],
    ) -> ModalStageResult {
        let mass = domain.assemble_mass_diagonal();
        // Influence vector per global translation direction: 1 on every free equation a node's
        // translation DOF maps to. (A rigid-diaphragm slave's translation is not exactly 1 on
        // its retained equation, so ratios in models using one are approximate.)
        let influence: Vec<Vec<f64>> = (0..NDIM)
            .map(|direction| {
                let mut r = vec![0.0; mass.len()];
                for &node in &self.node_ids {
                    if let Some(eq) = domain.equation_of(node, direction) {
                        r[eq] = 1.0;
                    }
                }
                r
            })
            .collect();
        let total_mass: Vec<f64> = influence
            .iter()
            .map(|r| r.iter().zip(mass.iter()).map(|(r, m)| r * r * m).sum())
            .collect();
        let modes = modes
            .iter()
            .map(|mode| {
                let sign = mode
                    .shape
                    .iter()
                    .copied()
                    .fold(0.0_f64, |best, v| if v.abs() > best.abs() { v } else { best })
                    .signum();
                let sign = if sign == 0.0 { 1.0 } else { sign };
                let mut shape = Vec::with_capacity(self.node_ids.len() * NDOF);
                for &node in &self.node_ids {
                    for dof in 0..NDOF {
                        shape.push(
                            domain
                                .equation_of(node, dof)
                                .map_or(0.0, |eq| sign * mode.shape[eq]),
                        );
                    }
                }
                let participation: Vec<f64> = influence
                    .iter()
                    .map(|r| {
                        sign * (0..mass.len()).map(|i| mode.shape[i] * mass[i] * r[i]).sum::<f64>()
                    })
                    .collect();
                let mass_ratio = participation
                    .iter()
                    .zip(&total_mass)
                    .map(|(g, total)| if *total > 0.0 { g * g / total } else { 0.0 })
                    .collect();
                ModeResult {
                    frequency: mode.frequency,
                    shape,
                    participation,
                    mass_ratio,
                }
            })
            .collect();
        ModalStageResult {
            stage_index: self.current_stage,
            stage_id: self.stages[self.current_stage].id.clone(),
            ndf: NDOF,
            modes,
            total_mass,
        }
    }

    fn record_sample(&mut self) {
        let Some(runner) = &self.runner else {
            return;
        };
        match runner {
            StageRunner::Static { analysis, .. } => {
                let progress = self.load_factor;
                for (recorder, batch) in self.recorders.iter().zip(self.current_batch.iter_mut()) {
                    let value = match *recorder {
                        ResolvedRecorder::NodeDisp { node, dof } => {
                            Some(analysis.domain().node(node).displacement[dof as usize])
                        }
                        ResolvedRecorder::ElementForce { element, component } => Some(
                            analysis.domain().element_end_force(element, progress)
                                [component as usize],
                        ),
                        ResolvedRecorder::ElementLoad { element, component } => {
                            Some(analysis.domain().element_load_component(
                                element,
                                progress,
                                component as usize,
                            ))
                        }
                        ResolvedRecorder::Reaction { node, dof } => {
                            Some(analysis.domain().reaction(node, dof as usize, progress))
                        }
                        ResolvedRecorder::Fiber {
                            element,
                            point,
                            fiber,
                            response,
                        } => fiber_value(analysis.domain(), element, point, fiber, response),
                        ResolvedRecorder::NodeVel { .. }
                        | ResolvedRecorder::NodeAccel { .. }
                        | ResolvedRecorder::ModeShape { .. } => None,
                    };
                    if let Some(value) = value {
                        batch.push((progress, value));
                    }
                }
            }
            StageRunner::Transient { analysis, .. } => {
                let progress = self.load_factor; // == analysis.time(), set by step_once
                for (recorder, batch) in self.recorders.iter().zip(self.current_batch.iter_mut()) {
                    let value = match *recorder {
                        ResolvedRecorder::NodeDisp { node, dof } => {
                            Some(analysis.domain().node(node).displacement[dof as usize])
                        }
                        ResolvedRecorder::NodeVel { node, dof } => {
                            Some(analysis.domain().node(node).velocity[dof as usize])
                        }
                        ResolvedRecorder::NodeAccel { node, dof } => {
                            Some(analysis.domain().node(node).acceleration[dof as usize])
                        }
                        ResolvedRecorder::ElementForce { element, component } => Some(
                            analysis.domain().element_end_force(element, progress)
                                [component as usize],
                        ),
                        ResolvedRecorder::ElementLoad { element, component } => {
                            Some(analysis.domain().element_load_component(
                                element,
                                progress,
                                component as usize,
                            ))
                        }
                        ResolvedRecorder::Reaction { node, dof } => {
                            Some(analysis.domain().reaction(node, dof as usize, progress))
                        }
                        ResolvedRecorder::Fiber {
                            element,
                            point,
                            fiber,
                            response,
                        } => fiber_value(analysis.domain(), element, point, fiber, response),
                        ResolvedRecorder::ModeShape { .. } => None,
                    };
                    if let Some(value) = value {
                        batch.push((progress, value));
                    }
                }
            }
            StageRunner::Reset { .. } => {}
            StageRunner::Modal { domain, modes, .. } => {
                for (recorder, batch) in self.recorders.iter().zip(self.current_batch.iter_mut()) {
                    let ResolvedRecorder::ModeShape { mode, node, dof } = *recorder else {
                        continue;
                    };
                    let Some(computed) = modes.get(mode as usize) else {
                        continue;
                    };
                    // A fixed DOF has no free-DOF equation number and so no
                    // entry in `Mode::shape` — its mode-shape component is
                    // trivially zero (a fixed DOF can't participate in any
                    // mode).
                    let value = domain
                        .equation_of(node, dof as usize)
                        .map_or(0.0, |eq| computed.shape[eq]);
                    batch.push((computed.frequency, value));
                }
            }
        }
    }
}
