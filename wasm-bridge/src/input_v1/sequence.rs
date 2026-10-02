//! Wire-format `AnalysisSequence` — mirrors pysees-handoff.md's
//! `AnalysisStage` discriminated union: `Static` (pushover-style), `Modal`
//! (eigenvalue), and `Transient` (Newmark time-history, optionally driven
//! by `GroundMotionSpec`).

use serde::{Deserialize, Serialize};

use super::tables::{ElementKind, TimeSeriesSpec};
use super::tables3::ElementKind3;

// `ArcLength` is much larger than the other variants; these are decoded a
// handful of times per input, so boxing would only complicate the wire type.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum IntegratorSpec {
    LoadControl {
        increment: f64,
    },
    /// `node`/`dof` index into `NodeTable`; the compiler is responsible for
    /// having registered a reference load with nonzero sensitivity on that
    /// DOF (pysees-handoff.md's "valid reference-load sensitivity for
    /// displacement control") — decode does not re-derive that check.
    DisplacementControl {
        node: u32,
        dof: u8,
        increment: f64,
    },
    /// Arc-length continuation (`core::ArcLength`, docs/arclength.md):
    /// traces the equilibrium path through load limits and snap-back. The
    /// stage's `steps` is the step cap; `stop` can end the stage earlier.
    /// Requires a Newton algorithm with the current tangent and no line
    /// search, and a `normUnbalance` or `combined` convergence test (a
    /// `combined` test is the default when `convergence` is omitted).
    ArcLength(ArcLengthSpec),
}

/// Wire form of `core::ArcLength`. Omitted `minRadius`/`maxRadius` default
/// to `initialRadius` (a fixed-radius run).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ArcLengthSpec {
    pub initial_radius: f64,
    #[tsify(optional)]
    pub min_radius: Option<f64>,
    #[tsify(optional)]
    pub max_radius: Option<f64>,
    #[tsify(optional)]
    pub target_iterations: Option<u32>,
    #[tsify(optional)]
    pub max_retries: Option<u32>,
    #[serde(default)]
    #[tsify(optional)]
    pub direction: ArcDirectionSpec,
    pub scales: ArcScalesSpec,
    #[serde(default)]
    #[tsify(optional)]
    pub predictor: ArcPredictorSpec,
    #[tsify(optional)]
    pub seed: Option<ArcSeedSpec>,
    #[tsify(optional)]
    pub arc_tolerance: Option<f64>,
    #[tsify(optional)]
    pub correction_tolerance: Option<f64>,
    #[tsify(optional)]
    pub backtracking: Option<BacktrackingSpec>,
    #[tsify(optional)]
    pub stop: Option<ArcStopSpec>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum ArcDirectionSpec {
    #[default]
    Increasing,
    Decreasing,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum ArcPredictorSpec {
    #[default]
    Secant,
    Tangent,
}

/// The continuation metric's scales (`core::ArcScales`): explicit
/// characteristic translation/rotation/load-factor sizes, or `auto`
/// (translation/rotation scales derived from the first elastic tangent).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    deny_unknown_fields
)]
pub enum ArcScalesSpec {
    Explicit {
        displacement: f64,
        #[tsify(optional)]
        rotation: Option<f64>,
        load: f64,
    },
    Auto {
        load: f64,
    },
}

/// One component of a seed direction: `node`/`dof` index into `NodeTable`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ArcSeedComponentSpec {
    pub node: u32,
    pub dof: u8,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ArcSeedSpec {
    pub components: Vec<ArcSeedComponentSpec>,
    #[serde(default)]
    #[tsify(optional)]
    pub load: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct BacktrackingSpec {
    pub armijo: f64,
    pub min_step: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct DisplacementTargetSpec {
    pub node: u32,
    pub dof: u8,
    pub value: f64,
    #[serde(default)]
    #[tsify(optional)]
    pub exact: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct LoadFactorTargetSpec {
    pub value: f64,
    #[serde(default)]
    #[tsify(optional)]
    pub exact: bool,
}

/// Optional criteria that end an arc-length stage before its step cap.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct ArcStopSpec {
    #[tsify(optional)]
    pub displacement: Option<DisplacementTargetSpec>,
    #[tsify(optional)]
    pub load_factor: Option<LoadFactorTargetSpec>,
    #[serde(default)]
    #[tsify(optional)]
    pub load_factor_zero_crossing: bool,
    #[tsify(optional)]
    pub max_chord_length: Option<f64>,
}

/// Accepts the original bare strings and configurable algorithm objects.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(untagged)]
pub enum AlgorithmSpec {
    Legacy(LegacyAlgorithmSpec),
    Config(AlgorithmConfigSpec),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum LegacyAlgorithmSpec {
    Linear,
    NewtonRaphson,
}

impl Default for AlgorithmSpec {
    fn default() -> Self {
        Self::Linear
    }
}

// Keep the existing Rust constructors as well as the existing JS strings.
#[allow(non_upper_case_globals)]
impl AlgorithmSpec {
    pub const Linear: Self = Self::Legacy(LegacyAlgorithmSpec::Linear);
    pub const NewtonRaphson: Self = Self::Legacy(LegacyAlgorithmSpec::NewtonRaphson);
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum TangentStrategySpec {
    #[default]
    Current,
    ReuseAtStepStart,
    Initial,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    deny_unknown_fields
)]
pub enum LineSearchSpec {
    Bisection {
        tol: f64,
        max_iter: u32,
        max_eta: f64,
    },
    RegulaFalsi {
        tol: f64,
        max_iter: u32,
        max_eta: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    deny_unknown_fields
)]
pub enum AlgorithmConfigSpec {
    Linear,
    Newton {
        #[serde(default)]
        tangent: TangentStrategySpec,
        #[tsify(optional)]
        line_search: Option<LineSearchSpec>,
    },
    KrylovNewton {
        tangent: TangentStrategySpec,
        max_dimension: u32,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ConvergenceSpec {
    NormUnbalance {
        tol: f64,
        max_iter: u32,
    },
    NormDispIncr {
        tol: f64,
        max_iter: u32,
    },
    EnergyIncr {
        tol: f64,
        max_iter: u32,
    },
    /// `core::ConvergenceTest::Combined`: per-equation force equilibrium
    /// (`|r_i| <= forceTol|momentTol + relativeTol * F_ref_i`), plus an
    /// optional displacement-correction bound; all must pass. `momentTol`
    /// defaults to `forceTol`, `relativeTol` to `1e-6`.
    Combined {
        force_tol: f64,
        #[tsify(optional)]
        moment_tol: Option<f64>,
        #[tsify(optional)]
        relative_tol: Option<f64>,
        #[tsify(optional)]
        displacement_tol: Option<f64>,
        max_iter: u32,
    },
}

impl ConvergenceSpec {
    /// The decoder's default when a stage omits `convergence` — the
    /// handoff's TS type marks the field optional but does not name a
    /// default, so this one is this decoder's own choice.
    pub const DEFAULT: ConvergenceSpec = ConvergenceSpec::NormUnbalance {
        tol: 1e-6,
        max_iter: 20,
    };

    /// The default for an arc-length stage that omits `convergence`: force
    /// equilibrium per equation (an arc-length step always also has to meet
    /// its constraint tolerance).
    pub const ARC_LENGTH_DEFAULT: ConvergenceSpec = ConvergenceSpec::Combined {
        force_tol: 1e-6,
        moment_tol: None,
        relative_tol: None,
        displacement_tol: None,
        max_iter: 30,
    };
}

/// Rayleigh damping (`C = alpha_m*M + beta_k*K`) — see
/// `core::RayleighDamping`'s doc comment. Use `alpha_m: 0.0, beta_k: 0.0`
/// for undamped (`core::RayleighDamping::NONE`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct DampingSpec {
    pub alpha_m: f64,
    pub beta_k: f64,
}

/// One `core::GroundMotion`: a uniform acceleration time history applied
/// along one global translational axis (`direction`: `0`=x, `1`=y, and for
/// the spatial profile `2`=z — never a rotational DOF). `series` is
/// typically `TimeSeriesSpec::Path` (the accelerogram itself), reusing the
/// same wire type a static `LoadPatternTable` entry uses — an accelerogram
/// is exactly a piecewise-linear series, no separate representation needed
/// (`core::GroundMotion`'s own doc comment).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct GroundMotionSpec {
    pub direction: u8,
    pub series: TimeSeriesSpec,
    pub scale_factor: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
#[allow(clippy::large_enum_variant)] // see `IntegratorSpec`
pub enum StageSpec {
    Static {
        id: String,
        steps: u32,
        integrator: IntegratorSpec,
        algorithm: AlgorithmSpec,
        #[tsify(optional)]
        convergence: Option<ConvergenceSpec>,
        /// Indices into `LoadPatternTable` to freeze (via
        /// `Domain::hold_pattern_constant`) at this stage's final load
        /// factor once it completes.
        hold_patterns_after: Vec<u32>,
    },
    /// A single eigenvalue solve (`core::modal_analysis`) for the domain's
    /// current state — not iterative like `Static`/`Transient`, so this
    /// stage always takes exactly one `advance()` step regardless of the
    /// caller's step budget (`decode.rs`'s `compile_stages` hardcodes
    /// `steps: 1` for it). `modes` is the number of lowest-frequency modes
    /// to compute (`1..=free_dof_count`).
    Modal { id: String, modes: u32 },
    /// Newmark-beta time-history analysis (`core::TransientAnalysis`,
    /// fixed at the unconditionally-stable "average acceleration"
    /// parameters — see its own doc comment). `steps` fixed-size `dt`
    /// increments; `ground_motions` (zero or more) each add a uniform
    /// inertial excitation along one axis, on top of whatever `LoadPattern`
    /// reference loads are already active (e.g. a preceding `Static`
    /// stage's frozen gravity).
    Transient {
        id: String,
        steps: u32,
        dt: f64,
        damping: DampingSpec,
        ground_motions: Vec<GroundMotionSpec>,
        /// Omitted algorithms preserve the original linear Newmark solve.
        #[serde(default)]
        algorithm: AlgorithmSpec,
        #[tsify(optional)]
        convergence: Option<ConvergenceSpec>,
    },
}

impl StageSpec {
    pub fn id(&self) -> &str {
        match self {
            StageSpec::Static { id, .. }
            | StageSpec::Modal { id, .. }
            | StageSpec::Transient { id, .. } => id,
        }
    }
}

/// One recorder: a single scalar channel's history against its stage's load
/// factor. `NodeDisp` was M10's first (and, until now, only) variant
/// (pysees-handoff.md's "selected node displacement/load-factor history for
/// the static pushover"); `ElementForce` (results-storage-indexeddb.md's
/// "several more types of recorders" plan) is the second, sharing the exact
/// same batching/storage machinery — see `session::ResolvedRecorder` and
/// `StepOutcome::recorder_batches`, neither of which needed to change shape
/// to add it. `ElementForce`'s `element_kind`/`element_index` pair mirrors
/// `ElementLoadTable`'s existing disambiguation between per-kind element
/// tables (`ElasticBeamColumnTable`, `ForceBeamColumnTable`, ...) — there is
/// no single flat element table to index into directly, unlike `NodeTable`.
/// Adding a future response kind (velocity, acceleration, ...) is one more
/// variant here plus one more match arm in `PlanarSession::record_sample`,
/// not a new parallel type or a new `Session`/`StepOutcome` field.
/// Which half of a fiber's `(strain, stress)` pair (`core::FiberSection::
/// fiber_responses`'s doc comment) a `RecorderSpec::Fiber`/`RecorderSpec3::
/// Fiber` reads — one scalar channel per recorder, same as every other
/// kind here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum FiberResponseKind {
    Strain,
    Stress,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "response",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum RecorderSpec {
    /// `node`/`dof` index into `NodeTable`.
    NodeDisp { node: u32, dof: u8 },
    /// Only meaningful during a `Transient` stage — `analysis.domain()`'s
    /// node velocity is otherwise always zero (a `Static` stage never
    /// touches it), so this recorder simply contributes no sample outside
    /// one.
    NodeVel { node: u32, dof: u8 },
    /// See `NodeVel`'s doc comment — same restriction, for acceleration.
    NodeAccel { node: u32, dof: u8 },
    /// `element_index` indexes into whichever per-kind table `element_kind`
    /// names (e.g. `ForceBeamColumnTable`), not a flat cross-kind element
    /// list. `component` indexes the element's local nodal force vector
    /// (`ElementOps::local_force` — axial/shear/moment at each end, in the
    /// element's own axis frame), width `2 * dofsPerNode` (6 for this
    /// planar decoder).
    ElementForce {
        element_kind: ElementKind,
        element_index: u32,
        component: u8,
    },
    /// Only meaningful during a `Modal` stage: `mode` indexes that stage's
    /// computed `Vec<Mode>` (ascending frequency, `0` = lowest), `node`/
    /// `dof` locate one entry of that mode's free-DOF-indexed shape vector
    /// (`core::Mode::shape`) via `Domain::equation_of` — a fixed DOF has no
    /// equation number and reads as `0.0` (trivially true: a fixed DOF
    /// can't participate in any mode shape). The recorded "time" component
    /// is the mode's frequency (rad/time), not a pseudo-time or load
    /// factor — there's no more natural scalar to pair a mode's shape
    /// component with.
    ModeShape { mode: u32, node: u32, dof: u8 },
    /// Support/equilibrium reaction at `node`/`dof` (`core::Domain::
    /// reaction`'s doc comment) — most useful at a fixed DOF, but not
    /// restricted to one. Valid during `Static` and `Transient` stages
    /// (paired with the stage's own load factor/time, same as `NodeDisp`);
    /// records nothing during a `Modal` stage, which has no applied-load
    /// state for `reaction` to evaluate against.
    Reaction { node: u32, dof: u8 },
    /// One fiber's strain or stress, at one integration point, of one
    /// `DispBeamColumn`/`ForceBeamColumn` element (`core::Domain::
    /// element_fiber_responses`'s doc comment) — `element_kind`/
    /// `element_index` mirror `ElementForce`'s own disambiguation.
    /// `point`/`fiber` index into that element's own integration-point/
    /// fiber ordering, exactly as originally given to its constructor
    /// (`point` therefore also matches `IntegrationSpec`'s point count,
    /// and `fiber` matches `FiberTable`'s per-section fiber range).
    /// Records nothing for any other element kind (no fiber section to
    /// read) or an out-of-range `point`/`fiber` — the same "unresolved by
    /// runtime state, not a decode-time error" handling `ModeShape` uses,
    /// since neither can be checked before the analysis actually runs.
    Fiber {
        element_kind: ElementKind,
        element_index: u32,
        point: u32,
        fiber: u32,
        quantity: FiberResponseKind,
    },
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct SequenceSpec {
    pub stages: Vec<StageSpec>,
    pub recorders: Vec<RecorderSpec>,
}

/// `RecorderSpec`'s spatial counterpart — same shape, just `element_kind`
/// naming a spatial [`ElementKind3`] and `component` indexing a width-12
/// (`2 * SPATIAL_NDF`) local nodal force vector instead of width-6. `Node`
/// disp is otherwise dimension-agnostic (`dof` just goes up to 5 instead of
/// 2), but a shared enum would need `RecorderSpec::ElementForce` to name a
/// type that's different per profile, so this stays its own enum rather
/// than a generic parameter.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "response",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum RecorderSpec3 {
    NodeDisp {
        node: u32,
        dof: u8,
    },
    NodeVel {
        node: u32,
        dof: u8,
    },
    NodeAccel {
        node: u32,
        dof: u8,
    },
    ElementForce {
        element_kind: ElementKind3,
        element_index: u32,
        component: u8,
    },
    ModeShape {
        mode: u32,
        node: u32,
        dof: u8,
    },
    Reaction {
        node: u32,
        dof: u8,
    },
    Fiber {
        element_kind: ElementKind3,
        element_index: u32,
        point: u32,
        fiber: u32,
        quantity: FiberResponseKind,
    },
}

/// `SequenceSpec`'s spatial counterpart. `stages: Vec<StageSpec>` is reused
/// verbatim — stage/integrator/algorithm/convergence compilation
/// (`decode.rs`'s `compile_stages`) never touches element physics, only
/// node indices and DOF numbers, so nothing about it is planar-specific.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct SequenceSpec3 {
    pub stages: Vec<StageSpec>,
    pub recorders: Vec<RecorderSpec3>,
}
