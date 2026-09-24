//! Wire-format `AnalysisSequence` — mirrors pysees-handoff.md's
//! `AnalysisStage` discriminated union. `Modal`/`Transient` variants exist
//! here from the start (per "Staged decode support, not a staged wire
//! format") even though [`decode`](super::decode::decode) rejects them
//! outright today: the point is that a newer wire producer's sequence
//! decodes structurally and fails with a named diagnostic, not a parse
//! error, when it uses a stage kind this decoder doesn't implement yet.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum AlgorithmSpec {
    Linear,
    NewtonRaphson,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ConvergenceSpec {
    NormUnbalance { tol: f64, max_iter: u32 },
    NormDispIncr { tol: f64, max_iter: u32 },
    EnergyIncr { tol: f64, max_iter: u32 },
}

impl ConvergenceSpec {
    /// The decoder's default when a stage omits `convergence` — the
    /// handoff's TS type marks the field optional but does not name a
    /// default, so this one is this decoder's own choice.
    pub const DEFAULT: ConvergenceSpec = ConvergenceSpec::NormUnbalance {
        tol: 1e-6,
        max_iter: 20,
    };
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum StageSpec {
    Static {
        id: String,
        steps: u32,
        integrator: IntegratorSpec,
        algorithm: AlgorithmSpec,
        convergence: Option<ConvergenceSpec>,
        /// Indices into `LoadPatternTable` to freeze (via
        /// `Domain::hold_pattern_constant`) at this stage's final load
        /// factor once it completes.
        hold_patterns_after: Vec<u32>,
    },
    Modal {
        id: String,
        modes: u32,
    },
    Transient {
        id: String,
    },
}

impl StageSpec {
    pub fn id(&self) -> &str {
        match self {
            StageSpec::Static { id, .. }
            | StageSpec::Modal { id, .. }
            | StageSpec::Transient { id } => id,
        }
    }
}

/// M10's first recorder: a single node/dof's displacement history against
/// its stage's load factor (pysees-handoff.md's "selected node
/// displacement/load-factor history for the static pushover"). Deliberately
/// narrower than the handoff's general `RecorderSpec` sketch (target
/// kind/tags, response kind, component layout, sampling spec) — those
/// belong to the worker/SQLite layer this Rust-side decoder doesn't own.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RecorderSpec {
    /// Index into `NodeTable`.
    pub node: u32,
    pub dof: u8,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SequenceSpec {
    pub stages: Vec<StageSpec>,
    pub recorders: Vec<RecorderSpec>,
}
