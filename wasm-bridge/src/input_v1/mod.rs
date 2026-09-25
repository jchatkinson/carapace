//! `CarapaceInputV1` — the pysees ↔ Carapace wire format and its decoder,
//! per docs/pysees-handoff.md. This module is the Rust-side half of M10
//! (implementation-plan.md): header/table DTOs, the hand-written per-table
//! decoder into a `carapace-core` `Domain`/`Analysis`, and the stepped
//! `Session`/`advance` API. Every type here derives `serde`'s
//! `Serialize`/`Deserialize` so `../boundary.rs` can hand them across
//! `wasm_bindgen` via `serde-wasm-bindgen` — but that boundary today
//! accepts plain JS objects/arrays, not the structured-clone-plus-
//! transferable-typed-array shapes `postMessage` will eventually carry.
//! That optimization is deferred until `pysees`'s compiler exists to
//! actually produce `CarapaceInputV1` bytes; see each table's doc comment
//! for the exact correspondence once it does.

pub mod error;
pub mod materials;
pub mod sequence;
pub mod tables;

mod decode;
mod session;

pub use decode::decode;
pub use error::DecodeError;
pub use materials::MaterialSpec;
pub use session::{AnalysisErrorDetail, RecorderBatch, Session, StepOutcome};

use sequence::SequenceSpec;
use serde::{Deserialize, Serialize};
use tables::{
    ElasticBeamColumnTable, ElementLoadTable, FiberBeamColumnTable, FiberTable, LoadPatternTable,
    NodalLoadTable, NodeTable, TrussTable, ZeroLengthTable,
};

/// Small structured-clone header fields — everything else in
/// `CarapaceInputV1` is a bulk table.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Header {
    /// The wire format's own version, independent of `engine_version`.
    pub schema_version: u32,
    /// `2` (planar, `NDM=2`/`NDF=3`) or `3` (spatial, `NDM=3`/`NDF=6`).
    /// [`decode`] rejects anything but `2` today.
    pub space: u8,
    /// `carapace-core`'s version, recorded for run provenance.
    pub engine_version: String,
}

/// The full wire payload: header plus one table per (profile, entity-kind)
/// pair, mirroring `core`'s closed-enum element/material catalog.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct CarapaceInputV1 {
    pub header: Header,
    pub nodes: NodeTable,
    pub materials: Vec<MaterialSpec>,
    pub fibers: FiberTable,
    pub trusses: TrussTable,
    pub elastic_beam_columns: ElasticBeamColumnTable,
    pub disp_beam_columns: FiberBeamColumnTable,
    pub force_beam_columns: FiberBeamColumnTable,
    pub zero_lengths: ZeroLengthTable,
    pub load_patterns: LoadPatternTable,
    pub nodal_loads: NodalLoadTable,
    pub element_loads: ElementLoadTable,
    pub sequence: SequenceSpec,
}
