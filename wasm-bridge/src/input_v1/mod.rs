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
pub mod tables3;

mod decode;
mod decode3;
mod session;

pub use decode::decode;
pub use error::DecodeError;
pub use materials::MaterialSpec;
pub use session::{AnalysisErrorDetail, RecorderBatch, Session, StepOutcome};

use sequence::{SequenceSpec, SequenceSpec3};
use serde::{Deserialize, Serialize};
use tables::{
    ElasticBeamColumnTable, ElementLoadTable, EqualDofTable, FiberBeamColumnTable, FiberTable,
    LoadPatternTable, NodalLoadTable, NodeTable, RigidDiaphragmTable, TrussTable,
    ZeroLengthSectionTable, ZeroLengthTable,
};
use tables3::{
    ElasticBeamColumnTable3, ElementLoadTable3, EqualDofTable3, FiberBeamColumnTable3, FiberTable3,
    NodeTable3, RigidDiaphragmTable3, TrussTable3, ZeroLengthSectionTable3, ZeroLengthTable3,
};

/// Small structured-clone header fields — everything else in
/// `CarapaceInputV1` is a bulk table.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Header {
    /// The wire format's own version, independent of `engine_version`.
    pub schema_version: u32,
    /// `2` (planar, `NDM=2`/`NDF=3`) or `3` (spatial, `NDM=3`/`NDF=6`).
    /// [`decode`] rejects any other value.
    pub space: u8,
    /// `carapace-core`'s version, recorded for run provenance.
    pub engine_version: String,
}

/// The full wire payload: header plus one table per (profile, entity-kind)
/// pair, mirroring `core`'s closed-enum element/material catalog.
///
/// Every table is always present, possibly empty (the handoff's per-table
/// "presence directory" isn't modeled yet — see `tables.rs`'s module doc
/// comment), including the `*3` spatial tables when `header.space == 2` and
/// vice versa: `decode` only ever reads the table set matching
/// `header.space`, ignoring the other profile's tables entirely.
/// `materials`/`load_patterns`/`nodal_loads` are dimension-agnostic and so
/// are shared by both profiles rather than duplicated as `materials3`/etc.
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
    pub zero_length_sections: ZeroLengthSectionTable,
    pub equal_dofs: EqualDofTable,
    pub rigid_diaphragms: RigidDiaphragmTable,
    pub load_patterns: LoadPatternTable,
    pub nodal_loads: NodalLoadTable,
    pub element_loads: ElementLoadTable,
    pub sequence: SequenceSpec,

    pub nodes3: NodeTable3,
    pub fibers3: FiberTable3,
    pub trusses3: TrussTable3,
    pub elastic_beam_columns3: ElasticBeamColumnTable3,
    pub disp_beam_columns3: FiberBeamColumnTable3,
    pub force_beam_columns3: FiberBeamColumnTable3,
    pub zero_lengths3: ZeroLengthTable3,
    pub zero_length_sections3: ZeroLengthSectionTable3,
    pub equal_dofs3: EqualDofTable3,
    pub rigid_diaphragms3: RigidDiaphragmTable3,
    /// Nodal loads for the spatial profile reuse [`NodalLoadTable`] as-is —
    /// `(pattern, node, dof, value)` needs nothing profile-specific, `dof`
    /// simply ranges up to 5 instead of 2.
    pub nodal_loads3: NodalLoadTable,
    pub element_loads3: ElementLoadTable3,
    pub sequence3: SequenceSpec3,
}
