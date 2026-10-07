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
pub use session::{
    AnalysisErrorDetail, ArcFailureDetail, ContinuationDetail, ContinuationStopDetail,
    ModalResultsReport, ModalStageResult, ModeResult, RecorderBatch, Session, Session2, Session3,
    StepOutcome, StopReasonDetail,
};

use sequence::SequenceSpec;
use serde::{Deserialize, Serialize};
use tables::{
    ElasticBeamColumn2dTable, ElasticBeamColumn3dTable, ElementLoadTable, EqualDofTable,
    FiberBeamColumn2dTable, FiberBeamColumn3dTable, FiberTable, LinearConstraintTable,
    LoadPatternTable, NodalLoadTable, NodeTable, PlaneMaterialSpec, QuadTable, RigidDiaphragmTable,
    RigidLinkTable, TriangleTable, TrussTable, ZeroLengthSectionTable, ZeroLengthTable,
};

/// Small structured-clone header fields — everything else in
/// `CarapaceInputV1` is a bulk table.
#[derive(Debug, Clone, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct Header {
    /// The wire format's own version, independent of `engine_version`.
    pub schema_version: u32,
    /// Model dimension: `2` (`NDM=2`/`NDF=3`) or `3` (`NDM=3`/`NDF=6`).
    /// [`decode`] rejects any other value.
    pub ndm: u8,
    /// `carapace-core`'s version, recorded for run provenance.
    pub engine_version: String,
    /// Record one sample of every supported recorder at the start of each `Static`/`Transient`
    /// stage, before its first step (the stage's initial conditions, at load factor/time `0`).
    /// Off by default: every stage's samples are then only the ones its steps produce.
    #[serde(default)]
    #[tsify(optional)]
    pub record_initial: bool,
}

/// The full wire payload: header plus one table per entity kind, mirroring
/// `core`'s closed-enum element/material catalog.
///
/// One format for both profiles (`header.ndm`). Dimension-agnostic entities
/// share a table; formulations that differ between 2D and 3D have a table per
/// formulation (`*2d`/`*3d`), and the other profile's must be empty or
/// omitted ([`DecodeError::TableNotInProfile`]). Every table may be omitted,
/// which is the same as empty.
#[derive(Debug, Clone, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct CarapaceInputV1 {
    pub header: Header,
    #[serde(default)]
    #[tsify(optional)]
    pub nodes: NodeTable,
    #[serde(default)]
    #[tsify(optional)]
    pub materials: Vec<MaterialSpec>,
    #[serde(default)]
    #[tsify(optional)]
    pub fibers: FiberTable,
    #[serde(default)]
    #[tsify(optional)]
    pub trusses: TrussTable,
    #[serde(default)]
    #[tsify(optional)]
    pub elastic_beam_columns_2d: ElasticBeamColumn2dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub elastic_beam_columns_3d: ElasticBeamColumn3dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub disp_beam_columns_2d: FiberBeamColumn2dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub disp_beam_columns_3d: FiberBeamColumn3dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub force_beam_columns_2d: FiberBeamColumn2dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub force_beam_columns_3d: FiberBeamColumn3dTable,
    #[serde(default)]
    #[tsify(optional)]
    pub zero_lengths: ZeroLengthTable,
    #[serde(default)]
    #[tsify(optional)]
    pub zero_length_sections: ZeroLengthSectionTable,
    #[serde(default)]
    #[tsify(optional)]
    pub equal_dofs: EqualDofTable,
    /// Arena of plane materials for `triangles`/`quads` (2D models only).
    #[serde(default)]
    #[tsify(optional)]
    pub plane_materials: Vec<PlaneMaterialSpec>,
    #[serde(default)]
    #[tsify(optional)]
    pub triangles: TriangleTable,
    #[serde(default)]
    #[tsify(optional)]
    pub quads: QuadTable,
    #[serde(default)]
    #[tsify(optional)]
    pub rigid_diaphragms: RigidDiaphragmTable,
    #[serde(default)]
    #[tsify(optional)]
    pub rigid_links: RigidLinkTable,
    #[serde(default)]
    #[tsify(optional)]
    pub linear_constraints: LinearConstraintTable,
    #[serde(default)]
    #[tsify(optional)]
    pub load_patterns: LoadPatternTable,
    #[serde(default)]
    #[tsify(optional)]
    pub nodal_loads: NodalLoadTable,
    #[serde(default)]
    #[tsify(optional)]
    pub element_loads: ElementLoadTable,
    #[serde(default)]
    #[tsify(optional)]
    pub sequence: SequenceSpec,
}
