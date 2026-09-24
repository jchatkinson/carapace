//! Bulk per-(entity-kind) tables from pysees-handoff.md's wire format.
//!
//! Node/element/pattern references inside these tables are already
//! resolved to dense 0-based row indices by the (not-yet-written)
//! `pysees` compiler — the "duplicate/missing tags and every reference
//! between entities" validation the handoff assigns to the compiler
//! happens upstream of decode, so decode only ever sees plain array
//! indices, never sparse OpenSees-style tags.
//!
//! Each field is a plain `Vec`, not yet the transferable `Float64Array`/
//! `Uint8Array` values the handoff's wire format ultimately specifies for
//! `postMessage`: `boundary.rs`'s `serde-wasm-bindgen` decoding today
//! accepts an ordinary JS array for each of these (or a typed array,
//! copied element-by-element) rather than transferring one. Once real
//! transferable-array support lands, each `Vec<f64>`/`Vec<u32>` here is
//! exactly the shape a `Float64Array`/`Uint32Array` copies into. Likewise,
//! the handoff's per-table `tables` presence directory (a table for an
//! unsupported variant is *absent*, not present-but-empty) isn't modeled
//! yet — every table below is always present, possibly empty.
//!
//! `#[serde(rename_all = "camelCase")]` throughout so the JS/TS shape
//! matches pysees-handoff.md's own naming (`nodeI`, not `node_i`).

use serde::{Deserialize, Serialize};

/// Node table: `coords` stride 2 (x, y); `fixed` one bitmask byte per node
/// (bit 0 = ux, bit 1 = uy, bit 2 = rz); mass is sparse, addressed via a
/// parallel node-index array since most nodes carry none.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NodeTable {
    pub coords: Vec<f64>,
    pub fixed: Vec<u8>,
    pub mass_node_index: Vec<u32>,
    /// Stride 3 (mass_x, mass_y, mass_rz), parallel to `mass_node_index`.
    pub mass: Vec<f64>,
}

impl NodeTable {
    pub fn node_count(&self) -> usize {
        self.coords.len() / 2
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum TransformSpec {
    Linear,
    PDelta,
    Corotational,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum IntegrationSpec {
    Legendre { points: u32 },
    Lobatto { points: u32 },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TrussTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub area: Vec<f64>,
    /// Index into the material arena.
    pub material: Vec<u32>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ElasticBeamColumnTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub e: Vec<f64>,
    pub a: Vec<f64>,
    pub iz: Vec<f64>,
    pub transform: Vec<TransformSpec>,
    pub density: Vec<f64>,
}

/// Shared shape for `DispBeamColumn` and `ForceBeamColumn` — both are one
/// prismatic fiber section (see [`FiberTable`]) replicated across
/// integration points by `core`'s own element constructors.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FiberBeamColumnTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Index into `FiberTable::section_offsets`.
    pub fiber_section: Vec<u32>,
    pub integration: Vec<IntegrationSpec>,
    pub corotational: Vec<bool>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ZeroLengthTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Sparse `(zero_length row index, dof, material arena index)` — most
    /// DOFs on a given `ZeroLength` carry no material at all.
    pub materials: Vec<(u32, u8, u32)>,
}

/// Fibers for every `DispBeamColumn`/`ForceBeamColumn` section, flattened
/// and offset-indexed: section `k` occupies
/// `section_offsets[k]..section_offsets[k + 1]` in `y`/`area`/`material`.
/// `section_offsets` therefore has `num_sections + 1` entries.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FiberTable {
    pub section_offsets: Vec<u32>,
    pub y: Vec<f64>,
    pub area: Vec<f64>,
    /// Index into the material arena, parallel to `y`/`area`.
    pub material: Vec<u32>,
}

impl FiberTable {
    pub fn section_count(&self) -> usize {
        self.section_offsets.len().saturating_sub(1)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum TimeSeriesSpec {
    Constant,
    Linear { slope: f64 },
    Path { times: Vec<f64>, factors: Vec<f64> },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LoadPatternTable {
    pub series: Vec<TimeSeriesSpec>,
    pub scale_factor: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NodalLoadTable {
    /// Index into `LoadPatternTable`.
    pub pattern: Vec<u32>,
    /// Index into `NodeTable`.
    pub node: Vec<u32>,
    pub dof: Vec<u8>,
    pub value: Vec<f64>,
    /// Index into `SequenceSpec::stages`: this load is only registered on
    /// the `Domain` when that stage starts, not at decode time. This
    /// matters whenever an earlier stage's `LoadControl` shares one
    /// pseudo-time with every unfrozen pattern (core/src/analysis/
    /// integrator.rs) — a pattern's reference load must not exist yet if
    /// an earlier stage isn't meant to ramp it too (the same reason
    /// core/tests/m8_force_beam_column.rs's native two-phase test adds its
    /// lateral pattern's load only between phases, not upfront).
    pub stage: Vec<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum ElementKind {
    Truss,
    ElasticBeamColumn,
    DispBeamColumn,
    ForceBeamColumn,
    ZeroLength,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ElementLoadSpec {
    /// Only `ElasticBeamColumn` currently honors this (`core`'s
    /// `Element::form_load_vector` falls through to zero for every other
    /// kind) — decode rejects any other `element_kind` explicitly rather
    /// than silently accepting a load that will never apply.
    ///
    /// A struct-like (not tuple) variant: serde's internally-tagged enum
    /// representation (`tag = "kind"`, needed for a JS-friendly
    /// discriminated union) can't tag a bare-scalar tuple variant.
    UniformTransverse { w: f64 },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ElementLoadTable {
    pub pattern: Vec<u32>,
    pub element_kind: Vec<ElementKind>,
    /// Row index into the table named by the parallel `element_kind` entry.
    pub element_index: Vec<u32>,
    pub load: Vec<ElementLoadSpec>,
    /// See [`NodalLoadTable::stage`] — same per-stage registration timing.
    pub stage: Vec<u32>,
}
