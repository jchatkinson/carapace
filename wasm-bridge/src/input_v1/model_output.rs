//! `ModelOutput`: what decoding resolved that the caller cannot recompute cheaply and needs to
//! interpret results (plan `docs/csfm-plan.md` §2.6, §2.7). Today that is the embedded-bar
//! segment table, which `embedded_bars` rows are discretized into, with the TCM/POM routing and
//! crack spacing each bar's `csfmRebar` material resolved to.

use serde::Serialize;

use super::tables::ElementKind;

/// Which average-strain model a `csfmRebar` bar's effective reinforcement ratio routed to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum RebarRouting {
    /// Tension Chord Model (`rhoEff >= rho_crit`).
    Tcm,
    /// Pull-Out Model.
    Pom,
}

/// One discretized piece of an embedded bar, in the undeformed configuration.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct EmbeddedBarSegment {
    /// `[x, y]` of the segment start, in the bar's drawing direction.
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub length: f64,
    /// The host element's kind (`tri3` or `quad4`).
    pub host_kind: ElementKind,
    /// The host's nodes as `nodes` table rows, in the host's own node order.
    pub host_nodes: Vec<u32>,
}

/// One `embeddedBars` row's segments. Results are indexed `(bar, segment)`: the segment `s` of
/// bar `b` is element `embeddedBar` index `firstSegment + s`.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct EmbeddedBarGroup {
    pub first_segment: u32,
    pub segments: Vec<EmbeddedBarSegment>,
    /// `Some` only for a `csfmRebar` material.
    #[tsify(optional)]
    pub routing: Option<RebarRouting>,
    /// Crack spacing `sr` (mm; `0` for POM). `Some` only for a `csfmRebar` material.
    #[tsify(optional)]
    pub sr: Option<f64>,
    /// Present when the bar has a `bond` spec: where each bar point went.
    #[tsify(optional)]
    pub bond: Option<BondGroup>,
}

/// The bond-slip discretization of one bar: one entry per bar point, in polyline order (segment `s` joins
/// points `s` and `s + 1`).
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct BondGroup {
    pub points: Vec<BondPointInfo>,
}

#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct BondPointInfo {
    /// Undeformed position `[x, y]` (mm).
    pub position: [f64; 2],
    /// Unit tangent used for the slip direction.
    pub tangent: [f64; 2],
    /// `true` at a slip-locked vertex (no bond spring there).
    pub locked: bool,
    /// Tributary length of the bond spring (mm); `0` at a locked vertex.
    pub tributary: f64,
    /// Row of the bond spring in the `zeroLengths` table (the table's own rows first, then the bond springs).
    #[tsify(optional)]
    pub spring: Option<u32>,
    /// Row of the anchorage spring, at a bar end that has one.
    #[tsify(optional)]
    pub end_spring: Option<u32>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ModelOutput {
    pub embedded_bars: Vec<EmbeddedBarGroup>,
}
