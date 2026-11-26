//! The building description: storeys of named walls, floors and lateral loads. The wall layout, anchorage and
//! load wire types are the stud-wall tables' own (`input_v1::walls`), re-exported here.

use serde::{Deserialize, Serialize};

pub use crate::input_v1::walls::*;

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BuildingSpec {
    pub storeys: Vec<StoreySpec>,
    pub foundation_anchorage: AnchorageSpec,
    pub diaphragm: DiaphragmKind,
    #[serde(default = "default_combination")]
    #[cfg_attr(feature = "tsify", tsify(optional))]
    pub combination: [f64; 8],
}

fn default_combination() -> [f64; 8] {
    [1.0; 8]
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct StoreySpec {
    pub wall_height: f64,
    pub floor_depth: f64,
    pub band_ga: f64,
    pub rim_ea: f64,
    pub anchorage: AnchorageSpec,
    pub bearing: f64,
    pub walls: Vec<WallSpec>,
    #[serde(default)]
    pub floors: Vec<FloorSpec>,
    #[serde(default)]
    pub lateral: Vec<LateralLoad>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum DiaphragmKind {
    Rigid,
    Flexible,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallSpec {
    pub name: String,
    pub start: [f64; 2],
    pub direction: [f64; 2],
    pub layout: WallLayoutSpec,
    #[serde(default)]
    pub loads: Vec<SlotWallLoad>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SlotWallLoad {
    pub slot: usize,
    pub load: WallLoadSpec,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FloorSpec {
    pub bounds: [f64; 4],
    pub supports: Vec<String>,
    pub span_axis: SpanAxis,
    pub loads: Vec<FloorLoadSpec>,
    pub continuity: FloorContinuitySpec,
}
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum SpanAxis {
    X,
    Z,
}
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    tag = "kind"
)]
pub enum FloorLoadSpec {
    Area {
        slot: usize,
        w: f64,
        region: Option<[f64; 4]>,
    },
    Line {
        slot: usize,
        w: f64,
        s: f64,
        t0: f64,
        t1: f64,
    },
    LineAlong {
        slot: usize,
        w: f64,
        t: f64,
        s0: f64,
        s1: f64,
    },
    Point {
        slot: usize,
        p: f64,
        s: f64,
        t: f64,
    },
}
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LateralLoad {
    pub slot: usize,
    pub fx: f64,
    pub fz: f64,
    pub at: [f64; 2],
}
