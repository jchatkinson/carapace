//! Serde-safe mirrors of the core wall and floor inputs used by [`super::build`].

use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BuildingSpec {
    pub storeys: Vec<StoreySpec>,
    pub foundation_anchorage: AnchorageSpec,
    pub diaphragm: DiaphragmKind,
}

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

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AnchorageSpec {
    pub tension: f64,
    pub slack: f64,
    pub bearing: AnchorageBearingSpec,
}

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "kind")]
pub enum AnchorageBearingSpec {
    Linear { k: f64 },
    O86Crush { qr: f64 },
}

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum DiaphragmKind {
    Rigid,
    Flexible,
}

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

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SlotWallLoad {
    pub slot: usize,
    pub load: WallLoadSpec,
}

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallLayoutSpec {
    pub length: f64,
    pub height: f64,
    pub studs: Vec<WallStudSpec>,
    pub top_ea: f64,
    pub top_ei: Option<f64>,
    pub top_distribution: ChordDistributionSpec,
    pub bottom_ea: f64,
    pub bottom_ei: f64,
    #[serde(default)]
    pub openings: Vec<WallOpeningSpec>,
    pub sheathing: Option<WallSheathingSpec>,
}

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallStudSpec {
    pub x: f64,
    pub ea: f64,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallOpeningSpec {
    pub x0: f64,
    pub x1: f64,
    pub head_height: f64,
    pub header_ei: f64,
    pub jack_ea: f64,
    pub door: bool,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum ChordDistributionSpec {
    Beam,
    Direct,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallSheathingSpec {
    pub law: WallShearLawSpec,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "kind")]
pub enum WallShearLawSpec {
    Linear {
        ga: f64,
    },
    O86 {
        bv: f64,
        nail_spacing: f64,
        nail_diameter: f64,
    },
}

#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "kind")]
pub enum WallLoadSpec {
    TopPoint {
        x: f64,
        fx: f64,
        fy: f64,
    },
    TopLine {
        x0: f64,
        x1: f64,
        wx0: f64,
        wx1: f64,
        wy0: f64,
        wy1: f64,
    },
}

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
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum SpanAxis {
    X,
    Z,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FloorContinuitySpec {
    Simple,
    Continuous,
}
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "kind")]
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
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LateralLoad {
    pub slot: usize,
    pub fx: f64,
    pub fz: f64,
    pub at: [f64; 2],
}
