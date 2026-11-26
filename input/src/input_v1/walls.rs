//! Stud-wall wire types: the `studWalls` table's layout, stations and bearing, `wallLoad` element loads,
//! `Anchorage` materials and `floorLoads` continuity, with their conversions to the core types. The building
//! description (`crate::building`) uses the same types.

use carapace_core::model::{
    AnchorageBearing, ChordDistribution, FloorContinuity, FloorZone, Material, StudWallLayout,
    WallBearing, WallChord, WallLoad, WallOpening, WallShearLaw, WallSheathing, WallStud,
};
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AnchorageSpec {
    pub tension: f64,
    pub slack: f64,
    pub bearing: AnchorageBearingSpec,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    tag = "kind"
)]
pub enum AnchorageBearingSpec {
    Linear { k: f64 },
    O86Crush { qr: f64 },
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum WallChordSpec {
    Top,
    Bottom,
    /// A node the wall stands on (with `layout.bearing`): on the bottom line, any `x`.
    Base,
}

/// The wall bears on its `base` stations through the plate sandwich: `kv` / `kh` per bearing position and the
/// pier-end hold-down law (stud-wall tables only; the building generates its own links).
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallBearingSpec {
    pub kv: f64,
    pub kh: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[cfg_attr(feature = "tsify", tsify(optional))]
    pub hold_down: Option<AnchorageSpec>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
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
    #[serde(default)]
    #[cfg_attr(feature = "tsify", tsify(optional))]
    pub foundation: bool,
    /// Bearing on base stations (stud-wall tables only; rejected in a building).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[cfg_attr(feature = "tsify", tsify(optional))]
    pub bearing: Option<WallBearingSpec>,
    /// Floor zone above the framing (stud-wall tables only; a building takes it from its storey).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    #[cfg_attr(feature = "tsify", tsify(optional))]
    pub floor: Option<FloorZoneSpec>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FloorZoneSpec {
    pub depth: f64,
    pub ga: f64,
    pub rim_ea: f64,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallStudSpec {
    pub x: f64,
    pub ea: f64,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
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

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum ChordDistributionSpec {
    Beam,
    Direct,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallSheathingSpec {
    pub law: WallShearLawSpec,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    tag = "kind"
)]
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

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    rename_all = "camelCase",
    rename_all_fields = "camelCase",
    tag = "kind"
)]
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

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FloorContinuitySpec {
    Simple,
    Continuous,
}

impl AnchorageBearingSpec {
    pub fn to_core(self) -> AnchorageBearing {
        match self {
            Self::Linear { k } => AnchorageBearing::Linear { k },
            Self::O86Crush { qr } => AnchorageBearing::O86Crush { qr },
        }
    }
}

impl WallChordSpec {
    pub fn to_core(self) -> WallChord {
        match self {
            Self::Top => WallChord::Top,
            Self::Bottom => WallChord::Bottom,
            Self::Base => WallChord::Base,
        }
    }
}

impl WallShearLawSpec {
    pub fn to_core(self) -> WallShearLaw {
        match self {
            Self::Linear { ga } => WallShearLaw::Linear { ga },
            Self::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            } => WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            },
        }
    }
}

impl WallLoadSpec {
    pub fn to_core(self) -> WallLoad {
        match self {
            Self::TopPoint { x, fx, fy } => WallLoad::TopPoint { x, fx, fy },
            Self::TopLine {
                x0,
                x1,
                wx0,
                wx1,
                wy0,
                wy1,
            } => WallLoad::TopLine {
                x0,
                x1,
                wx0,
                wx1,
                wy0,
                wy1,
            },
        }
    }
}

impl FloorContinuitySpec {
    pub fn to_core(self) -> FloorContinuity {
        match self {
            Self::Simple => FloorContinuity::Simple,
            Self::Continuous => FloorContinuity::Continuous,
        }
    }
}

impl WallLayoutSpec {
    /// The core layout (the bearing's hold-down becomes an `Anchorage` material). `StudWall::new` validates it.
    pub fn to_core(&self) -> Result<StudWallLayout, String> {
        let bearing = self
            .bearing
            .map(|b| {
                let hold_down = b
                    .hold_down
                    .map(|a| {
                        Material::anchorage(a.tension, a.slack, a.bearing.to_core())
                            .map_err(|e| format!("bearing.holdDown: {e}"))
                    })
                    .transpose()?;
                Ok::<_, String>(WallBearing {
                    kv: b.kv,
                    kh: b.kh,
                    hold_down,
                })
            })
            .transpose()?;
        Ok(StudWallLayout {
            length: self.length,
            height: self.height,
            studs: self
                .studs
                .iter()
                .map(|s| WallStud { x: s.x, ea: s.ea })
                .collect(),
            top_ea: self.top_ea,
            top_ei: self.top_ei,
            top_distribution: match self.top_distribution {
                ChordDistributionSpec::Beam => ChordDistribution::Beam,
                ChordDistributionSpec::Direct => ChordDistribution::Direct,
            },
            bottom_ea: self.bottom_ea,
            bottom_ei: self.bottom_ei,
            openings: self
                .openings
                .iter()
                .map(|o| WallOpening {
                    x0: o.x0,
                    x1: o.x1,
                    head_height: o.head_height,
                    header_ei: o.header_ei,
                    jack_ea: o.jack_ea,
                    door: o.door,
                })
                .collect(),
            sheathing: self.sheathing.as_ref().map(|s| WallSheathing {
                law: s.law.to_core(),
            }),
            foundation: self.foundation,
            bearing,
            floor: self.floor.map(|z| FloorZone {
                depth: z.depth,
                ga: z.ga,
                rim_ea: z.rim_ea,
            }),
        })
    }
}
