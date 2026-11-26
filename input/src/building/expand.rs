//! `expand`: a [`BuildingSpec`] becomes an ordinary explicit model (nodes, materials, `studWalls`, `zeroLengths`,
//! `equalDofs`, rigid diaphragms, load patterns and loads) plus an index map. See `docs/building.md`.

use std::collections::{BTreeMap, BTreeSet};

use carapace_core::model::Material;
use serde::{Deserialize, Serialize};

use super::build::{convert_layout, err, floor_transfer_for, BuildError, SLOTS};
use super::spec::*;
use crate::input_v1::materials::MaterialSpec;
use crate::input_v1::sequence::SequenceSpec;
use crate::input_v1::tables::*;
use crate::input_v1::{CarapaceInputV1, Header};

/// Which model axis is vertical. The spec's plan coordinates `(x, z)` map to model `(x, z)` for `Y` (the
/// default) and to model `(x, y)` for `Z`; elevations go along `+Y` or `+Z`.
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum UpAxis {
    #[default]
    Y,
    Z,
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", default)]
pub struct ExpandOptions {
    pub up: UpAxis,
}

/// The generated model and where its parts are. `input.sequence` is empty and every load row has `stage` 0:
/// append at least one stage (and recorders) before decoding.
#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Expansion {
    pub input: CarapaceInputV1,
    pub map: ExpansionMap,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ExpansionMap {
    pub storeys: Vec<StoreyMap>,
    /// Per load slot (`0..8`): its row in `loadPatterns`, or `null` when no load uses the slot.
    pub slot_patterns: Vec<Option<u32>>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct StoreyMap {
    /// In spec order.
    pub walls: Vec<WallMap>,
    /// `zeroLengths` rows joining this storey's walls to the walls below (empty for storey 0).
    pub links_below: Vec<u32>,
    /// The rigid diaphragm's retained node (`nodes` index); `null` for flexible diaphragms.
    pub master: Option<u32>,
    /// Every wall's exterior top node (`nodes` indices, ascending): the diaphragm nodes.
    pub nodes: Vec<u32>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WallMap {
    pub name: String,
    /// Row in `studWalls`.
    pub element: u32,
    /// `nodes` indices: bottom start, bottom end, top end, top start.
    pub corners: [u32; 4],
    /// Exterior top row (the floor zone's top when there is one): wall-local `x` and `nodes` index.
    pub top_nodes: Vec<StationNode>,
    /// Bottom stations: pier ends on the foundation, the supporting stud positions above it.
    pub bottom_nodes: Vec<StationNode>,
    /// The anchorage `zeroLengths` row at each pier end: to the ground (storey 0) or to the wall below.
    pub pier_links: Vec<PierLink>,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct StationNode {
    pub x: f64,
    pub node: u32,
}

#[cfg_attr(feature = "schema", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "tsify", derive(tsify::Tsify))]
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PierLink {
    pub x: f64,
    pub element: u32,
}

struct WallGen {
    corners: [u32; 4],
    top: BTreeMap<i64, u32>,
    bottom: BTreeMap<i64, u32>,
    pier_links: Vec<PierLink>,
    row: u32,
}

fn key(x: f64) -> i64 {
    (x * 1e6).round() as i64
}
fn stations(m: &BTreeMap<i64, u32>) -> Vec<StationNode> {
    m.iter()
        .map(|(k, &node)| StationNode {
            x: *k as f64 / 1e6,
            node,
        })
        .collect()
}

fn anchorage(a: &AnchorageSpec) -> Result<MaterialSpec, BuildError> {
    if a.tension < 0. || a.slack < 0. || !a.tension.is_finite() || !a.slack.is_finite() {
        return Err(err("invalid anchorage values"));
    }
    Material::anchorage(a.tension, a.slack, a.bearing.to_core()).map_err(|e| err(e.to_string()))?;
    Ok(MaterialSpec::Anchorage {
        tension: a.tension,
        slack: a.slack,
        bearing: a.bearing,
    })
}

struct Gen {
    up: UpAxis,
    coords: Vec<[f64; 3]>,
    fixed: Vec<u8>,
    materials: Vec<MaterialSpec>,
    zero_lengths: ZeroLengthTable,
}
impl Gen {
    /// Model DOFs of plan x and plan z.
    fn plan(&self) -> [u8; 2] {
        match self.up {
            UpAxis::Y => [0, 2],
            UpAxis::Z => [0, 1],
        }
    }
    fn vertical(&self) -> u8 {
        match self.up {
            UpAxis::Y => 1,
            UpAxis::Z => 2,
        }
    }
    fn node(&mut self, plan: [f64; 2], elevation: f64) -> u32 {
        self.coords.push(match self.up {
            UpAxis::Y => [plan[0], elevation, plan[1]],
            UpAxis::Z => [plan[0], plan[1], elevation],
        });
        self.fixed.push(0);
        self.coords.len() as u32 - 1
    }
    fn material(&mut self, m: MaterialSpec) -> u32 {
        let at = self
            .materials
            .iter()
            .position(|x| format!("{x:?}") == format!("{m:?}"));
        at.unwrap_or_else(|| {
            self.materials.push(m);
            self.materials.len() - 1
        }) as u32
    }
    /// A `ZeroLength3` row `i -> j` with a material on each listed DOF.
    fn zero_length(&mut self, i: u32, j: u32, materials: &[(u8, MaterialSpec)]) -> u32 {
        let t = &mut self.zero_lengths;
        t.node_i.push(i);
        t.node_j.push(j);
        let row = t.node_i.len() as u32 - 1;
        for (dof, m) in materials {
            let m = self.material(m.clone());
            self.zero_lengths.materials.push((row, *dof, m));
        }
        row
    }
}

/// Expand `spec` into an explicit model. See `docs/building.md` for what is generated, the index map and the
/// axis conventions.
pub fn expand(spec: &BuildingSpec, options: &ExpandOptions) -> Result<Expansion, BuildError> {
    if spec.storeys.is_empty() {
        return Err(err("building must have at least one storey"));
    }
    if spec.combination.iter().any(|v| !v.is_finite()) {
        return Err(err("combination factors must be finite"));
    }
    let foundation = anchorage(&spec.foundation_anchorage)?;
    let mut g = Gen {
        up: options.up,
        coords: Vec::new(),
        fixed: Vec::new(),
        materials: Vec::new(),
        zero_lengths: ZeroLengthTable::default(),
    };
    let [px, pz] = g.plan();
    let vdof = g.vertical();
    let elastic = MaterialSpec::Elastic { e: 1e12 };
    let mut gens: Vec<Vec<WallGen>> = Vec::new();
    let mut stud_walls: Vec<StudWallSpec> = Vec::new();
    for (si, s) in spec.storeys.iter().enumerate() {
        if !s.wall_height.is_finite()
            || s.wall_height <= 0.
            || !s.floor_depth.is_finite()
            || s.floor_depth < 0.
        {
            return Err(err(format!("storey {si}: wallHeight/floorDepth invalid")));
        }
        if !s.band_ga.is_finite()
            || s.band_ga < 0.
            || !s.rim_ea.is_finite()
            || s.rim_ea < 0.
            || !s.bearing.is_finite()
            || s.bearing <= 0.
        {
            return Err(err(format!(
                "storey {si}: band, rim, or bearing properties are invalid"
            )));
        }
        anchorage(&s.anchorage)?;
        if s.walls.is_empty() {
            return Err(err(format!("storey {si}: no walls")));
        }
        let y0: f64 = 0.
            + spec.storeys[..si]
                .iter()
                .map(|p| p.wall_height + p.floor_depth)
                .sum::<f64>();
        // The exterior top row is the floor zone's top: `H + depth` above the base.
        let y_top = y0 + s.wall_height + s.floor_depth;
        let mut level: Vec<WallGen> = Vec::new();
        for w in &s.walls {
            let len = w.layout.length;
            let [ex, ez] = w.direction;
            if !w
                .start
                .iter()
                .chain(w.direction.iter())
                .all(|v| v.is_finite())
                || (ex.hypot(ez) - 1.).abs() > 1e-9
            {
                return Err(err(format!(
                    "storey {si} wall {}: direction must be finite and unit",
                    w.name
                )));
            }
            if !len.is_finite()
                || len <= 0.
                || (w.layout.height - s.wall_height).abs() > 1e-9 * s.wall_height.max(1.)
            {
                return Err(err(format!(
                    "storey {si} wall {}: invalid length or height",
                    w.name
                )));
            }
            if s.walls.iter().filter(|o| o.name == w.name).count() > 1 {
                return Err(err(format!("storey {si}: duplicate wall {}", w.name)));
            }
            if w.layout.bearing.is_some() {
                return Err(err(format!(
                    "storey {si} wall {}: a building generates its own links; the wall layout takes no bearing",
                    w.name
                )));
            }
            if w.layout.floor.is_some() {
                return Err(err(format!(
                    "storey {si} wall {}: the floor zone comes from the storey (floorDepth, bandGa, rimEa), not the wall layout",
                    w.name
                )));
            }
            let mut wire = w.layout.clone();
            wire.foundation = si == 0;
            wire.floor = (s.floor_depth > 0.).then_some(FloorZoneSpec {
                depth: s.floor_depth,
                ga: s.band_ga,
                rim_ea: s.rim_ea,
            });
            let layout = convert_layout(&wire, si, &w.name)?;
            if si > 0 {
                let matching = spec.storeys[si - 1].walls.iter().any(|b| {
                    b.name == w.name
                        && b.start == w.start
                        && b.direction == w.direction
                        && (b.layout.length - len).abs() <= 1e-9 * len
                });
                if !matching {
                    return Err(err(format!(
                        "storey {si} wall {} does not stack on a wall below",
                        w.name
                    )));
                }
            }
            let at = |x: f64| [w.start[0] + ex * x, w.start[1] + ez * x];
            let corners = [
                g.node(at(0.), y0),
                g.node(at(len), y0),
                g.node(at(len), y_top),
                g.node(at(0.), y_top),
            ];
            let mut xs: Vec<f64> = layout.studs.iter().map(|x| x.x).collect();
            xs.extend(layout.pier_ends());
            xs.sort_by(f64::total_cmp);
            xs.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * len);
            let mut top = BTreeMap::new();
            for &x in &xs {
                let n = if x == 0. {
                    corners[3]
                } else if x == len {
                    corners[2]
                } else {
                    g.node(at(x), y_top)
                };
                top.insert(key(x), n);
            }
            // Above the first storey the bottom stations stand on the wall below: one at each of its top
            // positions (the corners and every stud), where the vertical links go. Pier ends (hold-downs)
            // must be among them. On the foundation only the pier ends have a bottom station.
            let mut bxs = layout.pier_ends();
            if si > 0 {
                let below = spec.storeys[si - 1]
                    .walls
                    .iter()
                    .find(|a| a.name == w.name)
                    .expect("stacking was checked above");
                let below = convert_layout(&below.layout, si - 1, &below.name)?;
                let supports: Vec<f64> = below.studs.iter().map(|x| x.x).collect();
                for x in layout.pier_ends() {
                    if !supports.iter().any(|b| (b - x).abs() < 1e-9 * len) {
                        return Err(err(format!(
                            "storey {si} wall {}: pier end at x = {x} has no stud in the wall below to hold down to",
                            w.name
                        )));
                    }
                }
                bxs = supports;
            }
            bxs.sort_by(f64::total_cmp);
            bxs.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * len);
            let mut bottom = BTreeMap::new();
            for &x in &bxs {
                let n = if x == 0. {
                    corners[0]
                } else if x == len {
                    corners[1]
                } else {
                    g.node(at(x), y0)
                };
                bottom.insert(key(x), n);
            }
            let mut st: Vec<WallStationSpec> = top
                .iter()
                .map(|(k, &node)| (*k as f64 / 1e6, node))
                .filter(|&(x, _)| x > 0. && x < len)
                .map(|(x, node)| WallStationSpec {
                    chord: WallChordSpec::Top,
                    x,
                    node,
                })
                .collect();
            let bottom_stations: Vec<(f64, u32)> = if si > 0 {
                bottom.iter().map(|(k, &n)| (*k as f64 / 1e6, n)).collect()
            } else {
                layout
                    .pier_ends()
                    .into_iter()
                    .map(|x| (x, bottom[&key(x)]))
                    .collect()
            };
            st.extend(
                bottom_stations
                    .into_iter()
                    .filter(|&(x, _)| x > 0. && x < len)
                    .map(|(x, node)| WallStationSpec {
                        chord: WallChordSpec::Bottom,
                        x,
                        node,
                    }),
            );
            let row = stud_walls.len() as u32;
            stud_walls.push(StudWallSpec {
                corners,
                stations: st,
                layout: wire,
            });
            // Foundation supports at the pier ends: fixed ground nodes joined by the anchorage law, the ground's
            // plan translations rigid. Vertical behavior is carried by the anchorage material.
            let mut pier_links = Vec::new();
            if si == 0 {
                for x in layout.pier_ends() {
                    let node = bottom[&key(x)];
                    let c = g.coords[node as usize];
                    g.coords.push(c);
                    g.fixed.push(0b111);
                    let ground = g.coords.len() as u32 - 1;
                    let mats: Vec<_> = (0..3u8)
                        .map(|d| {
                            (
                                d,
                                if d == vdof {
                                    foundation.clone()
                                } else {
                                    elastic.clone()
                                },
                            )
                        })
                        .collect();
                    let element = g.zero_length(ground, node, &mats);
                    pier_links.push(PierLink { x, element });
                }
            }
            level.push(WallGen {
                corners,
                top,
                bottom,
                pier_links,
                row,
            });
        }
        gens.push(level);
    }

    // Join each wall to the wall below: its bottom nodes sit on the lower wall's exterior top nodes (the floor
    // zone's top row), tied in the wall plane and linked vertically, with the hold-down anchorage at pier ends
    // and bearing (compression only) elsewhere.
    let mut equal = EqualDofTable::default();
    let mut links_below: Vec<Vec<u32>> = vec![Vec::new(); spec.storeys.len()];
    for si in 1..spec.storeys.len() {
        for (wi, w) in spec.storeys[si].walls.iter().enumerate() {
            let lower = spec.storeys[si - 1]
                .walls
                .iter()
                .position(|b| b.name == w.name)
                .expect("stacking was checked above");
            let low_top = stations(&gens[si - 1][lower].top);
            let up_bottom = stations(&gens[si][wi].bottom);
            let wire = &stud_walls[gens[si][wi].row as usize].layout;
            let pier = convert_layout(wire, si, &w.name)?.pier_ends();
            let dofs: Vec<u8> = [(px, w.direction[0]), (pz, w.direction[1])]
                .into_iter()
                .filter(|(_, c)| c.abs() > 1e-9)
                .map(|(d, _)| d)
                .collect();
            for StationNode { x, node: upper } in up_bottom {
                let lower = low_top
                    .iter()
                    .find(|a| (a.x - x).abs() < 1e-8)
                    .ok_or_else(|| {
                        err(format!(
                            "storey {si} wall {}: no node below at x = {x} to link to",
                            w.name
                        ))
                    })?
                    .node;
                let is_pier = pier.iter().any(|p| (*p - x).abs() < 1e-8);
                let mat = if is_pier {
                    anchorage(&spec.storeys[si].anchorage)?
                } else {
                    anchorage(&AnchorageSpec {
                        tension: 0.,
                        slack: 0.,
                        bearing: AnchorageBearingSpec::Linear {
                            k: spec.storeys[si].bearing,
                        },
                    })?
                };
                let element = g.zero_length(lower, upper, &[(vdof, mat)]);
                equal.retained.push(lower);
                equal.constrained.push(upper);
                let row = equal.retained.len() as u32 - 1;
                equal.dofs.extend(dofs.iter().map(|&d| (row, d)));
                links_below[si].push(element);
                if is_pier {
                    gens[si][wi].pier_links.push(PierLink { x, element });
                }
            }
        }
    }

    // Diaphragms and lateral loads.
    let mut rigid = RigidDiaphragmTable::default();
    let mut linear = LinearConstraintTable::default();
    let mut nodal: BTreeMap<(usize, u32, u8), f64> = BTreeMap::new();
    let mut levels: Vec<(Option<u32>, Vec<u32>)> = Vec::new();
    let rotation = match options.up {
        UpAxis::Y => 4,
        UpAxis::Z => 5,
    };
    for (si, s) in spec.storeys.iter().enumerate() {
        let mut nodes: Vec<u32> = gens[si]
            .iter()
            .flat_map(|w| w.top.values().copied())
            .collect();
        nodes.sort_unstable();
        nodes.dedup();
        let master = if spec.diaphragm == DiaphragmKind::Rigid {
            if nodes.is_empty() {
                return Err(err(format!("storey {si}: no diaphragm nodes")));
            }
            let xyz = nodes.iter().fold([0.; 3], |mut a, &n| {
                for (i, v) in g.coords[n as usize].iter().enumerate() {
                    a[i] += v / nodes.len() as f64;
                }
                a
            });
            g.coords.push(xyz);
            g.fixed.push(0);
            let m = g.coords.len() as u32 - 1;
            rigid.retained.push(m);
            rigid.normal.push(match options.up {
                UpAxis::Y => Axis3Spec::Y,
                UpAxis::Z => Axis3Spec::Z,
            });
            let row = rigid.retained.len() as u32 - 1;
            rigid.constrained.extend(nodes.iter().map(|&n| (row, n)));
            Some(m)
        } else {
            None
        };
        for l in &s.lateral {
            if l.slot >= SLOTS || ![l.fx, l.fz, l.at[0], l.at[1]].iter().all(|v| v.is_finite()) {
                return Err(err(format!("storey {si}: invalid lateral load")));
            }
            // Rigid load is attached to the retained node; flexible distribution is handled by floor helpers.
            if let Some(m) = master {
                *nodal.entry((l.slot, m, px)).or_default() += l.fx;
                *nodal.entry((l.slot, m, pz)).or_default() += l.fz;
                let c = g.coords[m as usize];
                // Moment about the vertical of the force at `at` about the master, `r x F` in the model's
                // right-handed axes: `fx (z - zm) - fz (x - xm)` about +Y, the negative of it about +Z (plan
                // `z` is model `y` there, so the plan frame is mirrored).
                let my = l.fx * (l.at[1] - c[pz as usize]) - l.fz * (l.at[0] - c[px as usize]);
                let sign = if options.up == UpAxis::Y { 1. } else { -1. };
                *nodal.entry((l.slot, m, rotation)).or_default() += sign * my;
            } else {
                distribute_flexible(&g, &mut nodal, &gens[si], &s.walls, l)?;
            }
        }
        if master.is_none() {
            for (w, wg) in s.walls.iter().zip(&gens[si]) {
                let [ex, ez] = w.direction;
                for &n in wg.top.values() {
                    if ez.abs() < 1e-10 {
                        g.fixed[n as usize] |= 1 << pz;
                    } else if ex.abs() < 1e-10 {
                        g.fixed[n as usize] |= 1 << px;
                    } else {
                        let (slave, term, c) = if ex.abs() >= ez.abs() {
                            (pz, px, ez / ex)
                        } else {
                            (px, pz, ex / ez)
                        };
                        if linear.term_offsets.is_empty() {
                            linear.term_offsets.push(0);
                        }
                        linear.slave_node.push(n);
                        linear.slave_dof.push(slave);
                        linear.term_node.push(n);
                        linear.term_dof.push(term);
                        linear.term_coeff.push(c);
                        linear.term_offsets.push(linear.term_node.len() as u32);
                    }
                }
            }
        }
        levels.push((master, nodes));
    }

    // Loads, by slot: wall loads, floor loads (expanded by the decoder into each support wall's top-plate loads)
    // and lateral loads. One pattern per used slot, scaled by its combination factor.
    let mut wall_loads: Vec<(usize, u32, WallLoadSpec)> = Vec::new();
    let mut floor_rows: Vec<(usize, FloorTransferSpec)> = Vec::new();
    for (si, s) in spec.storeys.iter().enumerate() {
        for (w, wg) in s.walls.iter().zip(&gens[si]) {
            for l in &w.loads {
                if l.slot >= SLOTS {
                    return Err(err(format!(
                        "storey {si} wall {}: load slot out of range",
                        w.name
                    )));
                }
                wall_loads.push((l.slot, wg.row, l.load));
            }
        }
        for (fi, f) in s.floors.iter().enumerate() {
            let (_, bounds, supports) = floor_transfer_for(f, &s.walls)
                .map_err(|e| err(format!("storey {si} floor {fi}: {e}")))?;
            let mut rows: BTreeMap<usize, Vec<FloorTransferLoadSpec>> = BTreeMap::new();
            for l in &f.loads {
                let (slot, load) = match *l {
                    FloorLoadSpec::Area { slot, w, region } => {
                        (slot, FloorTransferLoadSpec::Area { w, region })
                    }
                    FloorLoadSpec::Line { slot, w, s, t0, t1 } => {
                        (slot, FloorTransferLoadSpec::Line { w, s, t0, t1 })
                    }
                    FloorLoadSpec::LineAlong { slot, w, t, s0, s1 } => {
                        (slot, FloorTransferLoadSpec::LineAlong { w, t, s0, s1 })
                    }
                    FloorLoadSpec::Point { slot, p, s, t } => {
                        (slot, FloorTransferLoadSpec::Point { p, s, t })
                    }
                };
                if slot >= SLOTS {
                    return Err(err(format!(
                        "storey {si} floor {fi}: load slot out of range"
                    )));
                }
                rows.entry(slot).or_default().push(load);
            }
            let support_rows = f
                .supports
                .iter()
                .zip(&supports)
                .map(|(name, support)| {
                    let wall = s
                        .walls
                        .iter()
                        .position(|w| &w.name == name)
                        .ok_or_else(|| {
                            err(format!(
                                "storey {si} floor {fi}: unknown support wall {name}"
                            ))
                        })?;
                    let w = &s.walls[wall];
                    let (coord, sign) = match f.span_axis {
                        SpanAxis::X => (w.start[1], w.direction[1]),
                        SpanAxis::Z => (w.start[0], w.direction[0]),
                    };
                    Ok(FloorTransferSupportSpec {
                        s: support.s,
                        t0: support.t0,
                        t1: support.t1,
                        wall: WallTargetSpec::StudWall {
                            index: gens[si][wall].row,
                        },
                        x_at_t0: (bounds.t0 - coord) / sign,
                        reversed: sign < 0.,
                    })
                })
                .collect::<Result<Vec<_>, BuildError>>()?;
            for (slot, loads) in rows {
                floor_rows.push((
                    slot,
                    FloorTransferSpec {
                        bounds: [bounds.s0, bounds.s1, bounds.t0, bounds.t1],
                        continuity: f.continuity,
                        supports: support_rows.clone(),
                        pattern: 0,
                        loads,
                    },
                ));
            }
        }
    }
    let used: BTreeSet<usize> = wall_loads
        .iter()
        .map(|l| l.0)
        .chain(floor_rows.iter().map(|l| l.0))
        .chain(nodal.keys().map(|k| k.0))
        .collect();
    let mut slot_patterns = vec![None; SLOTS];
    let mut load_patterns = LoadPatternTable::default();
    for (row, &slot) in used.iter().enumerate() {
        slot_patterns[slot] = Some(row as u32);
        load_patterns
            .series
            .push(TimeSeriesSpec::Linear { slope: 1. });
        load_patterns.scale_factor.push(spec.combination[slot]);
    }
    let pattern = |slot: usize| slot_patterns[slot].expect("used slot");
    let mut nodal_loads = NodalLoadTable::default();
    for ((slot, node, dof), value) in nodal {
        if value != 0. {
            nodal_loads.pattern.push(pattern(slot));
            nodal_loads.node.push(node);
            nodal_loads.dof.push(dof);
            nodal_loads.value.push(value);
            nodal_loads.stage.push(0);
        }
    }
    let mut element_loads = ElementLoadTable::default();
    for (slot, wall, load) in wall_loads {
        element_loads.pattern.push(pattern(slot));
        element_loads.element_kind.push(ElementKind::StudWall);
        element_loads.element_index.push(wall);
        element_loads.load.push(ElementLoadSpec::WallLoad { load });
        element_loads.stage.push(0);
    }
    let floor_loads = floor_rows
        .into_iter()
        .map(|(slot, mut row)| {
            row.pattern = pattern(slot);
            row
        })
        .collect();

    let input = CarapaceInputV1 {
        header: Header {
            schema_version: 1,
            ndm: 3,
            engine_version: env!("CARGO_PKG_VERSION").into(),
            record_initial: false,
        },
        nodes: NodeTable {
            coords: g.coords.iter().flatten().copied().collect(),
            fixed: g.fixed,
            ..Default::default()
        },
        materials: g.materials,
        fibers: Default::default(),
        trusses: Default::default(),
        elastic_beam_columns_2d: Default::default(),
        elastic_beam_columns_3d: Default::default(),
        disp_beam_columns_2d: Default::default(),
        disp_beam_columns_3d: Default::default(),
        force_beam_columns_2d: Default::default(),
        force_beam_columns_3d: Default::default(),
        zero_lengths: g.zero_lengths,
        zero_length_sections: Default::default(),
        equal_dofs: equal,
        plane_materials: Vec::new(),
        triangles: Default::default(),
        quads: Default::default(),
        embedded_bars: Vec::new(),
        shell_sections: Vec::new(),
        shell3s: Default::default(),
        shell4s: Default::default(),
        stud_walls,
        diaphragms: Vec::new(),
        floor_loads,
        rigid_diaphragms: rigid,
        rigid_links: Default::default(),
        linear_constraints: linear,
        load_patterns,
        nodal_loads,
        element_loads,
        sequence: SequenceSpec::default(),
    };
    let storeys = spec
        .storeys
        .iter()
        .zip(&gens)
        .zip(levels)
        .zip(links_below)
        .map(|(((s, level), (master, nodes)), links_below)| StoreyMap {
            walls: s
                .walls
                .iter()
                .zip(level)
                .map(|(w, wg)| WallMap {
                    name: w.name.clone(),
                    element: wg.row,
                    corners: wg.corners,
                    top_nodes: stations(&wg.top),
                    bottom_nodes: stations(&wg.bottom),
                    pier_links: wg.pier_links.clone(),
                })
                .collect(),
            links_below,
            master,
            nodes,
        })
        .collect();
    Ok(Expansion {
        input,
        map: ExpansionMap {
            storeys,
            slot_patterns,
        },
    })
}

/// A flexible diaphragm's lateral load: the floor transfer statics between the wall lines parallel to the load,
/// each reaction lumped onto the wall's exterior top stations (summed into `nodal`).
fn distribute_flexible(
    g: &Gen,
    nodal: &mut BTreeMap<(usize, u32, u8), f64>,
    gens: &[WallGen],
    walls: &[WallSpec],
    load: &LateralLoad,
) -> Result<(), BuildError> {
    use carapace_core::model::{
        flexible_diaphragm_loads, lump_line_reaction, FloorBounds, FloorLoad, FloorSupport,
    };
    let [px, pz] = g.plan().map(usize::from);
    let plan = |n: u32| {
        let c = g.coords[n as usize];
        [c[px], c[pz]]
    };
    for (axis, force) in [(SpanAxis::X, load.fx), (SpanAxis::Z, load.fz)] {
        if force == 0. {
            continue;
        }
        let candidates: Vec<usize> = (0..walls.len())
            .filter(|&i| match axis {
                SpanAxis::X => walls[i].direction[0].abs() > 1. - 1e-9,
                SpanAxis::Z => walls[i].direction[1].abs() > 1. - 1e-9,
            })
            .collect();
        if candidates.is_empty() {
            return Err(err(
                "flexible diaphragm has no wall line perpendicular to lateral load",
            ));
        }
        let (a, b) = match axis {
            SpanAxis::X => (1, 0),
            SpanAxis::Z => (0, 1),
        };
        let pts: Vec<[f64; 2]> = candidates
            .iter()
            .flat_map(|&i| gens[i].top.values().map(|&n| plan(n)))
            .collect();
        let range = |k: usize| {
            (
                pts.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min),
                pts.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max),
            )
        };
        let ((s0, s1), (t0, t1)) = (range(a), range(b));
        if s1 <= s0 || t1 <= t0 {
            return Err(err(
                "flexible diaphragm plan bounds must have positive area",
            ));
        }
        let supports: Vec<_> = candidates
            .iter()
            .map(|&i| {
                let (p, q) = (plan(gens[i].corners[0]), plan(gens[i].corners[1]));
                FloorSupport {
                    s: (p[a] + q[a]) * 0.5,
                    t0,
                    t1,
                }
            })
            .collect();
        let reactions = flexible_diaphragm_loads(
            FloorBounds { s0, s1, t0, t1 },
            &supports,
            &[FloorLoad::Area {
                slot: load.slot,
                w: force / ((s1 - s0) * (t1 - t0)),
                region: None,
            }],
        )
        .map_err(|e| err(e.to_string()))?;
        let dof = match axis {
            SpanAxis::X => px as u8,
            SpanAxis::Z => pz as u8,
        };
        for (k, &i) in candidates.iter().enumerate() {
            let mut nodes: Vec<(f64, u32)> =
                gens[i].top.values().map(|&n| (plan(n)[b], n)).collect();
            nodes.sort_by(|p, q| p.0.total_cmp(&q.0));
            nodes.dedup_by(|p, q| (p.0 - q.0).abs() < 1e-9);
            for line in &reactions.per_support[k].lines {
                let vals =
                    lump_line_reaction(&nodes.iter().map(|x| x.0).collect::<Vec<_>>(), *line)
                        .map_err(|e| err(e.to_string()))?;
                for ((_, n), v) in nodes.iter().zip(vals) {
                    *nodal.entry((load.slot, *n, dof)).or_default() += v;
                }
            }
        }
    }
    Ok(())
}
