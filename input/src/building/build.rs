use std::collections::BTreeMap;

use super::spec::*;
use carapace_core::model::{
    flexible_diaphragm_loads, lump_line_reaction, transfer_floor, AnchorageBearing,
    ChordDistribution, Domain3, Element3, Element3Id, ElementLoad3, FloorBounds, FloorContinuity,
    FloorLoad, FloorReactions, FloorSupport, FloorZone, LoadPatternId, LoadSeries, Material, Node3,
    Node3Id, StudWall, StudWallLayout, WallChord, WallLoad, WallOpening, WallShearLaw,
    WallSheathing, WallStation, WallStud, ZeroLength3, WALL_SLOTS,
};
use slotmap::Key;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildError(pub String);
impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for BuildError {}
fn err(s: impl Into<String>) -> BuildError {
    BuildError(s.into())
}
type GeneratedWallNodes = (
    WallSpec,
    StudWallLayout,
    [Node3Id; 4],
    BTreeMap<i64, Node3Id>,
    BTreeMap<i64, Node3Id>,
);

#[derive(Debug, Clone)]
pub struct BuildingWall {
    pub start: [f64; 2],
    pub direction: [f64; 2],
    pub length: f64,
    pub element: Element3Id,
    pub corners: [Node3Id; 4],
    pub top_nodes: Vec<(f64, Node3Id)>,
    pub bottom_nodes: Vec<(f64, Node3Id)>,
    pub pier_links: Vec<(f64, Element3Id)>,
}
#[derive(Debug, Clone)]
pub struct LevelHandle {
    pub master: Option<Node3Id>,
    pub nodes: Vec<Node3Id>,
}
#[derive(Debug, Clone, PartialEq)]
pub enum StoreyDrift {
    Rigid { ux: f64, uz: f64, ry: f64 },
    Flexible { wall_drift: BTreeMap<String, f64> },
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoldDownForce {
    pub x: f64,
    pub pier_m_jd: f64,
    pub link_force: f64,
    pub link: Element3Id,
}
pub struct Building {
    pub domain: Domain3,
    pub walls: Vec<BTreeMap<String, BuildingWall>>,
    /// Vertical links (anchorage at pier ends, bearing elsewhere) joining each storey's walls to the walls below.
    pub links: Vec<Vec<Element3Id>>,
    pub levels: Vec<LevelHandle>,
    pub floor_reactions: Vec<FloorReactions>,
    patterns: [LoadPatternId; WALL_SLOTS],
}

impl Building {
    /// Apply slot combination factors to the builder's eight independent load patterns.
    /// Set the eight pattern factors before building an analysis.
    pub fn apply_combination(&mut self, factors: [f64; WALL_SLOTS]) -> Result<(), BuildError> {
        if factors.iter().any(|v| !v.is_finite()) {
            return Err(err("combination factors must be finite"));
        }
        for (pattern, factor) in self.patterns.iter().zip(factors) {
            self.domain.set_load_pattern_scale_factor(*pattern, factor);
        }
        Ok(())
    }
    pub fn wall_response(
        &self,
        domain: &Domain3,
        storey: usize,
        name: &str,
    ) -> Option<carapace_core::model::StudWallResponse> {
        let w = self.walls.get(storey)?.get(name)?;
        domain.element_stud_wall_response(w.element, 1.)
    }
    pub fn storey_drift(&self, domain: &Domain3, storey: usize) -> Option<StoreyDrift> {
        let level = self.levels.get(storey)?;
        if let Some(master) = level.master {
            let u = domain.node(master).displacement;
            let below = if storey == 0 {
                [0.; 6]
            } else {
                self.levels
                    .get(storey - 1)?
                    .master
                    .map(|n| domain.node(n).displacement)
                    .unwrap_or([0.; 6])
            };
            return Some(StoreyDrift::Rigid {
                ux: u[0] - below[0],
                uz: u[2] - below[2],
                ry: u[4] - below[4],
            });
        }
        let mut wall_drift = BTreeMap::new();
        for (name, w) in self.walls.get(storey)? {
            let [ex, ez] = w.direction;
            let projected = |nodes: &[Node3Id]| {
                let mut values = Vec::new();
                for &n in nodes {
                    let c = domain.node(n).coords;
                    let dx = c[0] - w.start[0];
                    let dz = c[2] - w.start[1];
                    let t = dx * ex + dz * ez;
                    if t >= -1e-7 && t <= w.length + 1e-7 && (dx * ez - dz * ex).abs() < 1e-7 {
                        values.push(
                            domain.node(n).displacement[0] * ex
                                + domain.node(n).displacement[2] * ez,
                        );
                    }
                }
                values
            };
            let top_values = projected(&level.nodes);
            let bottom_values = if storey == 0 {
                vec![
                    domain.node(w.corners[0]).displacement[0] * ex
                        + domain.node(w.corners[0]).displacement[2] * ez,
                    domain.node(w.corners[1]).displacement[0] * ex
                        + domain.node(w.corners[1]).displacement[2] * ez,
                ]
            } else {
                projected(&self.levels.get(storey - 1)?.nodes)
            };
            if top_values.is_empty() || bottom_values.is_empty() {
                continue;
            }
            let bottom = bottom_values.iter().sum::<f64>() / bottom_values.len() as f64;
            wall_drift.insert(
                name.clone(),
                top_values.iter().sum::<f64>() / top_values.len() as f64 - bottom,
            );
        }
        Some(StoreyDrift::Flexible { wall_drift })
    }
    pub fn holddown_forces(
        &self,
        domain: &Domain3,
        storey: usize,
        name: &str,
    ) -> Option<Vec<HoldDownForce>> {
        let r = self.wall_response(domain, storey, name)?;
        let links = &self.walls.get(storey)?.get(name)?.pier_links;
        let mut out = Vec::new();
        for p in &r.piers {
            for (x, value) in [(p.xa, p.chords[0]), (p.xb, p.chords[1])] {
                let link = links
                    .iter()
                    .find(|(at, _)| (*at - x).abs() < 1e-7)
                    .map(|x| x.1)?;
                let f = domain.element_local_force(link);
                out.push(HoldDownForce {
                    x,
                    pier_m_jd: value,
                    link_force: f.as_slice().get(1).copied().unwrap_or(0.).abs(),
                    link,
                });
            }
        }
        Some(out)
    }
}

pub fn build(spec: &BuildingSpec) -> Result<Building, BuildError> {
    if spec.storeys.is_empty() {
        return Err(err("building must have at least one storey"));
    }
    anchorage(&spec.foundation_anchorage)?;
    let mut domain = Domain3::new();
    let mut wall_nodes: Vec<BTreeMap<String, GeneratedWallNodes>> = Vec::new();
    let mut floor_reactions = Vec::new();
    let patterns =
        std::array::from_fn(|_| domain.add_load_pattern(LoadSeries::Linear { slope: 1. }));
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
        let mut level = BTreeMap::new();
        for (wi, w) in s.walls.iter().enumerate() {
            let len = w.layout.length;
            let (ex, ez) = (w.direction[0], w.direction[1]);
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
            if level.contains_key(&w.name) {
                return Err(err(format!("storey {si}: duplicate wall {}", w.name)));
            }
            if w.layout.floor.is_some() {
                return Err(err(format!(
                    "storey {si} wall {}: the floor zone comes from the storey (floorDepth, bandGa, rimEa), not the wall layout",
                    w.name
                )));
            }
            let mut layout = convert_layout(&w.layout, si, &w.name)?;
            if s.floor_depth > 0. {
                layout.floor = Some(FloorZone {
                    depth: s.floor_depth,
                    ga: s.band_ga,
                    rim_ea: s.rim_ea,
                });
            }
            if si > 0 {
                let below = &spec.storeys[si - 1].walls;
                let matching = below.iter().any(|b| {
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
            let y0: f64 = spec.storeys[..si]
                .iter()
                .map(|p| p.wall_height + p.floor_depth)
                .sum();
            // The exterior top row is the floor zone's top: `H + depth` above the base.
            let y_top = y0 + s.wall_height + s.floor_depth;
            let corners = [
                domain.add_node(Node3::new([w.start[0], y0, w.start[1]])),
                domain.add_node(Node3::new([
                    w.start[0] + ex * len,
                    y0,
                    w.start[1] + ez * len,
                ])),
                domain.add_node(Node3::new([
                    w.start[0] + ex * len,
                    y_top,
                    w.start[1] + ez * len,
                ])),
                domain.add_node(Node3::new([w.start[0], y_top, w.start[1]])),
            ];
            let mut xs: Vec<f64> = layout.studs.iter().map(|x| x.x).collect();
            xs.extend(layout.pier_ends());
            xs.sort_by(f64::total_cmp);
            xs.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * len);
            let top: BTreeMap<i64, Node3Id> = xs
                .iter()
                .map(|&x| {
                    (
                        (x * 1e6).round() as i64,
                        if x == 0. {
                            corners[3]
                        } else if x == len {
                            corners[2]
                        } else {
                            domain.add_node(Node3::new([
                                w.start[0] + ex * x,
                                y_top,
                                w.start[1] + ez * x,
                            ]))
                        },
                    )
                })
                .collect();
            // Above the first storey the bottom stations stand on the wall below: one at each of its top
            // positions (the corners and every stud), where the vertical links go. Pier ends (hold-downs) must
            // be among them.
            let mut bxs: Vec<f64> = layout.studs.iter().map(|x| x.x).collect();
            bxs.extend(layout.pier_ends());
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
            let bottom: BTreeMap<i64, Node3Id> = bxs
                .iter()
                .map(|&x| {
                    (
                        (x * 1e6).round() as i64,
                        if x == 0. {
                            corners[0]
                        } else if x == len {
                            corners[1]
                        } else {
                            domain.add_node(Node3::new([
                                w.start[0] + ex * x,
                                y0,
                                w.start[1] + ez * x,
                            ]))
                        },
                    )
                })
                .collect();
            let _ = wi;
            level.insert(w.name.clone(), (w.clone(), layout, corners, top, bottom));
        }
        wall_nodes.push(level);
    }
    // Resolve gravity supports by wall name and accumulate their reactions as local top-plate loads.
    for (si, s) in spec.storeys.iter().enumerate() {
        for (fi, f) in s.floors.iter().enumerate() {
            let (reactions, supports) = floor_transfer_for(f, &s.walls)
                .map_err(|e| err(format!("storey {si} floor {fi}: {e}")))?;
            for (i, name) in f.supports.iter().enumerate() {
                let wall = wall_nodes[si].get_mut(name).ok_or_else(|| {
                    err(format!(
                        "storey {si} floor {fi}: unknown support wall {name}"
                    ))
                })?;
                let ws = &wall.0;
                let t0 = match f.span_axis {
                    SpanAxis::X => f.bounds[2],
                    SpanAxis::Z => f.bounds[0],
                };
                let (coord, sign) = match f.span_axis {
                    SpanAxis::X => (ws.start[1], ws.direction[1]),
                    SpanAxis::Z => (ws.start[0], ws.direction[0]),
                };
                let x_at_t0 = (t0 - coord) / sign;
                let loads = reactions.per_support[i].wall_loads(x_at_t0, sign < 0.);
                for (slot, group) in loads.into_iter().enumerate() {
                    for load in group {
                        wall.0.loads.push(SlotWallLoad {
                            slot,
                            load: from_wall_load(load),
                        });
                    }
                }
                let _ = &supports[i];
            }
            floor_reactions.push(reactions);
        }
    }
    // Foundation supports at all pier ends. The ground's plan translations are fixed; vertical behavior is
    // carried by the anchorage material.
    let mut walls: Vec<BTreeMap<String, BuildingWall>> = vec![BTreeMap::new(); spec.storeys.len()];
    let mut links = vec![Vec::new(); spec.storeys.len()];
    for si in 0..spec.storeys.len() {
        for (name, (ws, layout, corners, top, bottom)) in wall_nodes[si].clone() {
            let mut load_slots: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
            for l in ws.loads {
                if l.slot >= WALL_SLOTS {
                    return Err(err(format!(
                        "storey {si} wall {name}: load slot out of range"
                    )));
                }
                load_slots[l.slot].push(convert_wall_load(l.load));
            }
            let mut sorted: Vec<f64> = top.keys().map(|k| *k as f64 / 1e6).collect();
            sorted.sort_by(f64::total_cmp);
            sorted.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * layout.length);
            let stations = sorted
                .iter()
                .filter(|&&x| x > 0. && x < layout.length)
                .map(|&x| WallStation {
                    chord: WallChord::Top,
                    x,
                    node: top[&key(x)],
                })
                .chain(if si > 0 {
                    bottom
                        .iter()
                        .map(|(k, &node)| (*k as f64 / 1e6, node))
                        .filter(|(x, _)| *x > 0. && *x < layout.length)
                        .map(|(x, node)| WallStation {
                            chord: WallChord::Bottom,
                            x,
                            node,
                        })
                        .collect::<Vec<_>>()
                } else {
                    layout
                        .pier_ends()
                        .into_iter()
                        .filter(|x| *x > 0. && *x < layout.length)
                        .map(|x| WallStation {
                            chord: WallChord::Bottom,
                            x,
                            node: bottom[&key(x)],
                        })
                        .collect::<Vec<_>>()
                })
                .collect();
            let wall = StudWall::new3(corners, stations, layout.clone(), load_slots)
                .map_err(|e| err(format!("storey {si} wall {name}: {e}")))?;
            let eid = domain.add_element(Element3::StudWall(wall));
            for (slot, pattern) in patterns.iter().enumerate() {
                let mut factors = [0.; WALL_SLOTS];
                factors[slot] = 1.;
                domain.add_element_load(
                    *pattern,
                    eid,
                    ElementLoad3 {
                        wall: factors,
                        ..Default::default()
                    },
                );
            }
            let mut links = Vec::new();
            if si == 0 {
                for x in layout.pier_ends() {
                    let k = key(x);
                    let node = bottom[&k];
                    let (tension, slack, bearing) = anchorage(&spec.foundation_anchorage)?;
                    let ground =
                        domain.add_node(Node3::new(domain.node(node).coords).fix(0).fix(1).fix(2));
                    let mat = Material::anchorage(tension, slack, bearing)
                        .map_err(|e| err(e.to_string()))?;
                    let z = ZeroLength3::new(ground, node)
                        .with_material(0, Material::Elastic { e: 1e12 })
                        .with_material(1, mat)
                        .with_material(2, Material::Elastic { e: 1e12 });
                    links.push((x, domain.add_element(Element3::ZeroLength3(z))));
                    let _ = k;
                }
            }
            walls[si].insert(
                name,
                BuildingWall {
                    start: ws.start,
                    direction: ws.direction,
                    length: layout.length,
                    element: eid,
                    corners,
                    top_nodes: top.iter().map(|(k, n)| (*k as f64 / 1e6, *n)).collect(),
                    bottom_nodes: bottom.iter().map(|(k, n)| (*k as f64 / 1e6, *n)).collect(),
                    pier_links: links,
                },
            );
        }
    }
    // Join each wall to the wall below: its bottom nodes sit on the lower wall's exterior top nodes (the floor
    // zone's top row), tied in the wall plane and linked vertically, with the hold-down anchorage at pier ends
    // and bearing (compression only) elsewhere.
    for si in 1..spec.storeys.len() {
        for (name, (ws, _, _, _, _)) in &wall_nodes[si] {
            let low_top = walls[si - 1]
                .get(name)
                .ok_or_else(|| {
                    err(format!(
                        "storey {} wall {name}: missing lower handle",
                        si - 1
                    ))
                })?
                .top_nodes
                .clone();
            let up = walls[si]
                .get(name)
                .ok_or_else(|| err(format!("storey {si} wall {name}: missing upper handle")))?
                .clone();
            let pier = wall_nodes[si][name].1.pier_ends();
            let dofs: Vec<usize> = [(0, ws.direction[0]), (2, ws.direction[1])]
                .into_iter()
                .filter(|(_, c)| c.abs() > 1e-9)
                .map(|(d, _)| d)
                .collect();
            for &(x, upper) in &up.bottom_nodes {
                let (_, lower) = low_top
                    .iter()
                    .find(|(a, _)| (*a - x).abs() < 1e-8)
                    .ok_or_else(|| {
                        err(format!(
                            "storey {si} wall {name}: no node below at x = {x} to link to"
                        ))
                    })?;
                let is_pier = pier.iter().any(|p| (*p - x).abs() < 1e-8);
                let mat = if is_pier {
                    let (t, sl, b) = anchorage(&spec.storeys[si].anchorage)?;
                    Material::anchorage(t, sl, b)
                } else {
                    Material::anchorage(
                        0.,
                        0.,
                        AnchorageBearing::Linear {
                            k: spec.storeys[si].bearing,
                        },
                    )
                }
                .map_err(|e| err(e.to_string()))?;
                let id = domain.add_element(Element3::ZeroLength3(
                    ZeroLength3::new(*lower, upper).with_material(1, mat),
                ));
                domain.equal_dof(*lower, upper, &dofs);
                links[si - 1].push(id);
                if is_pier {
                    walls[si].get_mut(name).unwrap().pier_links.push((x, id));
                }
            }
        }
    }
    let mut levels = Vec::new();
    let mut nodal = NodalLoads::new();
    for (si, s) in spec.storeys.iter().enumerate() {
        let mut nodes: Vec<Node3Id> = walls[si]
            .values()
            .flat_map(|w| w.top_nodes.iter().map(|x| x.1))
            .collect();
        nodes.sort_by_key(|n| n.data().as_ffi());
        nodes.dedup();
        let master = if spec.diaphragm == DiaphragmKind::Rigid {
            if nodes.is_empty() {
                return Err(err(format!("storey {si}: no diaphragm nodes")));
            }
            let xyz = nodes.iter().fold([0.; 3], |mut a, n| {
                for (i, v) in domain.node(*n).coords.iter().enumerate() {
                    a[i] += v / nodes.len() as f64;
                }
                a
            });
            let m = domain.add_node(Node3::new(xyz));
            domain.rigid_diaphragm(m, &nodes);
            Some(m)
        } else {
            None
        };
        for l in &s.lateral {
            if l.slot >= WALL_SLOTS || ![l.fx, l.fz, l.at[0], l.at[1]].iter().all(|v| v.is_finite())
            {
                return Err(err(format!("storey {si}: invalid lateral load")));
            }
            // Rigid load is attached to the retained node; flexible distribution is handled by floor helpers.
            if let Some(m) = master {
                let pattern = patterns[l.slot];
                *nodal.entry((pattern, m, 0)).or_default() += l.fx;
                *nodal.entry((pattern, m, 2)).or_default() += l.fz;
                let c = domain.node(m).coords;
                // Moment about +Y of the force at `at` about the master: r x F with r = (dx, 0, dz),
                // matching `rigid_diaphragm` (a rotation `θ` moves a point by `θ dz` in x and `-θ dx` in z).
                let my = l.fx * (l.at[1] - c[2]) - l.fz * (l.at[0] - c[0]);
                *nodal.entry((pattern, m, 4)).or_default() += my;
            } else {
                distribute_flexible(&domain, &mut nodal, &walls[si], l, patterns[l.slot])?;
            }
        }
        if master.is_none() {
            for w in walls[si].values() {
                for &n in &w.top_nodes.iter().map(|p| p.1).collect::<Vec<_>>() {
                    let (ex, ez) = (w.direction[0], w.direction[1]);
                    if ez.abs() < 1e-10 {
                        domain.fix_node_dof(n, 2);
                    } else if ex.abs() < 1e-10 {
                        domain.fix_node_dof(n, 0);
                    } else if ex.abs() >= ez.abs() {
                        domain.add_constraint((n, 2), &[(n, 0, ez / ex)]);
                    } else {
                        domain.add_constraint((n, 0), &[(n, 2, ex / ez)]);
                    }
                }
            }
        }
        let _ = s;
        levels.push(LevelHandle { master, nodes });
        let _ = si;
    }
    for ((pattern, node, dof), v) in nodal {
        domain.add_nodal_load(pattern, node, dof, v);
    }
    Ok(Building {
        domain,
        walls,
        links,
        levels,
        floor_reactions,
        patterns,
    })
}

fn key(x: f64) -> i64 {
    (x * 1e6).round() as i64
}
fn anchorage(a: &AnchorageSpec) -> Result<(f64, f64, AnchorageBearing), BuildError> {
    let b = match a.bearing {
        AnchorageBearingSpec::Linear { k } => AnchorageBearing::Linear { k },
        AnchorageBearingSpec::O86Crush { qr } => AnchorageBearing::O86Crush { qr },
    };
    if a.tension < 0. || a.slack < 0. || !a.tension.is_finite() || !a.slack.is_finite() {
        return Err(err("invalid anchorage values"));
    }
    Ok((a.tension, a.slack, b))
}
fn convert_layout(l: &WallLayoutSpec, si: usize, name: &str) -> Result<StudWallLayout, BuildError> {
    let sheathing = l.sheathing.as_ref().map(|s| WallSheathing {
        law: match s.law {
            WallShearLawSpec::Linear { ga } => WallShearLaw::Linear { ga },
            WallShearLawSpec::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            } => WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            },
        },
    });
    let x = StudWallLayout {
        length: l.length,
        height: l.height,
        studs: l
            .studs
            .iter()
            .map(|s| WallStud { x: s.x, ea: s.ea })
            .collect(),
        top_ea: l.top_ea,
        top_ei: l.top_ei,
        top_distribution: match l.top_distribution {
            ChordDistributionSpec::Beam => ChordDistribution::Beam,
            ChordDistributionSpec::Direct => ChordDistribution::Direct,
        },
        bottom_ea: l.bottom_ea,
        bottom_ei: l.bottom_ei,
        openings: l
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
        sheathing,
        foundation: si == 0,
        floor: None,
    };
    if [l.length, l.height, l.top_ea, l.bottom_ea, l.bottom_ei]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.)
    {
        return Err(err(format!(
            "storey {si} wall {name}: non-finite or negative layout values"
        )));
    }
    Ok(x)
}
fn convert_wall_load(l: WallLoadSpec) -> WallLoad {
    match l {
        WallLoadSpec::TopPoint { x, fx, fy } => WallLoad::TopPoint { x, fx, fy },
        WallLoadSpec::TopLine {
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
fn from_wall_load(l: WallLoad) -> WallLoadSpec {
    match l {
        WallLoad::TopPoint { x, fx, fy } => WallLoadSpec::TopPoint { x, fx, fy },
        WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => WallLoadSpec::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        },
    }
}

/// Nodal loads summed per (pattern, node, dof) before they go into the domain: `Domain::add_nodal_load`
/// replaces an existing value for the same node and DOF rather than adding to it.
type NodalLoads = BTreeMap<(LoadPatternId, Node3Id, usize), f64>;

fn distribute_flexible(
    domain: &Domain3,
    nodal: &mut NodalLoads,
    walls: &BTreeMap<String, BuildingWall>,
    load: &LateralLoad,
    pattern: LoadPatternId,
) -> Result<(), BuildError> {
    for (axis, force) in [(SpanAxis::X, load.fx), (SpanAxis::Z, load.fz)] {
        if force == 0. {
            continue;
        }
        let mut candidates: Vec<_> = walls
            .values()
            .filter(|w| match axis {
                SpanAxis::X => w.direction[0].abs() > 1. - 1e-9,
                SpanAxis::Z => w.direction[1].abs() > 1. - 1e-9,
            })
            .collect();
        if candidates.is_empty() {
            return Err(err(
                "flexible diaphragm has no wall line perpendicular to lateral load",
            ));
        }
        let pts: Vec<_> = candidates
            .iter()
            .flat_map(|w| w.top_nodes.iter().map(|(_, n)| domain.node(*n).coords))
            .collect();
        let (s0, s1, t0, t1) = match axis {
            SpanAxis::X => (
                pts.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min),
                pts.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max),
                pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max),
            ),
            SpanAxis::Z => (
                pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max),
                pts.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min),
                pts.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max),
            ),
        };
        if s1 <= s0 || t1 <= t0 {
            return Err(err(
                "flexible diaphragm plan bounds must have positive area",
            ));
        }
        let supports: Vec<_> = candidates
            .iter()
            .map(|w| {
                let a = domain.node(w.corners[0]).coords;
                let b = domain.node(w.corners[1]).coords;
                let s = match axis {
                    SpanAxis::X => (a[2] + b[2]) * 0.5,
                    SpanAxis::Z => (a[0] + b[0]) * 0.5,
                };
                FloorSupport { s, t0, t1 }
            })
            .collect();
        let bounds = FloorBounds { s0, s1, t0, t1 };
        let reactions = flexible_diaphragm_loads(
            bounds,
            &supports,
            &[FloorLoad::Area {
                slot: load.slot,
                w: force / ((s1 - s0) * (t1 - t0)),
                region: None,
            }],
        )
        .map_err(|e| err(e.to_string()))?;
        for (i, w) in candidates.drain(..).enumerate() {
            let mut nodes: Vec<_> = w
                .top_nodes
                .iter()
                .map(|(_, n)| {
                    let c = domain.node(*n).coords;
                    let t = match axis {
                        SpanAxis::X => c[0],
                        SpanAxis::Z => c[2],
                    };
                    (t, *n)
                })
                .collect();
            nodes.sort_by(|a, b| a.0.total_cmp(&b.0));
            nodes.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9);
            for line in &reactions.per_support[i].lines {
                let vals =
                    lump_line_reaction(&nodes.iter().map(|x| x.0).collect::<Vec<_>>(), *line)
                        .map_err(|e| err(e.to_string()))?;
                let dof = match axis {
                    SpanAxis::X => 0,
                    SpanAxis::Z => 2,
                };
                for ((_, n), v) in nodes.iter().zip(vals) {
                    *nodal.entry((pattern, *n, dof)).or_default() += v;
                }
            }
        }
    }
    Ok(())
}
fn floor_transfer_for(
    f: &FloorSpec,
    walls: &[WallSpec],
) -> Result<(FloorReactions, Vec<FloorSupport>), BuildError> {
    // `bounds` is the plan rectangle `[x0, x1, z0, z1]`; `s` runs along the span, `t` across it.
    let [x0, x1, z0, z1] = f.bounds;
    let b = match f.span_axis {
        SpanAxis::X => FloorBounds {
            s0: x0,
            s1: x1,
            t0: z0,
            t1: z1,
        },
        SpanAxis::Z => FloorBounds {
            s0: z0,
            s1: z1,
            t0: x0,
            t1: x1,
        },
    };
    let supports: Vec<_> = f
        .supports
        .iter()
        .enumerate()
        .map(|(i, n)| {
            let wall = walls
                .iter()
                .find(|w| &w.name == n)
                .ok_or_else(|| err(format!("floor support {n}: unknown wall")))?;
            let axis = match f.span_axis {
                SpanAxis::X => [1., 0.],
                SpanAxis::Z => [0., 1.],
            };
            if (wall.direction[0] * axis[0] + wall.direction[1] * axis[1]).abs() > 1e-9 {
                return Err(err(format!(
                    "floor support {n}: wall is not perpendicular to span"
                )));
            }
            let axis_s = match f.span_axis {
                SpanAxis::X => wall.start[0],
                SpanAxis::Z => wall.start[1],
            };
            let t0 = match f.span_axis {
                SpanAxis::X => f.bounds[2],
                SpanAxis::Z => f.bounds[0],
            };
            let t1 = match f.span_axis {
                SpanAxis::X => f.bounds[3],
                SpanAxis::Z => f.bounds[1],
            };
            let _ = i;
            Ok(FloorSupport { s: axis_s, t0, t1 })
        })
        .collect::<Result<_, BuildError>>()?;
    let loads = f
        .loads
        .iter()
        .map(|l| match *l {
            FloorLoadSpec::Area { slot, w, region } => FloorLoad::Area { slot, w, region },
            FloorLoadSpec::Line { slot, w, s, t0, t1 } => FloorLoad::Line { slot, w, s, t0, t1 },
            FloorLoadSpec::LineAlong { slot, w, t, s0, s1 } => {
                FloorLoad::LineAlong { slot, w, t, s0, s1 }
            }
            FloorLoadSpec::Point { slot, p, s, t } => FloorLoad::Point { slot, p, s, t },
        })
        .collect::<Vec<_>>();
    let reactions = transfer_floor(
        b,
        &supports,
        &loads,
        match f.continuity {
            FloorContinuitySpec::Simple => FloorContinuity::Simple,
            FloorContinuitySpec::Continuous => FloorContinuity::Continuous,
        },
    )
    .map_err(|e| err(e.to_string()))?;
    Ok((reactions, supports))
}
