use std::collections::{BTreeMap, HashMap};

use super::expand::{expand, ExpandOptions, Expansion, StationNode};
use super::spec::*;
use crate::input_v1::sequence::{AlgorithmSpec, IntegratorSpec, StageSpec};
use crate::input_v1::tables::ElementKind;
use crate::input_v1::{decode, Session};
use carapace_core::model::{
    transfer_floor, Domain3, Element3Id, FloorBounds, FloorLoad, FloorReactions, FloorSupport,
    LoadPatternId, Node3Id, StudWallLayout,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildError(pub String);
impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for BuildError {}
pub(super) fn err(s: impl Into<String>) -> BuildError {
    BuildError(s.into())
}

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
    /// One load pattern per slot that carries a load (the slot is the index).
    patterns: [Option<LoadPatternId>; SLOTS],
}

/// Number of load slots of a building (`BuildingSpec::combination`).
pub(super) const SLOTS: usize = 8;

impl Building {
    /// Apply slot combination factors: slot `k`'s load pattern (one per slot that carries a load; the others are
    /// ignored) is scaled by `factors[k]`. Set them before building an analysis.
    pub fn apply_combination(&mut self, factors: [f64; SLOTS]) -> Result<(), BuildError> {
        if factors.iter().any(|v| !v.is_finite()) {
            return Err(err("combination factors must be finite"));
        }
        for (pattern, factor) in self.patterns.iter().zip(factors) {
            if let Some(pattern) = pattern {
                self.domain.set_load_pattern_scale_factor(*pattern, factor);
            }
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
    let Expansion { mut input, map } = expand(spec, &ExpandOptions::default())?;
    // One static stage so the decoder accepts the (stage 0) load rows; `into_loaded` registers them.
    input.sequence.stages.push(StageSpec::Static {
        id: "build".into(),
        steps: 1,
        integrator: IntegratorSpec::LoadControl { increment: 1. },
        algorithm: AlgorithmSpec::Linear,
        convergence: None,
        hold_patterns_after: Vec::new(),
    });
    let Session::D3(session) = decode(input).map_err(|e| err(format!("{e:?}")))? else {
        unreachable!("expansions are 3D")
    };
    let (domain, nodes, located) = session.into_loaded();
    let node = |i: u32| nodes[i as usize];
    let rows = |kind: fn(&ElementKind) -> bool| -> HashMap<u32, Element3Id> {
        located
            .iter()
            .filter(|(_, (k, _))| kind(k))
            .map(|(id, (_, row))| (*row, *id))
            .collect()
    };
    let walls_by_row = rows(|k| matches!(k, ElementKind::StudWall));
    let links_by_row = rows(|k| matches!(k, ElementKind::ZeroLength));
    let pattern_ids = domain.load_pattern_ids();
    let mut patterns = [None; SLOTS];
    for (slot, row) in map.slot_patterns.iter().enumerate() {
        patterns[slot] = row.map(|r| pattern_ids[r as usize + 1]);
    }
    let mut walls = Vec::new();
    let mut links = vec![Vec::new(); spec.storeys.len()];
    let mut levels = Vec::new();
    let mut floor_reactions = Vec::new();
    for (si, (s, m)) in spec.storeys.iter().zip(&map.storeys).enumerate() {
        let mut level = BTreeMap::new();
        for (w, wm) in s.walls.iter().zip(&m.walls) {
            let station = |v: &[StationNode]| v.iter().map(|p| (p.x, node(p.node))).collect();
            level.insert(
                w.name.clone(),
                BuildingWall {
                    start: w.start,
                    direction: w.direction,
                    length: w.layout.length,
                    element: walls_by_row[&wm.element],
                    corners: wm.corners.map(node),
                    top_nodes: station(&wm.top_nodes),
                    bottom_nodes: station(&wm.bottom_nodes),
                    pier_links: wm
                        .pier_links
                        .iter()
                        .map(|p| (p.x, links_by_row[&p.element]))
                        .collect(),
                },
            );
        }
        walls.push(level);
        if si > 0 {
            links[si - 1] = m.links_below.iter().map(|r| links_by_row[r]).collect();
        }
        levels.push(LevelHandle {
            master: m.master.map(node),
            nodes: m.nodes.iter().map(|&n| node(n)).collect(),
        });
        for (fi, f) in s.floors.iter().enumerate() {
            floor_reactions.push(
                floor_transfer_for(f, &s.walls)
                    .map_err(|e| err(format!("storey {si} floor {fi}: {e}")))?
                    .0,
            );
        }
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

pub(super) fn convert_layout(
    l: &WallLayoutSpec,
    si: usize,
    name: &str,
) -> Result<StudWallLayout, BuildError> {
    let x = l
        .to_core()
        .map_err(|e| err(format!("storey {si} wall {name}: {e}")))?;
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

pub(super) fn floor_transfer_for(
    f: &FloorSpec,
    walls: &[WallSpec],
) -> Result<(FloorReactions, FloorBounds, Vec<FloorSupport>), BuildError> {
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
    let reactions = transfer_floor(b, &supports, &loads, f.continuity.to_core())
        .map_err(|e| err(e.to_string()))?;
    Ok((reactions, b, supports))
}
