//! Stud-wall macro-element: one storey of a 2D light-frame wall (studs, plates, openings with headers, jacks
//! and cripples, optional sheathing), condensed onto its four corners and optional plate stations.
//!
//! The framing is linear elastic. Each full-height pier carries one constant shear field whose law is linear
//! (`ga`), the CSA O86:19 shear-plus-nail-slip law, or any uniaxial [`Material`]. Linear walls are condensed
//! once at construction; nonlinear walls condense the bare frame and solve a small per-pier system at each
//! state. Exterior DOFs are `ux, uy` only. An optional [`FloorZone`] puts the floor depth above the framing
//! (exterior top nodes at `y = H + depth`) and condenses it with the rest.
//!
//! User guide (layout, laws, supports, loads, analysis, results, assumptions, limitations):
//! `docs/stud-wall.md`. Implementation spec and verification: `docs/stud-wall-v1.md`.
//!
//! This file holds the layout types and the element interface (DOF masks, assembly, commit, forces); `build`
//! assembles and condenses the internal frame, `validate` checks the layout, `loads` discretizes top-plate
//! loads, `nonlinear` solves the per-pier system of a nonlinear law and `response` recovers member forces. The
//! shear laws are shared with diaphragm bays (`super::shear_law`), the floor-zone bays are `super::floor_zone`.

use arrayvec::ArrayVec;
use nalgebra::{DMatrix, DVector};

use super::floor_zone::{bay_response, bay_stiffness};
use super::shear_law::{shear_law_response, WallShearLaw};
use super::wall_plane::PlaneMap;
use super::{DofMask, NodeView, TangentSink, VectorSink};
use crate::model::{ElementForce, ElementLoad, Material, Node3Id, NodeId, NodeKey, WallCaseLoad};

mod build;
mod loads;
mod nonlinear;
mod response;
mod validate;

use build::*;
use loads::*;
use nonlinear::*;
pub use response::{
    StudWallResponse, WallAxialResponse, WallBaseResponse, WallChordResponse, WallFloorBay,
    WallHeaderResponse, WallHoldDownResponse, WallMemberResponse, WallPierResponse,
};
use validate::*;

/// Prepared assembly supports up to 256 exterior nodes (512 translational DOFs).
const MAX_EXTERIOR_NODES: usize = 256;
const MAX_EXTERIOR_DOFS: usize = MAX_EXTERIOR_NODES * 3;
/// Maximum number of piers in a wall with a nonlinear sheathing law (`O86` or `Material`): the local pier
/// solve uses fixed-size stack buffers.
pub const MAX_WALL_PIERS: usize = 32;
thread_local! {
    /// Dense tangent of a nonlinear wall (`(m n)²`, up to `MAX_EXTERIOR_DOFS²`): too large for the stack, so each
    /// thread reuses one buffer, which allocates only when it first grows.
    static NONLINEAR_TANGENT: std::cell::RefCell<Vec<f64>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Which plate a [`WallStation`] sits on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WallChord {
    /// Top plate, at local `y = H`.
    Top,
    /// Bottom plate, at local `y = 0`.
    Bottom,
    /// A base node the wall stands on (with [`StudWallLayout::bearing`]): on the bottom line (`y = 0`) at any
    /// `x`, also outside `0..L`. Not part of the wall's framing; see [`WallBearing`].
    Base,
}
/// How the top plate distributes load to the studs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChordDistribution {
    /// Continuous Euler–Bernoulli beam over the studs (`top_ei` required): continuity reactions.
    Beam,
    /// Axial bar hinged at every stud and cripple: simple-span (lever-rule) tributary distribution, the usual
    /// hand takedown. Top stations must sit at a stud or cripple.
    Direct,
}
/// One stud row: an axial bar from the bottom to the top plate, pinned at both. A pack of plies is one row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallStud {
    /// Local position along the wall, `0..=L`. The first row is at 0 and the last at `L` (the end posts).
    pub x: f64,
    /// Axial stiffness `EA` of the row (sum of its plies).
    pub ea: f64,
}
/// A window or door. Jacks stand at `x0` and `x1`, the header spans between them, simply supported on the
/// jacks, and stud rows strictly inside become cripples above the header. Nothing is modelled below a window.
/// A stud row exactly at `x0` or `x1` is the king stud: king and jack share one bottom-plate node (one end
/// post), and the adjacent pier ends there. The top and bottom plates are hinged at `x0` and `x1` for the
/// segments over the opening, so piers on either side are tied axially only.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallOpening {
    /// Left edge (left jack position), `0 < x0 < x1`.
    pub x0: f64,
    /// Right edge (right jack position), `x1 < L`.
    pub x1: f64,
    /// Height of the header underside above the bottom plate, `0 < head_height < H`.
    pub head_height: f64,
    /// Header bending stiffness `EI`.
    pub header_ei: f64,
    /// Axial stiffness `EA` of each jack (all plies on one side).
    pub jack_ea: f64,
    /// A door cuts the bottom plate between `x0` and `x1`.
    pub door: bool,
}
/// Sheathing of the wall's full-height piers. Sheathing over and under openings is not modelled.
#[derive(Debug, Clone)]
pub struct WallSheathing {
    /// The per-pier law.
    pub law: WallShearLaw,
}
/// An extra exterior node on a plate (`ux, uy` only), where something connects between the corners: the floor
/// or wall below at each stud, a drag strut, a beam bearing. The node must lie at the plate point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallStation<K = NodeId> {
    /// Plate the station is on.
    pub chord: WallChord,
    /// Local position along the plate, strictly between the corners.
    pub x: f64,
    /// The exterior node; distinct from the corners and other stations.
    pub node: K,
}
/// The floor or roof depth carried by the wall top, between the framing top (`y = H`) and the exterior top row
/// (`y = H + depth`): constant-shear bays between consecutive top positions (the corners and every top station),
/// an axial rim bar along the exterior top row, and rigid vertical posts from each top position's framing node to
/// its exterior node (the vertical load path; crushing and hold-downs belong to the links between walls). It
/// reproduces the O86:19 Fig. 11.7 levers `H + depth` for the wall's overturning.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloorZone {
    /// Depth `d` of the floor system, `> 0`.
    pub depth: f64,
    /// Shear rigidity `GA` (force per length) of the floor diaphragm, `> 0`. Bays over an opening of this wall
    /// carry no shear field (segmented design); the rim spans them.
    pub ga: f64,
    /// Rim axial stiffness `EA` (all of it, `>= 0`), along the exterior top row, in addition to the top plate.
    pub rim_ea: f64,
}
/// The wall stands on the surface below through its bottom plate: the plates, subfloor, rim or blocking and the
/// top plate below (the bearing "sandwich"). The surface is the wall's [`WallChord::Base`] stations, interpolated
/// linearly between consecutive base positions (several nodes at one `x` share equally). Every bottom-plate
/// bearing position (each stud and jack, the corners, bottom stations) gets a vertical spring `kv` and a
/// horizontal spring `kh` to the interpolated point; at each pier end the `hold_down` law, if given, replaces
/// `kv`. Every bearing position must lie within the base positions' range. Design: `docs/stud-wall-bearing.md`.
#[derive(Debug, Clone)]
pub struct WallBearing {
    /// Vertical bearing stiffness per position (force/length), `> 0`: the sandwich's crushing under one stud.
    pub kv: f64,
    /// Horizontal stiffness per position (force/length), `> 0`: sill nailing (a penalty for a rigid tie).
    pub kh: f64,
    /// Pier-end anchorage (`Material::anchorage` or any uniaxial material with a positive initial tangent):
    /// deformation is the post's rise relative to the surface, tension positive. `None`: pier ends bear on `kv`.
    pub hold_down: Option<Material>,
}
/// Framing of one wall, in its local frame (`x` along the bottom plate from `c0` to `c1`, `y` up the studs) and
/// consistent units. Validated by [`StudWall::new`].
#[derive(Debug, Clone)]
pub struct StudWallLayout {
    /// Wall length `L` (corner to corner along the bottom plate).
    pub length: f64,
    /// Stud length `H` (plate centreline to centreline).
    pub height: f64,
    /// Stud rows, strictly increasing in `x`, first at 0 and last at `L`.
    pub studs: Vec<WallStud>,
    /// Top plate axial stiffness `EA` (all plies).
    pub top_ea: f64,
    /// Top plate `EI`: required with [`ChordDistribution::Beam`], ignored with `Direct`.
    pub top_ei: Option<f64>,
    /// Top plate load distribution.
    pub top_distribution: ChordDistribution,
    /// Bottom plate axial stiffness `EA`.
    pub bottom_ea: f64,
    /// Bottom plate `EI` (the bottom plate is always a beam).
    pub bottom_ei: f64,
    /// Windows and doors, non-overlapping, with a full-height stud between consecutive openings.
    pub openings: Vec<WallOpening>,
    /// Pier sheathing; `None` for a gravity-only wall (whose racking mode must then be restrained).
    pub sheathing: Option<WallSheathing>,
    /// The wall bears on a rigid foundation: every interior bottom-plate position (each stud and jack) is fixed
    /// internally in `u` (sill anchorage) and `v` (bearing), with reactions in
    /// [`StudWallResponse::base_nodes`]. The bottom corners stay exterior. A bottom station replaces the
    /// internal support at its position, so hold-downs can attach at interior pier ends
    /// ([`StudWallLayout::pier_ends`]).
    pub foundation: bool,
    /// The wall bears on its [`WallChord::Base`] stations ([`WallBearing`]); exclusive with `foundation`.
    pub bearing: Option<WallBearing>,
    /// Floor zone above the framing. With it the exterior top corners and top stations sit at local
    /// `y = H + depth`, top-slot loads act on that row, and every top station must coincide with a stud, jack or
    /// cripple position.
    pub floor: Option<FloorZone>,
}
impl StudWallLayout {
    /// Pier end-post positions (each pier's `xa` and `xb`, in `x` order): where segmented design needs
    /// overturning anchorage. `0` and `L` are the corners `c0`/`c1`; attach hold-downs at the others through
    /// bottom stations (with `foundation`, a station there replaces the internal support).
    pub fn pier_ends(&self) -> Vec<f64> {
        pier_bounds(self)
            .into_iter()
            .flat_map(|(a, b)| [a, b])
            .collect()
    }
}

/// A load on the top plate, in local axes (`fx` along the wall, `fy` along the studs, negative down). A wall owns a
/// list of load cases (lists of these, given at construction); a load pattern selects case `k` and scales it through
/// `ElementLoad::wall`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WallLoad {
    /// Point load at `x`.
    TopPoint { x: f64, fx: f64, fy: f64 },
    /// Line load from `x0` to `x1`, varying linearly from `(wx0, wy0)` to `(wx1, wy1)` per unit length.
    TopLine {
        x0: f64,
        x1: f64,
        wx0: f64,
        wx1: f64,
        wy0: f64,
        wy1: f64,
    },
}
/// Invalid layout, stations or loads, naming the offending item.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StudWallError(pub String);
impl std::fmt::Display for StudWallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for StudWallError {}

#[derive(Clone, Copy, Debug)]
struct ZoneBayRec {
    x0: f64,
    x1: f64,
    /// `[bottom i, bottom j, top i, top j]`, each `(u, v)`.
    dofs: [usize; 8],
    /// Shear rigidity; zero over an opening.
    ga: f64,
}
#[derive(Clone, Copy, Debug)]
struct BarRec {
    x: f64,
    length: f64,
    k: f64,
    a: usize,
    b: usize,
    stud: bool,
    jack: bool,
}
#[derive(Clone, Copy, Debug)]
struct ChordAxialRec {
    x0: f64,
    x1: f64,
    a: usize,
    b: usize,
    ea: f64,
    chord: WallChord,
}
#[derive(Clone, Copy, Debug)]
struct BeamRec {
    x0: f64,
    x1: f64,
    dofs: [usize; 6],
    ei: f64,
    chord: WallChord,
}
#[derive(Clone, Copy, Debug)]
struct HeaderRec {
    x0: f64,
    x1: f64,
    dofs: [usize; 4],
    ei: f64,
}
#[derive(Clone, Debug)]
struct PierRec {
    xa: f64,
    xb: f64,
    beta: DVector<f64>,
    c: f64,
    area: f64,
}
/// One bearing position of a [`WallBearing`] wall: the plate node's `(u, v)` DOFs and the interpolated base
/// point (exterior node index, weight).
#[derive(Clone, Debug)]
struct BearingRec {
    x: f64,
    u: usize,
    v: usize,
    base: Vec<(usize, f64)>,
    /// Nonlinear component index of the hold-down at a pier end (replacing `kv`).
    hold: Option<usize>,
}
#[derive(Debug, Clone)]
struct PreparedAssembly {
    stiffness: Vec<f64>,
    loads: Vec<Vec<f64>>,
}
#[derive(Debug, Clone)]
struct NonlinearCache {
    coupling: DMatrix<f64>,
    strain_map: DMatrix<f64>,
    compliance: DMatrix<f64>,
    load_projection: DMatrix<f64>,
    load_strain: Vec<DVector<f64>>,
    internal_recovery: DMatrix<f64>,
}
#[derive(Debug, Clone)]
struct PreparedNonlinear {
    coupling: Vec<f64>,
    strain_map: Vec<f64>,
}
type CondensedSystem = (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>);

#[derive(Debug, Clone)]
/// Stud-wall element: build with [`StudWall::new`], add as `Element::StudWall`, read results with
/// `Domain::element_stud_wall_response`. Exterior nodes are the corners `c0..c3` (counter-clockwise from
/// bottom-left) then the stations, each with `ux, uy`. Limits: 256 exterior nodes; with a nonlinear law, 32
/// exterior nodes and `MAX_WALL_PIERS` piers.
pub struct StudWall<K: NodeKey = NodeId> {
    corners: [K; 4],
    stations: Vec<WallStation<K>>,
    nodes: Vec<K>,
    layout: StudWallLayout,
    stiffness: DMatrix<f64>,
    load_slots: Vec<DVector<f64>>,
    beam_load_slots: Vec<Vec<[f64; 6]>>,
    axial_load_slots: Vec<Vec<[f64; 2]>>,
    recovery: DMatrix<f64>,
    load_recovery: Vec<DVector<f64>>,
    bars: Vec<BarRec>,
    axial_chords: Vec<ChordAxialRec>,
    beams: Vec<BeamRec>,
    headers: Vec<HeaderRec>,
    piers: Vec<PierRec>,
    top_nodes: Vec<(f64, [usize; 3])>,
    bottom_nodes: Vec<(f64, [usize; 3])>,
    zone_bays: Vec<ZoneBayRec>,
    prepared: Option<PreparedAssembly>,
    prepared_nonlinear: Option<PreparedNonlinear>,
    full_stiffness: Option<DMatrix<f64>>,
    full_load_slots: Vec<DVector<f64>>,
    foundation_supports: Vec<(f64, usize, usize)>,
    bearing: Vec<BearingRec>,
    /// A bearing wall's out-of-plane ties (3D): each bottom exterior node (corner or bottom station, by exterior index) and its
    /// interpolated base (exterior indices and weights), joined by `kh` along the wall's normal. The wall stiffens nothing else
    /// out of plane, so without them its bottom would not follow the floor.
    plan_ties: Vec<(usize, Vec<(usize, f64)>)>,
    /// Per exterior node: 1 when it is a tied bottom node (its normal slots are equations of this element), 2 when it is only a
    /// base node of a tie (its normal slots join only an equation something else makes, `conditional_dof_mask`), else 0.
    plan_tied: Vec<u8>,
    /// Nonlinear components: the piers when the sheathing law is nonlinear (the first `pier_comps`), then the
    /// hold-downs of a bearing wall. Each is a `PierRec` (`area` 1 and `c` 0 for a hold-down).
    comps: Vec<PierRec>,
    pier_comps: usize,
    nonlinear: Option<NonlinearCache>,
    /// Per component: its uniaxial material state (a `Material` sheathing law or a hold-down), else `None`.
    material_states: Vec<Option<Material>>,
    committed_gamma: ArrayVec<f64, MAX_WALL_PIERS>,
    committed_q: ArrayVec<f64, MAX_WALL_PIERS>,
}

impl StudWall<NodeId> {
    /// Builds and condenses the wall. `corners` are `[c0, c1, c2, c3]`, counter-clockwise from bottom-left;
    /// `cases[k]` are the top-plate loads of load case `k` (any number of cases). Fails with [`StudWallError`] for an invalid layout,
    /// station or load, a reused node, or a singular internal stiffness. Corner geometry is checked later,
    /// by `validate`, when the analysis is built.
    pub fn new(
        corners: [NodeId; 4],
        stations: Vec<WallStation>,
        layout: StudWallLayout,
        cases: Vec<Vec<WallLoad>>,
    ) -> Result<Self, StudWallError> {
        Self::build(corners, stations, layout, cases)
    }
}

impl<K: NodeKey> StudWall<K> {
    pub fn nodes(&self) -> &[K] {
        &self.nodes
    }
    pub fn layout(&self) -> &StudWallLayout {
        &self.layout
    }
    pub fn stations(&self) -> &[WallStation<K>] {
        &self.stations
    }
    /// The in-plane slots, and the normal slots of a 3D bearing wall's plan ties.
    fn masks<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> (DofMask, DofMask) {
        let p = self.corners.map(|id| nodes.get(id).coords);
        match PlaneMap::from_corners(p[0], p[1], p[3]) {
            Ok(map) => {
                let plane = map
                    .slots
                    .iter()
                    .fold(DofMask::none(), |m, s| m.with(*s as usize));
                let normal = if NDIM >= 3 && !self.plan_ties.is_empty() {
                    normal_slots(&map)
                        .iter()
                        .fold(DofMask::none(), |m, s| m.with(s.0 as usize))
                } else {
                    DofMask::none()
                };
                (plane, normal)
            }
            Err(_) => (
                (0..NDIM.min(3)).fold(DofMask::none(), |m, s| m.with(s)),
                DofMask::none(),
            ),
        }
    }
    pub fn dof_mask<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DofMask {
        let (plane, normal) = self.masks(nodes);
        plane.union(normal)
    }
    /// The in-plane slots at every node, plus the normal slots at the tied bottom nodes.
    pub fn node_dof_mask<const NDIM: usize, const NDOF: usize>(
        &self,
        k: usize,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DofMask {
        let (plane, normal) = self.masks(nodes);
        if self.plan_tied.get(k) == Some(&1) {
            plane.union(normal)
        } else {
            plane
        }
    }
    /// The normal slots at a tie's base nodes: the bottom follows the base out of plane only where something else (a
    /// diaphragm, a beam, a support) holds the base that way; elsewhere the base stays at zero, so no mechanism is made.
    pub fn conditional_dof_mask<const NDIM: usize, const NDOF: usize>(
        &self,
        k: usize,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DofMask {
        if self.plan_tied.get(k) == Some(&2) {
            self.masks(nodes).1
        } else {
            DofMask::none()
        }
    }
    pub fn accepts_load(&self, load: &ElementLoad) -> bool {
        !load.has_beam() && !load.has_body() && !load.has_edge() && self.has_case(load.wall)
    }
    /// Number of load cases this wall owns.
    pub fn case_count(&self) -> usize {
        self.load_slots.len()
    }
    /// Whether `load` selects no case or one of this wall's cases (a [`WallCaseLoad::MIXED`] never does).
    pub fn has_case(&self, load: Option<WallCaseLoad>) -> bool {
        load.is_none_or(|w| (w.case as usize) < self.load_slots.len())
    }
    pub fn validate<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> Result<(), &'static str> {
        let p = self.corners.map(|id| nodes.get(id).coords);
        let map = PlaneMap::from_corners(p[0], p[1], p[3])?;
        let l = self.layout.length;
        let h = self.layout.height + self.layout.floor.map_or(0., |z| z.depth);
        let tol = 1e-9 * self.layout.length.max(h).max(1.);
        let expected = [(0., 0.), (l, 0.), (l, h), (0., h)];
        for (i, point) in p.iter().enumerate() {
            let (x, y, off) = map.point(point);
            if off
                .abs()
                .max((x - expected[i].0).abs())
                .max((y - expected[i].1).abs())
                > tol
            {
                return Err(if i == 1 || i == 3 {
                    "corner lengths do not match layout L and H"
                } else {
                    "corner geometry is not a counter-clockwise rectangle"
                });
            }
        }
        for s in &self.stations {
            let q = nodes.get(s.node).coords;
            let (x, y, off) = map.point(&q);
            let target = if s.chord == WallChord::Top { h } else { 0. };
            if (x - s.x).abs().max((y - target).abs()).max(off.abs()) > tol {
                return Err("station does not lie at its declared chord location");
            }
        }
        Ok(())
    }
    pub fn assemble_tangent<const NDIM: usize, const NDOF: usize, S: TangentSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
        sink: &mut S,
    ) {
        let map = self.frame(nodes);
        let mut exterior_dofs = ArrayVec::<(K, u8), MAX_EXTERIOR_DOFS>::new();
        exterior_dofs.extend(
            self.nodes
                .iter()
                .flat_map(|id| map.slots.iter().map(move |slot| (*id, *slot))),
        );
        let mut u = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        for &id in &self.nodes {
            u.extend(map.local(&nodes.get(id).displacement));
        }
        let n = self.nodes.len() * map.slots.len();
        let p = self.comps.len();
        let owned = self.prepared.is_none().then(|| self.prepared_for(&map));
        let prepared = self
            .prepared
            .as_ref()
            .unwrap_or_else(|| &owned.as_ref().unwrap().0);
        let prepared_nonlinear = self
            .prepared_nonlinear
            .as_ref()
            .or_else(|| owned.as_ref().and_then(|data| data.1.as_ref()));
        let mut gd = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        for &id in &self.nodes {
            let d = nodes.get(id).displacement;
            gd.extend(map.slots.iter().map(|&slot| d[slot as usize]));
        }
        let trial = self.nonlinear_trial(u.as_slice(), cases);
        let mut q = ArrayVec::<f64, MAX_WALL_PIERS>::new();
        let mut m = ArrayVec::<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>::new();
        if let Some((_, tq, tm)) = trial {
            q = tq;
            m = tm;
        } else if self.nonlinear.is_some() {
            q.extend((0..p).map(|_| f64::NAN));
            m.extend((0..p * p).map(|_| f64::NAN));
        }
        let mut rg = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        for i in 0..n {
            let mut v = (0..n)
                .map(|j| prepared.stiffness[i * n + j] * gd[j])
                .sum::<f64>();
            if let Some(c) = prepared_nonlinear {
                for a in 0..p {
                    v += c.coupling[i * p + a] * q[a];
                }
            }
            rg.push(v);
        }
        if self.nonlinear.is_some() {
            let c = prepared_nonlinear.expect("prepared nonlinear wall data");
            NONLINEAR_TANGENT.with_borrow_mut(|kg| {
                kg.clear();
                for i in 0..n {
                    for j in 0..n {
                        let mut v = prepared.stiffness[i * n + j];
                        for a in 0..p {
                            for b in 0..p {
                                v += c.coupling[i * p + a] * m[a * p + b] * c.strain_map[b * n + j];
                            }
                        }
                        kg.push(v);
                    }
                }
                self.emit(nodes, &map, &exterior_dofs, kg, &rg, sink);
            });
        } else {
            self.emit(nodes, &map, &exterior_dofs, &prepared.stiffness, &rg, sink);
        }
    }
    /// Hands the condensed in-plane tangent to `sink`, with a 3D bearing wall's plan ties added (one call, so a sink that keeps
    /// only the last call sees the whole element). A tie is `kh (a·d) a`, `a` the bottom node's normal slots minus its base's.
    fn emit<const NDIM: usize, const NDOF: usize, S: TangentSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        map: &PlaneMap,
        dofs: &[(K, u8)],
        k: &[f64],
        r: &[f64],
        sink: &mut S,
    ) {
        let kh = self.layout.bearing.as_ref().map_or(0., |b| b.kh);
        if NDIM < 3 || self.plan_ties.is_empty() || kh == 0. {
            sink.add_dyn(dofs, k, r);
            return;
        }
        let ns = normal_slots(map);
        let m = map.slots.len();
        let mut all: Vec<(K, u8)> = dofs.to_vec();
        // Position of (exterior node, normal slot t) in `all`: its in-plane slot when the normal shares one, else appended.
        let mut pos = vec![ArrayVec::<usize, 3>::new(); self.nodes.len()];
        for (e, p) in pos
            .iter_mut()
            .enumerate()
            .filter(|(e, _)| self.plan_tied[*e] > 0)
        {
            for &(s, _) in &ns {
                p.push(match map.slots.iter().position(|&q| q == s) {
                    Some(i) => e * m + i,
                    None => {
                        all.push((self.nodes[e], s));
                        all.len() - 1
                    }
                });
            }
        }
        let (n0, n) = (dofs.len(), all.len());
        let mut kk = vec![0.; n * n];
        let mut rr = vec![0.; n];
        for i in 0..n0 {
            kk[i * n..i * n + n0].copy_from_slice(&k[i * n0..(i + 1) * n0]);
            rr[i] = r[i];
        }
        for (e, w) in &self.plan_ties {
            let mut a: Vec<(usize, f64)> = Vec::new();
            for (t, &(_, c)) in ns.iter().enumerate() {
                a.push((pos[*e][t], c));
                a.extend(w.iter().map(|&(b, wt)| (pos[b][t], -wt * c)));
            }
            let gap: f64 = a
                .iter()
                .map(|&(i, c)| c * nodes.get(all[i].0).displacement[all[i].1 as usize])
                .sum();
            for &(i, ci) in &a {
                rr[i] += kh * ci * gap;
                for &(j, cj) in &a {
                    kk[i * n + j] += kh * ci * cj;
                }
            }
        }
        sink.add_dyn(&all, &kk, &rr);
    }
    pub fn assemble_load<const NDIM: usize, const NDOF: usize, S: VectorSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
        sink: &mut S,
    ) {
        if cases.is_empty() {
            return;
        }
        let map = self.frame(nodes);
        let mut exterior_dofs = ArrayVec::<(K, u8), MAX_EXTERIOR_DOFS>::new();
        exterior_dofs.extend(
            self.nodes
                .iter()
                .flat_map(|id| map.slots.iter().map(move |slot| (*id, *slot))),
        );
        let owned = self.prepared.is_none().then(|| self.prepared_for(&map));
        let prepared = self
            .prepared
            .as_ref()
            .unwrap_or_else(|| &owned.as_ref().unwrap().0);
        let mut f = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        f.extend((0..self.nodes.len() * map.slots.len()).map(|_| 0.));
        for c in cases {
            for (i, v) in f.iter_mut().enumerate() {
                *v += c.factor * prepared.loads[c.case as usize][i];
            }
        }
        sink.add_dyn(&exterior_dofs, &f);
    }
    pub fn prepare<const NDIM: usize, const NDOF: usize>(
        &mut self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) {
        let map = self.frame(nodes);
        let (prepared, nonlinear) = self.prepared_for(&map);
        self.prepared = Some(prepared);
        self.prepared_nonlinear = nonlinear;
    }
    pub fn commit<const NDIM: usize, const NDOF: usize>(
        &mut self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
    ) {
        if self.nonlinear.is_none() {
            return;
        }
        let ue = self.local_displacement(nodes);
        if let Some((gamma, q, _)) = self.nonlinear_trial(ue.as_slice(), cases) {
            self.committed_gamma.clear();
            self.committed_gamma.extend(gamma.iter().copied());
            self.committed_q.clear();
            self.committed_q.extend(q.iter().copied());
            for (material, strain) in self.material_states.iter_mut().zip(gamma.iter()) {
                if let Some(m) = material {
                    *m = m.commit(*strain);
                }
            }
        } else {
            self.committed_gamma.clear();
            self.committed_gamma
                .extend((0..self.comps.len()).map(|_| f64::NAN));
            self.committed_q.clear();
            self.committed_q
                .extend((0..self.comps.len()).map(|_| f64::NAN));
        }
    }
    pub fn local_force_width(&self) -> usize {
        self.nodes.len() * 2
    }
    /// Nonlinear wall resistance uses the pier flows cached by `commit`, matching the load state
    /// passed to that commit; `element_end_force` subtracts slot loads separately.
    pub fn local_force<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> ElementForce {
        let u = self.local_displacement(nodes);
        if self.nonlinear.is_some() {
            let mut r = &self.stiffness * &u;
            if let Some(cache) = &self.nonlinear {
                for i in 0..r.len() {
                    for j in 0..self.committed_q.len() {
                        r[i] += cache.coupling[(i, j)] * self.committed_q[j];
                    }
                }
            }
            return ElementForce::from_vec(r.as_slice().to_vec());
        }
        ElementForce::from_vec((&self.stiffness * u).as_slice().to_vec())
    }
    pub fn local_load_force<const NDIM: usize, const NDOF: usize>(
        &self,
        _nodes: &NodeView<'_, NDIM, NDOF, K>,
        cases: &[WallCaseLoad],
    ) -> ElementForce {
        let n = self.nodes.len() * 2;
        let mut f = vec![0.; n];
        for c in cases {
            for (i, value) in f.iter_mut().enumerate() {
                *value += c.factor * self.load_slots[c.case as usize][i];
            }
        }
        ElementForce::from_vec(f)
    }
    /// Recovered internal response at the nodes' current state, with `cases` the effective load cases (use
    /// `Domain::element_stud_wall_response`, which supplies them). Nonlinear laws are evaluated from their
    /// committed state, so after a committed step this is the committed response.
    fn frame<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> PlaneMap {
        let a = nodes.get(self.corners[0]).coords;
        let b = nodes.get(self.corners[1]).coords;
        PlaneMap::from_corners(a, b, nodes.get(self.corners[3]).coords)
            .expect("validated wall geometry")
    }
    fn local_displacement<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> DVector<f64> {
        let map = self.frame(nodes);
        DVector::from_iterator(
            self.nodes.len() * 2,
            self.nodes.iter().flat_map(|id| {
                let d = nodes.get(*id).displacement;
                map.local(&d)
            }),
        )
    }

    fn prepared_for(&self, map: &PlaneMap) -> (PreparedAssembly, Option<PreparedNonlinear>) {
        let stiffness = map.to_global_matrix(&self.stiffness, self.nodes.len());
        let loads = self
            .load_slots
            .iter()
            .map(|f| map.to_global_vector(f.as_slice(), self.nodes.len()))
            .collect();
        let nonlinear = self.nonlinear.as_ref().map(|cache| {
            let coupling = map.to_global_cols(&cache.coupling, self.nodes.len());
            let strain_map = map.from_global_rows(&cache.strain_map, self.nodes.len());
            PreparedNonlinear {
                coupling,
                strain_map,
            }
        });
        (PreparedAssembly { stiffness, loads }, nonlinear)
    }
}

impl StudWall<Node3Id> {
    /// Builds a 3D wall; its local axes are derived from corner geometry.
    pub fn new3(
        corners: [Node3Id; 4],
        stations: Vec<WallStation<Node3Id>>,
        layout: StudWallLayout,
        cases: Vec<Vec<WallLoad>>,
    ) -> Result<Self, StudWallError> {
        Self::build(corners, stations, layout, cases)
    }
}
