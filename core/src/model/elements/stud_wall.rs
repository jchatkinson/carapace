//! Stud-wall macro-element: one storey of a 2D light-frame wall (studs, plates, openings with headers, jacks
//! and cripples, optional sheathing), condensed onto its four corners and optional plate stations.
//!
//! The framing is linear elastic. Each full-height pier carries one constant shear field whose law is linear
//! (`ga`), the CSA O86:19 shear-plus-nail-slip law, or any uniaxial [`Material`]. Linear walls are condensed
//! once at construction; nonlinear walls condense the bare frame and solve a small per-pier system at each
//! state. Exterior DOFs are `ux, uy` only.
//!
//! User guide (layout, laws, supports, loads, analysis, results, assumptions, limitations):
//! `docs/stud-wall.md`. Implementation spec and verification: `docs/stud-wall-v1.md`.

use arrayvec::ArrayVec;
use nalgebra::{DMatrix, DVector};

use super::wall_plane::{PlaneMap, WallPlane};
use super::{DofMask, NodeView, TangentSink, VectorSink};
use crate::model::{ElementForce, ElementLoad, Material, Node3Id, NodeId, NodeKey};

/// Number of load slots per wall; see [`WallLoad`] and `ElementLoad::wall`.
pub const WALL_SLOTS: usize = crate::model::WALL_SLOTS;
/// Prepared assembly supports up to 256 exterior nodes (512 translational DOFs).
const MAX_EXTERIOR_NODES: usize = 256;
const MAX_EXTERIOR_DOFS: usize = MAX_EXTERIOR_NODES * 3;
/// Maximum number of piers in a wall with a nonlinear sheathing law (`O86` or `Material`): the local pier
/// solve uses fixed-size stack buffers.
pub const MAX_WALL_PIERS: usize = 16;
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
/// Per-pier sheathing law: shear flow `q` (force per length) against sheathing shear strain `γ_sh`. It is the
/// whole sheathing assembly's response (panel shear and fastener slip, both sides); plate and stud axial
/// flexibility are added by the element. Every pier of the wall uses the same law with its own state.
#[derive(Debug, Clone)]
pub enum WallShearLaw {
    /// `q = ga·γ_sh`. Linear: load cases superpose and `Algorithm::Linear` is exact.
    Linear {
        /// Apparent shear stiffness (force per length), e.g. SDPWS `Ga` in consistent units.
        ga: f64,
    },
    /// CSA O86:19 11.7.1.2 shear and nail-slip terms per unit height, nonlinear elastic and odd in `q`:
    /// `γ_sh = q/Bv + 0.0025·e_n`, `e_n = (0.013·q·s/d_F²)²` (Annex A.11.7: nails in wood-based sheathing on
    /// dry lumber). **N and mm only.** Needs an iterative algorithm (Newton).
    O86 {
        /// Shear-through-thickness rigidity `Bv` (N/mm).
        bv: f64,
        /// Panel-edge nail spacing `s` (mm). Use `s·√2` for green lumber (doubles `e_n`).
        nail_spacing: f64,
        /// Nail diameter `d_F` (mm).
        nail_diameter: f64,
    },
    /// `q = σ(γ_sh)` from a uniaxial material with history, committed per pier at each converged step.
    /// The initial tangent must be positive. Needs an iterative algorithm (Newton).
    Material(Material),
}

/// Evaluate a shared wall or diaphragm bay shear law at one generalized strain.
pub(crate) fn shear_law_response(
    law: &WallShearLaw,
    material: Option<&Material>,
    gamma: f64,
) -> (f64, f64) {
    match law {
        WallShearLaw::Linear { ga } => (*ga * gamma, *ga),
        WallShearLaw::O86 {
            bv,
            nail_spacing,
            nail_diameter,
        } => {
            let a = 0.0025 * (0.013 * nail_spacing / nail_diameter.powi(2)).powi(2);
            let den = 1. / bv + (1. / (bv * bv) + 4. * a * gamma.abs()).sqrt();
            let q = 2. * gamma / den;
            (q, 1. / (1. / bv + 2. * a * q.abs()))
        }
        WallShearLaw::Material(_) => material
            .expect("material law has a state")
            .trial_stress_tangent(gamma),
    }
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

/// A load on the top plate, in local axes (`fx` along the wall, `fy` along the studs, negative down). Loads are
/// grouped into `WALL_SLOTS` slots at construction; a load pattern scales slot `k` through
/// `ElementLoad::wall[k]`.
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

/// Full-height stud axial force (tension positive). Inside a sheathed pier the two ends differ by the shear
/// lag `(q_left − q_right)·H`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallAxialResponse {
    pub x: f64,
    pub bottom: f64,
    pub top: f64,
}
/// Jack or cripple axial force (tension positive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallMemberResponse {
    pub x: f64,
    pub axial: f64,
}
/// Header moment (positive sagging) and shear at a header node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallHeaderResponse {
    pub x: f64,
    pub moment: f64,
    pub shear: f64,
}
/// Plate axial force (tension positive), moment (positive sagging) and shear at a plate node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallChordResponse {
    pub x: f64,
    pub axial: f64,
    pub moment: f64,
    pub shear: f64,
}
/// One full-height pier: its shear field and overturning.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallPierResponse {
    /// First full-height stud of the pier.
    pub xa: f64,
    /// Last full-height stud of the pier.
    pub xb: f64,
    /// Shear flow (force per length); positive for a +x load at the top.
    pub q: f64,
    /// Compatible pier shear strain (sheathing strain plus shear lag).
    pub gamma: f64,
    /// Moment of the pier's stud axial forces at their bottom ends about the pier mid-width, CCW positive: the
    /// overturning moment of everything above the base (negative for a +x load).
    pub moment_base: f64,
    /// As `moment_base`, at the stud top ends; `moment_top = moment_base + q·b·H`.
    pub moment_top: f64,
    /// End-post centroid spacing, `xb - xa`.
    pub jd: f64,
    /// Overturning share only, `[−moment_base/jd, +moment_base/jd]`, tension positive. This is the
    /// CSA O86 / SDPWS M/jd design force; gravity share comes from each post's stud force or gravity analysis.
    pub chords: [f64; 2],
}
/// Force on the wall at one supported bottom position (local axes; gravity gives positive `fy`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallBaseResponse {
    pub x: f64,
    pub fx: f64,
    pub fy: f64,
}
/// Recovered internal response, in local axes, tension positive. See `docs/stud-wall.md` §7 for design use.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StudWallResponse {
    /// Full-height studs, in `x` order.
    pub studs: Vec<WallAxialResponse>,
    /// Jacks, two per opening.
    pub jacks: Vec<WallMemberResponse>,
    /// Cripples above openings.
    pub cripples: Vec<WallMemberResponse>,
    /// Header end reactions `[left, right]` per opening (the jack forces).
    pub header_reactions: Vec<[f64; 2]>,
    /// Header moment and shear at each header node.
    pub header: Vec<WallHeaderResponse>,
    /// Top plate results at each top-plate node.
    pub top_chord: Vec<WallChordResponse>,
    /// Bottom plate results at each bottom-plate node.
    pub bottom_chord: Vec<WallChordResponse>,
    /// Full-height piers, in `x` order (empty without sheathing).
    pub piers: Vec<WallPierResponse>,
    /// Bottom exterior nodes (`c0`, `c1`, bottom stations) and, with `foundation`, every internal support, by `x`.
    pub base_nodes: Vec<WallBaseResponse>,
    /// `[Fx, Fy, M]`: sums over `base_nodes`, `M = Σ (x − L/2)·fy`.
    pub base: [f64; 3],
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
/// exterior nodes and `MAX_WALL_PIERS` (16) piers.
pub struct StudWall<K: NodeKey = NodeId> {
    corners: [K; 4],
    stations: Vec<WallStation<K>>,
    nodes: Vec<K>,
    exterior_dofs: Vec<(K, u8)>,
    layout: StudWallLayout,
    plane: WallPlane,
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
    prepared: Option<PreparedAssembly>,
    prepared_nonlinear: Option<PreparedNonlinear>,
    full_stiffness: Option<DMatrix<f64>>,
    full_load_slots: Vec<DVector<f64>>,
    foundation_supports: Vec<(f64, usize, usize)>,
    nonlinear: Option<NonlinearCache>,
    material_states: Vec<Material>,
    committed_gamma: ArrayVec<f64, MAX_WALL_PIERS>,
    committed_q: ArrayVec<f64, MAX_WALL_PIERS>,
}

impl StudWall<NodeId> {
    /// Builds and condenses the wall. `corners` are `[c0, c1, c2, c3]`, counter-clockwise from bottom-left;
    /// `loads[k]` are the top-plate loads of slot `k`. Fails with [`StudWallError`] for an invalid layout,
    /// station or load, a reused node, or a singular internal stiffness. Corner geometry is checked later,
    /// by `validate`, when the analysis is built.
    pub fn new(
        corners: [NodeId; 4],
        stations: Vec<WallStation>,
        layout: StudWallLayout,
        loads: [Vec<WallLoad>; WALL_SLOTS],
    ) -> Result<Self, StudWallError> {
        Self::build(corners, stations, layout, loads, WallPlane::Planar)
    }
}

impl<K: NodeKey> StudWall<K> {
    fn build(
        corners: [K; 4],
        stations: Vec<WallStation<K>>,
        layout: StudWallLayout,
        loads: [Vec<WallLoad>; WALL_SLOTS],
        plane: WallPlane,
    ) -> Result<Self, StudWallError> {
        validate_layout(&layout, &stations)?;
        let is_nonlinear = layout
            .sheathing
            .as_ref()
            .is_some_and(|s| !matches!(&s.law, WallShearLaw::Linear { .. }));
        let potential_piers = if is_nonlinear {
            pier_bounds(&layout).len()
        } else {
            0
        };
        if potential_piers > MAX_WALL_PIERS {
            return Err(StudWallError(format!(
                "pier count {potential_piers} exceeds nonlinear stud-wall limit of {MAX_WALL_PIERS}"
            )));
        }
        let mut nodes = corners.to_vec();
        nodes.extend(stations.iter().map(|s| s.node));
        if nodes
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>()
            .len()
            != nodes.len()
        {
            return Err(StudWallError(
                "stations and corners must use distinct exterior nodes".into(),
            ));
        }
        if nodes.len() > MAX_EXTERIOR_NODES {
            return Err(StudWallError(format!(
                "exterior node count {} exceeds stud-wall limit of {MAX_EXTERIOR_NODES}",
                nodes.len()
            )));
        }
        let ne = nodes.len() * 2;
        let mut model = FrameBuild::new(ne);
        let mut top_x: Vec<f64> = layout.studs.iter().map(|s| s.x).collect();
        for o in &layout.openings {
            for s in &layout.studs {
                if s.x > o.x0 && s.x < o.x1 {
                    top_x.push(s.x);
                }
            }
        }
        top_x.extend(
            stations
                .iter()
                .filter(|s| s.chord == WallChord::Top)
                .map(|s| s.x),
        );
        if layout.top_distribution == ChordDistribution::Beam {
            for o in &layout.openings {
                top_x.extend([o.x0, o.x1]);
            }
        }
        top_x.sort_by(f64::total_cmp);
        dedup_x(&mut top_x, layout.length);
        let bottom_studs: Vec<f64> = layout
            .studs
            .iter()
            .filter(|s| !layout.openings.iter().any(|o| s.x > o.x0 && s.x < o.x1))
            .map(|s| s.x)
            .collect();
        let mut bottom_x = bottom_studs.clone();
        for o in &layout.openings {
            bottom_x.extend([o.x0, o.x1]);
        }
        bottom_x.extend(
            stations
                .iter()
                .filter(|s| s.chord == WallChord::Bottom)
                .map(|s| s.x),
        );
        bottom_x.sort_by(f64::total_cmp);
        dedup_x(&mut bottom_x, layout.length);
        let top_dofs = if layout.top_distribution == ChordDistribution::Beam {
            3
        } else {
            2
        };
        let top_nodes = make_chord_nodes(
            &mut model,
            &top_x,
            layout.height,
            WallChord::Top,
            top_dofs,
            &corners,
            &stations,
            &nodes,
        );
        let bottom_nodes = make_chord_nodes(
            &mut model,
            &bottom_x,
            0.,
            WallChord::Bottom,
            3,
            &corners,
            &stations,
            &nodes,
        );
        let mut header_nodes = Vec::new();
        for o in &layout.openings {
            let mut xs = vec![o.x0, o.x1];
            xs.extend(
                layout
                    .studs
                    .iter()
                    .filter(|s| s.x > o.x0 && s.x < o.x1)
                    .map(|s| s.x),
            );
            xs.sort_by(f64::total_cmp);
            dedup_x(&mut xs, layout.length);
            let group: Vec<_> = xs
                .iter()
                .map(|&x| {
                    let d = model.add_internal_node(x, o.head_height, 2);
                    let r = [d[0].unwrap(), d[1].unwrap()];
                    header_nodes.push((x, r));
                    r
                })
                .collect();
            for i in 0..xs.len() - 1 {
                let d = xs[i + 1] - xs[i];
                model.add_beam4(
                    [group[i][0], group[i][1], group[i + 1][0], group[i + 1][1]],
                    o.header_ei,
                    d,
                );
            }
        }
        let mut bars = Vec::new();
        let mut beams = Vec::new();
        let mut headers = Vec::new();
        let mut axial_chords = Vec::new();
        // Chord elements, with rotations retained internally and exterior translations mapped directly.
        let top_segments = add_chord_elements(
            &mut model,
            &top_nodes,
            layout.top_ea,
            layout.top_ei.unwrap_or(0.),
            layout.top_distribution,
            WallChord::Top,
            &layout.openings,
            &mut beams,
            &mut axial_chords,
        );
        for pair in bottom_nodes.windows(2) {
            if layout
                .openings
                .iter()
                .any(|o| o.door && pair[0].0 >= o.x0 - 1e-12 && pair[1].0 <= o.x1 + 1e-12)
            {
                continue;
            }
            let [a, b] = released_segment(&mut model, pair, &layout.openings);
            add_beam3(
                &mut model,
                a.0,
                b.0,
                a.1,
                b.1,
                layout.bottom_ea,
                layout.bottom_ei,
                WallChord::Bottom,
                &mut beams,
            );
            axial_chords.push(ChordAxialRec {
                x0: pair[0].0,
                x1: pair[1].0,
                a: pair[0].1[0],
                b: pair[1].1[0],
                ea: layout.bottom_ea,
                chord: WallChord::Bottom,
            });
        }
        for (x, stud) in layout.studs.iter().map(|s| (s.x, s)) {
            if layout.openings.iter().any(|o| x > o.x0 && x < o.x1) {
                continue;
            }
            let a = find_chord(&bottom_nodes, x)
                .ok_or_else(|| StudWallError("missing bottom stud node".into()))?;
            let b = find_chord(&top_nodes, x)
                .ok_or_else(|| StudWallError("missing top stud node".into()))?;
            let k = stud.ea / layout.height;
            model.add_bar(a[1], b[1], k);
            bars.push(BarRec {
                x,
                length: layout.height,
                k,
                a: a[1],
                b: b[1],
                stud: true,
                jack: false,
            });
        }
        for o in &layout.openings {
            let hnode = header_nodes
                .iter()
                .find(|(x, _)| (*x - o.x0).abs() < 1e-12)
                .unwrap()
                .1;
            let hright = header_nodes
                .iter()
                .find(|(x, _)| (*x - o.x1).abs() < 1e-12)
                .unwrap()
                .1;
            for (x, h) in [(o.x0, hnode[0]), (o.x1, hright[0])] {
                let b = find_chord(&bottom_nodes, x).unwrap()[1];
                let k = o.jack_ea / o.head_height;
                model.add_bar(b, h, k);
                bars.push(BarRec {
                    x,
                    length: o.head_height,
                    k,
                    a: b,
                    b: h,
                    stud: false,
                    jack: true,
                });
            }
            for s in layout.studs.iter().filter(|s| s.x > o.x0 && s.x < o.x1) {
                let a = header_nodes
                    .iter()
                    .find(|(x, _)| (*x - s.x).abs() < 1e-12)
                    .unwrap()
                    .1[0];
                let b = find_chord(&top_nodes, s.x).unwrap()[1];
                let len = layout.height - o.head_height;
                let k = s.ea / len;
                model.add_bar(a, b, k);
                bars.push(BarRec {
                    x: s.x,
                    length: len,
                    k,
                    a,
                    b,
                    stud: false,
                    jack: false,
                });
            }
        }
        for o in &layout.openings {
            let mut xs: Vec<_> = header_nodes
                .iter()
                .filter(|(x, _)| *x >= o.x0 && *x <= o.x1)
                .copied()
                .collect();
            xs.sort_by(|a, b| a.0.total_cmp(&b.0));
            for w in xs.windows(2) {
                headers.push(HeaderRec {
                    x0: w[0].0,
                    x1: w[1].0,
                    dofs: [w[0].1[0], w[0].1[1], w[1].1[0], w[1].1[1]],
                    ei: o.header_ei,
                });
            }
        }
        let piers = add_sheathing(&mut model, &layout, &top_nodes, &bottom_nodes)?;
        let mut foundation_supports = Vec::new();
        if layout.foundation {
            // A bottom station replaces the internal support at its position (for a hold-down there).
            for &(x, d) in &bottom_nodes {
                if x > 0. && x < layout.length && d[1] >= ne {
                    model.fixed_internal.extend([d[0], d[1]]);
                    foundation_supports.push((x, d[0], d[1]));
                }
            }
        }
        let full_stiffness = layout.foundation.then(|| model.k.clone());
        let mut nonlinear = if is_nonlinear {
            Some(make_nonlinear_cache(&model, ne, &piers)?)
        } else {
            None
        };
        let (kc, recovery, _ki) = model.condense(ne)?;
        let mut load_slots = Vec::with_capacity(WALL_SLOTS);
        let mut load_recovery = Vec::with_capacity(WALL_SLOTS);
        let mut full_load_slots = Vec::with_capacity(WALL_SLOTS);
        for (slot, slot_loads) in loads.iter().enumerate() {
            let mut p = DVector::zeros(model.ndof);
            for l in slot_loads {
                discretize_load(l, slot, &layout, &top_nodes, &top_segments, &mut p)?;
            }
            if let Some(cache) = &mut nonlinear {
                cache
                    .load_strain
                    .push(pier_load_strain(&cache.load_projection, ne, &p));
            }
            let (f, w) = model.condense_load(ne, &p)?;
            if layout.foundation {
                full_load_slots.push(p);
            }
            load_slots.push(f);
            load_recovery.push(w);
        }
        let beam_load_slots = (0..WALL_SLOTS)
            .map(|slot| {
                beams
                    .iter()
                    .map(|b| {
                        if b.chord == WallChord::Top {
                            segment_load(&loads[slot], b.x0, b.x1, layout.length, true)
                        } else {
                            [0.; 6]
                        }
                    })
                    .collect()
            })
            .collect();
        let axial_load_slots = (0..WALL_SLOTS)
            .map(|slot| {
                axial_chords
                    .iter()
                    .map(|b| {
                        let v = if b.chord == WallChord::Top {
                            segment_load(
                                &loads[slot],
                                b.x0,
                                b.x1,
                                layout.length,
                                layout.top_distribution == ChordDistribution::Beam,
                            )
                        } else {
                            [0.; 6]
                        };
                        [v[0], v[3]]
                    })
                    .collect()
            })
            .collect();
        let material_states = match layout.sheathing.as_ref().map(|s| &s.law) {
            Some(WallShearLaw::Material(m)) => (0..piers.len()).map(|_| m.clone()).collect(),
            _ => Vec::new(),
        };
        Ok(Self {
            corners,
            stations,
            exterior_dofs: nodes
                .iter()
                .flat_map(|&id| plane.slots().into_iter().map(move |s| (id, s)))
                .collect(),
            nodes,
            layout,
            plane,
            stiffness: kc,
            load_slots,
            beam_load_slots,
            axial_load_slots,
            recovery,
            load_recovery,
            bars,
            axial_chords,
            beams,
            headers,
            piers,
            top_nodes,
            bottom_nodes,
            prepared: None,
            prepared_nonlinear: None,
            full_stiffness,
            full_load_slots,
            foundation_supports,
            nonlinear,
            material_states,
            committed_gamma: ArrayVec::new(),
            committed_q: ArrayVec::new(),
        })
    }
    pub fn nodes(&self) -> &[K] {
        &self.nodes
    }
    pub fn layout(&self) -> &StudWallLayout {
        &self.layout
    }
    pub fn stations(&self) -> &[WallStation<K>] {
        &self.stations
    }
    pub fn dof_mask(&self) -> DofMask {
        self.plane
            .slots()
            .into_iter()
            .fold(DofMask::none(), |m, s| m.with(s as usize))
    }
    pub fn accepts_load(&self, load: &ElementLoad) -> bool {
        !load.has_beam() && !load.has_body() && !load.has_edge()
    }
    pub fn validate<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> Result<(), &'static str> {
        let p = self.corners.map(|id| nodes.get(id).coords);
        let map = self.plane.map(p[0], p[1]);
        let l = self.layout.length;
        let h = self.layout.height;
        let tol = 1e-9 * self.layout.length.max(self.layout.height).max(1.);
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
        if NDIM == 3 {
            let dx = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
            if dx[1].abs() > tol || (dx[0] * map.ex[0] + dx[2] * map.ex[2] - l).abs() > tol {
                return Err("corner lengths do not match layout L and H");
            }
        }
        for s in &self.stations {
            let q = nodes.get(s.node).coords;
            let (x, y, off) = map.point(&q);
            let target = if s.chord == WallChord::Top {
                self.layout.height
            } else {
                0.
            };
            if (x - s.x).abs().max((y - target).abs()).max(off.abs()) > tol {
                return Err("station does not lie at its declared chord location");
            }
        }
        Ok(())
    }
    pub fn assemble_tangent<const NDIM: usize, const NDOF: usize, S: TangentSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        _load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let map = self.frame(nodes);
        let mut u = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        for &id in &self.nodes {
            u.extend(map.local(&nodes.get(id).displacement));
        }
        let n = self.nodes.len() * map.slots.len();
        let p = self.piers.len();
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
        let trial = self.nonlinear_trial(u.as_slice(), _load);
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
                sink.add_dyn(&self.exterior_dofs, kg, &rg);
            });
        } else {
            sink.add_dyn(&self.exterior_dofs, &prepared.stiffness, &rg);
        }
    }
    pub fn assemble_load<const NDIM: usize, const NDOF: usize, S: VectorSink<K>>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let Some(load) = load else { return };
        let map = self.frame(nodes);
        let owned = self.prepared.is_none().then(|| self.prepared_for(&map));
        let prepared = self
            .prepared
            .as_ref()
            .unwrap_or_else(|| &owned.as_ref().unwrap().0);
        let mut f = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        f.extend((0..self.nodes.len() * map.slots.len()).map(|_| 0.));
        for k in 0..WALL_SLOTS {
            for (i, v) in f.iter_mut().enumerate() {
                *v += load.wall[k] * prepared.loads[k][i];
            }
        }
        sink.add_dyn(&self.exterior_dofs, &f);
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
        load: Option<&ElementLoad>,
    ) {
        if self.nonlinear.is_none() {
            return;
        }
        let ue = self.local_displacement(nodes);
        if let Some((gamma, q, _)) = self.nonlinear_trial(ue.as_slice(), load) {
            self.committed_gamma.clear();
            self.committed_gamma.extend(gamma.iter().copied());
            self.committed_q.clear();
            self.committed_q.extend(q.iter().copied());
            for (material, strain) in self.material_states.iter_mut().zip(gamma.iter()) {
                *material = material.commit(*strain);
            }
        } else {
            self.committed_gamma.clear();
            self.committed_gamma
                .extend((0..self.piers.len()).map(|_| f64::NAN));
            self.committed_q.clear();
            self.committed_q
                .extend((0..self.piers.len()).map(|_| f64::NAN));
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
        load: Option<&ElementLoad>,
    ) -> ElementForce {
        let n = self.nodes.len() * 2;
        let mut f = vec![0.; n];
        if let Some(l) = load {
            for k in 0..WALL_SLOTS {
                for (i, value) in f.iter_mut().enumerate() {
                    *value += l.wall[k] * self.load_slots[k][i];
                }
            }
        }
        ElementForce::from_vec(f)
    }
    /// Recovered internal response at the nodes' current state, with `load` the effective slot factors (use
    /// `Domain::element_stud_wall_response`, which supplies them). Nonlinear laws are evaluated from their
    /// committed state, so after a committed step this is the committed response.
    pub fn response<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
        load: Option<&ElementLoad>,
    ) -> StudWallResponse {
        let ue = self.local_displacement(nodes);
        let mut ui = -(&self.recovery * &ue);
        if let Some(l) = load {
            for k in 0..WALL_SLOTS {
                ui += &self.load_recovery[k] * l.wall[k];
            }
        }
        let nonlinear_state = self
            .nonlinear
            .as_ref()
            .and_then(|_| self.nonlinear_trial(ue.as_slice(), load));
        let nonlinear_q = nonlinear_state.as_ref().map(|(_, q, _)| q.clone());
        if let (Some(cache), Some((_, q, _))) = (&self.nonlinear, &nonlinear_state) {
            for i in 0..ui.len() {
                for j in 0..q.len() {
                    ui[i] -= cache.internal_recovery[(i, j)] * q[j];
                }
            }
        }
        let mut ufull = ue.as_slice().to_vec();
        ufull.extend(ui.iter().copied());
        let q: Vec<f64> = if let Some(q) = nonlinear_q.as_ref() {
            q.iter().copied().collect()
        } else if self.nonlinear.is_some() {
            vec![f64::NAN; self.piers.len()]
        } else {
            self.piers
                .iter()
                .map(|p| {
                    let gamma = p.beta.dot(&DVector::from_column_slice(&ufull));
                    let WallShearLaw::Linear { ga } = &self.layout.sheathing.as_ref().unwrap().law
                    else {
                        unreachable!("nonlinear laws use the local pier solve")
                    };
                    ga * p.area * gamma / (p.area + ga * p.c)
                })
                .collect()
        };
        let disp = |d: usize| {
            if d < ue.len() {
                ue[d]
            } else {
                ui[d - ue.len()]
            }
        };
        let mut r = StudWallResponse::default();
        let factors = load.map(|l| l.wall).unwrap_or([0.; WALL_SLOTS]);
        for b in &self.bars {
            let n = b.k * (disp(b.b) - disp(b.a));
            if b.stud {
                let q_left = self.field_side(b.x, true, &q);
                let q_right = self.field_side(b.x, false, &q);
                let delta = q_left - q_right;
                r.studs.push(WallAxialResponse {
                    x: b.x,
                    bottom: n - delta * b.length / 2.,
                    top: n + delta * b.length / 2.,
                });
            } else if b.jack {
                r.jacks.push(WallMemberResponse { x: b.x, axial: n });
            } else {
                r.cripples.push(WallMemberResponse { x: b.x, axial: n });
            }
        }
        r.header_reactions = r.jacks.iter().map(|j| [j.x, j.axial]).collect();
        for h in &self.headers {
            let l = h.x1 - h.x0;
            let [v0, t0, v1, t1] = h.dofs;
            let d = [disp(v0), disp(t0), disp(v1), disp(t1)];
            let k = beam4(h.ei, l);
            let f: [f64; 4] = std::array::from_fn(|a| (0..4).map(|b| k[a][b] * d[b]).sum());
            for (x, m, shear) in [(h.x0, f[1], f[0]), (h.x1, f[3], -f[2])] {
                r.header.push(WallHeaderResponse {
                    x,
                    moment: m,
                    shear,
                });
            }
        }
        r.top_chord = self.chord_response(&self.top_nodes, &ue, &ui, WallChord::Top, &q, &factors);
        r.bottom_chord = self.chord_response(
            &self.bottom_nodes,
            &ue,
            &ui,
            WallChord::Bottom,
            &q,
            &factors,
        );
        for (i, p) in self.piers.iter().enumerate() {
            let gamma = p.beta.dot(&DVector::from_column_slice(&ufull));
            let mid = (p.xa + p.xb) * 0.5;
            let mut moment_base = 0.;
            let mut moment_top = 0.;
            for stud in r.studs.iter().filter(|s| s.x >= p.xa && s.x <= p.xb) {
                moment_base += (stud.x - mid) * stud.bottom;
                moment_top += (stud.x - mid) * stud.top;
            }
            let jd = p.xb - p.xa;
            r.piers.push(WallPierResponse {
                xa: p.xa,
                xb: p.xb,
                q: q[i],
                gamma,
                moment_base,
                moment_top,
                jd,
                chords: [-moment_base / jd, moment_base / jd],
            });
        }
        let mut resisting = &self.stiffness * &ue;
        if let (Some(cache), Some(q)) = (&self.nonlinear, &nonlinear_q) {
            for i in 0..resisting.len() {
                for j in 0..q.len() {
                    resisting[i] += cache.coupling[(i, j)] * q[j];
                }
            }
        }
        for (k, slot) in self.load_slots.iter().enumerate() {
            resisting -= slot * factors[k];
        }
        for (i, id) in self.nodes.iter().enumerate() {
            if i < 2
                || self
                    .stations
                    .iter()
                    .any(|s| s.node == *id && s.chord == WallChord::Bottom)
            {
                let fx = resisting[2 * i];
                let fy = resisting[2 * i + 1];
                let x = if i == 0 {
                    0.
                } else if i == 1 {
                    self.layout.length
                } else {
                    self.stations
                        .iter()
                        .find(|s| s.node == *id)
                        .map(|s| s.x)
                        .unwrap_or(0.)
                };
                r.base[0] += fx;
                r.base[1] += fy;
                r.base[2] += (x - self.layout.length / 2.) * fy;
                r.base_nodes.push(WallBaseResponse { x, fx, fy });
            }
        }
        if !self.foundation_supports.is_empty() {
            let mut residual =
                self.full_stiffness.as_ref().unwrap() * DVector::from_column_slice(&ufull);
            if let Some(q) = &nonlinear_q {
                for (i, pier) in self.piers.iter().enumerate() {
                    residual += &pier.beta * (pier.area * q[i]);
                }
            }
            for (k, p) in self.full_load_slots.iter().enumerate() {
                residual -= p * factors[k];
            }
            for &(x, u, v) in &self.foundation_supports {
                let fx = residual[u];
                let fy = residual[v];
                r.base[0] += fx;
                r.base[1] += fy;
                r.base[2] += (x - self.layout.length / 2.) * fy;
                r.base_nodes.push(WallBaseResponse { x, fx, fy });
            }
            r.base_nodes.sort_by(|a, b| a.x.total_cmp(&b.x));
        }
        r
    }
    fn chord_response(
        &self,
        _nodes: &[(f64, [usize; 3])],
        ue: &DVector<f64>,
        ui: &DVector<f64>,
        chord: WallChord,
        q: &[f64],
        factors: &[f64; WALL_SLOTS],
    ) -> Vec<WallChordResponse> {
        let get = |d: usize| {
            if d < ue.len() {
                ue[d]
            } else {
                ui[d - ue.len()]
            }
        };
        let mut out = Vec::with_capacity(2 * self.axial_chords.len());
        for (ai, b) in self
            .axial_chords
            .iter()
            .enumerate()
            .filter(|(_, b)| b.chord == chord)
        {
            let l = b.x1 - b.x0;
            let n = b.ea * (get(b.b) - get(b.a)) / l;
            let p0 = (0..WALL_SLOTS)
                .map(|k| factors[k] * self.axial_load_slots[k][ai][0])
                .sum::<f64>();
            let p1 = (0..WALL_SLOTS)
                .map(|k| factors[k] * self.axial_load_slots[k][ai][1])
                .sum::<f64>();
            let qseg = self
                .piers
                .iter()
                .enumerate()
                .find(|(_, p)| (b.x0 + b.x1) * 0.5 >= p.xa && (b.x0 + b.x1) * 0.5 <= p.xb)
                .map(|(i, _)| q[i])
                .unwrap_or(0.);
            let beam_i = self.beams.iter().position(|m| {
                m.chord == chord && (m.x0 - b.x0).abs() < 1e-10 && (m.x1 - b.x1).abs() < 1e-10
            });
            let (m0, v0, m1, v1) = if let Some(bi) = beam_i {
                let m = &self.beams[bi];
                let w = m.dofs;
                let d = [get(w[1]), get(w[2]), get(w[4]), get(w[5])];
                let k = beam4(m.ei, l);
                let f: [f64; 4] = std::array::from_fn(|a| (0..4).map(|j| k[a][j] * d[j]).sum());
                let p: [f64; 6] = std::array::from_fn(|j| {
                    (0..WALL_SLOTS)
                        .map(|s| factors[s] * self.beam_load_slots[s][bi][j])
                        .sum()
                });
                (f[1] - p[2], f[0] - p[1], f[3] - p[5], -(f[2] - p[4]))
            } else {
                (0., 0., 0., 0.)
            };
            for (x, force, moment, shear) in [(b.x0, n + p0, m0, v0), (b.x1, n - p1, m1, v1)] {
                let field = if chord == WallChord::Top {
                    qseg * (x - (b.x0 + b.x1) / 2.)
                } else {
                    -qseg * (x - (b.x0 + b.x1) / 2.)
                };
                out.push(WallChordResponse {
                    x,
                    axial: force + field,
                    moment,
                    shear,
                });
            }
        }
        out
    }
    fn field_side(&self, x: f64, left: bool, q: &[f64]) -> f64 {
        self.piers
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                if left {
                    (p.xb - x).abs() < 1e-10
                } else {
                    (p.xa - x).abs() < 1e-10
                }
            })
            .map(|(i, _)| q[i])
            .sum()
    }
    fn frame<const NDIM: usize, const NDOF: usize>(
        &self,
        nodes: &NodeView<'_, NDIM, NDOF, K>,
    ) -> PlaneMap {
        let a = nodes.get(self.corners[0]).coords;
        let b = nodes.get(self.corners[1]).coords;
        self.plane.map(a, b)
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

    fn nonlinear_trial(
        &self,
        ue: &[f64],
        load: Option<&ElementLoad>,
    ) -> Option<(
        ArrayVec<f64, MAX_WALL_PIERS>,
        ArrayVec<f64, MAX_WALL_PIERS>,
        ArrayVec<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>,
    )> {
        let cache = self.nonlinear.as_ref()?;
        let count = self.piers.len();
        let mut y = ArrayVec::<f64, MAX_WALL_PIERS>::new();
        let mut gamma = ArrayVec::<f64, MAX_WALL_PIERS>::new();
        let mut q = ArrayVec::<f64, MAX_WALL_PIERS>::new();
        let mut d = ArrayVec::<f64, MAX_WALL_PIERS>::new();
        let factors = load.map(|l| l.wall).unwrap_or([0.; WALL_SLOTS]);
        for i in 0..count {
            let mut yi = (0..ue.len())
                .map(|j| cache.strain_map[(i, j)] * ue[j])
                .sum::<f64>();
            for (k, factor) in factors.iter().enumerate() {
                yi += cache.load_strain[k][i] * factor;
            }
            y.push(yi);
            gamma.push(yi);
            q.push(0.);
            d.push(0.);
        }
        let mut jac = ArrayVec::<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>::new();
        jac.extend((0..count * count).map(|_| 0.));
        for _ in 0..40 {
            for i in 0..count {
                let law = &self.layout.sheathing.as_ref()?.law;
                let (stress, tangent) =
                    shear_law_response(law, self.material_states.get(i), gamma[i]);
                q[i] = stress;
                d[i] = tangent;
            }
            let mut norm = 0.0f64;
            for i in 0..count {
                let residual = gamma[i] - y[i]
                    + (0..count)
                        .map(|j| cache.compliance[(i, j)] * q[j])
                        .sum::<f64>();
                norm = norm.max(residual.abs());
                for j in 0..count {
                    jac[i * count + j] = f64::from(i == j) + cache.compliance[(i, j)] * d[j];
                }
            }
            let scale = 1. + y.iter().map(|v| v.abs()).fold(0.0, f64::max);
            if norm <= 1e-14 * scale {
                let mut dqdy = ArrayVec::<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>::new();
                dqdy.extend((0..count * count).map(|_| 0.));
                for col in 0..count {
                    let mut a = jac.clone();
                    let mut rhs = ArrayVec::<f64, MAX_WALL_PIERS>::new();
                    rhs.extend((0..count).map(|_| 0.));
                    rhs[col] = 1.;
                    if !small_lu_solve(&mut a, &mut rhs, count) {
                        return None;
                    }
                    for row in 0..count {
                        dqdy[row * count + col] = d[row] * rhs[row];
                    }
                }
                return Some((gamma, q, dqdy));
            }
            let mut a = jac.clone();
            let mut rhs = ArrayVec::<f64, MAX_WALL_PIERS>::new();
            for i in 0..count {
                rhs.push(
                    gamma[i] - y[i]
                        + (0..count)
                            .map(|j| cache.compliance[(i, j)] * q[j])
                            .sum::<f64>(),
                );
            }
            if !small_lu_solve(&mut a, &mut rhs, count) {
                return None;
            }
            for i in 0..count {
                gamma[i] -= rhs[i];
            }
            if !gamma.iter().all(|v| v.is_finite()) {
                return None;
            }
        }
        None
    }
}

impl StudWall<Node3Id> {
    /// Builds a 3D wall in a vertical plane with global Y up.
    pub fn new3(
        corners: [Node3Id; 4],
        stations: Vec<WallStation<Node3Id>>,
        layout: StudWallLayout,
        loads: [Vec<WallLoad>; WALL_SLOTS],
        plane: WallPlane,
    ) -> Result<Self, StudWallError> {
        let plane = match plane {
            WallPlane::Planar => {
                return Err(StudWallError(
                    "3D stud wall requires a vertical wall plane".into(),
                ))
            }
            WallPlane::Vertical { ex, ez } => WallPlane::vertical(ex, ez).map_err(StudWallError)?,
        };
        Self::build(corners, stations, layout, loads, plane)
    }
}

// Internal frame builder uses exterior DOFs first and internal DOFs after them.
struct FrameBuild {
    ndof: usize,
    k: DMatrix<f64>,
    fixed_internal: Vec<usize>,
}
impl FrameBuild {
    fn new(ne: usize) -> Self {
        Self {
            ndof: ne,
            k: DMatrix::zeros(ne, ne),
            fixed_internal: Vec::new(),
        }
    }
    fn grow(&mut self) {
        let n = self.ndof + 1;
        let mut k = DMatrix::zeros(n, n);
        k.view_mut((0, 0), (self.ndof, self.ndof))
            .copy_from(&self.k);
        self.k = k;
        self.ndof = n;
    }
    fn add_internal_node(&mut self, _x: f64, _y: f64, n: usize) -> [Option<usize>; 3] {
        let mut d = [None; 3];
        for slot in d.iter_mut().take(n) {
            self.grow();
            *slot = Some(self.ndof - 1);
        }
        d
    }
    fn add_bar(&mut self, a: usize, b: usize, k: f64) {
        self.k[(a, a)] += k;
        self.k[(b, b)] += k;
        self.k[(a, b)] -= k;
        self.k[(b, a)] -= k;
    }
    fn add_beam3(&mut self, a: [usize; 3], b: [usize; 3], ea: f64, ei: f64, l: f64) {
        let k = beam6(ea, ei, l);
        let d = [a[0], a[1], a[2], b[0], b[1], b[2]];
        for i in 0..6 {
            for j in 0..6 {
                self.k[(d[i], d[j])] += k[i][j];
            }
        }
    }
    fn add_beam4(&mut self, d: [usize; 4], ei: f64, l: f64) {
        let k = beam4(ei, l);
        for i in 0..4 {
            for j in 0..4 {
                self.k[(d[i], d[j])] += k[i][j];
            }
        }
    }
    fn condense(&self, ne: usize) -> Result<CondensedSystem, StudWallError> {
        let ni = self.ndof - ne;
        let kee = self.k.view((0, 0), (ne, ne)).into_owned();
        if ni == 0 {
            return Ok((kee, DMatrix::zeros(0, ne), DMatrix::zeros(0, 0)));
        }
        let free: Vec<_> = (ne..self.ndof)
            .filter(|d| !self.fixed_internal.contains(d))
            .collect();
        let mut recovery = DMatrix::zeros(ni, ne);
        if free.is_empty() {
            return Ok((kee, recovery, DMatrix::zeros(0, 0)));
        }
        let kei = DMatrix::from_fn(ne, free.len(), |i, j| self.k[(i, free[j])]);
        let kie = DMatrix::from_fn(free.len(), ne, |i, j| self.k[(free[i], j)]);
        let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| self.k[(free[i], free[j])]);
        let lu = kii.clone().lu();
        let Some(r) = lu.solve(&kie) else {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        };
        let kc = kee - &kei * &r;
        if !kc.iter().all(|x| x.is_finite()) {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        }
        for (j, &d) in free.iter().enumerate() {
            recovery.row_mut(d - ne).copy_from(&r.row(j));
        }
        Ok((kc, recovery, kii))
    }
    fn condense_load(
        &self,
        ne: usize,
        p: &DVector<f64>,
    ) -> Result<(DVector<f64>, DVector<f64>), StudWallError> {
        let pe = p.rows(0, ne).into_owned();
        let free: Vec<_> = (ne..self.ndof)
            .filter(|d| !self.fixed_internal.contains(d))
            .collect();
        if free.is_empty() {
            return Ok((pe, DVector::zeros(0)));
        }
        let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| self.k[(free[i], free[j])]);
        let pi = DVector::from_iterator(free.len(), free.iter().map(|&d| p[d]));
        let w_free = kii
            .lu()
            .solve(&pi)
            .ok_or_else(|| StudWallError("internal stiffness K_ii is singular".into()))?;
        let kei = DMatrix::from_fn(ne, free.len(), |i, j| self.k[(i, free[j])]);
        let mut w = DVector::zeros(self.ndof - ne);
        for (j, &d) in free.iter().enumerate() {
            w[d - ne] = w_free[j];
        }
        Ok((pe - kei * &w_free, w))
    }
}

fn make_nonlinear_cache(
    model: &FrameBuild,
    ne: usize,
    piers: &[PierRec],
) -> Result<NonlinearCache, StudWallError> {
    let free: Vec<_> = (ne..model.ndof)
        .filter(|d| !model.fixed_internal.contains(d))
        .collect();
    let p = piers.len();
    let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| model.k[(free[i], free[j])]);
    let lu = kii.lu();
    let beta_i = DMatrix::from_fn(free.len(), p, |i, j| piers[j].beta[free[i]]);
    let mut force_i = beta_i.clone();
    for (j, pier) in piers.iter().enumerate() {
        force_i.column_mut(j).scale_mut(pier.area);
    }
    let solved = if free.is_empty() {
        DMatrix::zeros(0, p)
    } else {
        lu.solve(&force_i)
            .ok_or_else(|| StudWallError("unsheathed internal stiffness K_ii is singular".into()))?
    };
    let mut internal_recovery = DMatrix::zeros(model.ndof - ne, p);
    for (i, &d) in free.iter().enumerate() {
        for j in 0..p {
            internal_recovery[(d - ne, j)] = solved[(i, j)];
        }
    }
    let kei = DMatrix::from_fn(ne, free.len(), |i, j| model.k[(i, free[j])]);
    let mut coupling = DMatrix::zeros(ne, p);
    for j in 0..p {
        for i in 0..ne {
            coupling[(i, j)] = piers[j].area * piers[j].beta[i]
                - (0..free.len())
                    .map(|k| kei[(i, k)] * solved[(k, j)])
                    .sum::<f64>();
        }
    }
    let strain_map = DMatrix::from_fn(p, ne, |i, j| coupling[(j, i)] / piers[i].area);
    let load_solved = if free.is_empty() {
        DMatrix::zeros(0, p)
    } else {
        lu.solve(&beta_i)
            .ok_or_else(|| StudWallError("unsheathed internal stiffness K_ii is singular".into()))?
    };
    let mut load_projection = DMatrix::zeros(p, model.ndof);
    for i in 0..p {
        for (k, &d) in free.iter().enumerate() {
            load_projection[(i, d)] = load_solved[(k, i)];
        }
    }
    let mut compliance: DMatrix<f64> = DMatrix::zeros(p, p);
    for i in 0..p {
        for j in 0..p {
            if !free.is_empty() {
                compliance[(i, j)] = (0..free.len())
                    .map(|k| beta_i[(k, i)] * solved[(k, j)])
                    .sum();
            }
            if i == j {
                compliance[(i, j)] += piers[i].c / piers[i].area;
            }
        }
    }
    if !compliance.iter().all(|x| x.is_finite()) {
        return Err(StudWallError(
            "nonlinear pier condensation is non-finite".into(),
        ));
    }
    Ok(NonlinearCache {
        coupling,
        strain_map,
        compliance,
        load_projection,
        load_strain: Vec::with_capacity(WALL_SLOTS),
        internal_recovery,
    })
}

fn pier_load_strain(projection: &DMatrix<f64>, ne: usize, load: &DVector<f64>) -> DVector<f64> {
    DVector::from_fn(projection.nrows(), |i, _| {
        (ne..load.len()).map(|d| projection[(i, d)] * load[d]).sum()
    })
}

fn small_lu_solve(
    a: &mut ArrayVec<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>,
    rhs: &mut ArrayVec<f64, MAX_WALL_PIERS>,
    n: usize,
) -> bool {
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        if !a[pivot * n + k].is_finite() || a[pivot * n + k].abs() <= f64::MIN_POSITIVE {
            return false;
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            rhs.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        let tail: f64 = (i + 1..n).map(|j| a[i * n + j] * rhs[j]).sum();
        rhs[i] = (rhs[i] - tail) / a[i * n + i];
    }
    rhs.iter().all(|v| v.is_finite())
}

fn dedup_x(xs: &mut Vec<f64>, l: f64) {
    xs.sort_by(f64::total_cmp);
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * l.max(1.));
}
fn validate_layout<K: NodeKey>(
    l: &StudWallLayout,
    stations: &[WallStation<K>],
) -> Result<(), StudWallError> {
    for (name, v) in [
        ("L", l.length),
        ("H", l.height),
        ("top.ea", l.top_ea),
        ("bottom.ea", l.bottom_ea),
        ("bottom.ei", l.bottom_ei),
    ] {
        if !v.is_finite() || v <= 0. {
            return Err(StudWallError(format!("{name} must be positive and finite")));
        }
    }
    if l.top_distribution == ChordDistribution::Beam
        && l.top_ei.is_none_or(|x| !x.is_finite() || x <= 0.)
    {
        return Err(StudWallError(
            "top.ei must be positive and finite for Beam".into(),
        ));
    }
    for (i, s) in l.studs.iter().enumerate() {
        if !s.x.is_finite() || !s.ea.is_finite() || s.ea <= 0. {
            return Err(StudWallError(format!(
                "studs[{i}] x and ea must be finite and ea positive"
            )));
        }
        if i > 0 && s.x <= l.studs[i - 1].x {
            return Err(StudWallError(format!(
                "studs[{i}] must be strictly increasing in x"
            )));
        }
    }
    if l.studs.first().is_none_or(|s| s.x != 0.)
        || l.studs
            .last()
            .is_none_or(|s| (s.x - l.length).abs() > 1e-10 * l.length)
    {
        return Err(StudWallError("studs must start at 0 and end at L".into()));
    }
    for (i, o) in l.openings.iter().enumerate() {
        if !(0. < o.x0
            && o.x0 < o.x1
            && o.x1 < l.length
            && 0. < o.head_height
            && o.head_height < l.height)
        {
            return Err(StudWallError(format!(
                "openings[{i}] has invalid bounds or head_height"
            )));
        }
        if !o.header_ei.is_finite()
            || o.header_ei <= 0.
            || !o.jack_ea.is_finite()
            || o.jack_ea <= 0.
        {
            return Err(StudWallError(format!(
                "openings[{i}] header_ei and jack_ea must be positive and finite"
            )));
        }
        // A stud row at `x0` or `x1` is a king stud: it shares the jack's bottom-plate node (one end post).
        if i > 0
            && (o.x0 <= l.openings[i - 1].x1
                || !l
                    .studs
                    .iter()
                    .any(|s| s.x >= l.openings[i - 1].x1 && s.x <= o.x0))
        {
            return Err(StudWallError(format!(
                "openings[{i}] overlaps or has no full-height stud separation"
            )));
        }
    }
    if let Some(s) = &l.sheathing {
        match &s.law {
            WallShearLaw::Linear { ga } if !ga.is_finite() || *ga <= 0. => {
                return Err(StudWallError(
                    "sheathing.ga must be positive and finite".into(),
                ));
            }
            WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            } if [bv, nail_spacing, nail_diameter]
                .iter()
                .any(|x| !x.is_finite() || **x <= 0.) =>
            {
                return Err(StudWallError(
                    "O86 bv, nail_spacing and nail_diameter must be positive and finite".into(),
                ));
            }
            WallShearLaw::Material(m)
                if !m.initial_tangent().is_finite() || m.initial_tangent() <= 0. =>
            {
                return Err(StudWallError(
                    "sheathing material must have a positive finite initial tangent".into(),
                ));
            }
            _ => {}
        }
    }
    for (i, s) in stations.iter().enumerate() {
        if !s.x.is_finite() || s.x < 0. || s.x > l.length {
            return Err(StudWallError(format!("stations[{i}].x must be in [0,L]")));
        }
        if s.x == 0. || s.x == l.length {
            return Err(StudWallError(format!(
                "stations[{i}] duplicates a corner location"
            )));
        }
        if stations[..i]
            .iter()
            .any(|p| p.chord == s.chord && (p.x - s.x).abs() < 1e-12 * l.length)
        {
            return Err(StudWallError(format!(
                "stations[{i}] duplicates another station location"
            )));
        }
        if s.chord == WallChord::Bottom
            && l.openings
                .iter()
                .any(|o| o.door && s.x > o.x0 && s.x < o.x1)
        {
            return Err(StudWallError(format!("stations[{i}] lies inside a door")));
        }
        if l.top_distribution == ChordDistribution::Direct
            && s.chord == WallChord::Top
            && !l.studs.iter().any(|w| w.x == s.x)
            && !l.openings.iter().any(|o| {
                l.studs
                    .iter()
                    .any(|w| w.x == s.x && w.x > o.x0 && w.x < o.x1)
            })
        {
            return Err(StudWallError(format!(
                "stations[{i}] direct top station must coincide with a stud or cripple"
            )));
        }
    }
    Ok(())
}
#[allow(clippy::too_many_arguments)]
fn make_chord_nodes<K: NodeKey>(
    model: &mut FrameBuild,
    xs: &[f64],
    y: f64,
    chord: WallChord,
    nslot: usize,
    corners: &[K; 4],
    stations: &[WallStation<K>],
    nodes: &[K],
) -> Vec<(f64, [usize; 3])> {
    xs.iter()
        .map(|&x| {
            // Each exterior node contributes translation slots; beam rotations remain internal.
            let external = if x.abs() < 1e-12 {
                Some(if chord == WallChord::Top {
                    (corners[3], 3)
                } else {
                    (corners[0], 0)
                })
            } else if (x - xs.last().copied().unwrap_or(x)).abs() < 1e-12 {
                Some(if chord == WallChord::Top {
                    (corners[2], 2)
                } else {
                    (corners[1], 1)
                })
            } else {
                stations
                    .iter()
                    .find(|s| s.chord == chord && (s.x - x).abs() < 1e-12)
                    .and_then(|s| nodes.iter().position(|n| *n == s.node).map(|i| (s.node, i)))
            };
            let mut d = [None; 3];
            if let Some((_, i)) = external {
                d[0] = Some(2 * i);
                d[1] = Some(2 * i + 1);
                if nslot == 3 {
                    model.grow();
                    d[2] = Some(model.ndof - 1);
                }
            } else {
                d = model.add_internal_node(x, y, nslot);
            }
            (
                x,
                [d[0].unwrap(), d[1].unwrap(), d[2].unwrap_or(usize::MAX)],
            )
        })
        .collect()
}
type ChordSegment = [(f64, [usize; 3]); 2];

/// A plate segment over an opening is hinged at the opening edges: its end there gets its own rotation DOF,
/// so the plate ties the piers on either side axially only (no frame action through plate bending).
fn released_segment(
    model: &mut FrameBuild,
    pair: &[(f64, [usize; 3])],
    openings: &[WallOpening],
) -> ChordSegment {
    let mut seg = [pair[0], pair[1]];
    let tol = 1e-12 * pair[1].0.abs().max(1.);
    for o in openings {
        if pair[0].0 < o.x0 - tol || pair[1].0 > o.x1 + tol || seg[0].1[2] == usize::MAX {
            continue;
        }
        for (end, edge) in [(0, o.x0), (1, o.x1)] {
            if (seg[end].0 - edge).abs() <= tol {
                model.grow();
                seg[end].1[2] = model.ndof - 1;
            }
        }
    }
    seg
}

#[allow(clippy::too_many_arguments)]
fn add_chord_elements(
    model: &mut FrameBuild,
    nodes: &[(f64, [usize; 3])],
    ea: f64,
    ei: f64,
    dist: ChordDistribution,
    chord: WallChord,
    openings: &[WallOpening],
    beams: &mut Vec<BeamRec>,
    axial: &mut Vec<ChordAxialRec>,
) -> Vec<ChordSegment> {
    let mut segments = Vec::new();
    for pair in nodes.windows(2) {
        let l = pair[1].0 - pair[0].0;
        if chord == WallChord::Top && dist == ChordDistribution::Direct {
            model.add_bar(pair[0].1[0], pair[1].1[0], ea / l);
            segments.push([pair[0], pair[1]]);
        } else {
            let [a, b] = released_segment(model, pair, openings);
            add_beam3(model, a.0, b.0, a.1, b.1, ea, ei, chord, beams);
            segments.push([a, b]);
        }
        axial.push(ChordAxialRec {
            x0: pair[0].0,
            x1: pair[1].0,
            a: pair[0].1[0],
            b: pair[1].1[0],
            ea,
            chord,
        });
    }
    segments
}
#[allow(clippy::too_many_arguments)]
fn add_beam3(
    model: &mut FrameBuild,
    x0: f64,
    x1: f64,
    a: [usize; 3],
    b: [usize; 3],
    ea: f64,
    ei: f64,
    chord: WallChord,
    beams: &mut Vec<BeamRec>,
) {
    let dofs = [a[0], a[1], a[2], b[0], b[1], b[2]];
    model.add_beam3(a, b, ea, ei, x1 - x0);
    beams.push(BeamRec {
        x0,
        x1,
        dofs,
        ei,
        chord,
    });
}
fn beam6(ea: f64, ei: f64, l: f64) -> [[f64; 6]; 6] {
    let mut k = [[0.; 6]; 6];
    let a = ea / l;
    for (i, j, v) in [(0, 0, a), (0, 3, -a), (3, 0, -a), (3, 3, a)] {
        k[i][j] += v;
    }
    let b = 12. * ei / l.powi(3);
    let c = 6. * ei / l.powi(2);
    let d = 4. * ei / l;
    let e = 2. * ei / l;
    let m = [[b, c, -b, c], [c, d, -c, e], [-b, -c, b, -c], [c, e, -c, d]];
    for i in 0..4 {
        for j in 0..4 {
            k[[1, 2, 4, 5][i]][[1, 2, 4, 5][j]] += m[i][j];
        }
    }
    k
}
fn beam4(ei: f64, l: f64) -> [[f64; 4]; 4] {
    let b = 12. * ei / l.powi(3);
    let c = 6. * ei / l.powi(2);
    let d = 4. * ei / l;
    let e = 2. * ei / l;
    [[b, c, -b, c], [c, d, -c, e], [-b, -c, b, -c], [c, e, -c, d]]
}
fn find_chord(nodes: &[(f64, [usize; 3])], x: f64) -> Option<[usize; 3]> {
    nodes
        .iter()
        .find(|(xx, _)| (*xx - x).abs() < 1e-10)
        .map(|(_, d)| *d)
}

fn add_sheathing(
    model: &mut FrameBuild,
    l: &StudWallLayout,
    top: &[(f64, [usize; 3])],
    bottom: &[(f64, [usize; 3])],
) -> Result<Vec<PierRec>, StudWallError> {
    let Some(sh) = &l.sheathing else {
        return Ok(Vec::new());
    };
    let ga = match &sh.law {
        WallShearLaw::Linear { ga } => Some(*ga),
        WallShearLaw::O86 { .. } | WallShearLaw::Material(_) => None,
    };
    let ranges = pier_bounds(l);
    let mut out = Vec::new();
    for (xa, xb) in ranges {
        let width = xb - xa;
        if width <= 0. {
            continue;
        }
        let mut beta = DVector::zeros(model.ndof);
        for (nodes, sign) in [(top, 1.), (bottom, -1.)] {
            for pair in nodes.windows(2) {
                let a = pair[0].0.max(xa);
                let b = pair[1].0.min(xb);
                if b <= a {
                    continue;
                }
                let len = pair[1].0 - pair[0].0;
                let t0 = (a - pair[0].0) / len;
                let t1 = (b - pair[0].0) / len;
                let w0 = (b - a) * (2. - t0 - t1) / 2.;
                let w1 = (b - a) * (t0 + t1) / 2.;
                beta[pair[0].1[0]] += sign * w0 / (width * l.height);
                beta[pair[1].1[0]] += sign * w1 / (width * l.height);
            }
        }
        for (x, sgn) in [(xa, -1.), (xb, 1.)] {
            let t = find_chord(top, x)
                .ok_or_else(|| StudWallError("pier boundary has no top stud node".into()))?;
            let b = find_chord(bottom, x)
                .ok_or_else(|| StudWallError("pier boundary has no bottom stud node".into()))?;
            beta[t[1]] += sgn * 0.5 / width;
            beta[b[1]] += sgn * 0.5 / width;
        }
        let mut c = l.height.powi(3) / 12. * (1. / stud_at(l, xa).ea + 1. / stud_at(l, xb).ea);
        for nodes in [top, bottom] {
            let ea = if std::ptr::eq(nodes, top) {
                l.top_ea
            } else {
                l.bottom_ea
            };
            for p in nodes.windows(2) {
                let d = (p[1].0.min(xb) - p[0].0.max(xa)).max(0.);
                c += d.powi(3) / (12. * ea);
            }
        }
        let area = width * l.height;
        if let Some(ga) = ga {
            let fac = area * area / (area / ga + c);
            for i in 0..model.ndof {
                for j in 0..model.ndof {
                    model.k[(i, j)] += fac * beta[i] * beta[j];
                }
            }
        }
        out.push(PierRec {
            xa,
            xb,
            beta,
            c,
            area,
        });
    }
    Ok(out)
}
fn stud_at(l: &StudWallLayout, x: f64) -> &WallStud {
    l.studs.iter().find(|s| (s.x - x).abs() < 1e-10).unwrap()
}
fn pier_bounds(l: &StudWallLayout) -> Vec<(f64, f64)> {
    let mut runs = Vec::new();
    let mut a = 0.;
    for o in &l.openings {
        runs.push((a, o.x0));
        a = o.x1;
    }
    runs.push((a, l.length));
    runs.into_iter()
        .filter_map(|(a, b)| {
            let xs: Vec<_> = l
                .studs
                .iter()
                .filter(|s| s.x >= a && s.x <= b)
                .map(|s| s.x)
                .collect();
            (xs.len() > 1).then(|| (*xs.first().unwrap(), *xs.last().unwrap()))
        })
        .collect()
}

fn discretize_load(
    load: &WallLoad,
    slot: usize,
    l: &StudWallLayout,
    nodes: &[(f64, [usize; 3])],
    segments: &[ChordSegment],
    p: &mut DVector<f64>,
) -> Result<(), StudWallError> {
    let mut point = |x: f64, fx: f64, fy: f64| -> Result<(), StudWallError> {
        if ![x, fx, fy].iter().all(|v| v.is_finite()) || x < 0. || x > l.length {
            return Err(StudWallError(format!(
                "loads[{slot}] contains invalid coordinates or intensities"
            )));
        }
        if let Some((_, d)) = nodes
            .iter()
            .find(|(xx, _)| (*xx - x).abs() < 1e-12 * l.length)
        {
            p[d[0]] += fx;
            p[d[1]] += fy;
            return Ok(());
        }
        let pair = segments
            .iter()
            .find(|w| x > w[0].0 && x < w[1].0)
            .ok_or_else(|| {
                StudWallError(format!("loads[{slot}] point is outside the top chord"))
            })?;
        let le = pair[1].0 - pair[0].0;
        let t = (x - pair[0].0) / le;
        if l.top_distribution == ChordDistribution::Direct {
            p[pair[0].1[0]] += (1. - t) * fx;
            p[pair[1].1[0]] += t * fx;
            p[pair[0].1[1]] += (1. - t) * fy;
            p[pair[1].1[1]] += t * fy;
        } else {
            let h1 = 1. - 3. * t * t + 2. * t * t * t;
            let h2 = le * (t - 2. * t * t + t * t * t);
            let h3 = 3. * t * t - 2. * t * t * t;
            let h4 = le * (-t * t + t * t * t);
            for (d, v) in [
                (pair[0].1[0], (1. - t) * fx),
                (pair[1].1[0], t * fx),
                (pair[0].1[1], h1 * fy),
                (pair[0].1[2], h2 * fy),
                (pair[1].1[1], h3 * fy),
                (pair[1].1[2], h4 * fy),
            ] {
                p[d] += v;
            }
        }
        Ok(())
    };
    match *load {
        WallLoad::TopPoint { x, fx, fy } => point(x, fx, fy)?,
        WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => {
            if ![x0, x1, wx0, wx1, wy0, wy1].iter().all(|v| v.is_finite())
                || x0 < 0.
                || x1 > l.length
                || x1 <= x0
            {
                return Err(StudWallError(format!(
                    "loads[{slot}] has invalid line bounds or intensities"
                )));
            }
            for w in segments {
                let a = x0.max(w[0].0);
                let b = x1.min(w[1].0);
                if b <= a {
                    continue;
                }
                let gauss = [
                    -0.906179845938664,
                    -0.538469310105683,
                    0.,
                    0.538469310105683,
                    0.906179845938664,
                ];
                let weights = [
                    0.236926885056189,
                    0.478628670499366,
                    0.568888888888889,
                    0.478628670499366,
                    0.236926885056189,
                ];
                let dx = (b - a) / 2.;
                let mid = (a + b) / 2.;
                for (z, gw) in gauss.iter().zip(weights) {
                    let x = mid + dx * z;
                    let t = (x - x0) / (x1 - x0);
                    let wx = wx0 + (wx1 - wx0) * t;
                    let wy = wy0 + (wy1 - wy0) * t;
                    let s = (x - w[0].0) / (w[1].0 - w[0].0);
                    let fac = dx * gw;
                    if l.top_distribution == ChordDistribution::Direct {
                        for (d, v) in [
                            (w[0].1[0], (1. - s) * wx),
                            (w[1].1[0], s * wx),
                            (w[0].1[1], (1. - s) * wy),
                            (w[1].1[1], s * wy),
                        ] {
                            p[d] += fac * v;
                        }
                    } else {
                        let le = w[1].0 - w[0].0;
                        for (d, v) in [
                            (w[0].1[0], (1. - s) * wx),
                            (w[1].1[0], s * wx),
                            (w[0].1[1], (1. - 3. * s * s + 2. * s.powi(3)) * wy),
                            (w[0].1[2], le * (s - 2. * s * s + s.powi(3)) * wy),
                            (w[1].1[1], (3. * s * s - 2. * s.powi(3)) * wy),
                            (w[1].1[2], le * (-s * s + s.powi(3)) * wy),
                        ] {
                            p[d] += fac * v;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn segment_load(loads: &[WallLoad], x0: f64, x1: f64, length: f64, beam: bool) -> [f64; 6] {
    let mut out = [0.; 6];
    let le = x1 - x0;
    let mut add = |x: f64, fx: f64, fy: f64, scale: f64| {
        if x <= x0 + 1e-12 * length.max(1.) || x >= x1 - 1e-12 * length.max(1.) {
            return;
        }
        let t = (x - x0) / le;
        out[0] += scale * (1. - t) * fx;
        out[3] += scale * t * fx;
        if beam {
            let n = [
                1. - 3. * t * t + 2. * t.powi(3),
                le * (t - 2. * t * t + t.powi(3)),
                3. * t * t - 2. * t.powi(3),
                le * (-t * t + t.powi(3)),
            ];
            for (j, d) in [1, 2, 4, 5].into_iter().enumerate() {
                out[d] += scale * n[j] * fy;
            }
        } else {
            out[1] += scale * (1. - t) * fy;
            out[4] += scale * t * fy;
        }
    };
    for load in loads {
        match *load {
            WallLoad::TopPoint { x, fx, fy } => add(x, fx, fy, 1.),
            WallLoad::TopLine {
                x0: a,
                x1: b,
                wx0,
                wx1,
                wy0,
                wy1,
            } => {
                let a = a.max(x0);
                let b = b.min(x1);
                if b <= a {
                    continue;
                }
                let z = [
                    -0.906179845938664,
                    -0.538469310105683,
                    0.,
                    0.538469310105683,
                    0.906179845938664,
                ];
                let w = [
                    0.236926885056189,
                    0.478628670499366,
                    0.568888888888889,
                    0.478628670499366,
                    0.236926885056189,
                ];
                let mid = (a + b) * 0.5;
                let half = (b - a) * 0.5;
                for i in 0..5 {
                    let x = mid + half * z[i];
                    let t =
                        (x - (match *load {
                            WallLoad::TopLine { x0, .. } => x0,
                            _ => 0.,
                        })) / (match *load {
                            WallLoad::TopLine { x0, x1, .. } => x1 - x0,
                            _ => 1.,
                        });
                    add(x, wx0 + (wx1 - wx0) * t, wy0 + (wy1 - wy0) * t, half * w[i]);
                }
            }
        }
    }
    out
}
