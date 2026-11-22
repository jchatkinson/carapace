//! Linear stud-wall macro-element with condensed internal framing DOFs.

use arrayvec::ArrayVec;
use nalgebra::{DMatrix, DVector};

use super::{DofMask, NodeView, TangentSink, VectorSink};
use crate::model::{ElementForce, ElementLoad, NodeId};

pub const WALL_SLOTS: usize = crate::model::WALL_SLOTS;
/// Prepared assembly supports up to 256 exterior nodes (512 translational DOFs).
const MAX_EXTERIOR_NODES: usize = 256;
const MAX_EXTERIOR_DOFS: usize = MAX_EXTERIOR_NODES * 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WallChord {
    Top,
    Bottom,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChordDistribution {
    Beam,
    Direct,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallStud {
    pub x: f64,
    pub ea: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallOpening {
    pub x0: f64,
    pub x1: f64,
    pub head_height: f64,
    pub header_ei: f64,
    pub jack_ea: f64,
    pub door: bool,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallSheathing {
    pub ga: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallStation {
    pub chord: WallChord,
    pub x: f64,
    pub node: NodeId,
}
#[derive(Debug, Clone, PartialEq)]
pub struct StudWallLayout {
    pub length: f64,
    pub height: f64,
    pub studs: Vec<WallStud>,
    pub top_ea: f64,
    pub top_ei: Option<f64>,
    pub top_distribution: ChordDistribution,
    pub bottom_ea: f64,
    pub bottom_ei: f64,
    pub openings: Vec<WallOpening>,
    pub sheathing: Option<WallSheathing>,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WallLoad {
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StudWallError(pub String);
impl std::fmt::Display for StudWallError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for StudWallError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallAxialResponse {
    pub x: f64,
    pub bottom: f64,
    pub top: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallMemberResponse {
    pub x: f64,
    pub axial: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallHeaderResponse {
    pub x: f64,
    pub moment: f64,
    pub shear: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallChordResponse {
    pub x: f64,
    pub axial: f64,
    pub moment: f64,
    pub shear: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallPierResponse {
    pub xa: f64,
    pub xb: f64,
    pub q: f64,
    pub gamma: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallBaseResponse {
    pub x: f64,
    pub fx: f64,
    pub fy: f64,
}
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StudWallResponse {
    pub studs: Vec<WallAxialResponse>,
    pub jacks: Vec<WallMemberResponse>,
    pub cripples: Vec<WallMemberResponse>,
    pub header_reactions: Vec<[f64; 2]>,
    pub header: Vec<WallHeaderResponse>,
    pub top_chord: Vec<WallChordResponse>,
    pub bottom_chord: Vec<WallChordResponse>,
    pub piers: Vec<WallPierResponse>,
    pub base_nodes: Vec<WallBaseResponse>,
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
type CondensedSystem = (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>);

#[derive(Debug, Clone)]
/// Linear condensed wall; exterior assembly is prepared for at most 256 nodes.
pub struct StudWall {
    corners: [NodeId; 4],
    stations: Vec<WallStation>,
    nodes: Vec<NodeId>,
    exterior_dofs: Vec<(NodeId, u8)>,
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
    prepared: Option<PreparedAssembly>,
}

impl StudWall {
    pub fn new(
        corners: [NodeId; 4],
        stations: Vec<WallStation>,
        layout: StudWallLayout,
        loads: [Vec<WallLoad>; WALL_SLOTS],
    ) -> Result<Self, StudWallError> {
        validate_layout(&layout, &stations)?;
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
        add_chord_elements(
            &mut model,
            &top_nodes,
            layout.top_ea,
            layout.top_ei.unwrap_or(0.),
            layout.top_distribution,
            WallChord::Top,
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
            add_beam3(
                &mut model,
                pair[0].0,
                pair[1].0,
                pair[0].1,
                pair[1].1,
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
        let (kc, recovery, _ki) = model.condense(ne)?;
        let mut load_slots = Vec::with_capacity(WALL_SLOTS);
        let mut load_recovery = Vec::with_capacity(WALL_SLOTS);
        for (slot, slot_loads) in loads.iter().enumerate() {
            let mut p = DVector::zeros(model.ndof);
            for l in slot_loads {
                discretize_load(l, slot, &layout, &top_nodes, &mut model, &mut p)?;
            }
            let (f, w) = model.condense_load(ne, &p)?;
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
        Ok(Self {
            corners,
            stations,
            exterior_dofs: nodes.iter().flat_map(|&id| [(id, 0), (id, 1)]).collect(),
            nodes,
            layout,
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
        })
    }
    pub fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }
    pub fn layout(&self) -> &StudWallLayout {
        &self.layout
    }
    pub fn stations(&self) -> &[WallStation] {
        &self.stations
    }
    pub fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }
    pub fn accepts_load(&self, load: &ElementLoad) -> bool {
        !load.has_beam() && !load.has_body() && !load.has_edge()
    }
    pub fn validate(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> Result<(), &'static str> {
        let p = self.corners.map(|id| nodes.get(id).coords);
        let dx = [p[1][0] - p[0][0], p[1][1] - p[0][1]];
        let dy = [p[3][0] - p[0][0], p[3][1] - p[0][1]];
        let l = dx[0].hypot(dx[1]);
        let h = dy[0].hypot(dy[1]);
        let tol = 1e-9 * self.layout.length.max(self.layout.height).max(1.);
        if (l - self.layout.length).abs() > tol || (h - self.layout.height).abs() > tol {
            return Err("corner lengths do not match layout L and H");
        }
        if (dx[0] * dy[0] + dx[1] * dy[1]).abs() > tol * l.max(h) {
            return Err("corner rectangle is not orthogonal");
        }
        if ((p[2][0] - p[1][0]) - dy[0]).hypot((p[2][1] - p[1][1]) - dy[1]) > tol
            || ((p[2][0] - p[3][0]) - dx[0]).hypot((p[2][1] - p[3][1]) - dx[1]) > tol
        {
            return Err("corner geometry is not a counter-clockwise rectangle");
        }
        let ux = [dx[0] / l, dx[1] / l];
        let uy = [-ux[1], ux[0]];
        for s in &self.stations {
            let q = nodes.get(s.node).coords;
            let x = (q[0] - p[0][0]) * ux[0] + (q[1] - p[0][1]) * ux[1];
            let y = (q[0] - p[0][0]) * uy[0] + (q[1] - p[0][1]) * uy[1];
            let target = if s.chord == WallChord::Top {
                self.layout.height
            } else {
                0.
            };
            if (x - s.x).abs() > tol || (y - target).abs() > tol {
                return Err("station does not lie at its declared chord location");
            }
        }
        Ok(())
    }
    pub fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, 2, 3, NodeId>,
        _load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let n = self.nodes.len() * 2;
        let mut ug = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        let mut r = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
        let (c, s) = self.frame(nodes);
        for &id in &self.nodes {
            let d = nodes.get(id).displacement;
            ug.push(c * d[0] + s * d[1]);
            ug.push(-s * d[0] + c * d[1]);
        }
        if let Some(prepared) = &self.prepared {
            for i in 0..n {
                r.push((0..n).map(|j| prepared.stiffness[i * n + j] * ug[j]).sum());
            }
            sink.add_dyn(&self.exterior_dofs, &prepared.stiffness, &r);
            return;
        }
        let mut kg = vec![0.; n * n];
        let t = [[c, s], [-s, c]];
        for ni in 0..self.nodes.len() {
            for nj in 0..self.nodes.len() {
                for gi in 0..2 {
                    for gj in 0..2 {
                        let mut v = 0.;
                        for li in 0..2 {
                            for lj in 0..2 {
                                v += t[li][gi]
                                    * self.stiffness[(2 * ni + li, 2 * nj + lj)]
                                    * t[lj][gj];
                            }
                        }
                        kg[(2 * ni + gi) * n + 2 * nj + gj] = v;
                    }
                }
            }
        }
        for i in 0..n {
            r.push((0..n).map(|j| kg[i * n + j] * ug[j]).sum());
        }
        sink.add_dyn(&self.exterior_dofs, &kg, &r);
    }
    pub fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, 2, 3, NodeId>,
        load: Option<&ElementLoad>,
        sink: &mut S,
    ) {
        let Some(load) = load else { return };
        let n = self.nodes.len() * 2;
        if let Some(prepared) = &self.prepared {
            let mut f = ArrayVec::<f64, MAX_EXTERIOR_DOFS>::new();
            for _ in 0..n {
                f.push(0.);
            }
            for k in 0..WALL_SLOTS {
                for (i, value) in f.iter_mut().enumerate() {
                    *value += load.wall[k] * prepared.loads[k][i];
                }
            }
            sink.add_dyn(&self.exterior_dofs, &f);
            return;
        }
        let (c, s) = self.frame(nodes);
        let mut f = vec![0.; n];
        for k in 0..WALL_SLOTS {
            for (i, value) in f.iter_mut().enumerate() {
                *value += load.wall[k] * self.load_slots[k][i];
            }
        }
        for i in 0..self.nodes.len() {
            let x = f[2 * i];
            let y = f[2 * i + 1];
            f[2 * i] = c * x - s * y;
            f[2 * i + 1] = s * x + c * y;
        }
        sink.add_dyn(&self.exterior_dofs, &f);
    }
    pub fn prepare(&mut self, nodes: &NodeView<'_, 2, 3, NodeId>) {
        let n = self.nodes.len() * 2;
        let (c, s) = self.frame(nodes);
        let t = [[c, s], [-s, c]];
        let mut stiffness = vec![0.; n * n];
        let mut loads = Vec::with_capacity(WALL_SLOTS);
        for ni in 0..self.nodes.len() {
            for nj in 0..self.nodes.len() {
                for gi in 0..2 {
                    for gj in 0..2 {
                        let mut value = 0.;
                        for li in 0..2 {
                            for lj in 0..2 {
                                value += t[li][gi]
                                    * self.stiffness[(2 * ni + li, 2 * nj + lj)]
                                    * t[lj][gj];
                            }
                        }
                        stiffness[(2 * ni + gi) * n + 2 * nj + gj] = value;
                    }
                }
            }
        }
        for local in &self.load_slots {
            let mut global = vec![0.; n];
            for i in 0..self.nodes.len() {
                let x = local[2 * i];
                let y = local[2 * i + 1];
                global[2 * i] = c * x - s * y;
                global[2 * i + 1] = s * x + c * y;
            }
            loads.push(global);
        }
        self.prepared = Some(PreparedAssembly { stiffness, loads });
    }
    pub fn commit(&mut self, _nodes: &NodeView<'_, 2, 3, NodeId>) {}
    pub fn local_force_width(&self) -> usize {
        self.nodes.len() * 2
    }
    pub fn local_force(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> ElementForce {
        let u = self.local_displacement(nodes);
        ElementForce::from_vec((&self.stiffness * u).as_slice().to_vec())
    }
    pub fn local_load_force(
        &self,
        _nodes: &NodeView<'_, 2, 3, NodeId>,
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
    pub fn response(
        &self,
        nodes: &NodeView<'_, 2, 3, NodeId>,
        load: Option<&ElementLoad>,
    ) -> StudWallResponse {
        let ue = self.local_displacement(nodes);
        let mut ui = -(&self.recovery * &ue);
        if let Some(l) = load {
            for k in 0..WALL_SLOTS {
                ui += &self.load_recovery[k] * l.wall[k];
            }
        }
        let mut ufull = ue.as_slice().to_vec();
        ufull.extend(ui.iter().copied());
        let q: Vec<f64> = self
            .piers
            .iter()
            .map(|p| {
                let gamma = p.beta.dot(&DVector::from_column_slice(&ufull));
                self.layout.sheathing.unwrap().ga * p.area * gamma
                    / (p.area + self.layout.sheathing.unwrap().ga * p.c)
            })
            .collect();
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
            r.piers.push(WallPierResponse {
                xa: p.xa,
                xb: p.xb,
                q: q[i],
                gamma,
            });
        }
        let resisting = &self.stiffness * &ue;
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
    fn frame(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> (f64, f64) {
        let a = nodes.get(self.corners[0]).coords;
        let b = nodes.get(self.corners[1]).coords;
        let l = (b[0] - a[0]).hypot(b[1] - a[1]);
        ((b[0] - a[0]) / l, (b[1] - a[1]) / l)
    }
    fn local_displacement(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> DVector<f64> {
        let (c, s) = self.frame(nodes);
        DVector::from_iterator(
            self.nodes.len() * 2,
            self.nodes.iter().flat_map(|id| {
                let d = nodes.get(*id).displacement;
                [c * d[0] + s * d[1], -s * d[0] + c * d[1]]
            }),
        )
    }
}

// Internal frame builder uses exterior DOFs first and internal DOFs after them.
struct FrameBuild {
    ndof: usize,
    k: DMatrix<f64>,
}
impl FrameBuild {
    fn new(ne: usize) -> Self {
        Self {
            ndof: ne,
            k: DMatrix::zeros(ne, ne),
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
        let kei = self.k.view((0, ne), (ne, ni)).into_owned();
        let kie = self.k.view((ne, 0), (ni, ne)).into_owned();
        let kii = self.k.view((ne, ne), (ni, ni)).into_owned();
        let lu = kii.clone().lu();
        let Some(r) = lu.solve(&kie) else {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        };
        let kc = kee - &kei * &r;
        if !kc.iter().all(|x| x.is_finite()) {
            return Err(StudWallError("internal stiffness K_ii is singular".into()));
        }
        Ok((kc, r, kii))
    }
    fn condense_load(
        &self,
        ne: usize,
        p: &DVector<f64>,
    ) -> Result<(DVector<f64>, DVector<f64>), StudWallError> {
        let pe = p.rows(0, ne).into_owned();
        let pi = p.rows(ne, self.ndof - ne).into_owned();
        if pi.is_empty() {
            return Ok((pe, DVector::zeros(0)));
        }
        let w = self
            .k
            .view((ne, ne), (self.ndof - ne, self.ndof - ne))
            .into_owned()
            .lu()
            .solve(&pi)
            .ok_or_else(|| StudWallError("internal stiffness K_ii is singular".into()))?;
        let kei = self.k.view((0, ne), (ne, self.ndof - ne));
        Ok((pe - kei * &w, w))
    }
}

fn dedup_x(xs: &mut Vec<f64>, l: f64) {
    xs.sort_by(f64::total_cmp);
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * l.max(1.));
}
fn validate_layout(l: &StudWallLayout, stations: &[WallStation]) -> Result<(), StudWallError> {
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
        if i > 0
            && (o.x0 < l.openings[i - 1].x1
                || !l
                    .studs
                    .iter()
                    .any(|s| s.x > l.openings[i - 1].x1 && s.x < o.x0))
        {
            return Err(StudWallError(format!(
                "openings[{i}] overlaps or has no full-height stud separation"
            )));
        }
        if l.studs.iter().any(|s| s.x == o.x0 || s.x == o.x1) {
            return Err(StudWallError(format!(
                "openings[{i}] x0/x1 coincides with a stud row"
            )));
        }
    }
    if let Some(s) = l.sheathing {
        if !s.ga.is_finite() || s.ga <= 0. {
            return Err(StudWallError(
                "sheathing.ga must be positive and finite".into(),
            ));
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
fn make_chord_nodes(
    model: &mut FrameBuild,
    xs: &[f64],
    y: f64,
    chord: WallChord,
    nslot: usize,
    corners: &[NodeId; 4],
    stations: &[WallStation],
    nodes: &[NodeId],
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
#[allow(clippy::too_many_arguments)]
fn add_chord_elements(
    model: &mut FrameBuild,
    nodes: &[(f64, [usize; 3])],
    ea: f64,
    ei: f64,
    dist: ChordDistribution,
    chord: WallChord,
    beams: &mut Vec<BeamRec>,
    axial: &mut Vec<ChordAxialRec>,
) {
    for pair in nodes.windows(2) {
        let l = pair[1].0 - pair[0].0;
        if chord == WallChord::Top && dist == ChordDistribution::Direct {
            model.add_bar(pair[0].1[0], pair[1].1[0], ea / l);
        } else {
            add_beam3(
                model, pair[0].0, pair[1].0, pair[0].1, pair[1].1, ea, ei, chord, beams,
            );
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
    let Some(sh) = l.sheathing else {
        return Ok(Vec::new());
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
        let fac = area * area / (area / sh.ga + c);
        for i in 0..model.ndof {
            for j in 0..model.ndof {
                model.k[(i, j)] += fac * beta[i] * beta[j];
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
    model: &mut FrameBuild,
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
        let pair = nodes
            .windows(2)
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
            for w in nodes.windows(2) {
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
    let _ = model;
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
