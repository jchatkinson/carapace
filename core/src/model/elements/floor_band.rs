//! Linear 2D floor band joining the top of one stud wall to the bottom of the next.
//!
//! Vertical bearing and hold-down connections are separate `ZeroLength` elements. Without those links the
//! band has no vertical stiffness and the stacked model is a mechanism.

use nalgebra::{SMatrix, SVector};

use super::{DofMask, DofRef, NodeView, TangentSink};
use crate::model::{ElementForce, NodeId};

/// Stations and properties for a floor band.
#[derive(Debug, Clone)]
pub struct FloorBandLayout {
    /// Length along the wall line.
    pub length: f64,
    /// Vertical distance between bottom and top rows.
    pub depth: f64,
    /// Sorted local station positions, including both ends.
    pub xs: Vec<f64>,
    /// Shear rigidity `GA` for each bay; `None` means no shear transfer.
    pub bays: Vec<Option<f64>>,
    /// Rim axial stiffness `EA` along each bay.
    pub rim_ea: f64,
}

/// Invalid floor-band layout or node mapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FloorBandError(pub String);
impl std::fmt::Display for FloorBandError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for FloorBandError {}

impl FloorBandLayout {
    /// Uses the union of both walls' stud and pier-end positions. A bay within an opening of `below` has no
    /// shear field, matching segmented shear-wall design.
    pub fn between(
        below: &super::StudWallLayout,
        above: &super::StudWallLayout,
        depth: f64,
        ga: f64,
        rim_ea: f64,
    ) -> Result<Self, FloorBandError> {
        if (below.length - above.length).abs() > 1e-9 * below.length.max(above.length) {
            return Err(FloorBandError("wall lengths must match".into()));
        }
        if !ga.is_finite() || ga <= 0.0 {
            return Err(FloorBandError("ga must be finite and positive".into()));
        }
        let mut xs: Vec<f64> = below
            .studs
            .iter()
            .chain(&above.studs)
            .map(|s| s.x)
            .collect();
        xs.extend(below.pier_ends());
        xs.extend(above.pier_ends());
        xs.sort_by(f64::total_cmp);
        xs.dedup_by(|a, b| (*a - *b).abs() <= 1e-9 * below.length);
        let bays = xs
            .windows(2)
            .map(|w| {
                if below.openings.iter().any(|o| w[0] >= o.x0 && w[1] <= o.x1) {
                    None
                } else {
                    Some(ga)
                }
            })
            .collect();
        let layout = Self {
            length: below.length,
            depth,
            xs,
            bays,
            rim_ea,
        };
        validate_layout(&layout)?;
        Ok(layout)
    }
}

/// One bay's local shear response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloorBandBay {
    /// Left station.
    pub x0: f64,
    /// Right station.
    pub x1: f64,
    /// Shear strain.
    pub gamma: f64,
    /// Shear flow, zero for an inactive bay.
    pub q: f64,
}
/// Floor-band bay shear and rim axial forces (tension positive).
#[derive(Debug, Clone, PartialEq)]
pub struct FloorBandResponse {
    /// Per-bay shear fields.
    pub bays: Vec<FloorBandBay>,
    /// Axial force per rim segment.
    pub rim: Vec<f64>,
}

/// 2D floor-depth strip. Node order is bottom row followed by top row.
#[derive(Debug, Clone)]
pub struct FloorBand {
    bottom: Vec<NodeId>,
    top: Vec<NodeId>,
    nodes: Vec<NodeId>,
    layout: FloorBandLayout,
}

impl FloorBand {
    /// Creates a band from bottom and top nodes in station order. Node geometry is checked by `validate`.
    pub fn new(
        bottom: Vec<NodeId>,
        top: Vec<NodeId>,
        layout: FloorBandLayout,
    ) -> Result<Self, FloorBandError> {
        validate_layout(&layout)?;
        let n = layout.xs.len();
        if bottom.len() != n {
            return Err(FloorBandError(format!(
                "bottom node count must equal station count {n}"
            )));
        }
        if top.len() != n {
            return Err(FloorBandError(format!(
                "top node count must equal station count {n}"
            )));
        }
        let mut nodes = bottom.clone();
        nodes.extend_from_slice(&top);
        let mut unique = std::collections::HashSet::new();
        if nodes.iter().any(|id| !unique.insert(*id)) {
            return Err(FloorBandError(
                "bottom and top nodes must be distinct".into(),
            ));
        }
        Ok(Self {
            bottom,
            top,
            nodes,
            layout,
        })
    }
    /// Connected nodes, bottom row then top row.
    pub fn nodes(&self) -> &[NodeId] {
        &self.nodes
    }
    /// Layout used to build the band.
    pub fn layout(&self) -> &FloorBandLayout {
        &self.layout
    }
    pub fn dof_mask(&self) -> DofMask {
        DofMask::none().with(0).with(1)
    }
    pub fn validate(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> Result<(), &'static str> {
        let p0 = nodes.get(self.bottom[0]).coords;
        let pn = nodes.get(*self.bottom.last().unwrap()).coords;
        let dx = pn[0] - p0[0];
        let dy = pn[1] - p0[1];
        let l = dx.hypot(dy);
        if !l.is_finite() || (l - self.layout.length).abs() > 1e-9 * self.layout.length {
            return Err("floor band bottom endpoints must match layout length");
        }
        let ex = [dx / l, dy / l];
        let ey = [-ex[1], ex[0]];
        for (i, &x) in self.layout.xs.iter().enumerate() {
            for (id, y) in [(self.bottom[i], 0.), (self.top[i], self.layout.depth)] {
                let p = nodes.get(id).coords;
                let lx = (p[0] - p0[0]) * ex[0] + (p[1] - p0[1]) * ex[1];
                let ly = (p[0] - p0[0]) * ey[0] + (p[1] - p0[1]) * ey[1];
                if (lx - x).abs().max((ly - y).abs()) > 1e-9 * self.layout.length {
                    return Err("floor band node does not lie at its local station");
                }
            }
        }
        Ok(())
    }
    pub fn prepare(&mut self, _nodes: &NodeView<'_, 2, 3, NodeId>) {}
    pub fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, 2, 3, NodeId>,
        sink: &mut S,
    ) {
        let ex = self.basis(nodes);
        for i in 0..self.layout.bays.len() {
            let bay = self.bay(nodes, ex, i);
            let dofs: [DofRef<NodeId>; 8] =
                std::array::from_fn(|j| (bay.ids[j / 2], (j % 2) as u8));
            sink.add(&dofs, &bay.k, &bay.r);
        }
    }
    pub fn accepts_load(&self, _load: &crate::model::ElementLoad) -> bool {
        false
    }
    pub fn assemble_load<S: super::VectorSink<NodeId>>(
        &self,
        _nodes: &NodeView<'_, 2, 3, NodeId>,
        _load: Option<&crate::model::ElementLoad>,
        _sink: &mut S,
    ) {
    }
    pub fn commit(
        &mut self,
        _nodes: &NodeView<'_, 2, 3, NodeId>,
        _load: Option<&crate::model::ElementLoad>,
    ) {
    }
    pub fn local_force_width(&self) -> usize {
        self.nodes.len() * 2
    }
    pub fn local_force(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> ElementForce {
        let n = self.layout.xs.len();
        let mut out = vec![0.; self.local_force_width()];
        let ex = self.basis(nodes);
        for i in 0..n - 1 {
            let r = self.bay(nodes, ex, i).r;
            for (j, station) in [i, i + 1, n + i, n + i + 1].into_iter().enumerate() {
                out[2 * station] += r[2 * j];
                out[2 * station + 1] += r[2 * j + 1];
            }
        }
        ElementForce::from_vec(out)
    }
    pub fn local_load_force(
        &self,
        _nodes: &NodeView<'_, 2, 3, NodeId>,
        _load: Option<&crate::model::ElementLoad>,
    ) -> ElementForce {
        ElementForce::from_vec(vec![0.; self.local_force_width()])
    }
    /// Bay `i`'s stiffness and resisting force over `[bottom i, bottom i+1, top i, top i+1]`, its shear
    /// strain and flow, and the rim force (tension positive).
    fn bay(&self, nodes: &NodeView<'_, 2, 3, NodeId>, ex: [f64; 2], i: usize) -> Bay {
        let a = self.layout.xs[i + 1] - self.layout.xs[i];
        let ids = [
            self.bottom[i],
            self.bottom[i + 1],
            self.top[i],
            self.top[i + 1],
        ];
        let u: [f64; 8] = std::array::from_fn(|j| nodes.get(ids[j / 2]).displacement[j % 2]);
        let b = SVector::<f64, 8>::from(self.shear_b(a, ex));
        let gamma = b.dot(&SVector::from(u));
        let ga = self.layout.bays[i].unwrap_or(0.);
        let scale = ga * a * self.layout.depth;
        let mut k = b * b.transpose() * scale;
        let mut r = b * (scale * gamma);
        // Rim bar between the two top nodes.
        let c = self.layout.rim_ea / a;
        let q = [-ex[0], -ex[1], ex[0], ex[1]];
        let strain = (u[6] - u[4]) * ex[0] + (u[7] - u[5]) * ex[1];
        for m in 0..4 {
            for n in 0..4 {
                k[(4 + m, 4 + n)] += c * q[m] * q[n];
            }
            r[4 + m] += c * q[m] * strain;
        }
        Bay {
            ids,
            k,
            r,
            gamma,
            q: ga * gamma,
            rim: c * strain,
        }
    }
    fn basis(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> [f64; 2] {
        let a = nodes.get(self.bottom[0]).coords;
        let b = nodes.get(*self.bottom.last().unwrap()).coords;
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let l = dx.hypot(dy);
        [dx / l, dy / l]
    }
    /// Shear strain operator over `[bottom i, bottom j, top i, top j]`:
    /// `gamma = (u_Ti + u_Tj - u_Bi - u_Bj)/(2d) + (v_Bj - v_Bi)/a`. The rotation term uses the bottom row only,
    /// so the vertical couple `q a d` of the bay lands on the wall below (where the rim bears), and a link
    /// across the band carries the moment at the band top: the base of the wall above, as in O86 Fig. 11.7
    /// (`R21 = F21 (H + h2)/L1`).
    fn shear_b(&self, a: f64, ex: [f64; 2]) -> [f64; 8] {
        let ey = [-ex[1], ex[0]];
        let d = self.layout.depth;
        let local = [
            (-1. / (2. * d), -1. / a),
            (-1. / (2. * d), 1. / a),
            (1. / (2. * d), 0.),
            (1. / (2. * d), 0.),
        ];
        std::array::from_fn(|j| {
            let (u, v) = local[j / 2];
            u * ex[j % 2] + v * ey[j % 2]
        })
    }
    pub fn response(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> FloorBandResponse {
        let ex = self.basis(nodes);
        let (bays, rim) = (0..self.layout.bays.len())
            .map(|i| {
                let bay = self.bay(nodes, ex, i);
                let (x0, x1) = (self.layout.xs[i], self.layout.xs[i + 1]);
                (
                    FloorBandBay {
                        x0,
                        x1,
                        gamma: bay.gamma,
                        q: bay.q,
                    },
                    bay.rim,
                )
            })
            .unzip();
        FloorBandResponse { bays, rim }
    }
}

struct Bay {
    ids: [NodeId; 4],
    k: SMatrix<f64, 8, 8>,
    r: SVector<f64, 8>,
    gamma: f64,
    q: f64,
    rim: f64,
}

fn validate_layout(l: &FloorBandLayout) -> Result<(), FloorBandError> {
    let n = l.xs.len();
    if n < 2 {
        return Err(FloorBandError(
            "at least two xs stations are required".into(),
        ));
    }
    if !l.length.is_finite() || l.length <= 0. {
        return Err(FloorBandError("length must be finite and positive".into()));
    }
    if !l.depth.is_finite() || l.depth <= 0. {
        return Err(FloorBandError("depth must be finite and positive".into()));
    }
    if l.xs[0] != 0.
        || *l.xs.last().unwrap() != l.length
        || l.xs.windows(2).any(|w| !w[0].is_finite() || w[1] <= w[0])
    {
        return Err(FloorBandError(
            "xs must be strictly increasing from 0 to length".into(),
        ));
    }
    if l.bays.len() != n - 1 {
        return Err(FloorBandError(
            "bays length must be station count minus one".into(),
        ));
    }
    if l.bays.iter().flatten().any(|g| !g.is_finite() || *g <= 0.) {
        return Err(FloorBandError("ga must be finite and positive".into()));
    }
    if !l.rim_ea.is_finite() || l.rim_ea < 0. {
        return Err(FloorBandError(
            "rim_ea must be finite and nonnegative".into(),
        ));
    }
    Ok(())
}
