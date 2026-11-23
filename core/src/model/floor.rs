//! Statically determinate one-way floor load transfer to bearing lines.
//!
//! This is a load-transfer calculation, not a floor element. See `docs/stud-wall.md`.

use nalgebra::{DMatrix, DVector};

use std::collections::BTreeMap;

use super::WallLoad;

/// Floor transfer input error, with the invalid item named.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FloorError(pub String);

impl std::fmt::Display for FloorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for FloorError {}

/// Rectangular floor bounds in span-aligned plan coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloorBounds {
    pub s0: f64,
    pub s1: f64,
    pub t0: f64,
    pub t1: f64,
}

/// Bearing line, indexed by its position in the input support vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloorSupport {
    pub s: f64,
    pub t0: f64,
    pub t1: f64,
}

/// Floor load in the floor's span-aligned coordinates; positive values act down. `slot` is a generic load-case
/// index: the transfer keeps each case's reactions apart (`SupportReaction::wall_loads`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FloorLoad {
    /// Pressure over the whole floor or a rectangular patch.
    Area {
        slot: usize,
        w: f64,
        region: Option<[f64; 4]>,
    },
    /// Line across joists, parallel to the supports.
    Line {
        slot: usize,
        w: f64,
        s: f64,
        t0: f64,
        t1: f64,
    },
    /// Line along one joist; its load is concentrated at `t`.
    LineAlong {
        slot: usize,
        w: f64,
        t: f64,
        s0: f64,
        s1: f64,
    },
    /// Point load from a post or other concentrated load.
    Point { slot: usize, p: f64, s: f64, t: f64 },
}

/// Joist continuity assumption.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FloorContinuity {
    /// Each support interval shares load by tributary statics.
    Simple,
    /// Constant-EI beam on rigid supports, including any overhangs.
    Continuous,
}

/// A support's transferred floor loads and tributary report.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SupportReaction {
    /// Clipped start coordinate of the corresponding support.
    pub support_t0: f64,
    pub lines: Vec<ReactionLine>,
    pub points: Vec<ReactionPoint>,
    /// `(t0, t1, left support index, right support index)` for each simple strip.
    pub trib: Vec<(f64, f64, usize, usize)>,
}

/// Uniform reaction over a support interval, positive down.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionLine {
    pub slot: usize,
    pub t0: f64,
    pub t1: f64,
    pub w: f64,
}

/// Concentrated reaction at a support coordinate, positive down.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionPoint {
    pub slot: usize,
    pub t: f64,
    pub p: f64,
}

/// Reactions in input support order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FloorReactions {
    pub per_support: Vec<SupportReaction>,
}

impl SupportReaction {
    /// Map floor reactions to the supporting wall's local coordinates. `x_at_t0` is the wall coordinate at
    /// the support's floor `t0`; `reversed` maps increasing `t` toward decreasing wall `x`. Floor reactions
    /// are positive down, while wall local `y` is up, so emitted `TopLine` and `TopPoint` forces are negative.
    ///
    /// The result maps each load case (the loads' `slot`) that has a reaction to its wall loads.
    pub fn wall_loads(&self, x_at_t0: f64, reversed: bool) -> BTreeMap<usize, Vec<WallLoad>> {
        let x = |t: f64| x_at_t0 + (if reversed { -1.0 } else { 1.0 }) * (t - self.support_t0);
        let mut loads: BTreeMap<usize, Vec<WallLoad>> = BTreeMap::new();
        for r in &self.lines {
            let (x0, x1) = (x(r.t0), x(r.t1));
            loads.entry(r.slot).or_default().push(WallLoad::TopLine {
                x0: x0.min(x1),
                x1: x0.max(x1),
                wx0: 0.0,
                wx1: 0.0,
                wy0: -r.w,
                wy1: -r.w,
            });
        }
        for r in &self.points {
            loads.entry(r.slot).or_default().push(WallLoad::TopPoint {
                x: x(r.t),
                fx: 0.0,
                fy: -r.p,
            });
        }
        loads
    }

    /// Return this slot's uniform beam-load pieces. Split the beam at each returned end before applying
    /// loads with `ElementLoad::beam_uniform` or `beam_taper`.
    pub fn beam_segments(&self, slot: usize) -> Vec<(f64, f64, f64)> {
        self.lines
            .iter()
            .filter(|r| r.slot == slot)
            .map(|r| (r.t0, r.t1, r.w))
            .collect()
    }
}

/// Transfer floor loads to bearing lines. Supports are clipped to the floor width; output preserves their
/// input order. All geometry and loads must be finite, and all load coordinates must lie within the floor.
pub fn transfer_floor(
    bounds: FloorBounds,
    supports: &[FloorSupport],
    loads: &[FloorLoad],
    continuity: FloorContinuity,
) -> Result<FloorReactions, FloorError> {
    validate_bounds(bounds)?;
    let mut supports = supports.to_vec();
    for (i, s) in supports.iter_mut().enumerate() {
        finite(&[s.s, s.t0, s.t1], &format!("support {i}"))?;
        if s.s < bounds.s0
            || s.s > bounds.s1
            || s.t0 >= s.t1
            || s.t1 <= bounds.t0
            || s.t0 >= bounds.t1
        {
            return Err(err(format!(
                "support {i} is outside the floor or has non-positive width"
            )));
        }
        s.t0 = s.t0.max(bounds.t0);
        s.t1 = s.t1.min(bounds.t1);
    }
    for (i, load) in loads.iter().enumerate() {
        validate_load(bounds, *load, i)?;
    }
    let mut tcuts = vec![bounds.t0, bounds.t1];
    for s in &supports {
        tcuts.extend([s.t0, s.t1]);
    }
    for l in loads {
        match *l {
            FloorLoad::Area {
                region: Some([_, _, a, b]),
                ..
            } => tcuts.extend([a, b]),
            FloorLoad::Line { t0, t1, .. } => tcuts.extend([t0, t1]),
            _ => {}
        }
    }
    sort_dedup(&mut tcuts);
    let mut out = FloorReactions {
        per_support: vec![SupportReaction::default(); supports.len()],
    };
    for (reaction, support) in out.per_support.iter_mut().zip(&supports) {
        reaction.support_t0 = support.t0;
    }
    for pair in tcuts.windows(2) {
        let (t0, t1) = (pair[0], pair[1]);
        if t1 - t0 <= tol(t0, t1) {
            continue;
        }
        let mid = (t0 + t1) * 0.5;
        let mut active: Vec<(usize, f64)> = supports
            .iter()
            .enumerate()
            .filter(|(_, s)| mid > s.t0 - tol(mid, s.t0) && mid < s.t1 + tol(mid, s.t1))
            .map(|(i, s)| (i, s.s))
            .collect();
        active.sort_by(|a, b| a.1.total_cmp(&b.1));
        if active
            .windows(2)
            .any(|w| (w[0].1 - w[1].1).abs() <= tol(w[0].1, w[1].1))
        {
            return Err(err("two supports coincide across a loaded strip"));
        }
        let mut dist = Vec::new();
        let mut concentrated = Vec::new();
        for load in loads {
            match *load {
                FloorLoad::Area { slot, w, region } => {
                    let [s0, s1, a, b] =
                        region.unwrap_or([bounds.s0, bounds.s1, bounds.t0, bounds.t1]);
                    if mid >= a && mid <= b {
                        dist.push((slot, w, s0, s1));
                    }
                }
                FloorLoad::Line {
                    slot,
                    w,
                    s,
                    t0: a,
                    t1: b,
                } if mid >= a && mid <= b => concentrated.push((slot, s, w)),
                _ => {}
            }
        }
        if dist.is_empty() && concentrated.is_empty() {
            continue;
        }
        if active.len() < 2 {
            return Err(err(format!(
                "loaded strip {t0}..{t1} has fewer than two supports"
            )));
        }
        if continuity == FloorContinuity::Simple {
            for &(slot, w, a, b) in &dist {
                distribute(&active, (a, b, w, slot), (t0, t1), &mut out)?;
            }
            for &(slot, s, p) in &concentrated {
                simple_point(&active, s, p, slot, t0, t1, &mut out)?;
            }
            for pair in active.windows(2) {
                for (i, _) in pair {
                    out.per_support[*i]
                        .trib
                        .push((t0, t1, pair[0].0, pair[1].0));
                }
            }
        } else {
            continuous(&active, bounds, &dist, &concentrated, t0, t1, &mut out)?;
        }
    }
    for load in loads {
        let (slot, t, a, z, p, is_line) = match *load {
            FloorLoad::LineAlong { slot, w, t, s0, s1 } => (slot, t, s0, s1, w, true),
            FloorLoad::Point { slot, p, s, t } => (slot, t, s, s, p, false),
            _ => continue,
        };
        let mut active: Vec<(usize, f64)> = supports
            .iter()
            .enumerate()
            .filter(|(_, support)| {
                t >= support.t0 - tol(t, support.t0) && t <= support.t1 + tol(t, support.t1)
            })
            .map(|(i, support)| (i, support.s))
            .collect();
        active.sort_by(|a, b| a.1.total_cmp(&b.1));
        if active.len() < 2 {
            return Err(err(format!(
                "concentrated load at t={t} has fewer than two supports"
            )));
        }
        if continuity == FloorContinuity::Simple {
            if is_line {
                distribute(&active, (a, z, p, slot), (t, t), &mut out)?;
            } else {
                simple_point(&active, a, p, slot, t, t, &mut out)?;
            }
        } else {
            let dist = if is_line {
                vec![(slot, p, a, z)]
            } else {
                Vec::new()
            };
            let points = if is_line {
                Vec::new()
            } else {
                vec![(slot, a, p)]
            };
            continuous(&active, bounds, &dist, &points, t, t, &mut out)?;
        }
    }
    for r in &mut out.per_support {
        merge_lines(&mut r.lines);
    }
    Ok(out)
}

/// Transfer a horizontal diaphragm area load to wall lines as simple spans. Supports are the wall lines
/// perpendicular to the load direction; the returned support reactions are positive in the input load's
/// direction. This uses the same statics as gravity `transfer_floor`, with the load interpreted as force per
/// floor area rather than weight.
pub fn flexible_diaphragm_loads(
    bounds: FloorBounds,
    supports: &[FloorSupport],
    area_loads: &[FloorLoad],
) -> Result<FloorReactions, FloorError> {
    transfer_floor(bounds, supports, area_loads, FloorContinuity::Simple)
}

/// Lump a piecewise-uniform line reaction onto ordered wall-top nodes using midpoint tributary boundaries.
/// The returned nodal forces sum exactly to the integrated line reaction.
pub fn lump_line_reaction(nodes: &[f64], line: ReactionLine) -> Result<Vec<f64>, FloorError> {
    if nodes.is_empty()
        || nodes.windows(2).any(|w| !w[0].is_finite() || w[1] <= w[0])
        || !line.t0.is_finite()
        || !line.t1.is_finite()
        || line.t0 >= line.t1
        || !line.w.is_finite()
        || nodes.first().is_some_and(|x| *x > line.t0)
        || nodes.last().is_some_and(|x| *x < line.t1)
    {
        return Err(err(
            "line reaction and node coordinates must be finite, ordered, and covered by nodes",
        ));
    }
    let mut out = vec![0.; nodes.len()];
    for i in 0..nodes.len() {
        let a = if i == 0 {
            line.t0
        } else {
            (nodes[i - 1] + nodes[i]) * 0.5
        }
        .max(line.t0);
        let b = if i + 1 == nodes.len() {
            line.t1
        } else {
            (nodes[i] + nodes[i + 1]) * 0.5
        }
        .min(line.t1);
        out[i] = (b - a).max(0.) * line.w;
    }
    Ok(out)
}

fn continuous(
    active: &[(usize, f64)],
    b: FloorBounds,
    dist: &[(usize, f64, f64, f64)],
    points: &[(usize, f64, f64)],
    t0: f64,
    t1: f64,
    out: &mut FloorReactions,
) -> Result<(), FloorError> {
    let slots: Vec<usize> = dist
        .iter()
        .map(|x| x.0)
        .chain(points.iter().map(|x| x.0))
        .collect();
    let mut unique_slots = slots;
    unique_slots.sort_unstable();
    unique_slots.dedup();
    for slot in unique_slots {
        let dslot: Vec<_> = dist.iter().copied().filter(|x| x.0 == slot).collect();
        let pslot: Vec<_> = points.iter().copied().filter(|x| x.0 == slot).collect();
        let mut xs = vec![b.s0, b.s1];
        xs.extend(active.iter().map(|x| x.1));
        for &(_, _, a, z) in &dslot {
            xs.extend([a, z]);
        }
        xs.extend(pslot.iter().map(|x| x.1));
        sort_dedup(&mut xs);
        let n = xs.len();
        let mut k = DMatrix::<f64>::zeros(2 * n, 2 * n);
        let mut f = DVector::<f64>::zeros(2 * n);
        for e in 0..n - 1 {
            let (a, z) = (xs[e], xs[e + 1]);
            let l = z - a;
            let ke = beam_k(l);
            for i in 0..4 {
                for j in 0..4 {
                    let gi = 2 * e + i;
                    let gj = 2 * e + j;
                    k[(gi, gj)] += ke[(i, j)];
                }
            }
            for &(_, w, q0, q1) in &dslot {
                let ov0 = a.max(q0);
                let ov1 = z.min(q1);
                if ov1 - ov0 > tol(ov0, ov1) {
                    let q = w;
                    let fe = [-q * l / 2., -q * l * l / 12., -q * l / 2., q * l * l / 12.];
                    for i in 0..4 {
                        f[2 * e + i] += fe[i];
                    }
                }
            }
        }
        for &(_, s, p) in &pslot {
            let idx = nearest(&xs, s);
            f[2 * idx] -= p;
        }
        let constrained: Vec<usize> = active.iter().map(|x| 2 * nearest(&xs, x.1)).collect();
        let free: Vec<usize> = (0..2 * n).filter(|i| !constrained.contains(i)).collect();
        let kff = DMatrix::from_fn(free.len(), free.len(), |i, j| k[(free[i], free[j])]);
        let ff = DVector::from_iterator(free.len(), free.iter().map(|&i| f[i]));
        let ufree = kff
            .lu()
            .solve(&ff)
            .ok_or_else(|| err("continuous beam stiffness is singular"))?;
        let mut u = DVector::<f64>::zeros(2 * n);
        for (i, &g) in free.iter().enumerate() {
            u[g] = ufree[i];
        }
        let residual = &k * u - f;
        for &(support_index, x) in active {
            let p = residual[2 * nearest(&xs, x)];
            if p.abs() > 1e-14 {
                if t1 > t0 {
                    out.per_support[support_index]
                        .lines
                        .push(ReactionLine { slot, t0, t1, w: p });
                } else {
                    out.per_support[support_index]
                        .points
                        .push(ReactionPoint { slot, t: t0, p });
                }
            }
        }
    }
    Ok(())
}

fn beam_k(l: f64) -> DMatrix<f64> {
    let l2 = l * l;
    let l3 = l2 * l;
    DMatrix::from_row_slice(
        4,
        4,
        &[
            12. / l3,
            6. / l2,
            -12. / l3,
            6. / l2,
            6. / l2,
            4. / l,
            -6. / l2,
            2. / l,
            -12. / l3,
            -6. / l2,
            12. / l3,
            -6. / l2,
            6. / l2,
            2. / l,
            -6. / l2,
            4. / l,
        ],
    )
}

fn simple_point(
    active: &[(usize, f64)],
    s: f64,
    p: f64,
    slot: usize,
    t0: f64,
    t1: f64,
    out: &mut FloorReactions,
) -> Result<(), FloorError> {
    if s < active[0].1 - tol(s, active[0].1) {
        let (i, x) = active[0];
        let (j, y) = active[1];
        add_reaction(out, i, slot, t0, t1, p * (y - s) / (y - x));
        add_reaction(out, j, slot, t0, t1, p * (s - x) / (y - x));
    } else if s > active.last().unwrap().1 + tol(s, active.last().unwrap().1) {
        let (i, x) = active[active.len() - 2];
        let (j, y) = active[active.len() - 1];
        add_reaction(out, i, slot, t0, t1, p * (y - s) / (y - x));
        add_reaction(out, j, slot, t0, t1, p * (s - x) / (y - x));
    } else {
        let i = active
            .partition_point(|x| x.1 < s)
            .saturating_sub(1)
            .min(active.len() - 2);
        let (li, x) = active[i];
        let (ri, y) = active[i + 1];
        add_reaction(out, li, slot, t0, t1, p * (y - s) / (y - x));
        add_reaction(out, ri, slot, t0, t1, p * (s - x) / (y - x));
    }
    Ok(())
}
fn distribute(
    active: &[(usize, f64)],
    (a, b, w, slot): (f64, f64, f64, usize),
    (t0, t1): (f64, f64),
    out: &mut FloorReactions,
) -> Result<(), FloorError> {
    let mut cuts = vec![a, b];
    cuts.extend(active.iter().map(|x| x.1).filter(|x| *x > a && *x < b));
    sort_dedup(&mut cuts);
    for c in cuts.windows(2) {
        let (x, y) = (c[0], c[1]);
        let q = w * (y - x);
        let mid = (x + y) * 0.5;
        let (li, lx, ri, rx) = if mid < active[0].1 {
            (active[0].0, active[0].1, active[1].0, active[1].1)
        } else if mid > active.last().unwrap().1 {
            (
                active[active.len() - 2].0,
                active[active.len() - 2].1,
                active.last().unwrap().0,
                active.last().unwrap().1,
            )
        } else {
            let j = active
                .partition_point(|p| p.1 <= mid)
                .saturating_sub(1)
                .min(active.len() - 2);
            (active[j].0, active[j].1, active[j + 1].0, active[j + 1].1)
        };
        let centroid = (x + y) * 0.5;
        let left = q * (rx - centroid) / (rx - lx);
        let right = q - left;
        add_line(out, li, slot, t0, t1, left);
        add_line(out, ri, slot, t0, t1, right);
    }
    Ok(())
}
fn add_line(out: &mut FloorReactions, i: usize, slot: usize, t0: f64, t1: f64, p: f64) {
    if t1 > t0 {
        out.per_support[i]
            .lines
            .push(ReactionLine { slot, t0, t1, w: p });
    } else {
        out.per_support[i]
            .points
            .push(ReactionPoint { slot, t: t0, p });
    }
}
fn add_point(out: &mut FloorReactions, i: usize, slot: usize, t0: f64, t1: f64, p: f64) {
    out.per_support[i].points.push(ReactionPoint {
        slot,
        t: (t0 + t1) * 0.5,
        p,
    });
}
fn add_reaction(out: &mut FloorReactions, i: usize, slot: usize, t0: f64, t1: f64, p: f64) {
    if t1 > t0 {
        add_line(out, i, slot, t0, t1, p);
    } else {
        add_point(out, i, slot, t0, t1, p);
    }
}
fn merge_lines(lines: &mut Vec<ReactionLine>) {
    lines.sort_by(|a, b| a.slot.cmp(&b.slot).then(a.t0.total_cmp(&b.t0)));
    let mut merged: Vec<ReactionLine> = Vec::new();
    for r in lines.drain(..) {
        if let Some(last) = merged.last_mut() {
            if last.slot == r.slot
                && (last.t1 - r.t0).abs() <= tol(last.t1, r.t0)
                && (last.w - r.w).abs() <= tol(last.w, r.w)
            {
                last.t1 = r.t1;
                continue;
            }
        }
        merged.push(r);
    }
    *lines = merged;
}
fn validate_bounds(b: FloorBounds) -> Result<(), FloorError> {
    finite(&[b.s0, b.s1, b.t0, b.t1], "floor bounds")?;
    if b.s0 >= b.s1 || b.t0 >= b.t1 {
        return Err(err("floor bounds must have positive width"));
    }
    Ok(())
}
fn validate_load(b: FloorBounds, l: FloorLoad, i: usize) -> Result<(), FloorError> {
    let vals = match l {
        FloorLoad::Area { w, region, .. } => {
            let r = region.unwrap_or([b.s0, b.s1, b.t0, b.t1]);
            if r[0] < b.s0
                || r[1] > b.s1
                || r[2] < b.t0
                || r[3] > b.t1
                || r[0] >= r[1]
                || r[2] >= r[3]
            {
                return Err(err(format!(
                    "area load {i} is outside the floor or has non-positive extent"
                )));
            }
            vec![w, r[0], r[1], r[2], r[3]]
        }
        FloorLoad::Line { w, s, t0, t1, .. } => {
            if s < b.s0 || s > b.s1 || t0 < b.t0 || t1 > b.t1 || t0 >= t1 {
                return Err(err(format!("line load {i} is outside the floor")));
            }
            vec![w, s, t0, t1]
        }
        FloorLoad::LineAlong { w, t, s0, s1, .. } => {
            if t < b.t0 || t > b.t1 || s0 < b.s0 || s1 > b.s1 || s0 >= s1 {
                return Err(err(format!("line-along load {i} is outside the floor")));
            }
            vec![w, t, s0, s1]
        }
        FloorLoad::Point { p, s, t, .. } => {
            if s < b.s0 || s > b.s1 || t < b.t0 || t > b.t1 {
                return Err(err(format!("point load {i} is outside the floor")));
            }
            vec![p, s, t]
        }
    };
    finite(&vals, &format!("load {i}"))
}
fn finite(v: &[f64], name: &str) -> Result<(), FloorError> {
    if v.iter().any(|x| !x.is_finite()) {
        Err(err(format!("{name} contains a non-finite value")))
    } else {
        Ok(())
    }
}
fn err(s: impl Into<String>) -> FloorError {
    FloorError(s.into())
}
fn tol(a: f64, b: f64) -> f64 {
    1e-9 * a.abs().max(b.abs()).max(1.0)
}
fn sort_dedup(v: &mut Vec<f64>) {
    v.sort_by(f64::total_cmp);
    let mut i = 1;
    while i < v.len() {
        if (v[i] - v[i - 1]).abs() <= tol(v[i], v[i - 1]) {
            v.remove(i);
        } else {
            i += 1;
        }
    }
}
fn nearest(xs: &[f64], x: f64) -> usize {
    xs.iter()
        .enumerate()
        .min_by(|a, b| (a.1 - x).abs().total_cmp(&(b.1 - x).abs()))
        .unwrap()
        .0
}
