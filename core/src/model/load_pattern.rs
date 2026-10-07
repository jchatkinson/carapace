use std::collections::HashMap;
use std::hash::Hash;

use slotmap::new_key_type;

use super::{ElementId, NodeId, NDF};

new_key_type! {
    /// Generational index into `Domain`'s load-pattern store (§2.2).
    pub struct LoadPatternId;
}

/// The pseudo-time-to-load-factor function a `LoadPattern` scales its
/// reference loads by — Xara/OpenSees's `TimeSeries` family, cut down to
/// exactly the three shapes Carapace's scope actually needs (`Trig`/
/// `Rectangular`/`Pulse` etc. are for shaping a *dynamic* excitation
/// pattern, not a quasi-static load/displacement protocol — out of scope
/// until something needs them, §3). Closed enum (§2.1).
#[derive(Debug, Clone)]
pub enum LoadSeries {
    /// `factor(t) = 1.0` for every `t` — the pattern is always fully
    /// applied, independent of pseudo-time.
    Constant,
    /// `factor(t) = slope * t`. `slope: 1.0` reproduces the single-pattern,
    /// `load_factor`-scales-everything behavior every `Analysis` had before
    /// multi-pattern support existed.
    Linear { slope: f64 },
    /// Piecewise-linear interpolation between `(times[i], factors[i])`
    /// control points — Xara's `PathSeries`/`PathTimeSeries`. This is the
    /// load-path / quasi-static cyclic-protocol vehicle: drive it with
    /// `Integrator::LoadControl` (stepping pseudo-time along `times`) for a
    /// load-controlled protocol, or read an accelerogram's acceleration
    /// values into `factors` for `GroundMotion`. Holds the last point's
    /// value beyond the last control point (Xara's `useLast` default),
    /// and the first point's value before the first. `times` must be
    /// sorted ascending and the same length as `factors` — not enforced by
    /// the type (§2.8 is about runtime failure; a malformed path is a
    /// model-construction error, same class as `AnalysisBuilder::build`'s
    /// `ConstraintHandler` assertion).
    Path { times: Vec<f64>, factors: Vec<f64> },
}

impl LoadSeries {
    pub(crate) fn factor(&self, pseudo_time: f64) -> f64 {
        match self {
            LoadSeries::Constant => 1.0,
            LoadSeries::Linear { slope } => slope * pseudo_time,
            LoadSeries::Path { times, factors } => path_interpolate(times, factors, pseudo_time),
        }
    }

    /// `d(factor)/d(pseudo_time)` at `pseudo_time` — what `Integrator::
    /// DisplacementControl`'s unit-load technique needs (it must know how
    /// the load *responds* to a unit pseudo-time increment, not its
    /// current value). Zero for `Constant` (doesn't respond to pseudo-time
    /// at all — same reason a `hold_pattern_constant`-frozen pattern is
    /// excluded from the sensitivity sum entirely, see `Domain::
    /// assemble_reference_load_sensitivity`), `slope` for `Linear`, the
    /// active segment's slope for `Path`.
    pub(crate) fn slope(&self, pseudo_time: f64) -> f64 {
        match self {
            LoadSeries::Constant => 0.0,
            LoadSeries::Linear { slope } => *slope,
            LoadSeries::Path { times, factors } => path_slope(times, factors, pseudo_time),
        }
    }
}

fn path_segment(times: &[f64], pseudo_time: f64) -> usize {
    // Index `i` such that `pseudo_time` falls in `[times[i], times[i+1])`,
    // clamped to the path's first/last segment outside its range (Xara's
    // `useLast`/pre-first-point behavior). `times` has at least 2 entries
    // here (the `times.len() == 1` case is handled by callers before this
    // is reached).
    let n = times.partition_point(|&t| t <= pseudo_time);
    n.saturating_sub(1).min(times.len() - 2)
}

fn path_interpolate(times: &[f64], factors: &[f64], pseudo_time: f64) -> f64 {
    if times.len() == 1 {
        return factors[0];
    }
    let i = path_segment(times, pseudo_time);
    let (t0, t1) = (times[i], times[i + 1]);
    let (f0, f1) = (factors[i], factors[i + 1]);
    if pseudo_time <= t0 {
        f0
    } else if pseudo_time >= t1 {
        f1
    } else {
        f0 + (f1 - f0) * (pseudo_time - t0) / (t1 - t0)
    }
}

fn path_slope(times: &[f64], factors: &[f64], pseudo_time: f64) -> f64 {
    if times.len() == 1 {
        return 0.0;
    }
    let i = path_segment(times, pseudo_time);
    let (t0, t1) = (times[i], times[i + 1]);
    if pseudo_time < t0 || pseudo_time > t1 {
        return 0.0; // held constant outside the path's range — no response to pseudo-time there
    }
    (factors[i + 1] - factors[i]) / (t1 - t0)
}

/// The load on a 2D element: an additive accumulator with one field group per
/// load kind, so that the sum of every pattern's load on an element (each
/// scaled by its factor) is just a field-wise sum. `Uniform + Body` is simply
/// both fields set. Which kinds an element accepts is decided by
/// `ElementOps::accepts_load`, checked by `Domain::validate` before any load
/// is accumulated.
///
/// Component frames differ by kind: `beam_uniform` is in the element's
/// *local* axes (local `x` is the member axis), as in OpenSees's
/// `-beamUniform`; `body`, `edge_traction` and `edge_pressure` are in *global*
/// axes (pressure along the inward normal).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ElementLoad {
    /// Beam-columns: uniform load (force/length) `[wx, wy]`, `wx` along the
    /// member axis and `wy` transverse in local +y.
    pub beam_uniform: [f64; 2],
    /// Continuum elements: body force per unit volume `[bx, by]`.
    pub body: [f64; 2],
    /// Continuum elements: traction `[tx, ty]` per unit edge length on edge
    /// `k`, which joins local node `k` to `k + 1`.
    pub edge_traction: [[f64; 2]; 4],
    /// Continuum elements: pressure on edge `k`, positive acting into the
    /// element along the inward normal.
    pub edge_pressure: [f64; 4],
}

impl ElementLoad {
    /// A uniform beam load (the former `ElementLoad::Uniform { wx, wy }`).
    pub fn uniform(wx: f64, wy: f64) -> Self {
        ElementLoad {
            beam_uniform: [wx, wy],
            ..Default::default()
        }
    }

    /// A body force per unit volume in global axes.
    pub fn body(bx: f64, by: f64) -> Self {
        ElementLoad {
            body: [bx, by],
            ..Default::default()
        }
    }

    /// A global-axes traction per unit length on edge `edge` (`0..4`).
    pub fn edge_traction(edge: usize, tx: f64, ty: f64) -> Self {
        let mut load = ElementLoad::default();
        load.edge_traction[edge] = [tx, ty];
        load
    }

    /// A pressure on edge `edge` (`0..4`), positive into the element.
    pub fn edge_pressure(edge: usize, pressure: f64) -> Self {
        let mut load = ElementLoad::default();
        load.edge_pressure[edge] = pressure;
        load
    }

    /// Whether any beam-load field is nonzero.
    pub fn has_beam(&self) -> bool {
        self.beam_uniform != [0.0; 2]
    }

    /// Whether the body force is nonzero.
    pub fn has_body(&self) -> bool {
        self.body != [0.0; 2]
    }

    /// Whether any edge `>= edges` carries a traction or pressure.
    pub fn has_edge_from(&self, edges: usize) -> bool {
        (edges..4).any(|k| self.edge_traction[k] != [0.0; 2] || self.edge_pressure[k] != 0.0)
    }

    /// Whether any edge carries a traction or pressure.
    pub fn has_edge(&self) -> bool {
        self.has_edge_from(0)
    }
}

impl std::ops::Add for ElementLoad {
    type Output = ElementLoad;

    /// Field-wise sum — several loads on one element accumulate, as in OpenSees.
    fn add(self, other: ElementLoad) -> ElementLoad {
        ElementLoad {
            beam_uniform: std::array::from_fn(|i| self.beam_uniform[i] + other.beam_uniform[i]),
            body: std::array::from_fn(|i| self.body[i] + other.body[i]),
            edge_traction: std::array::from_fn(|k| {
                std::array::from_fn(|i| self.edge_traction[k][i] + other.edge_traction[k][i])
            }),
            edge_pressure: std::array::from_fn(|k| self.edge_pressure[k] + other.edge_pressure[k]),
        }
    }
}

/// Scalar view of an element load's components (local axes), for recording
/// the load an element actually carries at some pseudo-time.
pub trait ElementLoadComponents {
    /// Number of components (`ElementLoad`: 2, `ElementLoad3`: 3); the bound a
    /// recorder's `component` index is checked against.
    const COUNT: usize;

    /// Component `index` (`ElementLoad`: `0` = `wx`, `1` = `wy`;
    /// `ElementLoad3`: `0..3` = `wx`, `wy`, `wz`); `0.0` past the last.
    fn component(&self, index: usize) -> f64;
}

/// Component layout: `[wx, wy, bx, by, t0x, t0y, p0, t1x, t1y, p1, t2x, t2y, p2, t3x, t3y, p3]`
/// (beam uniform, body force, then traction and pressure per edge).
impl ElementLoadComponents for ElementLoad {
    const COUNT: usize = 16;

    fn component(&self, index: usize) -> f64 {
        match index {
            0 | 1 => self.beam_uniform[index],
            2 | 3 => self.body[index - 2],
            4..=15 => {
                let (edge, slot) = ((index - 4) / 3, (index - 4) % 3);
                match slot {
                    0 | 1 => self.edge_traction[edge][slot],
                    _ => self.edge_pressure[edge],
                }
            }
            _ => 0.0,
        }
    }
}

impl ElementLoadComponents for ElementLoad3 {
    const COUNT: usize = 3;

    fn component(&self, index: usize) -> f64 {
        let ElementLoad3::Uniform { wx, wy, wz } = *self;
        [wx, wy, wz].get(index).copied().unwrap_or(0.0)
    }
}

impl std::ops::Mul<f64> for ElementLoad {
    type Output = ElementLoad;

    fn mul(self, factor: f64) -> ElementLoad {
        ElementLoad {
            beam_uniform: self.beam_uniform.map(|v| v * factor),
            body: self.body.map(|v| v * factor),
            edge_traction: self.edge_traction.map(|t| t.map(|v| v * factor)),
            edge_pressure: self.edge_pressure.map(|v| v * factor),
        }
    }
}

/// `ElementLoad`'s spatial counterpart — a uniform load on
/// `ElasticBeamColumn3`, in the member's local axes (`vec_xz` defines local
/// `y`/`z`): axial `wx` and biaxial transverse `wy`/`wz` (local `y` and `z`
/// both carry bending in a spatial member). Xara/OpenSees's
/// `Beam3dUniformLoad`. `ElementOps::Load` for the spatial profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementLoad3 {
    /// Uniform load (force/length): `wx` along the member axis, `wy`/`wz`
    /// transverse in local +y/+z.
    Uniform { wx: f64, wy: f64, wz: f64 },
}

impl std::ops::Add for ElementLoad3 {
    type Output = ElementLoad3;

    /// See `ElementLoad`'s `Add`.
    fn add(self, other: ElementLoad3) -> ElementLoad3 {
        let (
            ElementLoad3::Uniform {
                wx: ax,
                wy: ay,
                wz: az,
            },
            ElementLoad3::Uniform {
                wx: bx,
                wy: by,
                wz: bz,
            },
        ) = (self, other);
        ElementLoad3::Uniform {
            wx: ax + bx,
            wy: ay + by,
            wz: az + bz,
        }
    }
}

impl std::ops::Mul<f64> for ElementLoad3 {
    type Output = ElementLoad3;

    fn mul(self, factor: f64) -> ElementLoad3 {
        let ElementLoad3::Uniform { wx, wy, wz } = self;
        ElementLoad3::Uniform {
            wx: wx * factor,
            wy: wy * factor,
            wz: wz * factor,
        }
    }
}

/// A named collection of reference loads (nodal and element), scaled by a
/// `LoadSeries` function of pseudo-time — Xara/OpenSees's `LoadPattern`.
/// `Domain` owns any number of these; multiple patterns coexist and sum
/// (`Domain::assemble_reference_load`), each independently scaled, so e.g.
/// gravity can be ramped to full scale and then held constant
/// (`Domain::hold_pattern_constant`) while a separate lateral pattern
/// continues ramping in a later analysis phase.
///
/// Generic over `NDOF` (nodal-load array size), the node/element ID types,
/// and `EL` (the element-load kind — `ElementLoad` for the planar profile,
/// `ElementLoad3` for the spatial one, matching `ElementOps::Load`). `NDOF` defaults to the planar profile so bare `LoadPattern`
/// keeps working unchanged, the same trick `Node`'s defaults use.
#[derive(Debug, Clone)]
pub(crate) struct LoadPattern<
    const NDOF: usize = NDF,
    NId = NodeId,
    EId = ElementId,
    EL = ElementLoad,
> where
    NId: Eq + Hash + Copy,
    EId: Eq + Hash + Copy,
{
    series: LoadSeries,
    scale_factor: f64,
    /// `Some(frozen_factor)` once `Domain::hold_pattern_constant` has been
    /// called — computed once, at freeze time, from the exact pseudo-time
    /// the caller was at (not lazily cached on every assemble call, which
    /// would require `assemble_reference_load` to take `&mut self` — see
    /// implementation-plan discussion). `None` means "still driven by
    /// `series`".
    frozen_factor: Option<f64>,
    nodal_loads: HashMap<NId, [f64; NDOF]>,
    element_loads: HashMap<EId, EL>,
}

impl<const NDOF: usize, NId, EId, EL> LoadPattern<NDOF, NId, EId, EL>
where
    NId: Eq + Hash + Copy,
    EId: Eq + Hash + Copy,
    EL: Copy + std::ops::Add<Output = EL>,
{
    pub(crate) fn has_element_loads(&self) -> bool {
        !self.element_loads.is_empty()
    }

    pub(crate) fn new(series: LoadSeries) -> Self {
        LoadPattern {
            series,
            scale_factor: 1.0,
            frozen_factor: None,
            nodal_loads: HashMap::new(),
            element_loads: HashMap::new(),
        }
    }

    pub(crate) fn with_scale_factor(mut self, scale_factor: f64) -> Self {
        self.scale_factor = scale_factor;
        self
    }

    pub(crate) fn factor(&self, pseudo_time: f64) -> f64 {
        self.frozen_factor
            .unwrap_or_else(|| self.series.factor(pseudo_time) * self.scale_factor)
    }

    /// Zero once frozen (§`LoadSeries::slope`'s doc comment) — a frozen
    /// pattern doesn't respond to further pseudo-time change.
    pub(crate) fn sensitivity(&self, pseudo_time: f64) -> f64 {
        if self.frozen_factor.is_some() {
            0.0
        } else {
            self.series.slope(pseudo_time) * self.scale_factor
        }
    }

    /// Whether this pattern still follows a `LoadSeries::Path` — arc-length
    /// continuation rejects these, since a path's slope isn't valid across
    /// a load reversal or a breakpoint.
    pub(crate) fn is_unfrozen_path(&self) -> bool {
        self.frozen_factor.is_none() && matches!(self.series, LoadSeries::Path { .. })
    }

    pub(crate) fn hold_constant(&mut self, pseudo_time: f64) {
        self.frozen_factor = Some(self.factor(pseudo_time));
    }

    pub(crate) fn add_nodal_load(&mut self, node: NId, dof: usize, value: f64) {
        self.nodal_loads.entry(node).or_insert([0.0; NDOF])[dof] = value;
    }

    /// Adds `load` to whatever this pattern already applies to `element`
    /// (loads accumulate rather than replace).
    pub(crate) fn add_element_load(&mut self, element: EId, load: EL) {
        self.element_loads
            .entry(element)
            .and_modify(|existing| *existing = *existing + load)
            .or_insert(load);
    }

    pub(crate) fn nodal_load(&self, node: NId) -> Option<&[f64; NDOF]> {
        self.nodal_loads.get(&node)
    }

    pub(crate) fn element_loads(&self) -> impl Iterator<Item = (&EId, &EL)> {
        self.element_loads.iter()
    }

    pub(crate) fn element_load(&self, element: EId) -> Option<&EL> {
        self.element_loads.get(&element)
    }
}

/// Patterns that currently contribute element loads, each with its factor
/// at `pseudo_time` — computed once per assembly so the per-element lookup
/// below doesn't re-evaluate every pattern's time series per element.
pub(crate) fn active_element_patterns<'a, const NDOF: usize, NId, EId, EL>(
    patterns: impl Iterator<Item = &'a LoadPattern<NDOF, NId, EId, EL>>,
    pseudo_time: f64,
) -> Vec<(f64, &'a LoadPattern<NDOF, NId, EId, EL>)>
where
    NId: Eq + Hash + Copy + 'a,
    EId: Eq + Hash + Copy + 'a,
    EL: Copy + std::ops::Add<Output = EL> + 'a,
{
    patterns
        .filter(|pattern| pattern.has_element_loads())
        .map(|pattern| (pattern.factor(pseudo_time), pattern))
        .filter(|(factor, _)| *factor != 0.0)
        .collect()
}

/// An element's *effective* load: every active pattern's load on it, each
/// scaled by that pattern's factor, summed — `None` if no active pattern
/// loads it. What state-dependent elements (`ForceBeamColumn`) need inside
/// their own state determination.
pub(crate) fn effective_element_load<const NDOF: usize, NId, EId, EL>(
    active: &[(f64, &LoadPattern<NDOF, NId, EId, EL>)],
    element: EId,
) -> Option<EL>
where
    NId: Eq + Hash + Copy,
    EId: Eq + Hash + Copy,
    EL: Copy + std::ops::Add<Output = EL> + std::ops::Mul<f64, Output = EL>,
{
    let mut total: Option<EL> = None;
    for (factor, pattern) in active {
        if let Some(load) = pattern.element_load(element) {
            let scaled = *load * *factor;
            total = Some(match total {
                Some(sum) => sum + scaled,
                None => scaled,
            });
        }
    }
    total
}
