use std::collections::HashMap;

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

/// An element load kind — currently just `ElasticBeamColumn`'s uniform
/// transverse load, the one element-load case that exists (§3.4). Closed
/// enum (§2.1); grow it only when a second element-load type is actually
/// needed.
#[derive(Debug, Clone, Copy)]
pub enum ElementLoad {
    /// Uniform transverse load (force/length) in the element's local +y
    /// direction.
    UniformTransverse(f64),
}

/// A named collection of reference loads (nodal and element), scaled by a
/// `LoadSeries` function of pseudo-time — Xara/OpenSees's `LoadPattern`.
/// `Domain` owns any number of these; multiple patterns coexist and sum
/// (`Domain::assemble_reference_load`), each independently scaled, so e.g.
/// gravity can be ramped to full scale and then held constant
/// (`Domain::hold_pattern_constant`) while a separate lateral pattern
/// continues ramping in a later analysis phase.
#[derive(Debug, Clone)]
pub(crate) struct LoadPattern {
    series: LoadSeries,
    scale_factor: f64,
    /// `Some(frozen_factor)` once `Domain::hold_pattern_constant` has been
    /// called — computed once, at freeze time, from the exact pseudo-time
    /// the caller was at (not lazily cached on every assemble call, which
    /// would require `assemble_reference_load` to take `&mut self` — see
    /// implementation-plan discussion). `None` means "still driven by
    /// `series`".
    frozen_factor: Option<f64>,
    nodal_loads: HashMap<NodeId, [f64; NDF]>,
    element_loads: HashMap<ElementId, ElementLoad>,
}

impl LoadPattern {
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
        self.frozen_factor.unwrap_or_else(|| self.series.factor(pseudo_time) * self.scale_factor)
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

    pub(crate) fn hold_constant(&mut self, pseudo_time: f64) {
        self.frozen_factor = Some(self.factor(pseudo_time));
    }

    pub(crate) fn add_nodal_load(&mut self, node: NodeId, dof: usize, value: f64) {
        self.nodal_loads.entry(node).or_insert([0.0; NDF])[dof] = value;
    }

    pub(crate) fn add_element_load(&mut self, element: ElementId, load: ElementLoad) {
        self.element_loads.insert(element, load);
    }

    pub(crate) fn nodal_load(&self, node: NodeId) -> Option<&[f64; NDF]> {
        self.nodal_loads.get(&node)
    }

    pub(crate) fn element_load(&self, element: ElementId) -> Option<&ElementLoad> {
        self.element_loads.get(&element)
    }
}
