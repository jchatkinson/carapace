use crate::model::LoadSeries;

/// Uniform ground-motion excitation in one direction — Xara/OpenSees's
/// `UniformExcitation`. `TransientAnalysis` uses this to add the effective
/// inertial force `-M·ι·ag(t)` (standard formulation, e.g. Chopra,
/// "Dynamics of Structures") to each Newmark step's load, where `ι` is the
/// unit influence vector selecting every free DOF in `direction`
/// (`Domain::direction_incidence`) and `ag(t)` is this motion's
/// acceleration at time `t`. An accelerogram is exactly a `LoadSeries::
/// Path` (piecewise-linear interpolation over recorded points) — no
/// separate series type needed; this reuses the same `LoadSeries` a static
/// `LoadPattern` would.
#[derive(Debug, Clone)]
pub struct GroundMotion {
    pub direction: usize,
    pub series: LoadSeries,
    pub scale_factor: f64,
}

impl GroundMotion {
    /// `direction`: which DOF index (0 = x, 1 = y) the acceleration acts
    /// along, applied uniformly to every node's mass in that direction.
    pub fn new(direction: usize, series: LoadSeries) -> Self {
        GroundMotion {
            direction,
            series,
            scale_factor: 1.0,
        }
    }

    /// E.g. a PGA scale on top of `series`'s own values.
    pub fn with_scale_factor(mut self, scale_factor: f64) -> Self {
        self.scale_factor = scale_factor;
        self
    }

    pub(crate) fn acceleration(&self, time: f64) -> f64 {
        self.series.factor(time) * self.scale_factor
    }
}
