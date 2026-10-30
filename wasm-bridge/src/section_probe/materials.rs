//! Section-probe materials.
//!
//! The probe evaluates a fiber section at a given strain state, with no load history, so every material here is
//! a stateless `stress(strain)`. Two families:
//!
//! - `Core`: any Carapace/OpenSees-compatible material from the probe's `coreMaterials` arena
//!   (`MaterialSpec`), evaluated from its virgin state. For the monotone-from-origin states a section probe
//!   evaluates, that is the material's envelope.
//! - Probe-only curves (`Concrete`, `Steel`, `General`): reversible stress-strain laws used by CalcsApp's
//!   sectional analysis that have no OpenSees equivalent (Whitney block, user curves, strain cutoffs, post-crush
//!   descent, elastic-plastic concrete). They are deliberately **not** `MaterialSpec` variants, so Carapace's
//!   core, arena and wire format stay OpenSees-compatible.
//!
//! Sign convention is OpenSees': tension positive, compression negative. Curves are nominal; capacity-reduction
//! factors are applied by the caller (see the per-group sums on `SectionProbe`).

use carapace_core::model::Material;
use serde::Deserialize;

/// A probe material, referenced by fibers by its position in `SectionProbeConfig::materials`.
#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ProbeMaterialSpec {
    /// A material from the `coreMaterials` arena, by arena index.
    Core {
        index: u32,
    },
    /// `sum(factor * stress(term))` over earlier probe materials (by position, which must be smaller than this
    /// material's own). A bar fiber is `steel - host concrete`: `[[steel, 1], [concrete, -1]]`.
    Combine {
        terms: Vec<(u32, f64)>,
    },
    Concrete(ConcreteSpec),
    Steel(SteelSpec),
    General(GeneralSpec),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum ConcreteModel {
    Elastic,
    ElasticNoTension,
    ElasticPlastic,
    Parabolic,
    Whitney,
    Hognestad,
    UserCurve,
}

/// CalcsApp concrete laws (`materialModels.js`). Magnitudes are used; signs of `fc`, `ec`, `ecu` are ignored.
#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ConcreteSpec {
    pub model: ConcreteModel,
    pub fc: f64,
    pub ec: f64,
    pub ecu: f64,
    /// Strain at the end of the post-crush descent; defaults to `ecu`.
    #[serde(default)]
    pub eu: Option<f64>,
    /// Residual stress magnitude at `eu`.
    #[serde(default)]
    pub fu: f64,
    /// Whitney stress-block factor.
    #[serde(default = "default_alpha")]
    pub alpha: f64,
    /// Whitney beta1; the block starts at `ecu * (1 - beta1)`. Required for `whitney`.
    #[serde(default)]
    pub beta1: Option<f64>,
    #[serde(default)]
    pub tension_capacity: bool,
    #[serde(default)]
    pub ft: f64,
    #[serde(default)]
    pub ets: f64,
    /// `[strain, stress]` pairs for `userCurve`.
    #[serde(default)]
    pub curve_points: Vec<[f64; 2]>,
    #[serde(default)]
    pub drop_after_last_point: bool,
}

fn default_alpha() -> f64 {
    0.85
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum SteelModel {
    Elastic,
    ElasticNoCompression,
    ElasticPlastic,
    UserCurve,
}

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct SteelSpec {
    pub model: SteelModel,
    pub fy: f64,
    pub es: f64,
    /// Post-yield modulus as a fraction of `es`, clamped to `[0, 1]`.
    #[serde(default)]
    pub hardening_ratio: f64,
    /// Compressive strain magnitude beyond which stress drops to zero (0 = no limit).
    #[serde(default, rename = "epsUC")]
    pub eps_uc: f64,
    /// Tensile strain beyond which stress drops to zero (0 = no limit).
    #[serde(default, rename = "epsUT")]
    pub eps_ut: f64,
    #[serde(default)]
    pub curve_points: Vec<[f64; 2]>,
    #[serde(default)]
    pub drop_after_last_point: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum GeneralModel {
    Elastic,
    UserCurve,
}

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct GeneralSpec {
    pub model: GeneralModel,
    pub e: f64,
    #[serde(default, rename = "epsUC")]
    pub eps_uc: f64,
    #[serde(default, rename = "epsUT")]
    pub eps_ut: f64,
    #[serde(default)]
    pub curve_points: Vec<[f64; 2]>,
    #[serde(default)]
    pub drop_after_last_point: bool,
}

/// Why a probe material spec was rejected.
#[derive(Debug, Clone, PartialEq)]
pub enum ProbeMaterialError {
    CoreIndexOutOfRange {
        material: usize,
        index: u32,
    },
    CombineTermNotEarlier {
        material: usize,
        term: u32,
    },
    WhitneyNeedsBeta1 {
        material: usize,
    },
    NonFinite {
        material: usize,
        field: &'static str,
    },
}

/// Piecewise-linear stress-strain table, sorted by strain. Outside the table the end stress is held, or zero
/// when `drop_after_last_point` is set.
#[derive(Debug, Clone)]
pub(super) struct Curve {
    points: Vec<[f64; 2]>,
    drop_outside: bool,
}

impl Curve {
    fn new(points: &[[f64; 2]], drop_outside: bool) -> Self {
        let mut points: Vec<[f64; 2]> = points
            .iter()
            .copied()
            .filter(|p| p[0].is_finite() && p[1].is_finite())
            .collect();
        points.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
        Curve {
            points,
            drop_outside,
        }
    }

    fn stress(&self, e: f64) -> f64 {
        let pts = &self.points;
        match pts.len() {
            0 => return 0.0,
            1 => return pts[0][1],
            _ => {}
        }
        if !e.is_finite() {
            return 0.0;
        }
        if e <= pts[0][0] {
            return if self.drop_outside { 0.0 } else { pts[0][1] };
        }
        let last = pts[pts.len() - 1];
        if e >= last[0] {
            return if self.drop_outside { 0.0 } else { last[1] };
        }
        for w in pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            if e >= a[0] && e <= b[0] {
                let span = b[0] - a[0];
                let t = (e - a[0]) / if span == 0.0 { 1.0 } else { span };
                return a[1] + t * (b[1] - a[1]);
            }
        }
        0.0
    }
}

#[derive(Debug, Clone)]
pub(super) struct ConcreteCurve {
    model: ConcreteModel,
    fc: f64,
    ec: f64,
    ecu: f64,
    eu: f64,
    fu: f64,
    alpha: f64,
    transition: f64,
    e0: f64,
    tension: Option<(f64, f64)>,
    curve: Curve,
}

impl ConcreteCurve {
    fn new(spec: &ConcreteSpec, material: usize) -> Result<Self, ProbeMaterialError> {
        for (field, value) in [("fc", spec.fc), ("ec", spec.ec), ("ecu", spec.ecu)] {
            if !value.is_finite() {
                return Err(ProbeMaterialError::NonFinite { material, field });
            }
        }
        let (fc, ec, ecu) = (spec.fc.abs(), spec.ec.abs(), spec.ecu.abs());
        let ecu = if ecu == 0.0 { 0.003 } else { ecu };
        let transition = match (spec.model, spec.beta1) {
            (ConcreteModel::Whitney, None) => {
                return Err(ProbeMaterialError::WhitneyNeedsBeta1 { material })
            }
            (_, Some(beta1)) => ecu * (1.0 - beta1),
            (_, None) => 0.0,
        };
        let eu = spec.eu.map(f64::abs).filter(|v| *v != 0.0).unwrap_or(ecu);
        Ok(ConcreteCurve {
            model: spec.model,
            fc,
            ec,
            ecu,
            eu: eu.max(ecu),
            fu: spec.fu.abs(),
            alpha: spec.alpha.max(0.0),
            transition,
            e0: if ec > 0.0 { 2.0 * fc / ec } else { 0.0 },
            tension: spec
                .tension_capacity
                .then_some((spec.ft.abs(), spec.ets.abs())),
            curve: Curve::new(&spec.curve_points, spec.drop_after_last_point),
        })
    }

    /// Linear descent from the stress at `ecu` to `-fu` at `eu`; `None` while the strain is within `ecu`.
    fn post_crush(&self, e_comp: f64, stress_at_ecu: f64) -> Option<f64> {
        if e_comp <= self.ecu {
            return None;
        }
        let residual = -self.fu;
        if self.eu <= self.ecu {
            return Some(residual);
        }
        let t = ((e_comp - self.ecu) / (self.eu - self.ecu)).clamp(0.0, 1.0);
        Some(stress_at_ecu + t * (residual - stress_at_ecu))
    }

    fn tension_stress(&self, e: f64) -> f64 {
        let Some((ft, ets)) = self.tension else {
            return 0.0;
        };
        if e <= 0.0 || self.ec <= 0.0 {
            return 0.0;
        }
        let e_peak = ft / self.ec;
        if e <= e_peak || ets <= 0.0 {
            (self.ec * e).min(ft)
        } else {
            (ft - ets * (e - e_peak)).max(0.0)
        }
    }

    /// `-fc * (2r - r^2)`, the Hognestad/Parabolic ascending branch, with `r = strain / e0`.
    fn ascending(&self, e_comp: f64) -> f64 {
        let r = e_comp / self.e0;
        -self.fc * (2.0 * r - r * r)
    }

    fn hognestad_at(&self, strain: f64) -> f64 {
        if self.e0 <= 0.0 {
            return -self.fc;
        }
        if strain <= self.e0 {
            return self.ascending(strain);
        }
        if self.ecu <= self.e0 {
            return -self.fc;
        }
        let r = ((strain - self.e0) / (self.ecu - self.e0)).clamp(0.0, 1.0);
        -self.fc * (1.0 - 0.15 * r)
    }

    fn parabolic_at(&self, strain: f64) -> f64 {
        if self.e0 <= 0.0 {
            return -self.fc;
        }
        if strain <= self.e0 {
            return self.ascending(strain);
        }
        -self.fc
    }

    fn stress(&self, e: f64) -> f64 {
        use ConcreteModel::*;
        match self.model {
            UserCurve => self.curve.stress(e),
            Elastic => {
                if e > 0.0 {
                    return self.ec * e;
                }
                let e_comp = e.abs();
                self.post_crush(e_comp, -self.ec * self.ecu)
                    .unwrap_or(self.ec * e)
            }
            ElasticNoTension => {
                if e > 0.0 {
                    return 0.0;
                }
                let e_comp = e.abs();
                self.post_crush(e_comp, -self.ec * self.ecu)
                    .unwrap_or(self.ec * e)
            }
            ElasticPlastic => {
                if e > 0.0 {
                    return 0.0;
                }
                let e_comp = e.abs();
                self.post_crush(e_comp, -(self.ec * self.ecu).min(self.fc))
                    .unwrap_or((self.ec * e).max(-self.fc))
            }
            Whitney => {
                if e > 0.0 {
                    return 0.0;
                }
                let e_comp = e.abs();
                let at_ecu = if self.ecu >= self.transition {
                    -self.alpha * self.fc
                } else {
                    0.0
                };
                if let Some(post) = self.post_crush(e_comp, at_ecu) {
                    return post;
                }
                if e_comp >= self.transition {
                    -self.alpha * self.fc
                } else {
                    0.0
                }
            }
            Parabolic => {
                if e > 0.0 {
                    return self.tension_stress(e);
                }
                let e_comp = e.abs();
                if let Some(post) = self.post_crush(e_comp, self.parabolic_at(self.ecu)) {
                    return post;
                }
                self.parabolic_at(e_comp)
            }
            Hognestad => {
                if e > 0.0 {
                    return self.tension_stress(e);
                }
                let e_comp = e.abs();
                if let Some(post) = self.post_crush(e_comp, self.hognestad_at(self.ecu)) {
                    return post;
                }
                self.hognestad_at(e_comp)
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct SteelCurve {
    model: SteelModel,
    fy: f64,
    es: f64,
    eh: f64,
    eps_y: f64,
    eps_uc: f64,
    eps_ut: f64,
    curve: Curve,
}

fn cut_off(e: f64, eps_uc: f64, eps_ut: f64) -> bool {
    (eps_uc > 0.0 && e < -eps_uc) || (eps_ut > 0.0 && e > eps_ut)
}

impl SteelCurve {
    fn new(spec: &SteelSpec, material: usize) -> Result<Self, ProbeMaterialError> {
        for (field, value) in [("fy", spec.fy), ("es", spec.es)] {
            if !value.is_finite() {
                return Err(ProbeMaterialError::NonFinite { material, field });
            }
        }
        let (fy, es) = (spec.fy.abs(), spec.es.abs());
        Ok(SteelCurve {
            model: spec.model,
            fy,
            es,
            eh: spec.hardening_ratio.clamp(0.0, 1.0) * es,
            eps_y: if es > 0.0 { fy / es } else { 0.0 },
            eps_uc: spec.eps_uc.abs(),
            eps_ut: spec.eps_ut.abs(),
            curve: Curve::new(&spec.curve_points, spec.drop_after_last_point),
        })
    }

    fn stress(&self, e: f64) -> f64 {
        if cut_off(e, self.eps_uc, self.eps_ut) {
            return 0.0;
        }
        match self.model {
            SteelModel::Elastic => self.es * e,
            SteelModel::ElasticNoCompression => (self.es * e).max(0.0),
            SteelModel::UserCurve => self.curve.stress(e),
            SteelModel::ElasticPlastic => {
                let a = e.abs();
                if a <= self.eps_y {
                    self.es * e
                } else {
                    let sign = if e < 0.0 { -1.0 } else { 1.0 };
                    sign * (self.fy + self.eh * (a - self.eps_y))
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct GeneralCurve {
    model: GeneralModel,
    e: f64,
    eps_uc: f64,
    eps_ut: f64,
    curve: Curve,
}

impl GeneralCurve {
    fn new(spec: &GeneralSpec) -> Self {
        GeneralCurve {
            model: spec.model,
            e: spec.e.abs(),
            eps_uc: spec.eps_uc.abs(),
            eps_ut: spec.eps_ut.abs(),
            curve: Curve::new(&spec.curve_points, spec.drop_after_last_point),
        }
    }

    fn stress(&self, e: f64) -> f64 {
        match self.model {
            GeneralModel::UserCurve => self.curve.stress(e),
            GeneralModel::Elastic => {
                if !e.is_finite() || cut_off(e, self.eps_uc, self.eps_ut) {
                    0.0
                } else {
                    self.e * e
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(super) enum ProbeMaterial {
    Core(Material),
    Combine(Vec<(usize, f64)>),
    Concrete(ConcreteCurve),
    Steel(SteelCurve),
    General(GeneralCurve),
}

impl ProbeMaterial {
    pub(super) fn build(
        specs: &[ProbeMaterialSpec],
        core: &[Material],
    ) -> Result<Vec<ProbeMaterial>, ProbeMaterialError> {
        let mut out = Vec::with_capacity(specs.len());
        for (material, spec) in specs.iter().enumerate() {
            out.push(match spec {
                ProbeMaterialSpec::Core { index } => ProbeMaterial::Core(
                    core.get(*index as usize)
                        .ok_or(ProbeMaterialError::CoreIndexOutOfRange {
                            material,
                            index: *index,
                        })?
                        .clone(),
                ),
                ProbeMaterialSpec::Combine { terms } => {
                    for (term, factor) in terms {
                        if *term as usize >= material {
                            return Err(ProbeMaterialError::CombineTermNotEarlier {
                                material,
                                term: *term,
                            });
                        }
                        if !factor.is_finite() {
                            return Err(ProbeMaterialError::NonFinite {
                                material,
                                field: "factor",
                            });
                        }
                    }
                    ProbeMaterial::Combine(
                        terms
                            .iter()
                            .map(|(term, factor)| (*term as usize, *factor))
                            .collect(),
                    )
                }
                ProbeMaterialSpec::Concrete(spec) => {
                    ProbeMaterial::Concrete(ConcreteCurve::new(spec, material)?)
                }
                ProbeMaterialSpec::Steel(spec) => {
                    ProbeMaterial::Steel(SteelCurve::new(spec, material)?)
                }
                ProbeMaterialSpec::General(spec) => ProbeMaterial::General(GeneralCurve::new(spec)),
            });
        }
        Ok(out)
    }
}

/// Stress of `materials[index]` at `strain`.
pub(super) fn stress(materials: &[ProbeMaterial], index: usize, strain: f64) -> f64 {
    match &materials[index] {
        ProbeMaterial::Core(material) => material.trial_stress_tangent(strain).0,
        ProbeMaterial::Combine(terms) => terms
            .iter()
            .map(|(term, factor)| factor * stress(materials, *term, strain))
            .sum(),
        ProbeMaterial::Concrete(curve) => curve.stress(strain),
        ProbeMaterial::Steel(curve) => curve.stress(strain),
        ProbeMaterial::General(curve) => curve.stress(strain),
    }
}

/// Steps in the stress law of each probe material: the height of the largest one, and the strains where there is
/// one. Found by comparing the stress either side of every strain where a law changes character. Smooth and merely
/// piecewise-linear laws give none; a Whitney block, a cracking drop, a strain cutoff or a curve that drops after its
/// last point give the height of the step and where it is. Section evaluation averages fibers over their strain window
/// near such strains, and isoline solves use the height to know how close to equilibrium a stepped section can get.
pub(super) fn material_steps(
    specs: &[ProbeMaterialSpec],
    materials: &[ProbeMaterial],
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut jumps: Vec<f64> = Vec::with_capacity(specs.len());
    let mut points: Vec<Vec<f64>> = Vec::with_capacity(specs.len());
    for (index, spec) in specs.iter().enumerate() {
        let (jump, at) = match spec {
            ProbeMaterialSpec::Core { .. } => (0.0, Vec::new()),
            ProbeMaterialSpec::Combine { terms } => (
                terms
                    .iter()
                    .map(|(term, factor)| factor.abs() * jumps[*term as usize])
                    .sum(),
                terms
                    .iter()
                    .flat_map(|(term, _)| points[*term as usize].iter().copied())
                    .collect(),
            ),
            _ => {
                let mut candidates = vec![0.0];
                let mut add =
                    |values: &[f64]| candidates.extend(values.iter().flat_map(|v| [*v, -*v]));
                match spec {
                    ProbeMaterialSpec::Concrete(c) => {
                        add(&[c.ecu, c.eu.unwrap_or(c.ecu)]);
                        if let Some(beta1) = c.beta1 {
                            add(&[c.ecu * (1.0 - beta1)]);
                        }
                        if c.ec > 0.0 {
                            add(&[c.fc / c.ec, c.ft / c.ec]);
                        }
                        if c.ets > 0.0 {
                            add(&[c.ft / c.ets]);
                        }
                        add(&c.curve_points.iter().map(|p| p[0]).collect::<Vec<_>>());
                    }
                    ProbeMaterialSpec::Steel(st) => {
                        add(&[st.eps_uc, st.eps_ut]);
                        if st.es > 0.0 {
                            add(&[st.fy / st.es]);
                        }
                        add(&st.curve_points.iter().map(|p| p[0]).collect::<Vec<_>>());
                    }
                    ProbeMaterialSpec::General(g) => {
                        add(&[g.eps_uc, g.eps_ut]);
                        add(&g.curve_points.iter().map(|p| p[0]).collect::<Vec<_>>());
                    }
                    _ => {}
                }
                let mut at = Vec::new();
                let mut largest: f64 = 0.0;
                for b in candidates.into_iter().filter(|b| b.is_finite()) {
                    let d = 1e-9 * b.abs().max(1.0);
                    let step =
                        (stress(materials, index, b + d) - stress(materials, index, b - d)).abs();
                    // Elastic slope over the probe width is ~1e-9 of a stress; a real step is far larger.
                    if step
                        > 1e-6
                            * stress(materials, index, b)
                                .abs()
                                .max(stress(materials, index, b + 1e-3).abs())
                                .max(1e-12)
                    {
                        at.push(b);
                        largest = largest.max(step);
                    }
                }
                (largest, at)
            }
        };
        jumps.push(jump);
        points.push(at);
    }
    (jumps, points)
}

/// Samples across a fiber's strain window where a stepped law is averaged.
pub(super) const WINDOW_SAMPLES: usize = 32;

/// Stress of material `index` averaged over the strain range `[strain - half, strain + half]` that a fiber of finite
/// size sees across the section's strain gradient. A fiber is a point sample of a strain field, so a stepped law
/// (a Whitney block, a cracking drop, a strain cutoff) turns the resultants into a staircase in the strain state; the
/// window average integrates the step over the fiber, which makes the resultants continuous (the stress block is
/// integrated exactly in the limit of fine fibers). Only the few fibers whose window straddles a step cost more
/// than a plain evaluation.
pub(super) fn windowed_stress(
    materials: &[ProbeMaterial],
    index: usize,
    steps: &[f64],
    strain: f64,
    half: f64,
) -> f64 {
    // Only a window that contains a step differs from the plain evaluation (a kink averages out to the same value to
    // first order, and the engine's laws are piecewise linear or smooth away from the steps).
    if half <= 1e-12 || !steps.iter().any(|b| (b - strain).abs() < half) {
        return stress(materials, index, strain);
    }
    let width = 2.0 * half / WINDOW_SAMPLES as f64;
    (0..WINDOW_SAMPLES)
        .map(|k| stress(materials, index, strain - half + (k as f64 + 0.5) * width))
        .sum::<f64>()
        / WINDOW_SAMPLES as f64
}
