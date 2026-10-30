//! Fiber-section probe.
//!
//! Evaluates a fiber section at prescribed strain states, without a `Domain`, nodes or a load history, and
//! solves the equilibrium problems built on that: PMM state sweeps, direct fixed-axial-load /
//! fixed-moment-angle isolines, and moment-curvature at constant axial load. It is the section analogue of
//! `material_probe`, for sectional-analysis previews and for CalcsApp's PMM/moment-curvature engine.
//!
//! Section conventions are Carapace's 3D `FiberSection3` (`core::model::fiber_section`):
//!
//! - fiber `(y, z)` coordinates; strain `eps = eps0 - y * kz + z * ky`;
//! - `n = sum(f)`, `mz = -sum(y f)`, `my = sum(z f)`, with `f = stress * area`;
//! - tension positive.
//!
//! "Engineering" moments `m2 = sum(y f) = -mz` and `m3 = sum(z f) = my` are also reported; they are CalcsApp's
//! `M2`/`M3` when its `x` axis is the probe's `z`.
//!
//! Angles: `theta` fixes a direction `r = y cos(theta) + z sin(theta)` across the section. A strain profile
//! `(et, eb)` at angle `theta` has strain `et` at the largest `r` and `eb` at the smallest, linear between.
//!
//! Materials are stateless (see [`materials`]); nothing here commits or remembers history, so every
//! evaluation is path independent. Probe-only materials live in [`materials`] and are not part of Carapace's
//! `MaterialSpec`.

mod factoring;
mod isolines;
mod materials;
mod moment_curvature;

use carapace_core::model::Material;
use serde::{Deserialize, Serialize};

use crate::input_v1::materials::{resolve_materials, MaterialSpec};

pub use factoring::{FactoredResponse, Factoring};
pub use isolines::{IsolineOptions, IsolinePoint, IsolineRequest, IsolineResult, STATE_STRIDE};
pub use materials::{
    ConcreteModel, ConcreteSpec, GeneralModel, GeneralSpec, ProbeMaterialError, ProbeMaterialSpec,
    SteelModel, SteelSpec,
};
pub use moment_curvature::{MomentCurvaturePoint, MomentCurvatureRequest, MomentCurvatureResult};

use materials::{material_jumps, stress, windowed_stress, ProbeMaterial, WINDOW_SAMPLES};

/// Fibers as parallel arrays (`y`, `z`, `area`, `material`, optional `group`).
#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ProbeFibers {
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub area: Vec<f64>,
    /// Index into `SectionProbeConfig::materials`.
    pub material: Vec<u32>,
    /// Group id per fiber (for per-group force sums, e.g. concrete vs steel resistance factors). Defaults to
    /// group 0 for every fiber.
    #[serde(default)]
    pub group: Vec<u32>,
}

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct SectionProbeConfig {
    /// Carapace materials that `core` entries of `materials` refer to by index (may be empty).
    #[serde(default)]
    pub core_materials: Vec<MaterialSpec>,
    pub materials: Vec<ProbeMaterialSpec>,
    pub fibers: ProbeFibers,
    /// Average fibers over the strain window they span where a stress law has a step (default `true`; see
    /// `materials::windowed_stress`). `false` evaluates every fiber at its centroid strain, as CalcsApp's JS engine did.
    #[serde(default = "default_true")]
    pub step_smoothing: bool,
}

fn default_true() -> bool {
    true
}

#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum SectionProbeError {
    InvalidFibers {
        reason: String,
    },
    InvalidMaterial {
        reason: String,
    },
    /// The strain profile has no extent across the chosen angle (all fibers at the same `r`).
    DegenerateSection,
    InvalidRequest {
        reason: String,
    },
    NoGrid,
}

/// Section resultants at one strain state.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct SectionResponse {
    pub n: f64,
    pub mz: f64,
    pub my: f64,
    /// `-mz`
    pub m2: f64,
    /// `my`
    pub m3: f64,
}

impl SectionResponse {
    pub(super) fn new(n: f64, mz: f64, my: f64) -> Self {
        SectionResponse {
            n,
            mz,
            my,
            m2: -mz,
            m3: my,
        }
    }
}

/// Maximum compressive and tensile axial force over uniform strain (zero curvature), and the strain where each
/// occurs. The compressive one is the concentric capacity `P0` and is generally interior to the crushing-limited
/// PMM boundary; its moment is non-zero when the section is not symmetric about its centroid.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct AxialExtremes {
    pub compression: AxialExtreme,
    pub tension: AxialExtreme,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct AxialExtreme {
    pub strain: f64,
    pub response: SectionResponse,
}

pub struct SectionProbe {
    pub(super) y: Vec<f64>,
    pub(super) z: Vec<f64>,
    area: Vec<f64>,
    material: Vec<u32>,
    pub(super) group: Vec<u32>,
    group_count: usize,
    materials: Vec<ProbeMaterial>,
    /// Largest force step one fiber can take from a jump in its stress law, and the section's largest lever arm
    /// (see `isolines`).
    quantum: f64,
    radius: f64,
    /// Materials with a step in their law: their fibers are averaged over the strain window the fiber spans.
    windowed: Vec<bool>,
    /// Half the size of each fiber, `sqrt(area) / 2`.
    half_size: Vec<f64>,
    grid: Option<isolines::Grid>,
}

fn invalid_fibers(reason: impl Into<String>) -> SectionProbeError {
    SectionProbeError::InvalidFibers {
        reason: reason.into(),
    }
}

impl SectionProbe {
    pub fn new(config: &SectionProbeConfig) -> Result<Self, SectionProbeError> {
        let f = &config.fibers;
        let n = f.y.len();
        if n == 0 {
            return Err(invalid_fibers("a section needs at least one fiber"));
        }
        if f.z.len() != n || f.area.len() != n || f.material.len() != n {
            return Err(invalid_fibers(
                "y, z, area and material must have one entry per fiber",
            ));
        }
        if !f.group.is_empty() && f.group.len() != n {
            return Err(invalid_fibers("group must be empty or one entry per fiber"));
        }
        if f.y
            .iter()
            .chain(&f.z)
            .chain(&f.area)
            .any(|v| !v.is_finite())
        {
            return Err(invalid_fibers("fiber coordinates and areas must be finite"));
        }
        let core: Vec<Material> = if config.core_materials.is_empty() {
            Vec::new()
        } else {
            resolve_materials(&config.core_materials).map_err(|error| {
                SectionProbeError::InvalidMaterial {
                    reason: format!("{error:?}"),
                }
            })?
        };
        let materials = ProbeMaterial::build(&config.materials, &core).map_err(|error| {
            SectionProbeError::InvalidMaterial {
                reason: format!("{error:?}"),
            }
        })?;
        if let Some(bad) = f.material.iter().find(|m| **m as usize >= materials.len()) {
            return Err(invalid_fibers(format!(
                "fiber material {bad} is not defined"
            )));
        }
        let group = if f.group.is_empty() {
            vec![0; n]
        } else {
            f.group.clone()
        };
        let group_count = group.iter().copied().max().unwrap_or(0) as usize + 1;
        let jumps = material_jumps(&config.materials, &materials);
        let windowed: Vec<bool> = jumps
            .iter()
            .map(|j| config.step_smoothing && *j > 0.0)
            .collect();
        // Averaged over a fiber's window a step is spread over `WINDOW_SAMPLES` smaller ones.
        let spread = if config.step_smoothing {
            WINDOW_SAMPLES as f64
        } else {
            1.0
        };
        let quantum = (0..n)
            .map(|i| f.area[i].abs() * jumps[f.material[i] as usize] / spread)
            .fold(0.0, f64::max);
        let half_size: Vec<f64> = f.area.iter().map(|a| 0.5 * a.abs().sqrt()).collect();
        let radius = (0..n).map(|i| f.y[i].hypot(f.z[i])).fold(0.0, f64::max);
        Ok(SectionProbe {
            y: f.y.clone(),
            z: f.z.clone(),
            area: f.area.clone(),
            material: f.material.clone(),
            group,
            group_count,
            materials,
            quantum,
            radius,
            windowed,
            half_size,
            grid: None,
        })
    }

    pub fn fiber_count(&self) -> usize {
        self.y.len()
    }

    pub fn group_count(&self) -> usize {
        self.group_count
    }

    #[inline]
    pub(super) fn fiber_force(&self, i: usize, eps0: f64, kz: f64, ky: f64) -> f64 {
        let strain = eps0 - self.y[i] * kz + self.z[i] * ky;
        let material = self.material[i] as usize;
        let stress = if self.windowed[material] {
            windowed_stress(
                &self.materials,
                material,
                strain,
                self.half_size[i] * kz.hypot(ky),
            )
        } else {
            stress(&self.materials, material, strain)
        };
        stress * self.area[i]
    }

    /// Resultants at `(eps0, kz, ky)`.
    pub fn evaluate(&self, eps0: f64, kz: f64, ky: f64) -> SectionResponse {
        let (mut n, mut mz, mut my) = (0.0, 0.0, 0.0);
        for i in 0..self.y.len() {
            let force = self.fiber_force(i, eps0, kz, ky);
            n += force;
            mz -= self.y[i] * force;
            my += self.z[i] * force;
        }
        SectionResponse::new(n, mz, my)
    }

    /// Total and per-group resultants: `[n, mz, my]` for the total followed by `[n, mz, my]` for each group.
    pub fn evaluate_grouped(&self, eps0: f64, kz: f64, ky: f64, out: &mut [f64]) {
        out.fill(0.0);
        for i in 0..self.y.len() {
            let force = self.fiber_force(i, eps0, kz, ky);
            let (dn, dmz, dmy) = (force, -self.y[i] * force, self.z[i] * force);
            let g = 3 + 3 * self.group[i] as usize;
            out[0] += dn;
            out[1] += dmz;
            out[2] += dmy;
            out[g] += dn;
            out[g + 1] += dmz;
            out[g + 2] += dmy;
        }
    }

    /// Stress of probe material `index` at each of `strains` (a material preview, and the pure steel stress of a
    /// bar fiber whose own material is `steel - host`). Out-of-range indices give NaN.
    pub fn material_stress(&self, index: usize, strains: &[f64]) -> Vec<f64> {
        if index >= self.materials.len() {
            return vec![f64::NAN; strains.len()];
        }
        strains
            .iter()
            .map(|e| stress(&self.materials, index, *e))
            .collect()
    }

    /// `(strain, stress)` of every fiber at `(eps0, kz, ky)`, flattened: strain, stress, strain, stress, ...
    pub fn fiber_responses(&self, eps0: f64, kz: f64, ky: f64) -> Vec<f64> {
        let mut out = Vec::with_capacity(2 * self.y.len());
        for i in 0..self.y.len() {
            let strain = eps0 - self.y[i] * kz + self.z[i] * ky;
            out.push(strain);
            out.push(stress(&self.materials, self.material[i] as usize, strain));
        }
        out
    }

    /// Resultants for many states. `states` holds `(eps0, kz, ky)` triples. The output stride is 3
    /// (`n, mz, my`), or `3 * (1 + group_count)` with `with_groups` (total, then each group).
    pub fn evaluate_batch(&self, states: &[f64], with_groups: bool) -> Vec<f64> {
        let stride = if with_groups {
            3 * (1 + self.group_count)
        } else {
            3
        };
        let mut out = vec![0.0; states.len() / 3 * stride];
        for (k, s) in states.chunks_exact(3).enumerate() {
            let slot = &mut out[k * stride..(k + 1) * stride];
            if with_groups {
                self.evaluate_grouped(s[0], s[1], s[2], slot);
            } else {
                let r = self.evaluate(s[0], s[1], s[2]);
                slot.copy_from_slice(&[r.n, r.mz, r.my]);
            }
        }
        out
    }

    /// `(min, max)` of `r = y cos(theta) + z sin(theta)` over the fibers.
    pub fn extents(&self, theta: f64) -> (f64, f64) {
        let (s, c) = theta.sin_cos();
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in 0..self.y.len() {
            let r = self.y[i] * c + self.z[i] * s;
            lo = lo.min(r);
            hi = hi.max(r);
        }
        (lo, hi)
    }

    /// `(eps0, kz, ky)` of the strain profile with strain `et` at the largest `r` and `eb` at the smallest.
    pub fn profile_state(
        &self,
        theta: f64,
        et: f64,
        eb: f64,
    ) -> Result<(f64, f64, f64), SectionProbeError> {
        let (lo, hi) = self.extents(theta);
        let h = hi - lo;
        if !(h > 0.0) {
            return Err(SectionProbeError::DegenerateSection);
        }
        let g = (et - eb) / h;
        let (s, c) = theta.sin_cos();
        Ok((eb - g * lo, -g * c, g * s))
    }

    /// Strain states `(eps0, kz, ky)` for every `(angle, profile)` pair, angles outermost. `profiles` holds
    /// `(et, eb)` pairs.
    pub fn profile_states(
        &self,
        angles: &[f64],
        profiles: &[f64],
    ) -> Result<Vec<f64>, SectionProbeError> {
        let mut out = Vec::with_capacity(angles.len() * profiles.len() / 2 * 3);
        for &theta in angles {
            for p in profiles.chunks_exact(2) {
                let (e0, kz, ky) = self.profile_state(theta, p[0], p[1])?;
                out.extend_from_slice(&[e0, kz, ky]);
            }
        }
        Ok(out)
    }

    /// Axial-force extremes over uniform strain in `[-max_strain, max_strain]`: a scan, then golden-section
    /// refinement around the best sample of each.
    pub fn axial_extremes(&self, max_strain: f64) -> Result<AxialExtremes, SectionProbeError> {
        if !(max_strain.is_finite() && max_strain > 0.0) {
            return Err(SectionProbeError::InvalidRequest {
                reason: "maxStrain must be positive".into(),
            });
        }
        const SAMPLES: usize = 800;
        let n_at = |e: f64| self.evaluate(e, 0.0, 0.0).n;
        let step = 2.0 * max_strain / SAMPLES as f64;
        let samples: Vec<f64> = (0..=SAMPLES)
            .map(|i| n_at(-max_strain + i as f64 * step))
            .collect();
        let refine = |sign: f64| {
            let mut best = 0usize;
            for (i, n) in samples.iter().enumerate() {
                if sign * n < sign * samples[best] {
                    best = i;
                }
            }
            let centre = -max_strain + best as f64 * step;
            let (mut a, mut b) = (
                (centre - step).max(-max_strain),
                (centre + step).min(max_strain),
            );
            let phi = 0.618_033_988_749_894_9;
            for _ in 0..60 {
                let (c, d) = (b - phi * (b - a), a + phi * (b - a));
                if sign * n_at(c) < sign * n_at(d) {
                    b = d;
                } else {
                    a = c;
                }
            }
            let strain = 0.5 * (a + b);
            AxialExtreme {
                strain,
                response: self.evaluate(strain, 0.0, 0.0),
            }
        };
        Ok(AxialExtremes {
            compression: refine(1.0),
            tension: refine(-1.0),
        })
    }
}
