//! Moment-curvature at constant axial force.
//!
//! For each curvature `kappa` along the direction `theta` the axial strain `eps0` is solved so that the section
//! carries the axial force `P`: strain is `eps0 - kappa * r` with `r = y cos(theta) + z sin(theta)`, so positive
//! `kappa` compresses the large-`r` side. The curvature direction is fixed (the perpendicular moment is whatever
//! the section develops), which is what CalcsApp's moment-curvature reports.
//!
//! Each solve starts from the previous step's `eps0`, expands a bracket in the direction that reduces the
//! residual `N(eps0) - P`, then closes it with Illinois false position. If that fails (a non-monotone `N`, as
//! past crushing), a coarse scan of a wide window picks the sign change closest to the previous `eps0`.
//! Evaluation is path independent, so there is no unloading history to disagree with a reversible law.

use serde::{Deserialize, Serialize};

use super::{SectionProbe, SectionProbeError};

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct MomentCurvatureRequest {
    /// Axial force carried (tension positive).
    pub axial: f64,
    /// Curvature direction (see the module docs).
    pub angle: f64,
    /// Curvatures to solve at, in the order given (reuse of the previous `eps0` assumes they are ordered).
    pub curvatures: Vec<f64>,
    /// Relative force tolerance; default `1e-9`.
    #[serde(default = "default_tolerance")]
    pub tolerance: f64,
}

fn default_tolerance() -> f64 {
    1e-9
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct MomentCurvaturePoint {
    pub curvature: f64,
    pub eps0: f64,
    pub n: f64,
    pub mz: f64,
    pub my: f64,
    pub m2: f64,
    pub m3: f64,
    /// Resultant `hypot(m2, m3)`.
    pub moment: f64,
    pub min_strain: f64,
    pub max_strain: f64,
    pub converged: bool,
    pub iterations: u32,
    /// `N - P` at the returned `eps0`.
    pub residual: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct MomentCurvatureResult {
    pub points: Vec<MomentCurvaturePoint>,
    pub non_converged: u32,
}

impl SectionProbe {
    pub fn moment_curvature(
        &self,
        request: &MomentCurvatureRequest,
    ) -> Result<MomentCurvatureResult, SectionProbeError> {
        if !request.axial.is_finite() || !request.angle.is_finite() || !(request.tolerance > 0.0) {
            return Err(SectionProbeError::InvalidRequest {
                reason: "axial, angle and a positive tolerance are required".into(),
            });
        }
        let (lo, hi) = self.extents(request.angle);
        let height = hi - lo;
        if !(height > 0.0) {
            return Err(SectionProbeError::DegenerateSection);
        }
        let (s, c) = request.angle.sin_cos();
        let mut points = Vec::with_capacity(request.curvatures.len());
        let mut non_converged = 0;
        let mut previous = 0.0;
        for &kappa in &request.curvatures {
            let (kz, ky) = (kappa * c, -kappa * s);
            let force_scale = |e0: f64| {
                let mut sum = request.axial.abs();
                for i in 0..self.y.len() {
                    sum += self.fiber_force(i, e0, kz, ky).abs();
                }
                sum.max(1.0)
            };
            let residual = |e0: f64| self.evaluate(e0, kz, ky).n - request.axial;
            let solved = solve_axial(
                &residual,
                &force_scale,
                previous,
                kappa.abs() * height,
                request.tolerance,
            );
            previous = solved.eps0;
            let r = self.evaluate(solved.eps0, kz, ky);
            let (a, b) = (solved.eps0 - kappa * hi, solved.eps0 - kappa * lo);
            if !solved.converged {
                non_converged += 1;
            }
            points.push(MomentCurvaturePoint {
                curvature: kappa,
                eps0: solved.eps0,
                n: r.n,
                mz: r.mz,
                my: r.my,
                m2: r.m2,
                m3: r.m3,
                moment: r.m2.hypot(r.m3),
                min_strain: a.min(b),
                max_strain: a.max(b),
                converged: solved.converged,
                iterations: solved.iterations,
                residual: r.n - request.axial,
            });
        }
        Ok(MomentCurvatureResult {
            points,
            non_converged,
        })
    }
}

struct Solved {
    eps0: f64,
    converged: bool,
    iterations: u32,
}

fn solve_axial(
    residual: &dyn Fn(f64) -> f64,
    force_scale: &dyn Fn(f64) -> f64,
    start: f64,
    strain_span: f64,
    tolerance: f64,
) -> Solved {
    let f_start = residual(start);
    let tol = |e: f64| tolerance * force_scale(e);
    if f_start.abs() <= tol(start) {
        return Solved {
            eps0: start,
            converged: true,
            iterations: 0,
        };
    }
    // Expand a bracket from `start`; N increases with eps0 for a stiffening-positive section, so step against the sign.
    let base = (0.1 * strain_span).max(1e-6);
    for direction in [
        if f_start > 0.0 { -1.0 } else { 1.0 },
        if f_start > 0.0 { 1.0 } else { -1.0 },
    ] {
        let (mut a, mut fa) = (start, f_start);
        let mut step = base;
        for _ in 0..40 {
            let b = a + direction * step;
            let fb = residual(b);
            if fa * fb <= 0.0 {
                return illinois(residual, force_scale, a, fa, b, fb, tolerance);
            }
            a = b;
            fa = fb;
            step *= 2.0;
            if a.abs() > 10.0 {
                break;
            }
        }
    }
    // Fallback scan of a wide window for the sign change closest to `start`.
    let window = (5.0 * strain_span).max(0.05);
    let samples = 400;
    let mut best: Option<(f64, f64, f64, f64, f64)> = None; // distance, a, fa, b, fb
    let mut prev = (start - window, residual(start - window));
    for i in 1..=samples {
        let e = start - window + 2.0 * window * i as f64 / samples as f64;
        let f = residual(e);
        if prev.1 * f <= 0.0 {
            let distance = ((prev.0 + e) * 0.5 - start).abs();
            if best.map_or(true, |b| distance < b.0) {
                best = Some((distance, prev.0, prev.1, e, f));
            }
        }
        prev = (e, f);
    }
    match best {
        Some((_, a, fa, b, fb)) => illinois(residual, force_scale, a, fa, b, fb, tolerance),
        None => Solved {
            eps0: start,
            converged: false,
            iterations: 0,
        },
    }
}

fn illinois(
    residual: &dyn Fn(f64) -> f64,
    force_scale: &dyn Fn(f64) -> f64,
    mut a: f64,
    mut fa: f64,
    mut b: f64,
    mut fb: f64,
    tolerance: f64,
) -> Solved {
    // Invariant: fa and fb have opposite signs (or one is zero).
    for iteration in 1..=100u32 {
        let c = if fb != fa {
            (a * fb - b * fa) / (fb - fa)
        } else {
            0.5 * (a + b)
        };
        let fc = residual(c);
        if fc.abs() <= tolerance * force_scale(c) || (b - a).abs() < 1e-15 {
            return Solved {
                eps0: c,
                converged: true,
                iterations: iteration,
            };
        }
        if fc * fb > 0.0 {
            fa *= 0.5; // Illinois: halve the retained end so the bracket keeps closing
        } else {
            a = b;
            fa = fb;
        }
        b = c;
        fb = fc;
    }
    let mid = 0.5 * (a + b);
    Solved {
        eps0: mid,
        converged: residual(mid).abs() <= 1e3 * tolerance * force_scale(mid),
        iterations: 100,
    }
}
