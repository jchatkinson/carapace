//! `Material::CsfmRebar`: the CSFM embedded-rebar law (plan
//! `docs/csfm-plan.md` §1 point 3, §2.2), ported from CalcsApp's
//! `resolveModel`/`makeRodLaw`/`makeBareRodLaw`
//! (`../calcs.app/lib/fem/csfm-nonlinear.ts`,
//! `../calcs.app/lib/fem/nonlinear-material.ts`).
//!
//! The law is stateless and nonlinear-elastic: every resolved quantity
//! (routing, crack spacing, branch boundaries) is computed once at
//! construction, and `evaluate` is then a total, branch-dispatching
//! function of the trial strain alone — unloading retraces the loading
//! curve exactly, which is monotonic-only by construction and not an
//! approximation of hysteresis.

use super::{LimitKind, LimitStateConfig, Material, MaterialError};

/// Which average-strain model a path's effective reinforcement ratio
/// routes to (plan §1): the Tension Chord Model above the critical ratio
/// `rho_crit`, the Pull-Out Model below it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CsfmRebarRouting {
    Tcm,
    Pom,
}

/// The full field set behind `Material::CsfmRebar` (see that variant's doc
/// comment for why it's boxed). `Copy` despite the field count — every
/// field is `f64`/`bool`/`CsfmRebarRouting`/`Option<f64>`.
///
/// `tau_b0`/`tau_b1` are the bond shear stresses (`2 * fctm`, `fctm`); `sr`
/// is the crack spacing (`lambda * d * fctm * (1 - rho_eff) / (2 * tau_b0 *
/// rho_eff)` for TCM, `0.0` for POM — unused by POM's forward relation);
/// `sigma_min` is the TCM crack stress `tau_b0 * sr / d` (`0.0` for POM).
/// `eps_m_at_fy`/`sigma_at_branch2_upper`/`eps_m_at_branch2_upper`/
/// `eps_m_at_ft` are the TCM/POM forward relation's branch boundaries,
/// resolved once via the forward relation itself (`forward_strain`) —
/// exactly as CalcsApp's `resolveModel` computes `branchBoundaries`.
/// `sigma_at_branch2_upper`/`eps_m_at_branch2_upper` are `None` when the
/// nominal TCM branch-2 upper stress `fy + 2 * tau_b1 * sr / d` is at or
/// past `ft` (branch 2 then runs all the way to rupture, with no branch 3).
#[derive(Debug, Clone, Copy)]
pub struct CsfmRebarFields {
    pub es: f64,
    pub fy: f64,
    pub ft: f64,
    pub eps_u: f64,
    pub esh: f64,
    pub tension_stiffening: bool,
    pub routing: CsfmRebarRouting,
    pub diameter_mm: f64,
    pub tau_b0: f64,
    pub tau_b1: f64,
    pub sr: f64,
    pub sigma_min: f64,
    pub eps_m_at_fy: f64,
    pub sigma_at_branch2_upper: Option<f64>,
    pub eps_m_at_branch2_upper: Option<f64>,
    pub eps_m_at_ft: f64,
}

/// POM has a square-root inverse at the origin; strains at or below this
/// are returned as the exact zero-stress, `Es`-tangent origin point rather
/// than evaluating the (numerically fine, but needlessly evaluated) square
/// root at a vanishing argument. Matches CalcsApp's `ORIGIN_STRAIN`.
const ORIGIN_STRAIN: f64 = 1e-9;

impl Material {
    /// Builds a `CsfmRebar`. `es`/`fy`/`ft`/`eps_u` are the bare-steel
    /// bilinear properties; `fctm`/`ecm` are the host concrete's tensile
    /// strength and modulus (the rebar law takes them explicitly rather
    /// than referencing a `PlaneMaterial`, keeping `Material`
    /// self-contained — plan §6 decision 2); `diameter_mm`/`rho_eff` are
    /// the bar diameter and effective reinforcement ratio; `lambda` is the
    /// crack-spacing calibration factor; `tension_stiffening` selects
    /// between the TCM/POM average-strain relation in tension (`true`) and
    /// symmetric bare bilinear steel (`false`) — compression is always
    /// bare bilinear steel either way.
    ///
    /// # Errors
    ///
    /// [`MaterialError::InvalidSteelProperties`] if `es`, `fy`, `ft` or
    /// `eps_u` is not finite and positive, `ft <= fy`, or `eps_u <= fy /
    /// es`; [`MaterialError::InvalidReinforcementRatio`] if `rho_eff` is
    /// not in `(0, 1)`; [`MaterialError::InvalidConcreteInputs`] if
    /// `fctm`, `ecm`, `diameter_mm` or `lambda` is not finite and
    /// positive; [`MaterialError::TcmStartupExceedsYield`] if TCM routing
    /// gives a startup bridge whose crack stress `2 * sigma_min` exceeds
    /// `fy` (plan §1, §6 — CalcsApp rejects this rather than reproducing
    /// the discontinuous TS bridge it used to have).
    #[allow(clippy::too_many_arguments)]
    pub fn csfm_rebar(
        es: f64,
        fy: f64,
        ft: f64,
        eps_u: f64,
        fctm: f64,
        ecm: f64,
        diameter_mm: f64,
        rho_eff: f64,
        lambda: f64,
        tension_stiffening: bool,
    ) -> Result<Material, MaterialError> {
        let positive_finite = |v: f64| v.is_finite() && v > 0.0;
        if !positive_finite(es) || !positive_finite(fy) || !positive_finite(ft) || ft <= fy {
            return Err(MaterialError::InvalidSteelProperties);
        }
        let ey = fy / es;
        if !positive_finite(eps_u) || eps_u <= ey {
            return Err(MaterialError::InvalidSteelProperties);
        }
        if !(rho_eff.is_finite() && rho_eff > 0.0 && rho_eff < 1.0) {
            return Err(MaterialError::InvalidReinforcementRatio);
        }
        if !positive_finite(fctm)
            || !positive_finite(ecm)
            || !positive_finite(diameter_mm)
            || !positive_finite(lambda)
        {
            return Err(MaterialError::InvalidConcreteInputs);
        }

        let esh = (ft - fy) / (eps_u - ey);
        let tau_b0 = 2.0 * fctm;
        let tau_b1 = fctm;
        // rho_crit = fctm / (fy - (Es/Ecm - 1) * fctm) (plan §1, derived in
        // core/tests/csfm/reference/tension_chord.md from first-crack
        // equilibrium — not a CSFM-specific formula).
        let rho_crit = fctm / (fy - (es / ecm - 1.0) * fctm);
        let routing = if rho_eff >= rho_crit {
            CsfmRebarRouting::Tcm
        } else {
            CsfmRebarRouting::Pom
        };
        let sr = match routing {
            CsfmRebarRouting::Tcm => {
                lambda * diameter_mm * fctm * (1.0 - rho_eff) / (2.0 * tau_b0 * rho_eff)
            }
            CsfmRebarRouting::Pom => 0.0,
        };
        let sigma_min = match routing {
            CsfmRebarRouting::Tcm => tau_b0 * sr / diameter_mm,
            CsfmRebarRouting::Pom => 0.0,
        };
        if 2.0 * sigma_min > fy {
            return Err(MaterialError::TcmStartupExceedsYield);
        }

        // The forward relation only reads the fields set below (not the
        // boundaries being computed from it), so resolve them against a
        // base value with placeholder boundaries first.
        let base = CsfmRebarFields {
            es,
            fy,
            ft,
            eps_u,
            esh,
            tension_stiffening,
            routing,
            diameter_mm,
            tau_b0,
            tau_b1,
            sr,
            sigma_min,
            eps_m_at_fy: 0.0,
            sigma_at_branch2_upper: None,
            eps_m_at_branch2_upper: None,
            eps_m_at_ft: 0.0,
        };
        let eps_m_at_fy = forward_strain(fy, &base);
        let nominal_upper = match routing {
            CsfmRebarRouting::Tcm => Some(fy + 2.0 * tau_b1 * sr / diameter_mm),
            CsfmRebarRouting::Pom => None,
        };
        let sigma_at_branch2_upper = nominal_upper.filter(|&upper| upper < ft);
        let eps_m_at_branch2_upper =
            sigma_at_branch2_upper.map(|sigma| forward_strain(sigma, &base));
        let eps_m_at_ft = forward_strain(ft, &base);

        Ok(Material::CsfmRebar(Box::new(CsfmRebarFields {
            eps_m_at_fy,
            sigma_at_branch2_upper,
            eps_m_at_branch2_upper,
            eps_m_at_ft,
            ..base
        })))
    }
}

/// The TCM/POM forward average-strain relation `eps_m(sigma)`
/// (`forwardStrain` in `csfm-nonlinear.ts`), used only at construction to
/// resolve the branch boundaries. Requires `0 <= sigma <= ft`.
fn forward_strain(sigma: f64, f: &CsfmRebarFields) -> f64 {
    match f.routing {
        CsfmRebarRouting::Pom => {
            let ratio = f.tau_b1 / f.tau_b0;
            let denominator = f.ft + f.fy * (ratio - 1.0);
            if sigma <= f.fy {
                sigma * sigma * ratio / (2.0 * f.es * denominator)
            } else {
                ((f.fy / f.es) * (sigma + f.fy * (ratio / 2.0 - 1.0))
                    + (sigma - f.fy).powi(2) / (2.0 * f.esh))
                    / denominator
            }
        }
        CsfmRebarRouting::Tcm => {
            let offset = f.tau_b0 * f.sr / (f.es * f.diameter_mm);
            let upper = f.fy + 2.0 * f.tau_b1 * f.sr / f.diameter_mm;
            let ey = f.fy / f.es;
            if sigma <= f.fy {
                (sigma / f.es - offset).max(0.0)
            } else if sigma <= upper {
                (sigma - f.fy).powi(2) * f.diameter_mm / (4.0 * f.esh * f.tau_b1 * f.sr)
                    * (1.0 - f.esh * f.tau_b0 / (f.es * f.tau_b1))
                    + (sigma - f.fy) / f.es * (f.tau_b0 / f.tau_b1)
                    + ey
                    - offset
            } else {
                ey + (sigma - f.fy) / f.esh - f.tau_b1 * f.sr / (f.esh * f.diameter_mm)
            }
        }
    }
}

/// `d(eps_m)/d(sigma)` along the same forward relation, branch for branch
/// (`forwardStrainSlope`). The rod tangent is its reciprocal.
fn forward_strain_slope(sigma: f64, f: &CsfmRebarFields) -> f64 {
    match f.routing {
        CsfmRebarRouting::Pom => {
            let ratio = f.tau_b1 / f.tau_b0;
            let denominator = f.ft + f.fy * (ratio - 1.0);
            if sigma <= f.fy {
                sigma * ratio / (f.es * denominator)
            } else {
                (f.fy / f.es + (sigma - f.fy) / f.esh) / denominator
            }
        }
        CsfmRebarRouting::Tcm => {
            let upper = f.fy + 2.0 * f.tau_b1 * f.sr / f.diameter_mm;
            if sigma <= f.fy {
                1.0 / f.es
            } else if sigma <= upper {
                2.0 * (sigma - f.fy) * f.diameter_mm / (4.0 * f.esh * f.tau_b1 * f.sr)
                    * (1.0 - f.esh * f.tau_b0 / (f.es * f.tau_b1))
                    + (f.tau_b0 / f.tau_b1) / f.es
            } else {
                1.0 / f.esh
            }
        }
    }
}

/// Cancellation-free quadratic solve for `a*x^2 + b*x + c = 0` (or, if `a`
/// is exactly zero, the affine fallback `-c/b`): `q = -(b + sign(b) *
/// sqrt(disc)) / 2`, roots `q/a` and `c/q` (plan §2.2 — avoids the
/// catastrophic cancellation the textbook `(-b +- sqrt(disc))/(2a)` formula
/// suffers when `b` and `sqrt(disc)` are close in magnitude and sign).
/// Picks whichever root lies in `[min, max]` (with a small tolerance for
/// floating-point roundoff at the branch's own boundary); a
/// `debug_assert!` guards that one does, since construction is required to
/// make every branch monotone and well posed (plan §2.2) — a release build
/// instead clamps to the branch to avoid propagating a `NaN`/out-of-range
/// result into an otherwise-total function.
fn solve_quadratic_branch(a: f64, b: f64, c: f64, min: f64, max: f64) -> f64 {
    let tolerance = 1e-9 * (1.0 + min.abs().max(max.abs()));
    if a.abs() < 1e-15 {
        let root = -c / b;
        debug_assert!(
            root >= min - tolerance && root <= max + tolerance,
            "linear inverse root {root} outside [{min}, {max}] (a={a}, b={b}, c={c})"
        );
        return root.clamp(min, max);
    }
    let discriminant = (b * b - 4.0 * a * c).max(0.0);
    let q = -(b + b.signum() * discriminant.sqrt()) / 2.0;
    let candidates = if q.abs() > f64::MIN_POSITIVE {
        [q / a, c / q]
    } else {
        // b and sqrt(disc) cancelled to (near) zero: the quadratic's own
        // vertex is at x = -b/(2a), a double root.
        [-b / (2.0 * a), -b / (2.0 * a)]
    };
    for root in candidates {
        if root >= min - tolerance && root <= max + tolerance {
            return root.clamp(min, max);
        }
    }
    debug_assert!(
        false,
        "no quadratic inverse root in [{min}, {max}]: a={a}, b={b}, c={c}, roots={candidates:?}"
    );
    candidates
        .into_iter()
        .map(|root| root.clamp(min, max))
        .min_by(|x, y| {
            let dx = (x - candidates[0]).abs();
            let dy = (y - candidates[0]).abs();
            dx.partial_cmp(&dy).unwrap()
        })
        .unwrap_or(min)
}

/// The TCM/POM inverse `sigma(eps_m)` (`inverseStress`), branch for
/// branch, on the stabilized-cracking relation (not the startup bridge —
/// see [`tension_stiffened_stress`]). Requires `strain >= 0` and `strain
/// <= eps_m_at_ft` (checked by the caller).
fn inverse_stress(strain: f64, f: &CsfmRebarFields) -> f64 {
    match f.routing {
        CsfmRebarRouting::Pom => {
            let ratio = f.tau_b1 / f.tau_b0;
            let denominator = f.ft + f.fy * (ratio - 1.0);
            if strain <= f.eps_m_at_fy {
                (strain.max(0.0) * 2.0 * f.es * denominator / ratio).sqrt()
            } else {
                // Solved directly for sigma in [fy, ft] (unlike the TCM
                // branch below, this is not shifted to x = sigma - fy).
                let a = 1.0 / (2.0 * f.esh);
                let b = f.fy / f.es - f.fy / f.esh;
                let c = f.fy * f.fy / (2.0 * f.esh) + (f.fy / f.es) * f.fy * (ratio / 2.0 - 1.0)
                    - strain * denominator;
                solve_quadratic_branch(a, b, c, f.fy, f.ft)
            }
        }
        CsfmRebarRouting::Tcm => {
            let offset = f.tau_b0 * f.sr / (f.es * f.diameter_mm);
            if strain <= f.eps_m_at_fy {
                f.es * strain + f.sigma_min
            } else {
                let branch2 = |upper_minus_fy: f64| {
                    let a = f.diameter_mm / (4.0 * f.esh * f.tau_b1 * f.sr)
                        * (1.0 - f.esh * f.tau_b0 / (f.es * f.tau_b1));
                    let b = (f.tau_b0 / f.tau_b1) / f.es;
                    let c = f.fy / f.es - offset - strain;
                    f.fy + solve_quadratic_branch(a, b, c, 0.0, upper_minus_fy)
                };
                match (f.sigma_at_branch2_upper, f.eps_m_at_branch2_upper) {
                    (None, _) => branch2(f.ft - f.fy),
                    (Some(upper), Some(eps_upper)) if strain <= eps_upper => branch2(upper - f.fy),
                    _ => {
                        f.fy + f.esh
                            * (strain - f.fy / f.es + f.tau_b1 * f.sr / (f.esh * f.diameter_mm))
                    }
                }
            }
        }
    }
}

/// Bridges the TCM forward relation's non-unique origin interval (every
/// stress in `[0, sigma_min]` maps to average strain `0`) with a cubic
/// that is `C1` at `eps_startup = sigma_min / Es`, joining `Es * eps +
/// sigma_min` (the inverse's first branch) with matching value and slope
/// (`tensionStiffenedStress`). A finite-element Newton iteration needs a
/// continuous, differentiable response even at an infinitesimal trial
/// strain, which the bare interval is not.
fn tension_stiffened_stress(strain: f64, f: &CsfmRebarFields) -> f64 {
    let inverse = inverse_stress(strain, f);
    if f.routing != CsfmRebarRouting::Tcm {
        return inverse;
    }
    let eps_startup = f.sigma_min / f.es;
    if eps_startup <= 0.0 || strain >= eps_startup {
        return inverse;
    }
    let t = strain / eps_startup;
    let blend = t * t * (3.0 - 2.0 * t);
    f.es * strain + f.sigma_min * blend
}

/// Clamps a reciprocal tangent to `[Es / 1e4, 10 * Es]` (`boundedRodTangent`)
/// — a numerical-conditioning limit on the assembled global diagonal, not a
/// material property. The POM forward relation is quadratic through the
/// origin, so its slope vanishes there and the raw reciprocal diverges.
fn bounded_tangent(tangent: f64, f: &CsfmRebarFields) -> f64 {
    let maximum = 10.0 * f.es;
    let minimum = f.es / 1e4;
    if !tangent.is_finite() || tangent > maximum {
        maximum
    } else {
        tangent.max(minimum)
    }
}

/// Bilinear bare steel, symmetric in sign (`evaluateSteel`): elastic below
/// yield, linear hardening (slope `Esh`) above, rupture past `eps_u` with
/// stress `sign(strain) * ft` and tangent `ft / |strain|` (a small
/// positive secant, per plan §1 — strain stays finite and nonzero past
/// `eps_u` since `eps_u > 0`).
fn bare_steel(strain: f64, f: &CsfmRebarFields) -> (f64, f64) {
    let ey = f.fy / f.es;
    if strain.abs() > f.eps_u {
        let sign = if strain >= 0.0 { 1.0 } else { -1.0 };
        return (sign * f.ft, f.ft / strain.abs().max(1e-9));
    }
    if strain.abs() <= ey {
        (f.es * strain, f.es)
    } else {
        (
            strain.signum() * (f.fy + f.esh * (strain.abs() - ey)),
            f.esh,
        )
    }
}

/// `(stress, tangent)` at `strain` (`makeRodLaw`/`makeBareRodLaw`
/// dispatch): bare steel in compression always, and in tension too when
/// `tension_stiffening` is off; otherwise the TCM/POM inverse (with the
/// startup bridge and rupture) in tension.
pub(super) fn evaluate(f: &CsfmRebarFields, strain: f64) -> (f64, f64) {
    if strain <= 0.0 || !f.tension_stiffening {
        return bare_steel(strain, f);
    }
    if f.routing == CsfmRebarRouting::Pom && strain <= ORIGIN_STRAIN {
        return (0.0, f.es);
    }
    if strain > f.eps_m_at_ft {
        return (f.ft, f.ft / strain);
    }
    let stress = tension_stiffened_stress(strain, f);
    let eps_startup = f.sigma_min / f.es;
    if f.routing == CsfmRebarRouting::Tcm && eps_startup > 0.0 && strain < eps_startup {
        let t = strain / eps_startup;
        let tangent = f.es * (1.0 + 6.0 * t * (1.0 - t));
        return (stress, bounded_tangent(tangent, f));
    }
    (
        stress,
        bounded_tangent(1.0 / forward_strain_slope(stress, f), f),
    )
}

/// Whether `strain` is past rupture (plan §1's sign- and mode-dependent
/// thresholds), mirroring `evaluate`'s own branch choice exactly.
pub(super) fn ruptured(f: &CsfmRebarFields, strain: f64) -> bool {
    if strain <= 0.0 || !f.tension_stiffening {
        strain.abs() > f.eps_u
    } else {
        strain > f.eps_m_at_ft
    }
}

/// `(utilization, LimitKind)` at `strain` (plan §2.5 point 3): the rupture
/// utilization always (the same branch split as [`ruptured`] — tension with
/// tension stiffening uses `eps / eps_m_at_ft`; compression and bare steel
/// both use `|eps| / eps_u`, since `-eps` and `|eps|` agree when `eps <=
/// 0`, which is exactly when that branch applies), and the user strain
/// limit's `|eps| / limit` instead when `rebar_strain_limit` is usable and
/// gives the higher of the two. Utilization `1` is "landed on the limit";
/// unlike [`ruptured`]'s strict `>`, this is not clamped and crosses `1`
/// continuously, so a path can be checked just below and above a landing
/// to tell the two apart (plan §4 G).
pub(super) fn limit_utilization(
    f: &CsfmRebarFields,
    strain: f64,
    rebar_strain_limit: Option<f64>,
) -> (f64, LimitKind) {
    let rupture = if strain <= 0.0 || !f.tension_stiffening {
        strain.abs() / f.eps_u
    } else {
        strain / f.eps_m_at_ft
    };
    match LimitStateConfig::usable(rebar_strain_limit) {
        Some(limit) => {
            let user = strain.abs() / limit;
            if user > rupture {
                (user, LimitKind::RebarStrainLimit)
            } else {
                (rupture, LimitKind::RebarRupture)
            }
        }
        None => (rupture, LimitKind::RebarRupture),
    }
}
