//! `PlaneMaterial::CsfmConcrete`: the CSFM concrete law (plan
//! `docs/csfm-plan.md` §1 point 2, §2.1), ported from CalcsApp's
//! `evaluateConcrete`/`makeConcreteContinuumLaw`
//! (`../calcs.app/lib/fem/nonlinear-material.ts`,
//! `../calcs.app/lib/fem/csfm-nonlinear.ts`).
//!
//! A stateless rotating principal-strain law with `nu = 0`: stress is a
//! function of the current principal strains alone, and the material
//! direction is always taken at the current principal angle (no separate
//! "fixed crack" memory). The Cartesian tangent is the **analytic**
//! consistent tangent of that rotating map, `D = T(theta)^T * C_p *
//! T(theta)` (plan §2.1) — not a rotated secant, and not CalcsApp's
//! regime-aware finite difference, which this law replaces.

use super::{PlaneMaterial, PlaneMaterialError, PlaneMatrix, PlaneVector};

/// Which concrete compression envelope this law uses (plan §1 point 2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CsfmConcreteLaw {
    /// Eurocode-style parabola (0 to `eps_c2`) then flat plateau (`eps_c2`
    /// to `eps_cu2` and beyond).
    ParabolaRectangle,
    /// Linear elastic (`Ecm`) up to the strength cap, then flat (ideally
    /// plastic).
    ElasticIdeallyPlastic,
}

/// The resolved field set behind `PlaneMaterial::CsfmConcrete` (see that
/// variant's doc comment for why it's boxed). `fctm`/`eta_fc`/`e0` are
/// derived once at construction: `fctm = 0.3 * fck^(2/3)`, `eta_fc =
/// min(1, (30/fck)^(1/3))`, and `e0` is the uncracked compression tangent
/// (`Ecm` for elastic-ideally-plastic, `2 * eta_fc * fck / |eps_c2|` for
/// parabola-rectangle — the parabola's own slope at the origin) used both
/// as the exact-origin isotropic seed and as the compression branch's
/// strength-independent elastic slope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CsfmConcreteFields {
    pub fck: f64,
    pub ecm: f64,
    pub law: CsfmConcreteLaw,
    pub softening: bool,
    /// Tensile regularization factor `gamma` (plan §1): `None` gives the
    /// exact zero-tension law; `Some(gamma)` (`0 < gamma < 1`) adds
    /// `sigma_t = gamma * fctm * (1 - exp(-eps / eps_reg))` with `eps_reg
    /// = gamma * fctm / Ecm` derived (not configured), so the tension
    /// tangent at the origin is exactly `Ecm`.
    pub gamma: Option<f64>,
    pub fctm: f64,
    pub eta_fc: f64,
    pub e0: f64,
}

const EPS_C2: f64 = -0.002;
/// The exact-origin isotropic seed applies within this strain magnitude on
/// every Cartesian component (plan §1, §4 A — matches CalcsApp's
/// `ORIGIN_STRAIN`). It is a numerical policy, not the law's derivative
/// there, and is never checked against a finite difference.
const ORIGIN_STRAIN: f64 = 1e-9;

impl PlaneMaterial {
    /// Builds a `CsfmConcrete`. `ecm_override` replaces the derived `Ecm =
    /// 22000 * ((fck + 8) / 10)^0.3` when given (plan §1). There is no
    /// crack-transition band (plan §1, §6): both principal directions
    /// always use this one scalar law.
    ///
    /// # Errors
    ///
    /// [`PlaneMaterialError::InvalidFck`] if `fck` is not in `0 < fck <=
    /// 50` MPa, or `ecm_override` is given and not finite and positive;
    /// [`PlaneMaterialError::InvalidRegularizationGamma`] if `gamma` is
    /// given and not in `0 < gamma < 1`.
    pub fn csfm_concrete(
        fck: f64,
        ecm_override: Option<f64>,
        law: CsfmConcreteLaw,
        softening: bool,
        gamma: Option<f64>,
    ) -> Result<Self, PlaneMaterialError> {
        if !(fck > 0.0 && fck <= 50.0) {
            return Err(PlaneMaterialError::InvalidFck);
        }
        let ecm = match ecm_override {
            Some(ecm) if ecm.is_finite() && ecm > 0.0 => ecm,
            Some(_) => return Err(PlaneMaterialError::InvalidFck),
            None => 22000.0 * ((fck + 8.0) / 10.0).powf(0.3),
        };
        if let Some(gamma) = gamma {
            if !(gamma.is_finite() && gamma > 0.0 && gamma < 1.0) {
                return Err(PlaneMaterialError::InvalidRegularizationGamma);
            }
        }
        let fctm = 0.3 * fck.powf(2.0 / 3.0);
        let eta_fc = (30.0 / fck).powf(1.0 / 3.0).min(1.0);
        let e0 = match law {
            CsfmConcreteLaw::ElasticIdeallyPlastic => ecm,
            CsfmConcreteLaw::ParabolaRectangle => 2.0 * eta_fc * fck / EPS_C2.abs(),
        };
        Ok(PlaneMaterial::CsfmConcrete(Box::new(CsfmConcreteFields {
            fck,
            ecm,
            law,
            softening,
            gamma,
            fctm,
            eta_fc,
            e0,
        })))
    }
}

/// Principal strains (`eps1 >= eps2`) and the angle from global `x` to the
/// `eps1` direction, from the Cartesian engineering-strain vector. Same
/// convention as CalcsApp's `principal` (`csfm-nonlinear.ts`).
fn principal(strain: &PlaneVector) -> (f64, f64, f64) {
    let (ex, ey, gamma) = (strain[0], strain[1], strain[2]);
    let average = (ex + ey) / 2.0;
    let radius = ((ex - ey) / 2.0).hypot(gamma / 2.0);
    let theta = 0.5 * gamma.atan2(ex - ey);
    (average + radius, average - radius, theta)
}

/// The strain-transformation matrix (global -> principal axes at `theta`),
/// identical in convention to the `Orthotropic` arm's `t` in this file:
/// `D = T^T * D' * T`, stress transforms by `T`'s inverse transpose.
fn rotation(theta: f64) -> PlaneMatrix {
    let (s, c) = theta.sin_cos();
    PlaneMatrix::new(
        c * c,
        s * s,
        c * s, //
        s * s,
        c * c,
        -c * s, //
        -2.0 * c * s,
        2.0 * c * s,
        c * c - s * s,
    )
}

/// Kaufmann's compression-softening coefficient (plan §1): `1` whenever
/// softening is off or `eps1 <= 0` (the `max(0, eps1)` inside the formula
/// already gives exactly `1` there, so no separate branch is needed).
fn kc2_of(eps1: f64, softening: bool) -> f64 {
    if !softening {
        1.0
    } else {
        (1.0 / (0.8 + 140.0 * eps1.max(0.0))).min(1.0)
    }
}

/// `d(kc2)/d(eps1)`: zero while the `min(1, ...)` cap is active (`eps1 <=
/// 0.2 / 140`, which also covers every `eps1 <= 0`), or softening is off;
/// otherwise the derivative of `1 / (0.8 + 140 * eps1)`.
fn dkc2_deps1(eps1: f64, softening: bool) -> f64 {
    if !softening || eps1 <= 0.2 / 140.0 {
        0.0
    } else {
        -140.0 / (0.8 + 140.0 * eps1).powi(2)
    }
}

/// The compression branch (`compressionStress`) and its derivatives
/// `(stress, d(stress)/d(eps), d(stress)/d(fc))` — the second needed only
/// for the `kc2` coupling term, where `fc` itself varies with `eps1`.
/// Requires `eps <= 0`.
fn compression(eps: f64, fc: f64, law: CsfmConcreteLaw, ecm: f64) -> (f64, f64, f64) {
    match law {
        CsfmConcreteLaw::ParabolaRectangle => {
            if eps >= EPS_C2 {
                // r in [0, 1] from 0 at eps=0 to 1 at eps=eps_c2.
                let r = eps / EPS_C2;
                let shape = 2.0 * r - r * r;
                (-fc * shape, -2.0 * fc * (1.0 - r) / EPS_C2, -shape)
            } else {
                (-fc, 0.0, -1.0)
            }
        }
        CsfmConcreteLaw::ElasticIdeallyPlastic => {
            let raw = ecm * eps;
            if raw > -fc {
                (raw, ecm, 0.0)
            } else {
                (-fc, 0.0, -1.0)
            }
        }
    }
}

/// The (optionally regularized) tension branch (`residualTensileStress`)
/// and its derivative. Requires `eps > 0`; has no dependence on `fc`.
fn tension(eps: f64, fctm: f64, ecm: f64, gamma: Option<f64>) -> (f64, f64) {
    match gamma {
        None => (0.0, 0.0),
        Some(gamma) => {
            let eps_reg = gamma * fctm / ecm;
            let x = eps / eps_reg;
            // 1 - exp(-x) written as -expm1(-x) (plan §2.1).
            let stress = gamma * fctm * -(-x).exp_m1();
            let slope = gamma * fctm / eps_reg * (-x).exp();
            (stress, slope)
        }
    }
}

/// `(stress, d(stress)/d(eps))` at `eps`, dispatching compression/tension
/// with the given inclusive-at-zero convention for which branch owns
/// `eps == 0` exactly (the two principal directions use different
/// conventions here, matching CalcsApp's `epsMax <= 0` / `epsMin < 0` —
/// plan §2.1, "the branch at the kink `eps_1 = 0` follows the TS regime
/// rule").
fn branch(eps: f64, fc: f64, f: &CsfmConcreteFields, compression_at_zero: bool) -> (f64, f64) {
    let is_compression = if compression_at_zero {
        eps <= 0.0
    } else {
        eps < 0.0
    };
    if is_compression {
        let (stress, slope, _d_fc) = compression(eps, fc, f.law, f.ecm);
        (stress, slope)
    } else {
        tension(eps, f.fctm, f.ecm, f.gamma)
    }
}

/// Stress and the analytic Cartesian tangent at `strain` (plan §2.1).
pub(super) fn evaluate(f: &CsfmConcreteFields, strain: &PlaneVector) -> (PlaneVector, PlaneMatrix) {
    if strain.iter().all(|value| value.abs() < ORIGIN_STRAIN) {
        let d = PlaneMatrix::new(
            f.e0,
            0.0,
            0.0, //
            0.0,
            f.e0,
            0.0, //
            0.0,
            0.0,
            f.e0 / 2.0,
        );
        return (PlaneVector::zeros(), d);
    }

    let (eps1, eps2, theta) = principal(strain);
    let fc1 = f.eta_fc * f.fck;
    let kc2 = kc2_of(eps1, f.softening);
    let fc2 = f.eta_fc * kc2 * f.fck;

    // eps1 == 0 exactly counts as compression (CalcsApp's `epsMax <= 0`);
    // eps2 == 0 exactly counts as tension (`epsMin < 0`, strict).
    let (sigma1, dsigma1_deps1) = branch(eps1, fc1, f, true);
    let (sigma2, dsigma2_deps2) = branch(eps2, fc2, f, false);

    let dfc2_deps1 = f.eta_fc * f.fck * dkc2_deps1(eps1, f.softening);
    let dsigma2_deps1 = if eps2 < 0.0 && dfc2_deps1 != 0.0 {
        let (_, _, d_fc) = compression(eps2, fc2, f.law, f.ecm);
        d_fc * dfc2_deps1
    } else {
        0.0
    };

    let gap = eps1 - eps2;
    let scale = eps1.abs().max(eps2.abs()).max(1e-12);
    let g_star = if gap.abs() <= 1e-8 * scale {
        // lim_{gap->0} (sigma1 - sigma2) / (2 * gap), with eps1 = avg +
        // gap/2, eps2 = avg - gap/2: by L'Hopital the derivative of the
        // numerator wrt gap at gap=0 is (dsigma1_deps1 + dsigma2_deps2) /
        // 2 (chain rule through the 1/2 slopes), and the limit is that
        // divided by the denominator's own slope (2), i.e. divided by 4
        // overall.
        0.25 * (dsigma1_deps1 + dsigma2_deps2)
    } else {
        (sigma1 - sigma2) / (2.0 * gap)
    };

    let c_p = PlaneMatrix::new(
        dsigma1_deps1,
        0.0,
        0.0, //
        dsigma2_deps1,
        dsigma2_deps2,
        0.0, //
        0.0,
        0.0,
        g_star,
    );
    let t = rotation(theta);
    let d = t.transpose() * c_p * t;
    let stress = t.transpose() * PlaneVector::new(sigma1, sigma2, 0.0);
    (stress, d)
}
