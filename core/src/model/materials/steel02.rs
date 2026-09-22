use super::Material;

/// `Steel02`'s `kon` history flag — `Initial` (0) until the first strain
/// increment, then `Positive`/`Negative` (1/2) for which asymptote the
/// current Menegotto-Pinto curve is heading toward. OpenSees also has a
/// `kon == 3` "pinned to `sigini`" state, not ported since `Material`
/// doesn't support `Steel02`'s optional initial-stress offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Steel02Kon {
    Initial,
    Positive,
    Negative,
}

impl Material {
    /// `sigini` (OpenSees' optional initial-stress offset) is not
    /// supported — see the variant's doc comment. Ten parameters is
    /// `Steel02`'s actual OpenSees signature (minus `sigini`/`density`),
    /// not something to trim for its own sake.
    #[allow(clippy::too_many_arguments)]
    pub fn steel02(fy: f64, e0: f64, b: f64, r0: f64, cr1: f64, cr2: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> Self {
        let epsy = fy / e0;
        Material::Steel02 {
            fy,
            e0,
            b,
            r0,
            cr1,
            cr2,
            a1,
            a2,
            a3,
            a4,
            min_strain: -epsy,
            max_strain: epsy,
            pl_strain: 0.0,
            asymptote_strain: 0.0,
            asymptote_stress: 0.0,
            reversal_strain: 0.0,
            reversal_stress: 0.0,
            kon: Steel02Kon::Initial,
            strain: 0.0,
            stress: 0.0,
        }
    }
}

/// Ported from `Steel02::setTrialStrain` (`Steel02.cpp`). OpenSees'
/// `sigini`-pinned quick return (`kon == 0 || kon == 3` with `|dStrain|`
/// under a `DBL_EPSILON`-scale tolerance) is not ported — it exists to
/// let the material sit at a nonzero initial stress before any strain
/// increment, a feature `Material::steel02` doesn't expose. Letting a
/// zero (or tiny) `dstrain` fall through the general formula instead is
/// safe here for the same reason it is in `steel01::evaluate`/
/// `concrete01::evaluate`: plugging `strain == self.strain` back in is
/// idempotent, since the committed state was itself produced by this
/// same formula.
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Steel02 {
        fy,
        e0,
        b,
        r0,
        cr1,
        cr2,
        a1,
        a2,
        a3,
        a4,
        min_strain,
        max_strain,
        pl_strain,
        asymptote_strain,
        asymptote_stress,
        reversal_strain,
        reversal_stress,
        kon,
        strain: cstrain,
        stress: cstress,
    } = m
    else {
        unreachable!()
    };
    let (fy, e0, b, r0, cr1, cr2, a1, a2, a3, a4) = (*fy, *e0, *b, *r0, *cr1, *cr2, *a1, *a2, *a3, *a4);
    let (mut min_strain, mut max_strain, mut pl_strain) = (*min_strain, *max_strain, *pl_strain);
    let (mut asymptote_strain, mut asymptote_stress) = (*asymptote_strain, *asymptote_stress);
    let (mut reversal_strain, mut reversal_stress) = (*reversal_strain, *reversal_stress);
    let mut kon = *kon;
    let (cstrain, cstress) = (*cstrain, *cstress);

    let esh = b * e0;
    let epsy = fy / e0;
    let dstrain = strain - cstrain;

    if kon == Steel02Kon::Initial {
        max_strain = epsy;
        min_strain = -epsy;
        if dstrain < 0.0 {
            kon = Steel02Kon::Negative;
            asymptote_strain = min_strain;
            asymptote_stress = -fy;
            pl_strain = min_strain;
        } else {
            kon = Steel02Kon::Positive;
            asymptote_strain = max_strain;
            asymptote_stress = fy;
            pl_strain = max_strain;
        }
    }

    // Load reversal: recompute the intersection of the elastic and
    // strain-hardening asymptotes, shifted by the isotropic-hardening
    // factor (a1/a2 compression side, a3/a4 tension side) — see
    // `Steel02.cpp`'s comments on the corresponding branches.
    if kon == Steel02Kon::Negative && dstrain > 0.0 {
        kon = Steel02Kon::Positive;
        reversal_strain = cstrain;
        reversal_stress = cstress;
        if cstrain < min_strain {
            min_strain = cstrain;
        }
        let d1 = (max_strain - min_strain) / (2.0 * (a4 * epsy));
        let shift = 1.0 + a3 * d1.powf(0.8);
        asymptote_strain = (fy * shift - esh * epsy * shift - reversal_stress + e0 * reversal_strain) / (e0 - esh);
        asymptote_stress = fy * shift + esh * (asymptote_strain - epsy * shift);
        pl_strain = max_strain;
    } else if kon == Steel02Kon::Positive && dstrain < 0.0 {
        kon = Steel02Kon::Negative;
        reversal_strain = cstrain;
        reversal_stress = cstress;
        if cstrain > max_strain {
            max_strain = cstrain;
        }
        let d1 = (max_strain - min_strain) / (2.0 * (a2 * epsy));
        let shift = 1.0 + a1 * d1.powf(0.8);
        asymptote_strain = (-fy * shift + esh * epsy * shift - reversal_stress + e0 * reversal_strain) / (e0 - esh);
        asymptote_stress = -fy * shift + esh * (asymptote_strain + epsy * shift);
        pl_strain = min_strain;
    }

    // Menegotto-Pinto transition curve between the reversal point
    // `(reversal_strain, reversal_stress)` and the asymptote
    // intersection `(asymptote_strain, asymptote_stress)`; `xi` (the
    // normalized plastic excursion) softens the curvature `r` further
    // from `r0` the larger the excursion has been.
    let xi = ((pl_strain - asymptote_strain) / epsy).abs();
    let r = r0 * (1.0 - (cr1 * xi) / (cr2 + xi));
    let eps_ratio = (strain - reversal_strain) / (asymptote_strain - reversal_strain);
    let dum1 = 1.0 + eps_ratio.abs().powf(r);
    let dum2 = dum1.powf(1.0 / r);

    let stress = (b * eps_ratio + (1.0 - b) * eps_ratio / dum2) * (asymptote_stress - reversal_stress) + reversal_stress;
    let tangent = (b + (1.0 - b) / (dum1 * dum2)) * (asymptote_stress - reversal_stress) / (asymptote_strain - reversal_strain);

    (
        stress,
        tangent,
        Material::Steel02 {
            fy,
            e0,
            b,
            r0,
            cr1,
            cr2,
            a1,
            a2,
            a3,
            a4,
            min_strain,
            max_strain,
            pl_strain,
            asymptote_strain,
            asymptote_stress,
            reversal_strain,
            reversal_stress,
            kon,
            strain,
            stress,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel02_gives_exact_initial_tangent_at_zero_strain() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        // At strain=reversal_strain=0 the Menegotto-Pinto curve's
        // eps_ratio is exactly 0 regardless of R, so this must be exact,
        // not just approximately E0.
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 29000.0));
    }

    #[test]
    fn steel02_tangent_approaches_hardening_slope_far_past_yield() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        let esh = 0.01 * 29000.0;
        // Far enough past the yield asymptote (epsy ≈ 0.00207) that the
        // Menegotto-Pinto curve has converged onto the strain-hardening
        // line.
        let (_stress, tangent) = m.trial_stress_tangent(0.05);
        assert!((tangent - esh).abs() < 1.0, "expected tangent near Esh={esh}, got {tangent}");
    }

    #[test]
    fn steel02_reversal_starts_out_near_elastic_stiffness() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        // Push well past yield and commit, then reverse by a small strain
        // step — right after a reversal the curve is close to its elastic
        // (E0) asymptote, unlike deep into either the loading or
        // hardening branch.
        let m = m.commit(0.02);
        let (_stress, tangent) = m.trial_stress_tangent(0.0199);
        assert!(tangent > 0.9 * 29000.0, "expected near-elastic unloading tangent, got {tangent}");
    }

    #[test]
    fn steel02_tracks_peak_strain_across_a_reversal() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        let m = m.commit(0.02);
        let m = m.commit(-0.01);
        let Material::Steel02 { max_strain, .. } = m else { panic!() };
        // max_strain should have latched onto the 0.02 peak on reversal,
        // not stayed at the virgin epsy ≈ 0.00207.
        assert!((max_strain - 0.02).abs() < 1e-9, "got {max_strain}");
    }
}
