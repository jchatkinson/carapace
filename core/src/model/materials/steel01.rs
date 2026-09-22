use super::Material;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Steel01Loading {
    None,
    Loading,
    Unloading,
}

impl Material {
    pub fn steel01(fy: f64, e0: f64, b: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> Self {
        Material::Steel01 {
            fy,
            e0,
            b,
            a1,
            a2,
            a3,
            a4,
            min_strain: 0.0,
            max_strain: 0.0,
            shift_p: 1.0,
            shift_n: 1.0,
            loading: Steel01Loading::None,
            strain: 0.0,
            stress: 0.0,
        }
    }
}

/// Ported from `Steel01::determineTrialState` (`Steel01.cpp`). `m` stands in
/// for the `C*` fields; `strain` is the trial (`Tstrain`).
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Steel01 {
        fy,
        e0,
        b,
        a1,
        a2,
        a3,
        a4,
        min_strain,
        max_strain,
        shift_p,
        shift_n,
        loading,
        strain: cstrain,
        stress: cstress,
    } = m
    else {
        unreachable!()
    };
    let (fy, e0, b, a1, a2, a3, a4) = (*fy, *e0, *b, *a1, *a2, *a3, *a4);
    let (mut min_strain, mut max_strain, mut shift_p, mut shift_n) =
        (*min_strain, *max_strain, *shift_p, *shift_n);
    let mut loading = *loading;
    let (cstrain, cstress) = (*cstrain, *cstress);

    let dstrain = strain - cstrain;

    let fy_one_minus_b = fy * (1.0 - b);
    let esh = b * e0;
    let epsy = fy / e0;

    let c1 = esh * strain;
    let c2 = shift_n * fy_one_minus_b;
    let c3 = shift_p * fy_one_minus_b;
    let c = cstress + e0 * dstrain;

    // The three-way clamp `max(c1-c2, min(c1+c3, c))` — see the
    // commented-out simplification in `Steel01.cpp`, which the source
    // says only exists split into two `if`s as an optimizer
    // workaround; the single-expression form is used here.
    let stress = (c1 - c2).max((c1 + c3).min(c));
    let tangent = if (stress - c).abs() < f64::EPSILON { e0 } else { esh };

    if loading == Steel01Loading::None && dstrain != 0.0 {
        loading = if dstrain > 0.0 {
            Steel01Loading::Loading
        } else {
            Steel01Loading::Unloading
        };
    }

    // Transition from loading to unloading.
    if loading == Steel01Loading::Loading && dstrain < 0.0 {
        loading = Steel01Loading::Unloading;
        if cstrain > max_strain {
            max_strain = cstrain;
        }
        shift_n = 1.0 + a1 * ((max_strain - min_strain) / (2.0 * a2 * epsy)).powf(0.8);
    }

    // Transition from unloading to loading.
    if loading == Steel01Loading::Unloading && dstrain > 0.0 {
        loading = Steel01Loading::Loading;
        if cstrain < min_strain {
            min_strain = cstrain;
        }
        shift_p = 1.0 + a3 * ((max_strain - min_strain) / (2.0 * a4 * epsy)).powf(0.8);
    }

    (
        stress,
        tangent,
        Material::Steel01 {
            fy,
            e0,
            b,
            a1,
            a2,
            a3,
            a4,
            min_strain,
            max_strain,
            shift_p,
            shift_n,
            loading,
            strain,
            stress,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel01_matches_elastic_on_first_loading_within_yield() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        let (stress, tangent) = m.trial_stress_tangent(0.001);
        assert!((stress - 29000.0 * 0.001).abs() < 1e-9, "got {stress}");
        assert_eq!(tangent, 29000.0);
    }

    #[test]
    fn steel01_bilinear_hardening_past_yield() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        // epsy = fy/E0 = 60/29000 ≈ 0.002069. Push well past it.
        let (stress, tangent) = m.trial_stress_tangent(0.02);
        let esh = 0.01 * 29000.0;
        // On first loading (shift_p=shift_n=1, min=max=0) the upper
        // bounding line is c1+c3 = Esh*strain + fy*(1-b).
        let expected = esh * 0.02 + 60.0 * (1.0 - 0.01);
        assert!((stress - expected).abs() < 1e-6, "expected {expected}, got {stress}");
        assert!((tangent - esh).abs() < 1e-9);
    }

    #[test]
    fn steel01_isotropic_hardening_shifts_the_bound_on_reversal() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        // Load well into the hardening range, commit, then reverse twice —
        // the first reversal (loading -> unloading) shifts shift_n, the
        // second (unloading -> loading) shifts shift_p.
        let m = m.commit(0.02);
        let m = m.commit(-0.02);
        let Material::Steel01 { shift_n, .. } = m else { panic!() };
        assert!(shift_n > 1.0, "shift_n should have grown after the first reversal: {shift_n}");

        let m = m.commit(0.02);
        let Material::Steel01 { shift_p, .. } = m else { panic!() };
        assert!(shift_p > 1.0, "shift_p should have grown after the second reversal: {shift_p}");
    }
}
