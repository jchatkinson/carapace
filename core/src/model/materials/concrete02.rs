use super::Material;

impl Material {
    /// `fc`/`epsc0`/`fcu`/`epscu` are normalized to be negative
    /// (compression), matching OpenSees' `Concrete02` constructor. `rat` is
    /// the ratio between the unloading slope at `epscu` and the initial
    /// slope; `ft`/`ets` are the tensile strength and tension-softening
    /// modulus.
    pub fn concrete02(fc: f64, epsc0: f64, fcu: f64, epscu: f64, rat: f64, ft: f64, ets: f64) -> Self {
        let fc = -fc.abs();
        let epsc0 = -epsc0.abs();
        let fcu = -fcu.abs();
        let epscu = -epscu.abs();
        Material::Concrete02 {
            fc,
            epsc0,
            fcu,
            epscu,
            rat,
            ft,
            ets,
            min_strain: 0.0,
            tension_strain: 0.0,
            strain: 0.0,
            stress: 0.0,
        }
    }
}

/// Ported from `Concrete02::setTrialStrain` (`Concrete02.cpp`).
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Concrete02 {
        fc,
        epsc0,
        fcu,
        epscu,
        rat,
        ft,
        ets,
        min_strain,
        tension_strain,
        strain: cstrain,
        stress: cstress,
    } = m
    else {
        unreachable!()
    };
    let (fc, epsc0, fcu, epscu, rat, ft, ets) = (*fc, *epsc0, *fcu, *epscu, *rat, *ft, *ets);
    let (mut min_strain, mut tension_strain) = (*min_strain, *tension_strain);
    let (cstrain, cstress) = (*cstrain, *cstress);

    let ec0 = 2.0 * fc / epsc0;

    let (stress, tangent) = if strain < min_strain {
        // Further into compression than ever committed: monotonic
        // compression envelope, and a new low-water mark.
        let (stress, tangent) = concrete02_compr_envlp(strain, fc, epsc0, fcu, epscu);
        min_strain = strain;
        (stress, tangent)
    } else {
        // Unloading/reloading branch: point R (`reload_strain`,
        // `reload_stress`) is where the reloading slope, extended from
        // the previous minimum, would cross the *initial* tangent line
        // — the classic EERC-report construction (`Concrete02.cpp`'s
        // comments cite Eq. 2.31-2.36).
        let reload_strain = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
        let reload_stress = ec0 * reload_strain;
        let (min_stress, _) = concrete02_compr_envlp(min_strain, fc, epsc0, fcu, epscu);
        let reload_slope = (min_stress - reload_stress) / (min_strain - reload_strain);
        // Strain where the reloading line (through the previous
        // minimum) crosses zero stress.
        let zero_stress_strain = min_strain - min_stress / reload_slope;

        if strain <= zero_stress_strain {
            let stress_lower_bound = min_stress + reload_slope * (strain - min_strain);
            let stress_upper_bound = reload_slope * 0.5 * (strain - zero_stress_strain);
            let mut stress = cstress + ec0 * (strain - cstrain);
            let mut tangent = ec0;
            if stress <= stress_lower_bound {
                stress = stress_lower_bound;
                tangent = reload_slope;
            }
            if stress >= stress_upper_bound {
                stress = stress_upper_bound;
                tangent = 0.5 * reload_slope;
            }
            (stress, tangent)
        } else {
            // Tension side: `tension_peak_strain` is the strain at the
            // peak of the tensile stress-strain relation reached so
            // far, shifted by `zero_stress_strain` (Eq. 2.42-2.43).
            let tension_peak_strain = zero_stress_strain + tension_strain;
            if strain <= tension_peak_strain {
                let (peak_stress, _) = concrete02_tens_envlp(tension_strain, ft, ec0, ets);
                let tangent = if tension_strain != 0.0 { peak_stress / tension_strain } else { ec0 };
                (tangent * (strain - zero_stress_strain), tangent)
            } else {
                tension_strain = strain - zero_stress_strain;
                concrete02_tens_envlp(tension_strain, ft, ec0, ets)
            }
        }
    };

    (
        stress,
        tangent,
        Material::Concrete02 {
            fc,
            epsc0,
            fcu,
            epscu,
            rat,
            ft,
            ets,
            min_strain,
            tension_strain,
            strain,
            stress,
        },
    )
}

/// `Concrete02::Compr_Envlp`: same Kent-Scott-Park compression envelope as
/// `Concrete01`'s (see `concrete01::concrete01_envelope`), but note the flat
/// crushed branch's tangent floor here is `1e-10` (matching `Concrete02.cpp`),
/// not the exact `0.0` `concrete01_envelope` uses for the same branch — a
/// genuine difference between the two ports, not a simplification to
/// reconcile away.
fn concrete02_compr_envlp(strain: f64, fc: f64, epsc0: f64, fcu: f64, epscu: f64) -> (f64, f64) {
    let ec0 = 2.0 * fc / epsc0;
    let ratio = strain / epsc0;
    if strain >= epsc0 {
        (fc * ratio * (2.0 - ratio), ec0 * (1.0 - ratio))
    } else if strain > epscu {
        let tangent = (fcu - fc) / (epscu - epsc0);
        (tangent * (strain - epsc0) + fc, tangent)
    } else {
        (fcu, 1.0e-10)
    }
}

/// `Concrete02::Tens_Envlp`: linear elastic up to `ft`, then linear
/// softening at slope `-ets` down to zero, then a flat zero branch (with
/// the same `1e-10` tangent floor as the compression envelope above, for
/// the same singular-tangent-avoidance reason).
fn concrete02_tens_envlp(strain: f64, ft: f64, ec0: f64, ets: f64) -> (f64, f64) {
    let eps0 = ft / ec0;
    let epsu = ft * (1.0 / ets + 1.0 / ec0);
    if strain <= eps0 {
        (strain * ec0, ec0)
    } else if strain <= epsu {
        (ft - ets * (strain - eps0), -ets)
    } else {
        (0.0, 1.0e-10)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concrete02_matches_compression_envelope_on_first_loading() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let strain = -0.001;
        let (stress, tangent) = m.trial_stress_tangent(strain);
        let eta = strain / -0.002;
        let expected_stress = -4.0 * (2.0 * eta - eta * eta);
        assert!((stress - expected_stress).abs() < 1e-9, "expected {expected_stress}, got {stress}");
        assert!((tangent - 2000.0).abs() < 1e-9, "got {tangent}");
    }

    #[test]
    fn concrete02_is_linear_elastic_in_tension_below_ft() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let ec0 = 2.0 * 4.0 / 0.002;
        // ft/Ec0 = 0.4/4000 = 0.0001, so 0.00005 is within the linear
        // range — unlike Concrete01, Concrete02 carries real tension.
        let (stress, tangent) = m.trial_stress_tangent(0.00005);
        assert!((stress - ec0 * 0.00005).abs() < 1e-9, "got {stress}");
        assert!((tangent - ec0).abs() < 1e-9, "got {tangent}");
    }

    #[test]
    fn concrete02_softens_in_tension_past_ft() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        // Past eps0=ft/Ec0=0.0001: linear softening at slope -Ets from ft.
        let (stress, tangent) = m.trial_stress_tangent(0.0005);
        let expected_stress = 0.4 - 50.0 * (0.0005 - 0.0001);
        assert!((stress - expected_stress).abs() < 1e-9, "expected {expected_stress}, got {stress}");
        assert_eq!(tangent, -50.0);
    }

    #[test]
    fn concrete02_unloads_with_degrading_stiffness_below_initial() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let ec0 = 2.0 * 4.0 / 0.002;

        let m = m.commit(-0.003);
        let (_stress, unload_tangent) = m.trial_stress_tangent(-0.002);
        assert!(
            unload_tangent < ec0 && unload_tangent > 0.0,
            "reload slope should be positive but degraded below Ec0={ec0}, got {unload_tangent}"
        );
    }
}
