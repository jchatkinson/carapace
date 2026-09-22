use super::Material;

impl Material {
    /// `fpc`/`epsc0`/`fpcu`/`epscu` are normalized to be negative
    /// (compression), matching OpenSees' `Concrete01` constructor.
    pub fn concrete01(fpc: f64, epsc0: f64, fpcu: f64, epscu: f64) -> Self {
        let fpc = -fpc.abs();
        let epsc0 = -epsc0.abs();
        let fpcu = -fpcu.abs();
        let epscu = -epscu.abs();
        let ec0 = 2.0 * fpc / epsc0;
        Material::Concrete01 {
            fpc,
            epsc0,
            fpcu,
            epscu,
            min_strain: 0.0,
            end_strain: 0.0,
            unload_slope: ec0,
            strain: 0.0,
            stress: 0.0,
        }
    }
}

/// Ported from `Concrete01::setTrialStrain` (`Concrete01.cpp`), the
/// version actually driving the reload/envelope/unload state machine
/// (`determineTrialState` is dead code in the original — never called
/// from `setTrialStrain`).
pub(super) fn evaluate(m: &Material, strain: f64) -> (f64, f64, Material) {
    let Material::Concrete01 {
        fpc,
        epsc0,
        fpcu,
        epscu,
        min_strain,
        end_strain,
        unload_slope,
        strain: cstrain,
        stress: cstress,
    } = m
    else {
        unreachable!()
    };
    let (fpc, epsc0, fpcu, epscu) = (*fpc, *epsc0, *fpcu, *epscu);
    let (min_strain, end_strain, unload_slope) = (*min_strain, *end_strain, *unload_slope);
    let (cstrain, cstress) = (*cstrain, *cstress);

    // Quick return: tension, zero stiffness, no compression-envelope
    // history update — but strain/stress still commit to the new
    // point (point 2 of the M7 stage-2 porting recipe).
    if strain > 0.0 {
        return (
            0.0,
            0.0,
            Material::Concrete01 {
                fpc,
                epsc0,
                fpcu,
                epscu,
                min_strain,
                end_strain,
                unload_slope,
                strain,
                stress: 0.0,
            },
        );
    }

    let temp_stress = cstress + unload_slope * strain - unload_slope * cstrain;

    let (stress, tangent, min_strain, end_strain, unload_slope) = if strain < cstrain {
        // Further into compression: `reload()`.
        let (mut stress, mut tangent, min_strain, end_strain, unload_slope) = if strain <= min_strain {
            let new_min_strain = strain;
            let (stress, tangent) = concrete01_envelope(new_min_strain, fpc, epsc0, fpcu, epscu);
            let (end_strain, unload_slope) = concrete01_unload(new_min_strain, stress, fpc, epsc0, epscu);
            (stress, tangent, new_min_strain, end_strain, unload_slope)
        } else if strain <= end_strain {
            let tangent = unload_slope;
            let stress = tangent * (strain - end_strain);
            (stress, tangent, min_strain, end_strain, unload_slope)
        } else {
            (0.0, 0.0, min_strain, end_strain, unload_slope)
        };

        if temp_stress > stress {
            stress = temp_stress;
            tangent = unload_slope;
        }

        (stress, tangent, min_strain, end_strain, unload_slope)
    } else if temp_stress <= 0.0 {
        (temp_stress, unload_slope, min_strain, end_strain, unload_slope)
    } else {
        (0.0, 0.0, min_strain, end_strain, unload_slope)
    };

    (
        stress,
        tangent,
        Material::Concrete01 {
            fpc,
            epsc0,
            fpcu,
            epscu,
            min_strain,
            end_strain,
            unload_slope,
            strain,
            stress,
        },
    )
}

/// `Concrete01::envelope()`: the Kent-Scott-Park backbone, evaluated at a
/// strain at or past the previous minimum (i.e. loading further into
/// compression along the virgin curve).
fn concrete01_envelope(strain: f64, fpc: f64, epsc0: f64, fpcu: f64, epscu: f64) -> (f64, f64) {
    if strain > epsc0 {
        let eta = strain / epsc0;
        let stress = fpc * (2.0 * eta - eta * eta);
        let ec0 = 2.0 * fpc / epsc0;
        let tangent = ec0 * (1.0 - eta);
        (stress, tangent)
    } else if strain > epscu {
        let tangent = (fpc - fpcu) / (epsc0 - epscu);
        let stress = fpc + tangent * (strain - epsc0);
        (stress, tangent)
    } else {
        (fpcu, 0.0)
    }
}

/// `Concrete01::unload()`: the Karsan-Jirsa degrading unload-slope formula,
/// given the new `min_strain` and the envelope stress just computed there.
fn concrete01_unload(min_strain: f64, stress_at_min: f64, fpc: f64, epsc0: f64, epscu: f64) -> (f64, f64) {
    let mut temp_strain = min_strain;
    if temp_strain < epscu {
        temp_strain = epscu;
    }
    let eta = temp_strain / epsc0;

    let mut ratio = 0.707 * (eta - 2.0) + 0.834;
    if eta < 2.0 {
        ratio = 0.145 * eta * eta + 0.13 * eta;
    }

    let mut end_strain = ratio * epsc0;
    let temp1 = min_strain - end_strain;
    let ec0 = 2.0 * fpc / epsc0;
    let temp2 = stress_at_min / ec0;

    let unload_slope;
    if temp1 > -f64::EPSILON {
        // Should always be negative in practice.
        unload_slope = ec0;
    } else if temp1 <= temp2 {
        end_strain = min_strain - temp1;
        unload_slope = stress_at_min / temp1;
    } else {
        end_strain = min_strain - temp2;
        unload_slope = ec0;
    }

    (end_strain, unload_slope)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concrete01_matches_parabolic_envelope_on_first_loading() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        // Envelope is symmetric-normalized negative; strain -0.001 is
        // within the parabola (|strain| < |epsc0|).
        let strain = -0.001;
        let (stress, _tangent) = m.trial_stress_tangent(strain);
        let eta = strain / -0.002;
        let expected = -4.0 * (2.0 * eta - eta * eta);
        assert!((stress - expected).abs() < 1e-9, "expected {expected}, got {stress}");
    }

    #[test]
    fn concrete01_carries_no_tension() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        assert_eq!(m.trial_stress_tangent(0.001), (0.0, 0.0));
    }

    #[test]
    fn concrete01_unloads_with_degrading_stiffness_below_initial() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        let ec0 = 2.0 * 4.0 / 0.002;

        // Push past peak stress (epsc0=-0.002) then commit, then partially
        // unload — but stay on the unload line (between min_strain and
        // end_strain), not far enough to cross back through zero stress.
        let m = m.commit(-0.003);
        let (_stress, unload_tangent) = m.trial_stress_tangent(-0.002);
        assert!(
            unload_tangent < ec0 && unload_tangent > 0.0,
            "Karsan-Jirsa unload slope should be positive but degraded below Ec0={ec0}, got {unload_tangent}"
        );
    }
}
