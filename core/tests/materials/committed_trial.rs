//! A trial at the strain a material was just committed at must report the
//! same stress and tangent the committed trial did. OpenSees materials return
//! their committed `(stress, tangent)` for a trial strain change under
//! `DBL_EPSILON`; displacement-control predictors, element-state
//! determination and arc-length seeds all evaluate there.

use carapace_core::model::{Material, Pinching4DmgCyc};

fn path(amplitude: f64) -> Vec<f64> {
    let peaks = [
        0.3, -0.3, 0.6, -0.6, 1.0, -1.0, 1.6, -0.4, 0.8, -1.8, 0.2, 2.0, -0.1, 0.0,
    ];
    let mut out = vec![];
    let mut last = 0.0;
    for p in peaks {
        let target = p * amplitude;
        for i in 1..=24 {
            out.push(last + (target - last) * i as f64 / 24.0);
        }
        last = target;
    }
    out
}

fn check(name: &str, material: Material, amplitude: f64) -> bool {
    let mut bad = vec![];
    for sign in [1.0, -1.0] {
        let mut m = material.clone();
        for (i, e) in path(amplitude).into_iter().map(|e| e * sign).enumerate() {
            let (s, t) = m.trial_stress_tangent(e);
            let next = m.commit(e);
            let (s2, t2) = next.trial_stress_tangent(e);
            let tol_s = 1e-9 * (1.0 + s.abs());
            let tol_t = 1e-9 * (1.0 + t.abs());
            if (s - s2).abs() > tol_s || (t - t2).abs() > tol_t {
                bad.push(format!("sign {sign} step {i} e={e:.6}: stress {s:.6}->{s2:.6}, tangent {t:.6}->{t2:.6}"));
            }
            m = next;
        }
    }
    if bad.is_empty() {
        println!("{name}: ok");
    } else {
        println!("{name}: {} mismatches, first: {}", bad.len(), bad[0]);
    }
    bad.is_empty()
}

#[test]
fn committed_trial_is_idempotent() {
    let results = all_materials();
    let failing: Vec<_> = results
        .iter()
        .filter(|(name, ok)| !ok && !KNOWN_FAILING.contains(name))
        .collect();
    assert!(failing.is_empty(), "new idempotence failures: {failing:?}");
}

/// Materials whose OpenSees counterpart returns its previous/committed
/// tangent for a trial at the committed strain, but whose port recomputes.
/// `elastic_pp` is a deliberate divergence: OpenSees treats a trial on the
/// yield surface as plastic (tangent 0), and Newton then oscillates when
/// unloading from yield (verified in OpenSees on `materials/state.rs`'s model).
/// Remove a name once fixed; the test above then guards it.
const KNOWN_FAILING: &[&str] = &["elastic_pp"];

#[test]
fn known_failing_materials_still_fail() {
    let results = all_materials();
    let fixed: Vec<_> = results
        .iter()
        .filter(|(name, ok)| *ok && KNOWN_FAILING.contains(name))
        .collect();
    assert!(
        fixed.is_empty(),
        "now passing, remove from KNOWN_FAILING: {fixed:?}"
    );
}

fn all_materials() -> Vec<(&'static str, bool)> {
    let steel01 = Material::steel01(355.0, 200000.0, 0.02, 0.0, 1.0, 0.0, 1.0);
    let steel01_iso = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
    let steel02 = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
    let concrete01 = Material::concrete01(4.0, 0.002, 3.0, 0.006);
    let concrete02 = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
    let hyst = Material::hysteretic(
        10.0, 0.01, 6.0, 0.02, 2.0, 0.03, -10.0, -0.01, -6.0, -0.02, -2.0, -0.03, 0.5, 0.5, 0.1,
        0.1, 0.1,
    );
    let hyst_plain = Material::hysteretic(
        10.0, 0.01, 6.0, 0.02, 2.0, 0.03, -10.0, -0.01, -6.0, -0.02, -2.0, -0.03, 1.0, 1.0, 0.0,
        0.0, 0.0,
    );
    let pinch = Material::pinching4(
        10.0,
        0.01,
        15.0,
        0.02,
        17.0,
        0.03,
        10.0,
        0.04,
        -10.0,
        -0.01,
        -15.0,
        -0.02,
        -17.0,
        -0.03,
        -10.0,
        -0.04,
        0.5,
        0.25,
        0.05,
        0.5,
        0.25,
        0.05,
        [0.0; 4],
        0.0,
        [0.0; 4],
        0.0,
        [0.0; 4],
        0.0,
        10.0,
        Pinching4DmgCyc::EnergyBased,
    );
    let steel = || Material::steel01(60.0, 29000.0, 0.01, 0.0, 1.0, 0.0, 1.0);
    let cases: Vec<(&'static str, Material, f64)> = vec![
        ("elastic", Material::Elastic { e: 100.0 }, 0.01),
        ("elastic_pp", Material::elastic_pp(100.0, 0.005), 0.02),
        (
            "gap",
            Material::Gap {
                e: 100.0,
                gap: -0.002,
            },
            0.01,
        ),
        ("ent", Material::Ent { e: 100.0 }, 0.01),
        ("steel01", steel01, 0.02),
        ("steel01_iso", steel01_iso, 0.02),
        ("steel02", steel02, 0.02),
        ("concrete01", concrete01, 0.004),
        ("concrete02", concrete02, 0.004),
        ("hysteretic", hyst, 0.02),
        ("hysteretic_plain", hyst_plain, 0.02),
        ("pinching4", pinch, 0.02),
        (
            "parallel",
            Material::parallel(vec![steel(), Material::concrete01(4.0, 0.002, 3.0, 0.006)]),
            0.004,
        ),
        (
            "series",
            Material::series(vec![Material::Elastic { e: 500.0 }, steel()]),
            0.02,
        ),
        ("min_max", Material::min_max(steel(), -0.01, 0.015), 0.02),
    ];
    cases
        .into_iter()
        .map(|(name, m, amplitude)| (name, check(name, m, amplitude)))
        .collect()
}
