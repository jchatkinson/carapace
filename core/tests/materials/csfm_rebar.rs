//! `Material::CsfmRebar` (`docs/csfm-plan.md` §4 B): stress/tangent parity
//! with the CalcsApp golden fixtures, TCM/POM routing against the
//! `rho_crit` closed form, forward-then-inverse round trips, the analytic
//! tangent against the central-difference checker on smooth branches,
//! the TCM startup bridge's continuity, the numerical policies (tangent
//! clamp, POM origin guard, rupture) tested on their own, bilinear
//! bare-steel behavior, construction errors, and (plan §2.5 point 3) the
//! `limit_utilization` query on every rupture sub-case (tension with
//! tension stiffening, compression, bare steel) plus the optional
//! caller-supplied strain limit, each checked just below and above its
//! utilization-1 landing and against `csfm_rebar_ruptured`'s strict `>`.

use carapace_core::model::{LimitKind, Material, MaterialError};
use carapace_core::testkit::{assert_material_tangent_matches_fd, FdConfig};
use serde_json::Value;

fn fixture() -> Value {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/csfm/fixtures/rebar-points.json"
    );
    let text = std::fs::read_to_string(path).expect("rebar-points.json fixture");
    serde_json::from_str(&text).expect("valid JSON")
}

fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

fn material_for_case(case: &Value) -> (Material, f64) {
    let resolved = &case["resolved"];
    let path = &case["path"];
    let es = resolved["Es"].as_f64().unwrap();
    let fy = resolved["fy"].as_f64().unwrap();
    let ft = resolved["ft"].as_f64().unwrap();
    let eps_u = path["ultimateStrain"].as_f64().unwrap();
    let fctm = resolved["fctm"].as_f64().unwrap();
    let ecm = case["ecm"].as_f64().unwrap();
    let diameter_mm = resolved["diameterMm"].as_f64().unwrap();
    let area_mm2 = resolved["areaMm2"].as_f64().unwrap();
    let effective_area_mm2 = path["effectiveAreaMm2"].as_f64().unwrap();
    let rho_eff = area_mm2 / effective_area_mm2;
    let lambda = case["lambda"].as_f64().unwrap();
    let tension_stiffening = case["tensionStiffening"].as_bool().unwrap();
    let material = Material::csfm_rebar(
        es,
        fy,
        ft,
        eps_u,
        fctm,
        ecm,
        diameter_mm,
        rho_eff,
        lambda,
        tension_stiffening,
    )
    .unwrap();
    (material, rho_eff)
}

/// Stress and tangent match every fixture strain to `rel` 1e-9 (the
/// cancellation-free quadratic roots may differ from TS by rounding; any
/// larger difference is a bug, per the plan).
#[test]
fn stress_and_tangent_match_the_fixtures_at_every_strain() {
    let fixture = fixture();
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let (material, _) = material_for_case(case);
        for point in case["points"].as_array().unwrap() {
            let strain = point["strain"].as_f64().unwrap();
            let expected_stress = point["stress"].as_f64().unwrap();
            let expected_tangent = point["tangent"].as_f64().unwrap();
            let (stress, tangent) = material.trial_stress_tangent(strain);
            assert!(
                close(stress, expected_stress, 1e-9, 1e-9),
                "case {}: strain {strain}: stress {stress} vs fixture {expected_stress}",
                case["label"]
            );
            assert!(
                close(tangent, expected_tangent, 1e-9, 1e-9),
                "case {}: strain {strain}: tangent {tangent} vs fixture {expected_tangent}",
                case["label"]
            );
            checked += 1;
        }
    }
    assert!(checked > 50, "expected many fixture points, got {checked}");
}

/// Routing: `rho_crit = fctm / (fy - (Es/Ecm - 1) * fctm)`, checked on
/// both sides (the TCM fixtures route to TCM, the POM fixture to POM).
#[test]
fn routing_matches_rho_crit_on_both_sides() {
    let fixture = fixture();
    for case in fixture["cases"].as_array().unwrap() {
        let resolved = &case["resolved"];
        let (_, rho_eff) = material_for_case(case);
        let es = resolved["Es"].as_f64().unwrap();
        let fy = resolved["fy"].as_f64().unwrap();
        let fctm = resolved["fctm"].as_f64().unwrap();
        let ecm = case["ecm"].as_f64().unwrap();
        let rho_crit = fctm / (fy - (es / ecm - 1.0) * fctm);
        let expected_routing = resolved["routing"].as_str().unwrap();
        let is_tcm = rho_eff >= rho_crit;
        assert_eq!(
            is_tcm,
            expected_routing == "TCM",
            "case {}: rho_eff={rho_eff} rho_crit={rho_crit}",
            case["label"]
        );
    }
}

/// Forward-then-inverse round trips: TS's forward `eps_m(sigma)` (the
/// fixture's `forward` entries, which already exclude the TCM `[0,
/// sigma_min]` interval and the startup bridge) fed through this law's
/// `evaluate` at that strain must reproduce the original stress to `rel`
/// 1e-9.
#[test]
fn forward_then_inverse_round_trips() {
    let fixture = fixture();
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let (material, _) = material_for_case(case);
        for pair in case["forward"].as_array().unwrap() {
            let expected_stress = pair["stress"].as_f64().unwrap();
            let strain = pair["strain"].as_f64().unwrap();
            let (stress, _) = material.trial_stress_tangent(strain);
            assert!(
                close(stress, expected_stress, 1e-6, 1e-9),
                "case {}: strain {strain}: round trip stress {stress} vs {expected_stress}",
                case["label"]
            );
            checked += 1;
        }
    }
    assert!(
        checked > 5,
        "expected several forward round trips, got {checked}"
    );
}

/// The analytic tangent matches the central-difference checker to `rel`
/// 1e-7 on smooth branches (every fixture point not marked `nearBoundary`,
/// and away from the TCM startup bridge and POM origin guard, both tested
/// on their own below).
#[test]
fn tangent_matches_the_central_difference_checker_on_smooth_branches() {
    let fixture = fixture();
    let config = FdConfig::default();
    let mut checked = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let (material, _) = material_for_case(case);
        let Material::CsfmRebar(fields) = &material else {
            panic!()
        };
        let eps_startup = fields.sigma_min / fields.es;
        for point in case["points"].as_array().unwrap() {
            if point["nearBoundary"].as_bool().unwrap_or(false) {
                continue;
            }
            let strain = point["strain"].as_f64().unwrap();
            if point["status"].as_str() == Some("ruptured") {
                continue;
            }
            // The startup bridge, the POM origin guard, and the tangent
            // clamp are deliberate numerical policies with their own
            // tests elsewhere, not smooth branches of the underlying law
            // (the POM forward relation is quadratic through the origin,
            // so its slope vanishes there and the reciprocal tangent is
            // clamped well before `clampActive` itself turns off).
            if eps_startup > 0.0 && strain > 0.0 && strain < eps_startup * 1.5 {
                continue;
            }
            if point["clampActive"].as_bool().unwrap_or(false) {
                continue;
            }
            if strain.abs() < 1e-7 {
                continue;
            }
            let ran = assert_material_tangent_matches_fd(&material, strain, &config, |_| false);
            assert!(
                ran,
                "checker refused strain {strain} in case {}",
                case["label"]
            );
            checked += 1;
        }
    }
    assert!(
        checked > 50,
        "expected many smooth-branch points, got {checked}"
    );
}

/// Numerical policy, tested on its own: the TCM startup bridge is `C1` at
/// `eps_startup = sigma_min / Es` (value and slope match the first
/// inverse branch there), the stress at the startup strain is `2 *
/// sigma_min`, and the stress just above zero strain carries no
/// `sigma_min` jump.
#[test]
fn tcm_startup_bridge_is_c1_and_has_no_jump_at_the_origin() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "tcm-lambda-1")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    assert!(2.0 * fields.sigma_min <= fields.fy);
    let eps_startup = fields.sigma_min / fields.es;

    let (stress_at_startup, _) = material.trial_stress_tangent(eps_startup);
    assert!(close(stress_at_startup, 2.0 * fields.sigma_min, 1e-6, 1e-9));

    // The join is C1, not C0 at a flat point: stress has slope ~Es through
    // eps_startup, so a naive stress_just_below ~= stress_just_above check
    // would reject a perfectly continuous, non-flat function. Check value
    // continuity by extrapolating the linear branch's own closed form
    // (valid on both sides, since the cubic's value matches it exactly at
    // the boundary) across the step, and check slope continuity directly.
    let h = 1e-5 * eps_startup;
    let (stress_just_below, tangent_just_below) = material.trial_stress_tangent(eps_startup - h);
    let (stress_just_above, tangent_just_above) = material.trial_stress_tangent(eps_startup + h);
    let linear_at = |strain: f64| fields.es * strain + fields.sigma_min;
    assert!(
        close(stress_just_below, linear_at(eps_startup - h), 1e-6, 1e-6),
        "cubic bridge should match the linear branch's value to O(h^2) just below eps_startup: {stress_just_below} vs {}",
        linear_at(eps_startup - h)
    );
    assert!(
        close(stress_just_above, linear_at(eps_startup + h), 1e-6, 1e-9),
        "stress just above eps_startup should be exactly the linear branch: {stress_just_above} vs {}",
        linear_at(eps_startup + h)
    );
    assert!(
        close(tangent_just_below, tangent_just_above, 1e-2, 1e-4),
        "tangent jump at eps_startup: {tangent_just_below} vs {tangent_just_above}"
    );

    let (stress_near_zero, _) = material.trial_stress_tangent(1e-9);
    assert!(
        stress_near_zero.abs() < fields.sigma_min,
        "stress near zero strain should not jump to sigma_min: {stress_near_zero} vs sigma_min {}",
        fields.sigma_min
    );
}

/// Numerical policy: the tangent clamp `[Es / 1e4, 10 * Es]`.
#[test]
fn tangent_clamp_bounds_the_pom_tangent_near_the_origin() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "pom")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    // Just past the origin guard, the POM forward relation's reciprocal
    // tangent is largest (it diverges at the origin itself).
    let (_, tangent) = material.trial_stress_tangent(2e-9);
    assert!(tangent <= 10.0 * fields.es + 1e-6);
    assert!(tangent >= fields.es / 1e4 - 1e-6);
}

/// Numerical policy: the POM origin guard (stress 0, tangent `Es`, strain
/// at or below `ORIGIN_STRAIN`).
#[test]
fn pom_origin_guard_gives_zero_stress_and_es_tangent() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "pom")
        .unwrap();
    let (material, _) = material_for_case(case);
    let (stress, tangent) = material.trial_stress_tangent(5e-10);
    assert_eq!(stress, 0.0);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    assert_eq!(tangent, fields.es);
}

/// Numerical policy: rupture gives stress `±ft` and a small positive
/// secant tangent, and `csfm_rebar_ruptured` agrees, checked just past
/// each threshold in both tension (with tension stiffening) and
/// compression/bare steel.
#[test]
fn rupture_gives_sign_ft_and_a_positive_secant_just_past_each_threshold() {
    let fixture = fixture();
    for case in fixture["cases"].as_array().unwrap() {
        let (material, _) = material_for_case(case);
        let Material::CsfmRebar(fields) = &material else {
            panic!()
        };

        let past_tension = if fields.tension_stiffening {
            fields.eps_m_at_ft * 1.001
        } else {
            fields.eps_u * 1.001
        };
        let (stress, tangent) = material.trial_stress_tangent(past_tension);
        assert!(
            close(stress, fields.ft, 1e-6, 1e-9),
            "case {}",
            case["label"]
        );
        assert!(tangent > 0.0);
        assert!(material.csfm_rebar_ruptured(past_tension));

        let past_compression = -fields.eps_u * 1.001;
        let (stress, tangent) = material.trial_stress_tangent(past_compression);
        assert!(
            close(stress, -fields.ft, 1e-6, 1e-9),
            "case {}",
            case["label"]
        );
        assert!(tangent > 0.0);
        assert!(material.csfm_rebar_ruptured(past_compression));

        let just_below_tension = if fields.tension_stiffening {
            fields.eps_m_at_ft * 0.999
        } else {
            fields.eps_u * 0.999
        };
        assert!(!material.csfm_rebar_ruptured(just_below_tension));
    }
}

/// Compression is bilinear (elastic then hardening), and bare mode
/// (tension stiffening off) is bilinear in both directions.
#[test]
fn compression_and_bare_mode_are_bilinear() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "bare")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    let ey = fields.fy / fields.es;
    for strain in [-0.5 * ey, 0.5 * ey, -1.5 * ey, 1.5 * ey] {
        let (stress, tangent) = material.trial_stress_tangent(strain);
        if strain.abs() <= ey {
            assert!(close(stress, fields.es * strain, 1e-6, 1e-9));
            assert_eq!(tangent, fields.es);
        } else {
            let expected = strain.signum() * (fields.fy + fields.esh * (strain.abs() - ey));
            assert!(close(stress, expected, 1e-6, 1e-9));
            assert_eq!(tangent, fields.esh);
        }
    }
}

/// `rho_eff` out of `(0, 1)` and `ft <= fy` are errors.
#[test]
fn invalid_reinforcement_ratio_and_steel_properties_are_errors() {
    let valid = (
        200_000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 20.0, 0.02, 0.67, true,
    );
    assert!(Material::csfm_rebar(
        valid.0, valid.1, valid.2, valid.3, valid.4, valid.5, valid.6, valid.7, valid.8, valid.9
    )
    .is_ok());
    assert!(matches!(
        Material::csfm_rebar(
            valid.0, valid.1, valid.2, valid.3, valid.4, valid.5, valid.6, 0.0, valid.8, valid.9
        ),
        Err(MaterialError::InvalidReinforcementRatio)
    ));
    assert!(matches!(
        Material::csfm_rebar(
            valid.0, valid.1, valid.2, valid.3, valid.4, valid.5, valid.6, 1.0, valid.8, valid.9
        ),
        Err(MaterialError::InvalidReinforcementRatio)
    ));
    assert!(matches!(
        Material::csfm_rebar(
            valid.0, valid.1, 400.0, valid.3, valid.4, valid.5, valid.6, valid.7, valid.8, valid.9
        ),
        Err(MaterialError::InvalidSteelProperties)
    ));
}

/// TS's discontinuous case (`fck = 30, fy = 500, ft = 550, rho = 0.01,
/// lambda = 2`): `2 * sigma_min > fy`, which CalcsApp now rejects rather
/// than reproducing the overshoot (plan §1, §6). `fctm`/`ecm` here match
/// `resolveConcreteProperties(30)`.
#[test]
fn tcm_startup_exceeding_yield_is_rejected() {
    let fctm = 0.3 * 30.0_f64.powf(2.0 / 3.0);
    let ecm = 22000.0 * (3.8_f64).powf(0.3);
    let result = Material::csfm_rebar(
        200_000.0, 500.0, 550.0, 0.05, fctm, ecm, 20.0, 0.01, 2.0, true,
    );
    assert!(matches!(result, Err(MaterialError::TcmStartupExceedsYield)));
}

/// Tension with tension stiffening: `limit_utilization` is `eps /
/// eps_m_at_ft`, landing at exactly `eps_m_at_ft` — the same strain
/// `csfm_rebar_ruptured` flips at, but with its strict `>` (plan §2.5:
/// "the query reports both ... tests use states just below and just above
/// the limit, not exactly at it").
#[test]
fn tension_rupture_utilization_lands_at_eps_m_at_ft_distinct_from_ruptured() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "tcm-lambda-0.67")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    let eps_m_at_ft = fields.eps_m_at_ft;

    for (strain, expect_above) in [
        (eps_m_at_ft * 0.999, false),
        (eps_m_at_ft, false),
        (eps_m_at_ft * 1.001, true),
    ] {
        let (utilization, kind) = material.csfm_rebar_limit_utilization(strain, None);
        assert_eq!(kind, LimitKind::RebarRupture);
        let expected = strain / eps_m_at_ft;
        assert!(close(utilization, expected, 1e-12, 1e-9), "{utilization}");
        assert_eq!(material.csfm_rebar_ruptured(strain), expect_above);
    }
    assert!(
        (material.csfm_rebar_limit_utilization(eps_m_at_ft, None).0 - 1.0).abs() < 1e-9,
        "landing must read utilization 1 while ruptured() ({}) is still false",
        material.csfm_rebar_ruptured(eps_m_at_ft)
    );
}

/// Compression (`-eps / eps_u`) and bare steel in either direction (`|eps|
/// / eps_u`) are the same formula — both apply exactly when the branch
/// split in `csfm_rebar::limit_utilization`/`ruptured` says so (`strain <=
/// 0` or tension stiffening off) — landing at `eps_u` just like
/// `csfm_rebar_ruptured`'s own threshold there.
#[test]
fn compression_and_bare_steel_utilization_both_land_at_eps_u() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "bare")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    let eps_u = fields.eps_u;
    assert!(!fields.tension_stiffening, "the \"bare\" case must be bare-mode (tension stiffening off) for the bare-steel sub-case below");

    // Compression: tension-stiffened rebar's own `-eps / eps_u` branch
    // (strain <= 0) coincides with bare mode's `|eps| / eps_u` here since
    // this fixture case has tension stiffening off either way, so build a
    // tension-stiffened sibling material to isolate the compression
    // sub-case explicitly.
    let tcm_case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "tcm-lambda-0.67")
        .unwrap();
    let (tcm_material, _) = material_for_case(tcm_case);
    let Material::CsfmRebar(tcm_fields) = &tcm_material else {
        panic!()
    };
    assert!(tcm_fields.tension_stiffening);
    for (strain, expect_above) in [
        (-tcm_fields.eps_u * 0.999, false),
        (-tcm_fields.eps_u, false),
        (-tcm_fields.eps_u * 1.001, true),
    ] {
        let (utilization, kind) = tcm_material.csfm_rebar_limit_utilization(strain, None);
        assert_eq!(kind, LimitKind::RebarRupture);
        assert!(
            close(utilization, strain.abs() / tcm_fields.eps_u, 1e-12, 1e-9),
            "{utilization}"
        );
        assert_eq!(tcm_material.csfm_rebar_ruptured(strain), expect_above);
    }

    // Bare steel, both directions.
    for sign in [1.0, -1.0] {
        for (strain, expect_above) in [
            (sign * eps_u * 0.999, false),
            (sign * eps_u, false),
            (sign * eps_u * 1.001, true),
        ] {
            let (utilization, kind) = material.csfm_rebar_limit_utilization(strain, None);
            assert_eq!(kind, LimitKind::RebarRupture);
            assert!(
                close(utilization, strain.abs() / eps_u, 1e-12, 1e-9),
                "{utilization}"
            );
            assert_eq!(material.csfm_rebar_ruptured(strain), expect_above);
        }
    }
}

/// A caller-supplied `rebar_strain_limit` is opt-in (absent by default,
/// unlike crushing's always-on default on the concrete side) and becomes
/// the controlling kind only once it exceeds the rupture utilization.
#[test]
fn rebar_strain_limit_is_opt_in_and_can_control() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["label"] == "bare")
        .unwrap();
    let (material, _) = material_for_case(case);
    let Material::CsfmRebar(fields) = &material else {
        panic!()
    };
    // Narrower than eps_u, so it controls before rupture does.
    let limit = fields.eps_u * 0.5;
    let (default_utilization, default_kind) = material.csfm_rebar_limit_utilization(limit, None);
    assert_eq!(default_kind, LimitKind::RebarRupture);
    assert!(default_utilization < 1.0, "{default_utilization}");

    let (below, _) = material.csfm_rebar_limit_utilization(limit * 0.999, Some(limit));
    let (at_limit, kind_at) = material.csfm_rebar_limit_utilization(limit, Some(limit));
    let (above, kind_above) = material.csfm_rebar_limit_utilization(limit * 1.001, Some(limit));
    assert_eq!(kind_at, LimitKind::RebarStrainLimit);
    assert_eq!(kind_above, LimitKind::RebarStrainLimit);
    assert!(below < 1.0, "{below}");
    assert!((at_limit - 1.0).abs() < 1e-9, "{at_limit}");
    assert!(above > 1.0, "{above}");

    // A non-finite/non-positive limit is treated as absent.
    for bad in [f64::NAN, 0.0, -1.0] {
        let (u, kind) = material.csfm_rebar_limit_utilization(limit, Some(bad));
        assert_eq!(kind, LimitKind::RebarRupture);
        assert!((u - default_utilization).abs() < 1e-12, "{u}");
    }
}
