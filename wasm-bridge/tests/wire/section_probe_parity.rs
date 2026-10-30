//! Section probe against CalcsApp's JS sectional-analysis engine. The fixtures are generated, never hand-edited:
//! `cd ../calcs.app && npm run export:section-probe-golden` (see `tests/fixtures/section_probe/README.md`).

use carapace_wasm::section_probe::*;
use serde::Deserialize;
use serde_json::Value;

fn fixture(name: &str) -> Value {
    let path = format!(
        "{}/tests/fixtures/section_probe/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}")))
        .unwrap()
}

fn section() -> SectionProbe {
    SectionProbe::new(&SectionProbeConfig::deserialize(fixture("section.json")).unwrap()).unwrap()
}

fn f64s(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

/// Every probe-only material law matches the engine's stress function at every sampled strain.
#[test]
fn material_laws_match_the_engine() {
    let cases = fixture("materials.json");
    let mut checked = 0;
    for case in cases["cases"].as_array().unwrap() {
        let spec = ProbeMaterialSpec::deserialize(&case["spec"]).unwrap();
        let probe = SectionProbe::new(&SectionProbeConfig {
            core_materials: vec![],
            materials: vec![spec],
            fibers: ProbeFibers {
                y: vec![0.0],
                z: vec![0.0],
                area: vec![1.0],
                material: vec![0],
                group: vec![],
            },
        })
        .unwrap();
        let (strains, stress) = (f64s(&case["strains"]), f64s(&case["stress"]));
        for (e, want) in strains.iter().zip(&stress) {
            let got = probe.evaluate(*e, 0.0, 0.0).n;
            assert!(
                (got - want).abs() <= 1e-9 * want.abs().max(1.0),
                "{}: strain {e}: probe {got} vs engine {want}",
                case["name"]
            );
            checked += 1;
        }
    }
    assert!(checked > 1000, "{checked}");
}

/// PMM boundary states: the probe's resultants at the engine's strain profiles.
#[test]
fn pmm_states_match_the_engine() {
    let probe = section();
    let states = fixture("states.json");
    let mut worst: f64 = 0.0;
    let list = states["states"].as_array().unwrap();
    let m_max = list
        .iter()
        .map(|s| s["m2"].as_f64().unwrap().hypot(s["m3"].as_f64().unwrap()))
        .fold(1.0, f64::max);
    let n_max = list
        .iter()
        .map(|s| s["n"].as_f64().unwrap().abs())
        .fold(1.0, f64::max);
    for s in list {
        let (theta, et, eb) = (
            s["theta"].as_f64().unwrap(),
            s["et"].as_f64().unwrap(),
            s["eb"].as_f64().unwrap(),
        );
        let (e0, kz, ky) = probe.profile_state(theta, et, eb).unwrap();
        let r = probe.evaluate(e0, kz, ky);
        worst = worst
            .max((r.n - s["n"].as_f64().unwrap()).abs() / n_max)
            .max((r.m2 - s["m2"].as_f64().unwrap()).abs() / m_max)
            .max((r.m3 - s["m3"].as_f64().unwrap()).abs() / m_max);
    }
    assert!(
        worst < 1e-9,
        "worst relative difference {worst:e} over {} states",
        list.len()
    );
}

/// Moment-curvature at constant axial force, including a tension case and a rotated direction.
#[test]
fn moment_curvature_matches_the_engine() {
    let probe = section();
    let fx = fixture("mphi.json");
    for case in fx["cases"].as_array().unwrap() {
        let pts = case["points"].as_array().unwrap();
        let peak = pts
            .iter()
            .map(|p| p["m2"].as_f64().unwrap().hypot(p["m3"].as_f64().unwrap()))
            .fold(1.0, f64::max);
        let curvatures: Vec<f64> = pts
            .iter()
            .map(|p| p["curvature"].as_f64().unwrap())
            .collect();
        let res = probe
            .moment_curvature(&MomentCurvatureRequest {
                axial: case["axialForce"].as_f64().unwrap(),
                angle: case["angleRad"].as_f64().unwrap(),
                curvatures,
                tolerance: 1e-9,
            })
            .unwrap();
        let failed: Vec<_> = res
            .points
            .iter()
            .enumerate()
            .filter(|(_, p)| !p.converged)
            .map(|(i, p)| (i, p.curvature, p.eps0, p.residual))
            .take(5)
            .collect();
        assert_eq!(
            res.non_converged, 0,
            "{}: first non-converged (index, curvature, eps0, residual): {failed:?}",
            case["name"]
        );
        let (mut worst, mut worst_strain): (f64, f64) = (0.0, 0.0);
        for (got, want) in res.points.iter().zip(pts).skip(1) {
            if !want["converged"].as_bool().unwrap() {
                continue;
            }
            let dm = (got.m2 - want["m2"].as_f64().unwrap())
                .abs()
                .max((got.m3 - want["m3"].as_f64().unwrap()).abs())
                / peak;
            worst = worst.max(dm);
            // the engine reports extremes clipped at zero
            worst_strain = worst_strain
                .max((got.min_strain.min(0.0) - want["minStrain"].as_f64().unwrap()).abs())
                .max((got.max_strain.max(0.0) - want["maxStrain"].as_f64().unwrap()).abs());
        }
        println!("{}: worst moment difference {worst:e} of peak, worst extreme-strain difference {worst_strain:e}", case["name"]);
        assert!(
            worst < 1e-5,
            "{}: worst moment difference {worst:e} of peak",
            case["name"]
        );
        assert!(
            worst_strain < 1e-7,
            "{}: worst extreme-strain difference {worst_strain:e}",
            case["name"]
        );
    }
}

/// Direct fixed-P / fixed-moment-angle solve against the engine's roots on the same sampled boundary.
#[test]
fn isolines_match_the_engine() {
    let mut probe = section();
    let fx = fixture("isolines.json");
    let (angles, profiles) = (f64s(&fx["angles"]), f64s(&fx["profiles"]));
    let (levels, moment_angles) = (f64s(&fx["levels"]), f64s(&fx["momentAnglesRad"]));
    probe.build_grid(&angles, &profiles, None).unwrap();
    let res = probe
        .solve_isolines(&IsolineRequest {
            levels: levels.clone(),
            moment_angles,
            options: {
                let mut o = IsolineOptions::default();
                if let Ok(v) = std::env::var("RS") {
                    o.residual_seeds = v.parse().unwrap();
                }
                if let Ok(v) = std::env::var("MS") {
                    o.moment_seeds = v.parse().unwrap();
                }
                if let Ok(v) = std::env::var("ST") {
                    o.seed_angle_stride = v.parse().unwrap();
                }
                o
            },
            factored: false,
        })
        .unwrap();
    let roots = fx["roots"].as_array().unwrap();
    let m_max = roots
        .iter()
        .filter_map(|r| r["moment"].as_f64())
        .fold(1.0, f64::max);
    let (mut expected, mut missing, mut wrong_root, mut worst) = (0, vec![], vec![], 0.0f64);
    for root in roots {
        let (ai, li) = (
            root["momentAngleIndex"].as_u64().unwrap() as usize,
            root["levelIndex"].as_u64().unwrap() as usize,
        );
        let got = &res[ai * levels.len() + li];
        let Some(want) = root["moment"].as_f64() else {
            continue;
        };
        expected += 1;
        if !got.found {
            missing.push((ai, li));
            continue;
        }
        let d = (got.moment - want).abs() / m_max;
        if d > 1e-4 {
            wrong_root.push((ai, li, got.moment, want));
        } else {
            worst = worst.max(d);
        }
    }
    println!(
        "{expected} engine roots; {} not found; {} on another root; worst agreement {worst:e}",
        missing.len(),
        wrong_root.len()
    );
    let evals: u32 = res.iter().map(|p| p.evals).sum();
    println!(
        "{:.0} section evaluations per point",
        evals as f64 / res.len() as f64
    );
    let mut counts: Vec<u32> = res.iter().map(|p| p.evals).collect();
    counts.sort();
    println!(
        "evals per point: min {} p25 {} median {} p75 {} p90 {} max {}",
        counts[0],
        counts[counts.len() / 4],
        counts[counts.len() / 2],
        counts[counts.len() * 3 / 4],
        counts[counts.len() * 9 / 10],
        counts[counts.len() - 1]
    );
    assert!(
        missing.is_empty(),
        "engine roots the probe did not find: {missing:?}"
    );
    // Near the axial apex the inner and outer roots merge at a fold and the outer basin is narrow (module docs of
    // `isolines`); only the highest-compression levels may pick the inner root.
    assert!(wrong_root.len() <= 3, "{wrong_root:?}");
    assert!(
        wrong_root.iter().all(|(_, li, _, _)| *li <= 1),
        "{wrong_root:?}"
    );
    assert!(worst < 2e-5, "{worst:e}"); // both sides stop at a normalized residual of 2e-6
}
