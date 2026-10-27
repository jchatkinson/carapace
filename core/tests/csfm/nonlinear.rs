//! Native end-to-end runs of the CSFM golden model fixtures (plan `docs/csfm-plan.md` §4 F),
//! against CalcsApp's TS `solveLoadPath` (`core/tests/csfm/fixtures/models/*.json`). Builds a real
//! `Domain` (`CsfmConcrete` + `CsfmRebar` + `EmbeddedBar` via `embed_polyline`, `model.rs`), solves
//! with the parity settings (plan §2.5 point 2), and compares every increment's displacements and
//! bar stress to TS.
//!
//! None of the current fixtures ever needed a load-control cutback (every `retryHistory` in every
//! fixture has exactly one entry), so this uses plain `Integrator::LoadControl` per schedule entry
//! — one `step()` call per TS `increments[]` entry, 1:1. `LoadControlCutbacks`'s own retry/revert
//! behavior is exercised directly in `core/tests/analysis/load_control_cutbacks.rs` instead.

use std::path::Path;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use serde_json::Value;

use crate::model::{build, BuiltModel};

fn fixture(name: &str) -> Value {
    let path = format!(
        "{}/tests/csfm/fixtures/models/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: invalid JSON: {e}"))
}

fn model_fixture_files() -> Vec<String> {
    let dir = format!("{}/tests/csfm/fixtures/models", env!("CARGO_MANIFEST_DIR"));
    let mut names: Vec<String> = std::fs::read_dir(Path::new(&dir))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.ends_with(".json"))
        .collect();
    names.sort();
    names
}

/// `|a - b| <= abs + rel * max(|a|, |b|)` (plan §4's componentwise tolerance convention).
fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

/// Runs one model fixture's full schedule natively and checks every increment against TS.
/// `displacement_rel`/`bar_stress_rel` are the plan §4 F tolerances (`1e-6`/`1e-5`); passed in
/// (rather than hardcoded) so a future fixture that's known to be more sensitive (e.g. nearer a
/// limit point) can still be run with a documented, explicit widening instead of silently reusing
/// the default.
fn run_and_compare(name: &str, displacement_rel: f64, bar_stress_rel: f64) {
    let fixture = fixture(name);
    let model = &fixture["model"];
    let built = build(model);
    let BuiltModel {
        domain,
        nodes,
        rod_paths,
    } = built;

    let max_load = model["loads"]
        .as_array()
        .unwrap()
        .iter()
        .map(|l| l["value"].as_f64().unwrap().abs())
        .fold(0.0_f64, f64::max);
    let solver = &fixture["solver"];
    let max_iter = solver["maxIterationsPerIncrement"].as_u64().unwrap() as usize;

    let schedule: Vec<f64> = fixture["schedule"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let increments = fixture["increments"].as_array().unwrap();
    assert_eq!(
        schedule.len(),
        increments.len(),
        "{name}: no fixture increment used a load-control cutback (see module doc), so one \
         schedule entry must be one recorded TS increment"
    );

    let locations = fixture["locations"].as_array().unwrap();

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl {
            increment: schedule[0],
        })
        .algorithm(Algorithm::Newton {
            tangent: carapace_core::analysis::TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: (1e-10 * max_load).max(1e-10),
            moment_tol: (1e-10 * max_load).max(1e-10),
            relative_tol: 1e-10,
            displacement_tol: None,
            max_iter,
        })
        .try_build(domain)
        .unwrap_or_else(|e| panic!("{name}: invalid model: {e:?}"));

    for (i, (&nominal, inc)) in schedule.iter().zip(increments.iter()).enumerate() {
        analysis.set_integrator(Integrator::LoadControl { increment: nominal });
        let result = analysis
            .step()
            .unwrap_or_else(|e| panic!("{name}: increment {i}: expected convergence, got {e:?}"));

        let expected_fraction = inc["acceptedLoadFraction"].as_f64().unwrap();
        assert!(
            close(result.load_factor, expected_fraction, 1e-12, 1e-9),
            "{name}: increment {i}: load factor {} vs TS's accepted fraction {}",
            result.load_factor,
            expected_fraction
        );

        let expected_u = inc["displacements"].as_array().unwrap();
        let max_u = expected_u
            .iter()
            .map(|v| v.as_f64().unwrap().abs())
            .fold(0.0_f64, f64::max);
        for (n, &node) in nodes.iter().enumerate() {
            for dof in 0..2 {
                let got = analysis.domain().node(node).displacement[dof];
                let want = expected_u[2 * n + dof].as_f64().unwrap();
                assert!(
                    close(got, want, 1e-12 * max_u.max(1.0), displacement_rel),
                    "{name}: increment {i}: node {n} dof {dof}: got {got}, TS {want}"
                );
            }
        }

        for path in &rod_paths {
            let fixture_rod_count = locations
                .iter()
                .filter(|loc| loc["kind"] == "rod" && loc["pathId"] == path.path_id.as_str())
                .count();
            assert_eq!(
                fixture_rod_count,
                path.bars.len(),
                "{name}: path {}: embed_polyline produced {} segments, TS fixture has {}",
                path.path_id,
                path.bars.len(),
                fixture_rod_count
            );
            for (rod_id, &bar) in path.bars.iter().enumerate() {
                let location_index = locations
                    .iter()
                    .position(|loc| {
                        loc["kind"] == "rod"
                            && loc["pathId"] == path.path_id.as_str()
                            && loc["rodId"].as_u64() == Some(rod_id as u64)
                    })
                    .unwrap_or_else(|| {
                        panic!(
                            "{name}: no location for path {} rodId {rod_id}",
                            path.path_id
                        )
                    });
                let want = inc["stress"][location_index][0].as_f64().unwrap();
                let got = domain_element_stress(analysis.domain(), bar);
                assert!(
                    close(got, want, 1e-9, bar_stress_rel),
                    "{name}: increment {i}: path {} rod {rod_id}: got {got} MPa, TS {want} MPa",
                    path.path_id
                );
            }
        }
    }
}

fn domain_element_stress(
    domain: &carapace_core::model::Domain,
    id: carapace_core::model::ElementId,
) -> f64 {
    domain
        .element_gauss_responses(id)
        .expect("EmbeddedBar reports one gauss-shaped response")[0]
        .stress[0]
}

/// §4 F: every CSFM golden model fixture, run natively end-to-end, must match TS's converged state
/// at every increment (displacements `rel` 1e-6, bar stress `rel` 1e-5).
#[test]
fn every_model_fixture_matches_ts_increment_by_increment() {
    let mut failures = Vec::new();
    for name in model_fixture_files() {
        let result = std::panic::catch_unwind(|| run_and_compare(&name, 1e-6, 1e-5));
        if let Err(payload) = result {
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "panicked with a non-string payload".to_string());
            failures.push(format!("{name}: {message}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} model fixtures failed:\n{}",
        failures.len(),
        model_fixture_files().len(),
        failures.join("\n\n")
    );
}

/// A load increment far smaller than the schedule's puts most Gauss points inside the origin-seed
/// zone (strains under 1e-9). Newton must still converge there: the seeded stress used to be zero
/// against an `E0` tangent, so the residual could not fall below `E0 * 1e-9 * volume` per point.
#[test]
fn a_tiny_load_increment_converges_with_gauss_points_in_the_origin_seed_zone() {
    let fixture = fixture("deep-beam.json");
    let BuiltModel { domain, nodes, .. } = build(&fixture["model"]);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0005 })
        .algorithm(Algorithm::Newton {
            tangent: carapace_core::analysis::TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-3,
            max_iter: 60,
        })
        .try_build(domain)
        .unwrap();
    let result = analysis.step().expect("a tiny increment converges");
    assert!((result.load_factor - 0.0005).abs() < 1e-15);
    assert!(
        nodes
            .iter()
            .any(|&node| analysis.domain().node(node).displacement[1] != 0.0),
        "the beam deflected"
    );
}
