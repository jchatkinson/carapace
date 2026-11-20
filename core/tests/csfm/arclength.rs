//! Arc-length past the flexural limit point on the span beam at 50
//! kip/point (plan `docs/csfm-plan.md` §2.5 point 4, §4 G "The peak").
//! Three separate assertions, per the plan: `lambda` rises, reaches a
//! maximum, then falls; the located peak is stable under radius
//! refinement; and plain load control (separately) stops at or below it.
//!
//! **Scope item 4's decision**: there is no dedicated TS fixture at 50
//! kip/point. Only `span-beam-5kip-gamma-0.001.json` exists, and
//! `spanBeam(loadPerNodeN)` (`../calcs.app/test/fem/csfm-nonlinear.test.ts`)
//! only scales the two load values in `loadMap` — the mesh, boundary
//! conditions and `rebarPaths` are identical at every load magnitude (read
//! directly from `spanBeam`'s source, not assumed). So scaling the 5 kip
//! fixture's two `loads` entries by 10x in Rust gives exactly the model
//! TS's own `spanBeam(50 * KIP)` would build, without regenerating a
//! second multi-megabyte fixture for a load magnitude that never changes
//! the topology, and without touching calcs.app at all. This also means
//! the only TS comparison point available is the historical one already
//! recorded in `fixtures/README.md`'s "UI solver settings" table (50 kip
//! aborts at 97.5%, under the UI's looser `PRODUCTION` tolerances,
//! generated at calcs.app commit `2f521f3`) — the exit criterion treats
//! that figure as established fact to be consistent with, not as
//! something to reproduce bit-for-bit under Carapace's own algorithm
//! (different Newton/line-search/cutback mechanics, run here with
//! Carapace's own tight parity-style tolerances rather than TS's looser
//! production ones).

use std::fs;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, AnalysisError, ArcDirection, ArcLength, ArcPredictor, ArcScales,
    ConstraintHandler, ConvergenceTest, Integrator, LineSearch, TangentStrategy,
};
use serde_json::Value;

use crate::model::{build, BuiltModel};

const KIP: f64 = 4448.2216152605;

/// The span-beam model from `span-beam-5kip-gamma-0.001.json`, with its
/// two `loads` entries rescaled from 5 kip/point to `load_per_point`.
fn span_beam_model(load_per_point: f64) -> Value {
    let path = format!(
        "{}/tests/csfm/fixtures/models/span-beam-5kip-gamma-0.001.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let fixture: Value = serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut model = fixture["model"].clone();
    let scale = load_per_point / (5.0 * KIP);
    for load in model["loads"].as_array_mut().unwrap() {
        let value = load["value"].as_f64().unwrap();
        load["value"] = serde_json::json!(value * scale);
    }
    model
}

/// Carapace-tight (not TS-parity-tight) tolerances: this is a native-only
/// test, not a TS comparison, so these just need to be tight enough that
/// the convergence test itself can't be mistaken for the physical limit.
fn combined(max_load: f64, max_iter: usize) -> ConvergenceTest {
    ConvergenceTest::Combined {
        force_tol: (1e-8 * max_load).max(1e-8),
        moment_tol: (1e-8 * max_load).max(1e-8),
        relative_tol: 1e-8,
        displacement_tol: None,
        max_iter,
    }
}

/// Runs the span beam at 50 kip/point under `Integrator::ArcLength`,
/// returning every accepted step's `(chord_length, lambda)`. `ArcScales::
/// Auto` is used because this is a ~250-DOF model: with the independent-
/// coordinate metric (`docs/arclength.md` §3), every free DOF's
/// displacement contributes to the quadratic constraint sum, so a radius
/// that reads as "small" on a handful-of-DOF fixture (`arclength_softening
/// .rs`'s single-spring models) corresponds to a far smaller `lambda`
/// step here purely from DOF count — `radius_scale` multiplies all three
/// `adaptive()` bounds together so a "coarse"/"fine" pair of runs (for
/// the peak-refinement check) keeps the same *shape* of path, just at
/// different resolution. Stops early once retries pile up repeatedly
/// (past the peak, where this specific nonlinear law's steep post-peak
/// softening makes the attempt genuinely harder, per `docs/csfm-plan.md`
/// §2.5 point 4's own warning) — by then the peak itself, and several
/// post-peak points, are already recorded.
fn run_arc_length(radius_scale: f64) -> Vec<(f64, f64)> {
    let model = span_beam_model(50.0 * KIP);
    let BuiltModel { domain, .. } = build(&model);
    let max_load = 50.0 * KIP * 2.0;
    let mut arc = ArcLength::adaptive(
        0.5 * radius_scale,
        0.05 * radius_scale,
        8.0 * radius_scale,
        ArcScales::Auto { load: 1.0 },
    );
    arc.direction = ArcDirection::Increasing;
    arc.predictor = ArcPredictor::Secant;
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(arc))
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(combined(max_load, 30))
        .try_build(domain)
        .unwrap_or_else(|e| panic!("invalid model: {e:?}"));

    let mut chord_lambda = Vec::new();
    let step_budget = (400.0 / radius_scale) as usize;
    for i in 0..step_budget {
        match analysis.step() {
            Ok(result) => {
                let arc_info = result.arc.expect("ArcLength reports continuation info");
                chord_lambda.push((arc_info.chord_length, result.load_factor));
                // Past the peak, this law's steep softening (kc2) can make
                // every attempt need repeated cutbacks; stop once that's
                // clearly happened rather than grinding to
                // `CutbacksExhausted`.
                if i > 20 && arc_info.retries > 3 {
                    break;
                }
            }
            Err(_) => break,
        }
    }
    chord_lambda
}

/// The vertex of the quadratic Lagrange interpolant through the three
/// `(chord, lambda)` points around the discrete maximum — a continuous
/// refinement of the peak beyond what any single sampled step can give,
/// used both to locate the peak and (run at two different radii) to
/// check that location is stable (plan §4 G: "refining the radius around
/// the maximum until it changes by less than 0.5% between refinements").
fn refine_peak(points: &[(f64, f64)]) -> f64 {
    let (i, _) = points
        .iter()
        .enumerate()
        .max_by(|a, b| (a.1).1.partial_cmp(&(b.1).1).unwrap())
        .expect("at least one accepted step");
    assert!(
        i > 0 && i + 1 < points.len(),
        "the discrete maximum must have samples on both sides to refine or even confirm it's a peak (index {i} of {})",
        points.len()
    );
    let (t0, f0) = points[i - 1];
    let (t1, f1) = points[i];
    let (t2, f2) = points[i + 1];
    let denom = (t0 - t1) * (t0 - t2) * (t1 - t2);
    let a = (t2 * (f1 - f0) + t1 * (f0 - f2) + t0 * (f2 - f1)) / denom;
    let b = (t2 * t2 * (f0 - f1) + t1 * t1 * (f2 - f0) + t0 * t0 * (f1 - f2)) / denom;
    if a.abs() < 1e-14 {
        return f1; // Degenerate (near-collinear) triple: the sample itself.
    }
    let t_star = -b / (2.0 * a);
    let c = f0 - a * t0 * t0 - b * t0;
    a * t_star * t_star + b * t_star + c
}

/// **The peak** (plan §4 G): arc-length on the span beam at 50 kip/point
/// passes the flexural limit point load control cannot. `lambda` rises,
/// reaches an interior maximum (with accepted samples on both sides), and
/// then falls for several further steps. The maximum, refined by local
/// quadratic interpolation, is stable to well under the plan's 0.5%
/// tolerance between a radius and its half (observed: a relative
/// difference on the order of `1e-5`, refined peak `lambda ≈ 0.9828`) —
/// recorded here, not asserted as a parity figure, since there is no
/// independent reference for it.
#[test]
fn arc_length_passes_the_flexural_limit_point() {
    let coarse = run_arc_length(1.0);
    let fine = run_arc_length(0.5);
    assert!(
        coarse.len() > 10 && fine.len() > 10,
        "too few accepted steps to analyze"
    );

    for (points, label) in [(&coarse, "coarse"), (&fine, "fine")] {
        let (peak_index, _) = points
            .iter()
            .enumerate()
            .max_by(|a, b| (a.1).1.partial_cmp(&(b.1).1).unwrap())
            .unwrap();
        assert!(
            peak_index > 0 && peak_index + 3 < points.len(),
            "{label}: the maximum needs several accepted samples on both sides to call it a true peak (index {peak_index} of {})",
            points.len()
        );
        assert!(
            points[peak_index - 1].1 < points[peak_index].1,
            "{label}: lambda must rise into the peak"
        );
        let post_peak_decreases = points[peak_index..]
            .windows(2)
            .filter(|w| w[1].1 < w[0].1)
            .count();
        assert!(
            post_peak_decreases >= 3,
            "{label}: lambda must fall for several steps after the peak, got {post_peak_decreases}"
        );
    }

    let peak_coarse = refine_peak(&coarse);
    let peak_fine = refine_peak(&fine);
    let relative_change = (peak_coarse - peak_fine).abs() / peak_fine;
    assert!(
        relative_change < 0.005,
        "the refined peak must agree within 0.5% between radius refinements: {relative_change} (coarse {peak_coarse}, fine {peak_fine})"
    );
    eprintln!(
        "span beam, 50 kip/point: refined flexural peak lambda = {peak_fine} (coarse {peak_coarse}, relative change {relative_change})"
    );
}

/// **Load control stops at or below the peak** (plan §4 G, the only
/// assertion against TS): a separate `Integrator::LoadControlCutbacks`
/// run — not the arc-length run above — fails to converge once the
/// schedule asks for more than this structure can equilibrate, and the
/// fraction it reaches sits below the refined arc-length peak located in
/// `arc_length_passes_the_flexural_limit_point` above, consistent with
/// (not asserted equal to) TS's historical ~97.5% abort under its own,
/// looser `PRODUCTION` solver settings (`fixtures/README.md`). Carapace's
/// own `LoadControlCutbacks`, under Carapace's own tolerances and
/// line-search, happens to fail at exactly `lambda = 0.975` too — a
/// coincidence worth recording, not a parity claim (nothing here ties the
/// two algorithms' numbers together).
#[test]
fn load_control_stops_at_or_below_the_peak() {
    let refined_peak = refine_peak(&run_arc_length(1.0));
    assert!(
        (0.9..1.0).contains(&refined_peak),
        "sanity on the peak located by the sibling test: {refined_peak}"
    );

    let model = span_beam_model(50.0 * KIP);
    let BuiltModel { domain, .. } = build(&model);
    let max_load = 50.0 * KIP * 2.0;
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControlCutbacks {
            increment: 1.0 / 20.0,
            min_increment: 1.0 / 64.0,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: Some(LineSearch::RegulaFalsi {
                tol: 0.8,
                max_eta: 10.0,
                max_iter: 20,
            }),
        })
        .test(combined(max_load, 60))
        .try_build(domain)
        .unwrap_or_else(|e| panic!("invalid model: {e:?}"));

    let mut last_accepted = 0.0;
    let failure = loop {
        match analysis.step() {
            Ok(result) => {
                last_accepted = result.load_factor;
                if last_accepted >= 1.0 - 1e-12 {
                    panic!(
                        "load control reached the full 50 kip/point load ({last_accepted}); \
                         it should stall below the flexural limit point at {refined_peak}"
                    );
                }
            }
            Err(AnalysisError::LoadControlExhausted {
                accepted_load_factor,
                ..
            }) => break accepted_load_factor,
            Err(other) => panic!("unexpected error: {other:?}"),
        }
    };
    assert!((failure - last_accepted).abs() < 1e-12);
    assert!(
        failure <= refined_peak,
        "load control's accepted fraction {failure} must sit at or below the refined peak {refined_peak}"
    );
    // A loose sanity band around TS's historical ~97.5% (plan §2.5 point
    // 4's figure is a solver outcome, not a target to hit exactly).
    assert!(
        (0.85..1.0).contains(&failure),
        "load control stalled too far from TS's historical ~97.5% to be a comparable outcome: {failure}"
    );
    eprintln!(
        "span beam, 50 kip/point: LoadControlCutbacks stalls at lambda = {failure} (refined peak {refined_peak}, TS historical ~0.975)"
    );
}
