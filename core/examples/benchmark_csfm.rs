//! Times a CSFM model on native Carapace (plan `docs/csfm-plan.md` phase 6).
//!
//! Reads a model in the golden fixture format (`core/tests/csfm/fixtures/models/*.json`: `model`,
//! `schedule`, `solver`), solves its load path with the fixture's solver settings and prints one JSON
//! object with the timings and the increment and iteration counts. Only `Analysis::step` is timed
//! (assembly, factorization and solves); model building is excluded, and one untimed warm-up run
//! precedes the timed runs. If the file has a `ts.finalDisplacements` array, the difference of the
//! final displacements from it is reported too.
//!
//! Usage:
//!   cargo run --release -p carapace-core --example benchmark_csfm -- <model.json> [reps]

use std::time::Instant;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use serde_json::{json, Value};

#[allow(dead_code)]
#[path = "../tests/csfm/model.rs"]
mod model;

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: benchmark_csfm <model.json> [reps]");
    let reps: usize = args.next().map_or(5, |r| r.parse().expect("reps"));
    let fixture: Value = serde_json::from_str(
        &std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}")),
    )
    .expect("valid JSON");
    let model = &fixture["model"];
    let schedule: Vec<f64> = fixture["schedule"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    let max_load = model["loads"]
        .as_array()
        .unwrap()
        .iter()
        .map(|l| l["value"].as_f64().unwrap().abs())
        .fold(0.0_f64, f64::max);
    let tolerance = (1e-10 * max_load).max(1e-10);
    let max_iter = fixture["solver"]["maxIterationsPerIncrement"]
        .as_u64()
        .unwrap() as usize;

    let mut times = Vec::new();
    let mut iterations = 0;
    let mut factorizations = 0;
    let mut displacements = Vec::new();
    // One untimed warm-up run, then `reps` timed runs.
    for rep in 0..=reps {
        let built = model::build(model);
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl {
                increment: schedule[0],
            })
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::Combined {
                force_tol: tolerance,
                moment_tol: tolerance,
                relative_tol: 1e-10,
                displacement_tol: None,
                max_iter,
            })
            .try_build(built.domain)
            .expect("valid model");
        let start = Instant::now();
        let (mut run_iterations, mut run_factorizations) = (0, 0);
        for &nominal in &schedule {
            analysis.set_integrator(Integrator::LoadControl { increment: nominal });
            let result = analysis.step().expect("the benchmark model converges");
            run_iterations += result.iterations;
            run_factorizations += result.factorizations;
        }
        let elapsed = start.elapsed().as_secs_f64() * 1e3;
        if rep > 0 {
            times.push(elapsed);
        }
        iterations = run_iterations;
        factorizations = run_factorizations;
        displacements = built
            .nodes
            .iter()
            .flat_map(|&node| analysis.domain().node(node).displacement[..2].to_vec())
            .collect::<Vec<f64>>();
    }

    let difference = fixture["ts"]["finalDisplacements"]
        .as_array()
        .map(|reference| {
            let reference: Vec<f64> = reference.iter().map(|v| v.as_f64().unwrap()).collect();
            let scale = reference.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            displacements
                .iter()
                .zip(&reference)
                .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs() / scale))
        });
    let mut sorted = times.clone();
    sorted.sort_by(f64::total_cmp);
    println!(
        "{}",
        json!({
            "engine": "native", "reps": reps, "timesMs": times, "medianMs": sorted[sorted.len() / 2],
            "dof": displacements.len(), "increments": schedule.len(),
            "iterations": iterations, "factorizations": factorizations,
            "maxDisplacementDifferenceFromReference": difference,
        })
    );
}
