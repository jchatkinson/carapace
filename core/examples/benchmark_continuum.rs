//! Continuum benchmark: an N x N Quad4 panel (clamped left edge, shear on the right edge).
//! Reports setup, one full static step, and the phases of one assembly + solve.
//! Run with `cargo run --release -p carapace-core --example benchmark_continuum -- [N] [reps] [enhanced]`.

use std::time::Instant;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Node, PlaneMaterial, Quad4, Quad4Formulation};

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let n: usize = args
        .first()
        .map_or(70, |a| a.parse().expect("N must be an integer"));
    let reps: usize = args
        .get(1)
        .map_or(5, |a| a.parse().expect("reps must be an integer"));
    let enhanced = args.get(2).is_some_and(|a| a == "enhanced");

    let setup = Instant::now();
    let mut domain = Domain::new();
    let mut grid = Vec::new();
    for j in 0..=n {
        let mut row = Vec::new();
        for i in 0..=n {
            let mut node = Node::new([i as f64, j as f64]);
            if i == 0 {
                node = node.fix(0).fix(1);
            }
            row.push(domain.add_node(node));
        }
        grid.push(row);
    }
    let material = PlaneMaterial::plane_stress(200_000.0, 0.3).unwrap();
    let formulation = if enhanced {
        Quad4Formulation::Enhanced
    } else {
        Quad4Formulation::Full
    };
    for j in 0..n {
        for i in 0..n {
            let nodes = [
                grid[j][i],
                grid[j][i + 1],
                grid[j + 1][i + 1],
                grid[j + 1][i],
            ];
            domain.add_element(Element::Quad4(
                Quad4::new(nodes, 1.0, material.clone()).with_formulation(formulation),
            ));
        }
    }
    for row in &grid {
        domain.load_node(row[n], 1, -1.0);
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    let setup_ms = setup.elapsed().as_secs_f64() * 1e3;

    let step = Instant::now();
    let result = analysis.step().expect("panel converges");
    let step_ms = step.elapsed().as_secs_f64() * 1e3;
    let tip = analysis.domain().node(grid[n][n]).displacement[1];

    let mut phases = vec![Vec::new(); 5];
    let mut triplets = 0;
    for _ in 0..reps {
        let (times, count) = analysis.domain().bench_assembly(1.0);
        triplets = count;
        for (phase, t) in phases.iter_mut().zip(times) {
            phase.push(t * 1e3);
        }
    }
    let m: Vec<f64> = phases.into_iter().map(median).collect();
    println!(
        "{n}x{n} Quad4 ({}), {} DOF, {triplets} triplets, {} Newton iterations, tip uy {tip:.6e}",
        if enhanced { "enhanced" } else { "full" },
        2 * n * (n + 1),
        result.iterations
    );
    println!("setup {setup_ms:.1} ms, full step {step_ms:.1} ms");
    println!("per assembly+solve (median of {reps}): elements {:.1} ms, matrix {:.1} ms, first factor {:.1} ms, repeat factor {:.1} ms, solve {:.1} ms", m[0], m[1], m[2], m[3], m[4]);
}
