//! Linear-elastic Quad4 panel benchmark: the same model as `comparison/opensees_panel.py`.
//!
//! An N x N mesh of unit square plane-stress `Quad4` elements (E = 200000, nu = 0.3, t = 1), clamped on the
//! left edge with a unit downward load on every right-edge node, solved in one Newton load step. The timed
//! region is that step (assembly, factorisation and solve).
//!
//! Usage: cargo run --release -p carapace-core --example benchmark_panel -- --n 80 [--json]

use std::env;
use std::time::Instant;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Node, PlaneMaterial, Quad4};

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args
        .iter()
        .position(|a| a == "--n")
        .map_or(50, |i| args[i + 1].parse().expect("--n must be an integer"));
    let json = args.iter().any(|a| a == "--json");
    assert!(n >= 1);

    let mut domain = Domain::new();
    let mut grid = Vec::new();
    for j in 0..=n {
        grid.push(
            (0..=n)
                .map(|i| {
                    let node = Node::new([i as f64, j as f64]);
                    domain.add_node(if i == 0 { node.fix(0).fix(1) } else { node })
                })
                .collect::<Vec<_>>(),
        );
    }
    let material = PlaneMaterial::plane_stress(200_000.0, 0.3).unwrap();
    for j in 0..n {
        for i in 0..n {
            let nodes = [
                grid[j][i],
                grid[j][i + 1],
                grid[j + 1][i + 1],
                grid[j + 1][i],
            ];
            domain.add_element(Element::Quad4(Quad4::new(nodes, 1.0, material.clone())));
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

    let start = Instant::now();
    let result = analysis.step().expect("panel converges");
    let solve_ms = start.elapsed().as_secs_f64() * 1000.0;

    let domain = analysis.domain();
    let tip = domain.node(grid[n][n]).displacement[1];
    let norm = grid
        .iter()
        .flatten()
        .map(|&id| {
            domain
                .node(id)
                .displacement
                .iter()
                .map(|d| d * d)
                .sum::<f64>()
        })
        .sum::<f64>()
        .sqrt();
    let free_dofs = 2 * n * (n + 1);
    if json {
        println!("{{\"engine\": \"carapace\", \"n\": {n}, \"elements\": {}, \"free_dofs\": {free_dofs}, \"solve_ms\": {solve_ms:.4}, \"iterations\": {}, \"tip_uy\": {tip:.12e}, \"disp_norm\": {norm:.12e}}}", n * n, result.iterations);
    } else {
        println!("{n}x{n} Quad4: {free_dofs} free DOFs, solve {solve_ms:.1} ms, {} iterations, tip uy {tip:.6e}", result.iterations);
    }
}
