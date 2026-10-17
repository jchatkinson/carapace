//! Nonlinear frame scaling benchmark: the same model as `comparison/opensees_scaling.py`.
//!
//! A fixed-base 2D moment frame of `dispBeamColumn` members (Legendre points, many-fibre
//! rectangular Steel02 sections) under inverted-triangle lateral loads on the left column line.
//! Load control scales `--load-base * story` N per story over `--steps`; displacement control
//! drives the roof x-displacement by `--displacement-increment` per step.
//!
//! Usage: cargo run --release -p carapace-core --example benchmark_scaling -- \
//!   --stories 20 --bays 4 [--steps 50] [--fibers 20] [--points 5] [--load-base 1.0] \
//!   [--displacement-increment 30.0] [--json]

use std::collections::HashMap;
use std::env;
use std::time::Instant;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, Element, Fiber, Material, Node,
};

const STORY_HEIGHT: f64 = 3000.0;
const BAY_WIDTH: f64 = 6000.0;
const E: f64 = 200_000.0;
const COLUMN: (f64, f64) = (8000.0, 2.0e8);
const BEAM: (f64, f64) = (6000.0, 3.0e8);

/// Rectangular section of `n` equal-area Steel02 fibres with the given area and inertia.
fn section(area: f64, inertia: f64, n: usize) -> Vec<Fiber> {
    let depth = (12.0 * inertia / area).sqrt();
    (0..n)
        .map(|i| {
            let y = -depth / 2.0 + (i as f64 + 0.5) * depth / n as f64;
            Fiber::new(
                y,
                area / n as f64,
                Material::steel02(355.0, E, 0.02, 18.5, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0),
            )
        })
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let get = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .map(|i| args[i + 1].clone())
    };
    let num = |name: &str, default: f64| {
        get(name).map_or(default, |v| {
            v.parse()
                .unwrap_or_else(|_| panic!("{name} must be a number"))
        })
    };
    let stories = num("--stories", 10.0) as usize;
    let bays = num("--bays", 3.0) as usize;
    let steps = num("--steps", 50.0) as usize;
    let fibers = num("--fibers", 20.0) as usize;
    let points = num("--points", 5.0) as usize;
    let load_base = num("--load-base", 1.0);
    let displacement_increment = get("--displacement-increment")
        .map(|v| v.parse::<f64>().expect("increment must be a number"));
    let json = args.iter().any(|a| a == "--json");
    assert!(stories >= 1 && bays >= 1 && steps >= 1 && fibers >= 2 && points >= 1);

    let setup_start = Instant::now();
    let mut domain = Domain::new();
    let mut nodes = HashMap::new();
    for s in 0..=stories {
        for b in 0..=bays {
            let mut node = Node::new([b as f64 * BAY_WIDTH, s as f64 * STORY_HEIGHT]);
            if s == 0 {
                node = node.fix(0).fix(1).fix(2);
            }
            nodes.insert((s, b), domain.add_node(node));
        }
    }
    let (columns, beams) = (
        section(COLUMN.0, COLUMN.1, fibers),
        section(BEAM.0, BEAM.1, fibers),
    );
    let mut element_count = 0;
    let mut add = |domain: &mut Domain, i, j, fibers: &Vec<Fiber>| {
        domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
            i,
            j,
            fibers.clone(),
            BeamIntegration::Legendre { points },
        )));
        element_count += 1;
    };
    for s in 0..stories {
        for b in 0..=bays {
            add(&mut domain, nodes[&(s, b)], nodes[&(s + 1, b)], &columns);
        }
    }
    for s in 1..=stories {
        for b in 0..bays {
            add(&mut domain, nodes[&(s, b)], nodes[&(s, b + 1)], &beams);
        }
    }
    for s in 1..=stories {
        domain.load_node(nodes[&(s, 0)], 0, load_base * s as f64);
    }

    let roof = nodes[&(stories, 0)];
    let integrator = match displacement_increment {
        Some(increment) => Integrator::DisplacementControl {
            node: roof,
            dof: 0,
            increment,
        },
        None => Integrator::LoadControl {
            increment: 1.0 / steps as f64,
        },
    };
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(integrator)
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-6,
            max_iter: 30,
        })
        .build(domain);
    let setup_ms = setup_start.elapsed().as_secs_f64() * 1000.0;

    let (mut roof_history, mut load_history, mut iterations) =
        (Vec::with_capacity(steps), Vec::with_capacity(steps), 0);
    let solve_start = Instant::now();
    for _ in 0..steps {
        let result = analysis.step().expect("step should converge");
        load_history.push(result.load_factor);
        iterations += result.iterations;
        roof_history.push(analysis.domain().node(roof).displacement[0]);
    }
    let solve_ms = solve_start.elapsed().as_secs_f64() * 1000.0;

    let list = |v: &[f64]| {
        v.iter()
            .map(|x| format!("{x:.10}"))
            .collect::<Vec<_>>()
            .join(", ")
    };
    let free_dofs = stories * (bays + 1) * 3;
    if json {
        println!(
            "{{\"engine\": \"carapace\", \"stories\": {stories}, \"bays\": {bays}, \"elements\": {element_count}, \"free_dofs\": {free_dofs}, \"steps\": {steps}, \"setup_ms\": {setup_ms:.4}, \"solve_ms\": {solve_ms:.4}, \"iterations\": {iterations}, \"final_roof_disp\": {:.10}, \"roof_disp_history\": [{}], \"load_factor_history\": [{}]}}",
            roof_history[steps - 1], list(&roof_history), list(&load_history),
        );
    } else {
        println!("{stories}x{bays} frame: {element_count} elements, {free_dofs} free DOFs, {steps} steps");
        println!("setup {setup_ms:.1} ms, solve {solve_ms:.1} ms, {iterations} iterations, final roof disp {:.4} mm", roof_history[steps - 1]);
    }
}
