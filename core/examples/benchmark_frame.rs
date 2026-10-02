//! Multi-story multi-bay structural frame benchmark for Carapace.
//!
//! Evaluates model assembly time, Newton-Raphson iteration solve time,
//! and tracks roof displacement and base shear across load steps.
//!
//! Usage:
//!   cargo run --release -p carapace-core --example benchmark_frame -- [OPTIONS]
//!
//! Options:
//!   --stories <usize>   Number of stories (default: 10)
//!   --bays <usize>      Number of bays (default: 3)
//!   --steps <usize>     Number of load steps (default: 50)
//!   --model <str>       Model type: 'elastic' or 'fiber' (default: 'elastic')
//!   --json              Output machine-readable JSON

use std::collections::HashMap;
use std::env;
use std::time::Instant;

use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, ElasticBeamColumn, Element, Fiber, GeomTransf, Material, Node,
};

struct BenchmarkConfig {
    stories: usize,
    bays: usize,
    steps: usize,
    model: String,
    json: bool,
}

fn parse_args() -> BenchmarkConfig {
    let args: Vec<String> = env::args().collect();
    let mut config = BenchmarkConfig {
        stories: 10,
        bays: 3,
        steps: 50,
        model: "elastic".to_string(),
        json: false,
    };

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--stories" => {
                if i + 1 < args.len() {
                    config.stories = args[i + 1].parse().unwrap_or(10);
                    i += 1;
                }
            }
            "--bays" => {
                if i + 1 < args.len() {
                    config.bays = args[i + 1].parse().unwrap_or(3);
                    i += 1;
                }
            }
            "--steps" => {
                if i + 1 < args.len() {
                    config.steps = args[i + 1].parse().unwrap_or(50);
                    i += 1;
                }
            }
            "--model" => {
                if i + 1 < args.len() {
                    config.model = args[i + 1].clone();
                    i += 1;
                }
            }
            "--json" => {
                config.json = true;
            }
            _ => {}
        }
        i += 1;
    }
    config
}

fn main() {
    let config = parse_args();
    assert!(config.stories >= 1, "stories must be >= 1");
    assert!(config.bays >= 1, "bays must be >= 1");
    assert!(config.steps >= 1, "steps must be >= 1");

    let h_story: f64 = 3000.0; // mm
    let l_bay: f64 = 6000.0;   // mm

    let col_area: f64 = 8000.0; // mm^2
    let col_iz: f64 = 2.0e8;    // mm^4
    let beam_area: f64 = 6000.0; // mm^2
    let beam_iz: f64 = 3.0e8;   // mm^4
    let e_modulus: f64 = 200000.0; // MPa

    let p_base: f64 = if config.model == "fiber" { 25000.0 } else { 10000.0 }; // N

    let setup_start = Instant::now();
    let mut domain = Domain::new();
    let mut nodes = HashMap::new();

    // 1. Create nodes: (stories + 1) rows, (bays + 1) columns
    for s in 0..=config.stories {
        for b in 0..=config.bays {
            let x = b as f64 * l_bay;
            let y = s as f64 * h_story;
            let mut node = Node::new([x, y]);
            if s == 0 {
                // Fixed base
                node = node.fix(0).fix(1).fix(2);
            }
            let id = domain.add_node(node);
            nodes.insert((s, b), id);
        }
    }

    // Fiber definitions if needed
    let col_h = (col_iz / col_area).sqrt();
    let col_fibers = vec![
        Fiber::new(col_h, col_area / 2.0, Material::steel01(355.0, e_modulus, 0.02, 0.0, 1.0, 0.0, 1.0)),
        Fiber::new(-col_h, col_area / 2.0, Material::steel01(355.0, e_modulus, 0.02, 0.0, 1.0, 0.0, 1.0)),
    ];

    let beam_h = (beam_iz / beam_area).sqrt();
    let beam_fibers = vec![
        Fiber::new(beam_h, beam_area / 2.0, Material::steel01(355.0, e_modulus, 0.02, 0.0, 1.0, 0.0, 1.0)),
        Fiber::new(-beam_h, beam_area / 2.0, Material::steel01(355.0, e_modulus, 0.02, 0.0, 1.0, 0.0, 1.0)),
    ];

    let mut element_count = 0;

    // 2. Add column elements: stories * (bays + 1)
    for s in 0..config.stories {
        for b in 0..=config.bays {
            let n1 = nodes[&(s, b)];
            let n2 = nodes[&((s + 1), b)];
            if config.model == "fiber" {
                domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
                    n1,
                    n2,
                    col_fibers.clone(),
                    BeamIntegration::Legendre { points: 3 },
                )));
            } else {
                domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                    n1,
                    n2,
                    e_modulus,
                    col_area,
                    col_iz,
                    GeomTransf::Linear,
                )));
            }
            element_count += 1;
        }
    }

    // 3. Add beam elements: stories * bays
    for s in 1..=config.stories {
        for b in 0..config.bays {
            let n1 = nodes[&(s, b)];
            let n2 = nodes[&(s, b + 1)];
            if config.model == "fiber" {
                domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
                    n1,
                    n2,
                    beam_fibers.clone(),
                    BeamIntegration::Legendre { points: 3 },
                )));
            } else {
                domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                    n1,
                    n2,
                    e_modulus,
                    beam_area,
                    beam_iz,
                    GeomTransf::Linear,
                )));
            }
            element_count += 1;
        }
    }

    // 4. Apply lateral floor loads (inverted triangular distribution)
    for s in 1..=config.stories {
        let load_mag = p_base * s as f64;
        let left_node = nodes[&(s, 0)];
        domain.load_node(left_node, 0, load_mag);
    }

    let increment = 1.0 / config.steps as f64;
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-6,
            max_iter: 30,
        })
        .build(domain);

    let setup_elapsed = setup_start.elapsed();

    // 5. Solve loop
    let roof_node = nodes[&(config.stories, 0)];
    let mut roof_history = Vec::with_capacity(config.steps);

    let solve_start = Instant::now();
    for _ in 0..config.steps {
        analysis.step().expect("step should converge successfully");
        let disp = analysis.domain().node(roof_node).displacement[0];
        roof_history.push(disp);
    }
    let solve_elapsed = solve_start.elapsed();

    let total_elapsed = setup_elapsed + solve_elapsed;
    let total_nodes = (config.stories + 1) * (config.bays + 1);
    let free_dofs = config.stories * (config.bays + 1) * 3;
    let final_roof_disp = *roof_history.last().unwrap_or(&0.0);

    let setup_ms = setup_elapsed.as_secs_f64() * 1000.0;
    let solve_ms = solve_elapsed.as_secs_f64() * 1000.0;
    let total_ms = total_elapsed.as_secs_f64() * 1000.0;

    if config.json {
        let history_str: Vec<String> = roof_history.iter().map(|d| format!("{:.10}", d)).collect();
        println!(
            "{{\n  \"engine\": \"carapace\",\n  \"model\": \"{}\",\n  \"stories\": {},\n  \"bays\": {},\n  \"nodes\": {},\n  \"elements\": {},\n  \"free_dofs\": {},\n  \"steps\": {},\n  \"setup_ms\": {:.4},\n  \"solve_ms\": {:.4},\n  \"total_ms\": {:.4},\n  \"final_roof_disp\": {:.10},\n  \"roof_disp_history\": [{}]\n}}",
            config.model,
            config.stories,
            config.bays,
            total_nodes,
            element_count,
            free_dofs,
            config.steps,
            setup_ms,
            solve_ms,
            total_ms,
            final_roof_disp,
            history_str.join(", ")
        );
    } else {
        println!("============================================================");
        println!("  Carapace Frame Benchmark");
        println!("============================================================");
        println!("  Model type:      {}", config.model);
        println!("  Grid dimensions: {} stories x {} bays", config.stories, config.bays);
        println!("  Mesh statistics: {} nodes, {} elements, {} free DOFs", total_nodes, element_count, free_dofs);
        println!("  Load steps:      {}", config.steps);
        println!("------------------------------------------------------------");
        println!("  Setup time:      {:.3} ms", setup_ms);
        println!("  Solve time:      {:.3} ms ({:.3} ms/step)", solve_ms, solve_ms / config.steps as f64);
        println!("  Total time:      {:.3} ms", total_ms);
        println!("  Final roof disp: {:.6} mm", final_roof_disp);
        println!("============================================================");
    }
}
