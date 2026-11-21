//! Large elastic shell benchmark: an N x N `Shell4` plate, simply supported edges, uniform pressure.
//! Times one linear step (assembly + factor + solve) and the phases of one assembly + solve.
//! Run with `cargo run --release -p carapace-core --features parallel --example benchmark_shell -- [N] [reps] [threads]`.
//! `threads = 1` is truly serial (the default of a `parallel` build is serial too).

use std::time::Instant;

use carapace_core::analysis::{
    set_parallelism, Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
    Parallelism,
};
use carapace_core::model::{Domain3, Element3, ElementLoad3, Node3, Shell4, ShellSection};

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |i: usize, d: usize| {
        args.get(i)
            .map_or(d, |a| a.parse().expect("integer argument"))
    };
    let (n, reps, threads) = (arg(0, 150), arg(1, 3), arg(2, 1));
    if threads != 1 {
        assert!(
            set_parallelism(Parallelism::Threads(threads)),
            "threads needs --features parallel"
        );
    }
    let (a, h, e, nu, p) = (10.0, 0.2, 2.0e5, 0.3, 0.01);

    let setup = Instant::now();
    let mut domain = Domain3::new();
    let mut at = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (x, y) = (a * i as f64 / n as f64, a * j as f64 / n as f64);
            let mut node = Node3::new([x, y, 0.0]);
            if i == 0 || j == 0 || i == n || j == n {
                node = node.fix(2);
            }
            if i == 0 && j == 0 {
                node = node.fix(0).fix(1);
            }
            if i == n && j == 0 {
                node = node.fix(1);
            }
            at.push(domain.add_node(node));
        }
    }
    let section = ShellSection::elastic_membrane_plate(e, nu, h, 0.0).unwrap();
    let pattern = domain.default_pattern();
    for j in 0..n {
        for i in 0..n {
            let ids = [
                at[j * (n + 1) + i],
                at[j * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i],
            ];
            let shell = domain.add_element(Element3::Shell4(Shell4::new(ids, section.clone())));
            domain.add_element_load(pattern, shell, ElementLoad3::pressure(p));
        }
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    let setup_ms = setup.elapsed().as_secs_f64() * 1e3;

    let step = Instant::now();
    analysis.step().expect("linear solve");
    let step_ms = step.elapsed().as_secs_f64() * 1e3;
    let centre = analysis
        .domain()
        .node(at[(n / 2) * (n + 1) + n / 2])
        .displacement[2];

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
        "{n}x{n} Shell4, {} elements, {} DOF, {triplets} triplets, threads {threads}, centre w {centre:.12e}",
        n * n,
        6 * (n + 1) * (n + 1)
    );
    println!("setup {setup_ms:.1} ms, full step {step_ms:.1} ms");
    println!("per assembly+solve (median of {reps}): elements {:.1} ms, matrix {:.1} ms, first factor {:.1} ms, repeat factor {:.1} ms, solve {:.1} ms", m[0], m[1], m[2], m[3], m[4]);
}
