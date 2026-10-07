//! NLTHA benchmark: multi-story RC frame, `DispBeamColumn` + RC fiber sections
//! (Concrete02 cover/core + Steel02 rebar), Newmark + Rayleigh + uniform ground
//! motion, Newton corrector.
//!
//! cargo run --release -p carapace-core --example benchmark_nltha -- [stories bays steps pga_g points tangent]
//!   tangent: current | reuse | initial   (default current)
//! Units: N, mm, s, tonne.

use std::time::Instant;

use carapace_core::analysis::{
    modal_analysis, Algorithm, ConvergenceTest, GroundMotion, RayleighDamping, TangentStrategy,
    TransientAnalysis,
};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, Element, Fiber, LoadSeries, Material, Node,
};

fn conc(fc: f64, e0: f64, fcu: f64, ecu: f64) -> Material {
    match std::env::var("CONC").as_deref() {
        Ok("el") => Material::Elastic { e: 30000.0 },
        Ok("01") => Material::concrete01(fc, e0, fcu, ecu),
        _ => Material::concrete02(fc, e0, fcu, ecu, 0.1, 2.5, std::env::var("ETS").ok().and_then(|s| s.parse().ok()).unwrap_or(3000.0)),
    }
}

fn rc_fibers(b: f64, h: f64, cover: f64, as_layer: f64, n_core: usize) -> Vec<Fiber> {
    let cover_mat = || conc(30.0, 0.002, 6.0, 0.006);
    let core_mat = || conc(38.0, 0.004, 15.0, 0.014);
    let steel = || if std::env::var("STEEL").as_deref() == Ok("01") { Material::steel01(420.0, 200000.0, 0.01, 0.0, 1.0, 0.0, 1.0) } else if std::env::var("STEEL").as_deref() == Ok("el") { Material::Elastic { e: 200000.0 } } else { Material::steel02(420.0, 200000.0, 0.01, 18.5, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0) };
    let mut f = Vec::new();
    // cover layers (top & bottom), 2 sub-layers each
    for s in 0..2 {
        let t = cover / 2.0;
        let y = h / 2.0 - t * (s as f64 + 0.5);
        f.push(Fiber::new(y, b * t, cover_mat()));
        f.push(Fiber::new(-y, b * t, cover_mat()));
    }
    // core layers
    let hc = h - 2.0 * cover;
    let dy = hc / n_core as f64;
    for i in 0..n_core {
        let y = hc / 2.0 - dy * (i as f64 + 0.5);
        f.push(Fiber::new(y, (b - 2.0 * cover) * dy, core_mat()));
    }
    // rebar: top/bottom + two intermediate layers
    let ys = h / 2.0 - cover;
    for y in [ys, -ys] {
        f.push(Fiber::new(y, as_layer, steel()));
    }
    for y in [ys / 3.0, -ys / 3.0] {
        f.push(Fiber::new(y, as_layer / 2.0, steel()));
    }
    f
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let g = |i: usize, d: f64| a.get(i).map(|s| s.parse::<f64>().unwrap()).unwrap_or(d);
    let (stories, bays, steps, pga_g, npts) = (g(0, 6.0) as usize, g(1, 3.0) as usize, g(2, 1000.0) as usize, g(3, 0.6), g(4, 5.0) as usize);
    let tangent = match a.get(5).map(String::as_str) {
        Some("reuse") => TangentStrategy::ReuseAtStepStart,
        Some("initial") => TangentStrategy::Initial,
        _ => TangentStrategy::Current,
    };
    let dt = 0.01;

    let (h, l) = (3500.0, 6000.0);
    let col = rc_fibers(500.0, 500.0, 40.0, 1900.0, 18);
    let beam = rc_fibers(300.0, 600.0, 40.0, 1200.0, 18);
    let mut domain = Domain::new();
    let mut nodes = Vec::new();
    for s in 0..=stories {
        let mut row = Vec::new();
        for b in 0..=bays {
            let mut n = Node::new([b as f64 * l, s as f64 * h]);
            n = if s == 0 { n.fix(0).fix(1).fix(2) } else { n.with_mass(0, 25.0).with_mass(1, 25.0).with_mass(2, 1.0e6) };
            row.push(domain.add_node(n));
        }
        nodes.push(row);
    }
    let integ = BeamIntegration::Lobatto { points: npts };
    for s in 1..=stories {
        for b in 0..=bays {
            domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(nodes[s - 1][b], nodes[s][b], col.clone(), integ)));
        }
        for b in 0..bays {
            domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(nodes[s][b], nodes[s][b + 1], beam.clone(), integ)));
        }
    }
    let nel = stories * (bays + 1) + stories * bays;
    let ndof = stories * (bays + 1) * 3;

    // Rayleigh from first two modes, 5%
    let modes = modal_analysis(&mut domain, 2).expect("modal");
    let (w1, w2) = (modes[0].frequency, modes[1].frequency);
    let _ = (w1, w2);
    // `frequency` is Hz
    let (w1, w2) = (2.0 * std::f64::consts::PI * w1, 2.0 * std::f64::consts::PI * w2);
    let xi = 0.05;
    let alpha = 2.0 * xi * w1 * w2 / (w1 + w2);
    let beta = 2.0 * xi / (w1 + w2);

    // synthetic accelerogram: modulated multi-sine
    let n_gm = steps + 1;
    let times: Vec<f64> = (0..n_gm).map(|i| i as f64 * dt).collect();
    let tot = steps as f64 * dt;
    let factors: Vec<f64> = times
        .iter()
        .map(|&t| {
            let env = (t / (0.15 * tot)).min(1.0) * (-(t - 0.5 * tot).max(0.0) / (0.4 * tot)).exp();
            let s = (2.0 * 3.14159 * 0.7 * t).sin() + 0.8 * (2.0 * 3.14159 * 1.3 * t + 1.0).sin() + 0.6 * (2.0 * 3.14159 * 2.9 * t + 2.0).sin() + 0.4 * (2.0 * 3.14159 * 5.3 * t).sin();
            9810.0 * pga_g * env * s / 1.8
        })
        .collect();

    let mut ta = TransientAnalysis::new(domain, RayleighDamping::new(alpha, beta), dt)
        .expect("transient")
        .with_ground_motion(GroundMotion::new(0, LoadSeries::Path { times, factors }))
        .with_algorithm(
            Algorithm::Newton { tangent, line_search: None },
            match std::env::var("TEST").as_deref() {
                Ok("du") => ConvergenceTest::NormDispIncr { tol: std::env::var("TOL").ok().and_then(|s| s.parse().ok()).unwrap_or(1e-5), max_iter: 25 },
                Ok("energy") => ConvergenceTest::EnergyIncr { tol: 1e-6, max_iter: 50 },
                _ => ConvergenceTest::NormUnbalance { tol: 1e-2, max_iter: 50 },
            },
        );
    let t0 = Instant::now();
    let (mut iters, mut facts, mut peak) = (0usize, 0usize, 0.0f64);
    let roof = nodes[stories][0];
    for _ in 0..steps {
        let r = ta.step().expect("step failed");
        iters += r.iterations;
        facts += r.factorizations;
        peak = peak.max(ta.domain().node(roof).displacement[0].abs());
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    println!(
        "{{\"stories\":{stories},\"bays\":{bays},\"elements\":{nel},\"free_dofs\":{ndof},\"fibers_per_section\":{},\"points\":{npts},\"steps\":{steps},\"T1\":{:.3},\"tangent\":\"{:?}\",\"iters\":{iters},\"factorizations\":{facts},\"peak_roof_mm\":{peak:.2},\"total_ms\":{ms:.1},\"ms_per_step\":{:.3}}}",
        col.len(), 1.0 / modes[0].frequency, tangent, ms / steps as f64
    );
}
