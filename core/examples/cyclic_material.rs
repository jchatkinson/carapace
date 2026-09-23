//! Runs one material through a reverse-cyclic `DisplacementControl`
//! protocol on a single-DOF `ZeroLength` spring and writes the recorded
//! (strain, stress) trace to CSV — the carapace side of the
//! carapace-vs-openseespy material comparison harness in `comparison/`.
//!
//! Usage: `cargo run -p carapace-core --example cyclic_material -- <material> <out.csv>`
//!
//! Material param sets here must match `comparison/materials.py` exactly —
//! that's what makes the two engines' output comparable.

use std::env;
use std::fs;
use std::io::Write;

use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest};
use carapace_core::model::{Domain, Element, Material, Node, Pinching4DmgCyc, ZeroLength};
use carapace_core::testkit::run_cyclic_protocol;

/// Parses `comparison/protocol.txt`'s "amplitude cycles" lines and its
/// `step` line — see that file's own doc comment for the format.
fn load_protocol(path: &str) -> (Vec<(f64, usize)>, f64) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    let mut step = None;
    let mut tiers = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts[0] == "step" {
            step = Some(parts[1].parse::<f64>().unwrap());
        } else {
            let amplitude: f64 = parts[0].parse().unwrap();
            let cycles: usize = parts[1].parse().unwrap();
            tiers.push((amplitude, cycles));
        }
    }
    (tiers, step.expect("protocol.txt must set `step`"))
}

/// Material param sets — keep in lockstep with `comparison/materials.py`.
fn material(name: &str) -> Material {
    match name {
        "elastic" => Material::Elastic { e: 200.0 },
        "elastic_pp" => Material::elastic_pp(200.0, 0.01),
        "steel01" => Material::steel01(2.0, 200.0, 0.02, 0.0, 1.0, 0.0, 1.0),
        "concrete01" => Material::concrete01(4.0, 0.002, 3.4, 0.006),
        "steel02" => Material::steel02(2.0, 200.0, 0.02, 18.0, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0),
        "concrete02" => Material::concrete02(4.0, 0.002, 3.4, 0.006, 0.1, 0.4, 200.0),
        "hysteretic" => Material::hysteretic(
            1.5, 0.005, 2.0, 0.02, 2.2, 0.04, //
            -1.5, -0.005, -2.0, -0.02, -2.2, -0.04, //
            0.5, 0.3, 0.0, 0.0, 0.0,
        ),
        "pinching4" => Material::pinching4(
            1.0, 0.005, 1.8, 0.015, 2.0, 0.025, 1.0, 0.04, //
            -1.0, -0.005, -1.8, -0.015, -2.0, -0.025, -1.0, -0.04, //
            0.5, 0.25, 0.05, 0.5, 0.25, 0.05, //
            [0.0, 0.0, 0.0, 0.0], 0.0, //
            [0.0, 0.0, 0.0, 0.0], 0.0, //
            [0.0, 0.0, 0.0, 0.0], 0.0, //
            10.0, Pinching4DmgCyc::EnergyBased,
        ),
        "parallel" => Material::parallel(vec![
            Material::steel01(1.0, 150.0, 0.02, 0.0, 1.0, 0.0, 1.0),
            Material::Elastic { e: 50.0 },
        ]),
        "series" => Material::series(vec![
            Material::Elastic { e: 300.0 },
            Material::steel01(1.0, 100.0, 0.02, 0.0, 1.0, 0.0, 1.0),
        ]),
        "min_max" => Material::min_max(
            Material::steel01(2.0, 200.0, 0.02, 0.0, 1.0, 0.0, 1.0),
            -0.025,
            0.025,
        ),
        other => panic!("unknown material {other:?}"),
    }
}

/// A single bare `ZeroLength` spring is only as stiff as the tested
/// material's *current* tangent — for a perfectly-plastic material fully
/// yielded (`elastic_pp`) or a material with no tensile stiffness at all
/// (`concrete01`, in tension), that tangent is exactly zero, leaving a
/// singular system under `DisplacementControl` (nothing to invert to
/// solve for the correcting load factor). A small parallel linear spring
/// — `SAFETY_FRACTION` of the material's own initial tangent, so it's
/// negligible next to the material's real response everywhere else — keeps
/// the combined system's tangent bounded away from zero. Must match
/// `comparison/materials.py`'s `SAFETY_FRACTION` and `INITIAL_TANGENT`
/// exactly so both engines record the same (material-only, safety spring
/// subtracted back out) stress trace.
const SAFETY_FRACTION: f64 = 0.005;

fn material_with_safety(name: &str) -> (Material, f64) {
    let inner = material(name);
    let safety_k = SAFETY_FRACTION * inner.initial_tangent();
    (Material::parallel(vec![inner, Material::Elastic { e: safety_k }]), safety_k)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let name = args.get(1).expect("usage: cyclic_material <material> <out.csv>");
    let out_path = args.get(2).expect("usage: cyclic_material <material> <out.csv>");

    let (combined_material, safety_k) = material_with_safety(name);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.load_node(node_j, 0, 1.0); // unit reference load: load_factor reads as spring force

    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, combined_material),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(carapace_core::analysis::Integrator::DisplacementControl {
            node: node_j,
            dof: 0,
            increment: 0.0, // placeholder — run_cyclic_protocol calls set_integrator before every step
        })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-10, max_iter: 30 })
        .build(domain);

    let protocol_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../comparison/protocol.txt");
    let (tiers, step) = load_protocol(protocol_path);

    let mut file = fs::File::create(out_path).unwrap_or_else(|e| panic!("failed to create {out_path}: {e}"));
    writeln!(file, "strain,stress").unwrap();
    writeln!(file, "0,0").unwrap(); // starting point, before any leg

    run_cyclic_protocol(&mut analysis, node_j, 0, &tiers, step, |strain, total_stress| {
        let stress = total_stress - safety_k * strain; // subtract out the safety spring's contribution
        writeln!(file, "{strain},{stress}").unwrap();
    });

    eprintln!("wrote {out_path}");
}
