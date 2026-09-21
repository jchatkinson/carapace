use wasm_bindgen::prelude::*;

use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, Material, Node, Truss};

#[wasm_bindgen]
pub fn axial_displacement(load: f64, length: f64, area: f64, modulus: f64) -> f64 {
    carapace_core::axial_displacement(load, length, area, modulus)
}

/// Same 2-node truss case as `axial_displacement`, but computed through the
/// real `Domain`/`Element::Truss`/`Analysis` architecture (M1) rather than
/// the closed-form placeholder — see implementation-plan.md M1 acceptance
/// criteria. Kept alongside the M0 export for wasm/Node verification.
#[wasm_bindgen]
pub fn axial_displacement_via_analysis(load: f64, length: f64, area: f64, modulus: f64) -> f64 {
    let mut domain = Domain::new();

    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).with_load(0, load));

    domain.add_element(Element::Truss(Truss::new(
        node_i,
        node_j,
        area,
        Material::Elastic { e: modulus },
    )));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    analysis.step().expect("linear elastic truss should solve");
    analysis.domain().node(node_j).displacement[0]
}
