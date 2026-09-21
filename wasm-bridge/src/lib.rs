use wasm_bindgen::prelude::*;

use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Material, Node, Truss, ZeroLength};

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

    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2).with_load(0, load));

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

/// M2 wiring check: a `ZeroLength` + `Material::Ent` ("no tension")
/// connector under a compressive load, computed through the same real
/// architecture — proves the `Element`/`Material` enum dispatch generalizes
/// beyond `Truss`/`Elastic` on wasm32 + Node, not just natively (see
/// `core/tests/m2_zero_length.rs` for the native-side equivalent and the
/// single-linear-regime caveat this test shares with it).
#[wasm_bindgen]
pub fn zero_length_ent_displacement(load: f64, modulus: f64) -> f64 {
    let mut domain = Domain::new();

    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2).with_load(0, load));

    domain.add_element(Element::ZeroLength(
        ZeroLength::new(node_i, node_j).with_material(0, Material::Ent { e: modulus }),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    analysis.step().expect("compressive Ent response should solve");
    analysis.domain().node(node_j).displacement[0]
}

/// M3 wiring check: a simply-supported `ElasticBeamColumn` under a uniform
/// transverse element load, returning the node_i end rotation — closed form
/// `theta = w*L^3/(24*E*I)`. See `core/tests/m3_beam.rs` for the native
/// equivalent and why a single element is exact here.
#[wasm_bindgen]
pub fn simply_supported_beam_end_rotation(e: f64, iz: f64, area: f64, length: f64, w: f64) -> f64 {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1));

    domain.add_element(Element::ElasticBeamColumn(
        ElasticBeamColumn::new(node_i, node_j, e, area, iz, GeomTransf::Linear).with_uniform_load(w),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);

    analysis.step().expect("simply supported beam under UDL should solve");
    analysis.domain().node(node_i).displacement[2]
}
