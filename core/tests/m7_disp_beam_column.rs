use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{
    BeamIntegration, Domain, ElasticBeamColumn, Element, Fiber, GeomTransf, Material, Node,
};

/// M7 stage 1 acceptance (implementation-plan §6): `DispBeamColumn`
/// (fiber-discretized, displacement-based) must exactly reproduce
/// `ElasticBeamColumn`'s closed-form linear-elastic response — not
/// approximately, exactly, to solver precision. Two symmetric elastic
/// fibers at `y = ±sqrt(Iz/A)`, area `A/2` each, reproduce `EA = E*A` and
/// `EI = E*Iz` exactly (`FiberSection`'s
/// `two_symmetric_elastic_fibers_reproduce_ea_ei_exactly` unit test checks
/// this piece in isolation). The element-level strain-displacement field
/// (`DispBeamColumn::strain_displacement`) uses the same cubic-Hermite/
/// linear-axial shape functions as `ElasticBeamColumn`, so for constant
/// `EA`/`EI` any >=2-point Gauss rule integrates the (at most quadratic)
/// integrand exactly — the two elements' stiffness matrices must therefore
/// be identical, not just close.
#[test]
fn two_fiber_elastic_section_matches_elastic_beam_column_exactly() {
    let (e, area, iz, length): (f64, f64, f64, f64) = (30000.0, 2.0, 1000.0, 100.0);
    let tip_load = -10.0;

    let disp_beam_tip = {
        let h = (iz / area).sqrt();
        let fibers = vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ];

        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]));
        domain.load_node(node_j, 1, tip_load);
        domain.add_element(Element::DispBeamColumn(carapace_core::model::DispBeamColumn::new(
            node_i,
            node_j,
            fibers,
            BeamIntegration::Legendre { points: 3 },
        )));

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever DispBeamColumn should solve");
        analysis.domain().node(node_j).displacement[1]
    };

    let elastic_beam_tip = {
        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]));
        domain.load_node(node_j, 1, tip_load);
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            node_i,
            node_j,
            e,
            area,
            iz,
            GeomTransf::Linear,
        )));

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever ElasticBeamColumn should solve");
        analysis.domain().node(node_j).displacement[1]
    };

    assert!(
        (disp_beam_tip - elastic_beam_tip).abs() < 1e-9,
        "DispBeamColumn={disp_beam_tip}, ElasticBeamColumn={elastic_beam_tip}, should match exactly"
    );

    // Also confirm against the closed-form cantilever tip deflection
    // directly, so this isn't just "the two happen to agree."
    let expected = tip_load * length.powi(3) / (3.0 * e * iz);
    assert!((disp_beam_tip - expected).abs() < 1e-9, "expected {expected}, got {disp_beam_tip}");
}

/// Sanity check that the exact match above isn't an artifact of a
/// specific point count or quadrature family — Lobatto (which includes
/// the endpoints, unlike Legendre) must give the identical exact answer
/// too, since the underlying integrand's polynomial degree doesn't care
/// which family of points samples it.
#[test]
fn lobatto_integration_gives_the_same_exact_result_as_legendre() {
    let (e, area, iz, length): (f64, f64, f64, f64) = (30000.0, 2.0, 1000.0, 100.0);
    let tip_load = -10.0;

    let run = |integration: BeamIntegration| {
        let h = (iz / area).sqrt();
        let fibers = vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ];
        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]));
        domain.load_node(node_j, 1, tip_load);
        domain.add_element(Element::DispBeamColumn(carapace_core::model::DispBeamColumn::new(
            node_i, node_j, fibers, integration,
        )));
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().unwrap();
        analysis.domain().node(node_j).displacement[1]
    };

    let legendre = run(BeamIntegration::Legendre { points: 2 });
    let lobatto = run(BeamIntegration::Lobatto { points: 4 });
    assert!((legendre - lobatto).abs() < 1e-9, "Legendre={legendre}, Lobatto={lobatto}");
}
