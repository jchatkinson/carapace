//! `DispBeamColumn` / `DispBeamColumn3` reproduce the elastic beam exactly for an elastic fiber section.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn3, Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3,
    Element, Element3, Fiber, Fiber3, GeomTransf, GeomTransf3, Material, Node, Node3, SpatialDof,
};

/// `DispBeamColumn`
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
        domain.add_element(Element::DispBeamColumn(
            carapace_core::model::DispBeamColumn::new(
                node_i,
                node_j,
                fibers,
                BeamIntegration::Legendre { points: 3 },
            ),
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
        analysis
            .step()
            .expect("cantilever DispBeamColumn should solve");
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
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);
        analysis
            .step()
            .expect("cantilever ElasticBeamColumn should solve");
        analysis.domain().node(node_j).displacement[1]
    };

    assert!(
        (disp_beam_tip - elastic_beam_tip).abs() < 1e-9,
        "DispBeamColumn={disp_beam_tip}, ElasticBeamColumn={elastic_beam_tip}, should match exactly"
    );

    // Also confirm against the closed-form cantilever tip deflection
    // directly, so this isn't just "the two happen to agree."
    let expected = tip_load * length.powi(3) / (3.0 * e * iz);
    assert!(
        (disp_beam_tip - expected).abs() < 1e-9,
        "expected {expected}, got {disp_beam_tip}"
    );
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
        domain.add_element(Element::DispBeamColumn(
            carapace_core::model::DispBeamColumn::new(node_i, node_j, fibers, integration),
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
        analysis.step().unwrap();
        analysis.domain().node(node_j).displacement[1]
    };

    let legendre = run(BeamIntegration::Legendre { points: 2 });
    let lobatto = run(BeamIntegration::Lobatto { points: 4 });
    assert!(
        (legendre - lobatto).abs() < 1e-9,
        "Legendre={legendre}, Lobatto={lobatto}"
    );
}

/// `DispBeamColumn3` (biaxial
/// fiber-discretized, displacement-based) must exactly reproduce
/// `ElasticBeamColumn3`'s closed-form linear-elastic response in *both*
/// bending planes plus axial and (decoupled) torsion — not approximately,
/// exactly, to solver precision. Four corner fibers at `(±hy, ±hz)`
/// reproduce `EA`, `EIy`, `EIz` exactly with zero product-of-inertia
/// coupling (`FiberSection3`'s `four_corner_elastic_fibers_reproduce_ea_
/// eiy_eiz_with_zero_cross_coupling` unit test checks this piece in
/// isolation). `DispBeamColumn3::strain_displacement`'s shape functions are
/// the same cubic-Hermite/linear-axial fields `ElasticBeamColumn3` is
/// built from, so for constant `EA`/`EIy`/`EIz` any >=2-point Gauss rule
/// integrates the (at most quadratic) integrand exactly.
#[test]
fn four_corner_fiber_section_matches_elastic_beam_column3_exactly_in_both_planes() {
    let (e, g, area, iy, iz, j, length): (f64, f64, f64, f64, f64, f64, f64) =
        (30_000.0, 12_000.0, 4.0, 500.0, 2000.0, 50.0, 100.0);
    let (fy, fz, torque, axial) = (-10.0, -6.0, 20.0, 500.0);

    let hz = (iy / area).sqrt();
    let hy = (iz / area).sqrt();
    let a4 = area / 4.0;
    let fibers = || {
        vec![
            Fiber3::new(hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(hy, -hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, hz, a4, Material::Elastic { e }),
            Fiber3::new(-hy, -hz, a4, Material::Elastic { e }),
        ]
    };

    let disp_tip = {
        let mut domain = Domain3::new();
        let node_i = domain.add_node(
            Node3::new([0.0, 0.0, 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        let node_j = domain.add_node(Node3::new([length, 0.0, 0.0]));
        domain.load_node(node_j, SpatialDof::Ux as usize, axial);
        domain.load_node(node_j, SpatialDof::Uy as usize, fy);
        domain.load_node(node_j, SpatialDof::Uz as usize, fz);
        domain.load_node(node_j, SpatialDof::Rx as usize, torque);
        domain.add_element(Element3::DispBeamColumn3(DispBeamColumn3::new(
            node_i,
            node_j,
            g,
            j,
            [0.0, 0.0, 1.0],
            fibers(),
            BeamIntegration::Legendre { points: 3 },
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
        analysis
            .step()
            .expect("cantilever DispBeamColumn3 should solve");
        analysis.domain().node(node_j).displacement
    };

    let elastic_tip = {
        let mut domain = Domain3::new();
        let node_i = domain.add_node(
            Node3::new([0.0, 0.0, 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        let node_j = domain.add_node(Node3::new([length, 0.0, 0.0]));
        domain.load_node(node_j, SpatialDof::Ux as usize, axial);
        domain.load_node(node_j, SpatialDof::Uy as usize, fy);
        domain.load_node(node_j, SpatialDof::Uz as usize, fz);
        domain.load_node(node_j, SpatialDof::Rx as usize, torque);
        domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
            node_i,
            node_j,
            e,
            g,
            area,
            j,
            iy,
            iz,
            GeomTransf3::linear([0.0, 0.0, 1.0]),
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
        analysis
            .step()
            .expect("cantilever ElasticBeamColumn3 should solve");
        analysis.domain().node(node_j).displacement
    };

    for dof in 0..6 {
        assert!(
            (disp_tip[dof] - elastic_tip[dof]).abs() < 1e-9,
            "dof {dof}: DispBeamColumn3={}, ElasticBeamColumn3={}, should match exactly",
            disp_tip[dof],
            elastic_tip[dof]
        );
    }

    // Also confirm against closed-form cantilever tip values directly, so
    // this isn't just "the two happen to agree."
    assert!((disp_tip[0] - axial * length / (e * area)).abs() < 1e-9);
    assert!((disp_tip[1] - fy * length.powi(3) / (3.0 * e * iz)).abs() < 1e-9);
    assert!((disp_tip[2] - fz * length.powi(3) / (3.0 * e * iy)).abs() < 1e-9);
    assert!((disp_tip[3] - torque * length / (g * j)).abs() < 1e-9);
}
