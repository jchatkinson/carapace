use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{
    BeamIntegration, Domain3, DispBeamColumn3, ElasticBeamColumn3, Element3, Fiber3, GeomTransf3, Material, Node3, SpatialDof,
};

/// M17 acceptance (spatial-architecture plan): `DispBeamColumn3` (biaxial
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
    let (e, g, area, iy, iz, j, length): (f64, f64, f64, f64, f64, f64, f64) = (30_000.0, 12_000.0, 4.0, 500.0, 2000.0, 50.0, 100.0);
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
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever DispBeamColumn3 should solve");
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
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever ElasticBeamColumn3 should solve");
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
