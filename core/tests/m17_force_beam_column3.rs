use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{
    BeamIntegration, Domain3, ElasticBeamColumn3, Element3, Fiber3, ForceBeamColumn3, GeomTransf3, Material, Node3, SpatialDof,
};

/// M17 acceptance (spatial-architecture plan): for a prismatic elastic
/// member, `ForceBeamColumn3` is exact regardless of point count — same
/// property `m8_force_beam_column.rs` checks planarly, now exercising both
/// bending planes plus axial and (decoupled) torsion. This is the actual
/// verification for `ForceBeamColumn3`'s basic system (`v4`/`v5`'s
/// rotation-coefficient sign flip relative to `v2`/`v3`) — see that type's
/// doc comment for why this test, not inspection, is the source of truth.
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

    let force_tip = {
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
        domain.add_element(Element3::ForceBeamColumn3(ForceBeamColumn3::new(
            node_i,
            node_j,
            g,
            j,
            [0.0, 0.0, 1.0],
            fibers(),
            BeamIntegration::Lobatto { points: 3 },
        )));

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);
        analysis.step().expect("cantilever ForceBeamColumn3 should solve");
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
            (force_tip[dof] - elastic_tip[dof]).abs() < 1e-6,
            "dof {dof}: ForceBeamColumn3={}, ElasticBeamColumn3={}, should match exactly",
            force_tip[dof],
            elastic_tip[dof]
        );
    }

    assert!((force_tip[0] - axial * length / (e * area)).abs() < 1e-6);
    assert!((force_tip[1] - fy * length.powi(3) / (3.0 * e * iz)).abs() < 1e-6);
    assert!((force_tip[2] - fz * length.powi(3) / (3.0 * e * iy)).abs() < 1e-6);
    assert!((force_tip[3] - torque * length / (g * j)).abs() < 1e-6);
}

/// M17 acceptance: `ForceBeamColumn3`'s state-determination Newton loop
/// must engage material nonlinearity in *both* bending planes independently
/// (not just reproduce the elastic case, and not just one plane) — the
/// spatial analogue of `m8_force_beam_column.rs`'s
/// `elastic_pp_section_softens_past_yield_and_shows_permanent_set_on_
/// unload`. A doubly-symmetric eight-fiber section (four fibers on each of
/// the `+y`/`-y`/`+z`/`-z` faces) yields independently under a `y`-only or
/// `z`-only tip load, each past its own `p_yield`.
#[test]
fn elastic_pp_section_softens_past_yield_in_both_planes_and_shows_permanent_set_on_unload() {
    let (e, g, j, length) = (30_000.0_f64, 12_000.0, 50.0, 100.0);
    let area = 0.25;
    let ys = [20.0, 15.0, 10.0, 5.0];
    let iz = area * ys.iter().map(|y| y * y * 2.0).sum::<f64>();
    let iy = iz; // doubly-symmetric: same layout mirrored onto z.
    let y_outer = ys[0];
    let eyp = 0.001;
    let fy_yield = e * eyp;
    let my = fy_yield * iz / y_outer;
    let p_yield = my / length;
    assert!(p_yield > 0.0);

    let fibers = || {
        ys.iter()
            .flat_map(|&y| [(y, 0.0), (-y, 0.0), (0.0, y), (0.0, -y)])
            .map(|(y, z)| Fiber3::new(y, z, area, Material::elastic_pp(e, eyp)))
            .collect::<Vec<_>>()
    };

    let build = |tip_load_dof: usize, tip_load: f64| {
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
        domain.load_node(node_j, tip_load_dof, tip_load);
        domain.add_element(Element3::ForceBeamColumn3(ForceBeamColumn3::new(
            node_i,
            node_j,
            g,
            j,
            [0.0, 0.0, 1.0],
            fibers(),
            BeamIntegration::Lobatto { points: 3 },
        )));
        (domain, node_j)
    };

    for (tip_load_dof, disp_dof, ei) in [(SpatialDof::Uy as usize, SpatialDof::Uy as usize, iz), (SpatialDof::Uz as usize, SpatialDof::Uz as usize, iy)] {
        // Elastic range: must match the closed form exactly.
        let elastic_tip_load = -0.5 * p_yield;
        let (domain1, node_j1) = build(tip_load_dof, elastic_tip_load);
        let mut analysis1 = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::NewtonRaphson)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-10, max_iter: 30 })
            .build(domain1);
        analysis1.step().expect("elastic-range step should converge");
        let u1 = analysis1.domain().node(node_j1).displacement[disp_dof];
        let expected_elastic = elastic_tip_load * length.powi(3) / (3.0 * e * ei);
        assert!(
            (u1 - expected_elastic).abs() < 1e-9,
            "dof {disp_dof}: expected {expected_elastic}, got {u1}"
        );

        // Past yield: measurably softer than the elastic prediction.
        let plastic_tip_load = -1.1 * p_yield;
        let (domain2, node_j2) = build(tip_load_dof, plastic_tip_load);
        let mut analysis2 = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::NewtonRaphson)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-10, max_iter: 30 })
            .build(domain2);
        analysis2.step().expect("past-yield step should converge");
        let u2 = analysis2.domain().node(node_j2).displacement[disp_dof];
        let expected_elastic_scaled = plastic_tip_load * length.powi(3) / (3.0 * e * ei);
        assert!(
            u2.abs() > expected_elastic_scaled.abs() * 1.01,
            "dof {disp_dof}: past-yield response should be measurably softer: got {u2}, elastic would be {expected_elastic_scaled}"
        );

        // Unload to zero: a permanent set must remain.
        let mut analysis3 = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 0.0 })
            .algorithm(Algorithm::NewtonRaphson)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-10, max_iter: 30 })
            .build(analysis2.domain().clone());
        analysis3.step().expect("unload-to-zero step should converge");
        let u3 = analysis3.domain().node(node_j2).displacement[disp_dof];
        assert!(u3.abs() > 1e-6, "dof {disp_dof}: unloading to zero load should leave a permanent set, got {u3}");
        assert!(u3.abs() < u2.abs(), "dof {disp_dof}: unloaded displacement should be smaller in magnitude than at peak load");
    }
}
