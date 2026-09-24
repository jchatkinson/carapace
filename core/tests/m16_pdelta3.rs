use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain3, Element3, ElasticBeamColumn3, GeomTransf3, Node3, SpatialDof};

/// `GeomTransf3::PDelta3` sanity check, the spatial analogue of
/// `m3_beam.rs`'s planar `pdelta_matches_linear_at_zero_axial_force_and_
/// softens_under_compression` (see its doc comment for why each case runs
/// as two `LoadControl` increments rather than one): at zero axial force
/// `PDelta3` must reduce exactly to `Linear3`, and under compression it
/// must soften the bending stiffness relative to `Linear3`.
///
/// Run in *both* bending planes (`Uy`, governed by `Iz`, and `Uz`, governed
/// by `Iy`) — not just one — because `ElasticBeamColumn3::geometric_
/// stiffness`'s `y`-plane block is a hand-derived sign flip of the
/// `z`-plane one (`ry = -dw/dx` vs `rz = +dv/dx`), not copied from a
/// reference; a bug there would only show up by actually exercising the
/// `Iy` plane, not just `Iz`.
#[test]
fn pdelta3_matches_linear3_at_zero_axial_force_and_softens_under_compression_in_both_planes() {
    let (e, g, a, j, iy, iz) = (30_000.0, 12_000.0, 100.0, 100.0, 1000.0, 2000.0);
    let length = 100.0;

    let build = |transform: GeomTransf3, axial_load: f64, transverse_dof: usize| {
        let mut domain = Domain3::new();
        let fixed = domain.add_node(
            Node3::new([0.0, 0.0, 0.0])
                .fix(SpatialDof::Ux as usize)
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Uz as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        let free = domain.add_node(Node3::new([length, 0.0, 0.0]));
        domain.load_node(free, SpatialDof::Ux as usize, axial_load);
        domain.load_node(free, transverse_dof, -1.0);
        domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
            fixed, free, e, g, a, j, iy, iz, transform,
        )));

        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 0.5 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
            .build(domain);

        analysis.step().expect("cantilever beam should solve (increment 1/2)");
        analysis.step().expect("cantilever beam should solve (increment 2/2)");
        analysis.domain().node(free).displacement[transverse_dof]
    };

    for transverse_dof in [SpatialDof::Uy as usize, SpatialDof::Uz as usize] {
        let v_linear_no_axial = build(GeomTransf3::linear([0.0, 0.0, 1.0]), 0.0, transverse_dof);
        let v_pdelta_no_axial = build(GeomTransf3::p_delta([0.0, 0.0, 1.0]), 0.0, transverse_dof);
        assert!(
            (v_linear_no_axial - v_pdelta_no_axial).abs() < 1e-9,
            "dof {transverse_dof}: PDelta3 with zero axial force should match Linear3 exactly"
        );

        let v_linear_compression = build(GeomTransf3::linear([0.0, 0.0, 1.0]), -1000.0, transverse_dof);
        let v_pdelta_compression = build(GeomTransf3::p_delta([0.0, 0.0, 1.0]), -1000.0, transverse_dof);
        assert!(
            v_pdelta_compression.abs() > v_linear_compression.abs(),
            "dof {transverse_dof}: compressive PDelta3 should soften the beam relative to Linear3: got {v_pdelta_compression} vs {v_linear_compression}"
        );
    }
}
