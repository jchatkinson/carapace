use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Node};

/// M3 acceptance (implementation-plan §7): a single `ElasticBeamColumn`
/// spanning a simply-supported beam, under a uniform transverse element
/// load, matches the classical closed-form end rotation
/// `theta = w*L^3 / (24*E*I)`.
///
/// A single 2-node beam element with consistent (virtual-work) equivalent
/// nodal loads exactly reproduces the nodal DOFs of the continuous beam
/// solution for element loads with no interior point loads — the same
/// identity underlying the classical slope-deflection / moment-distribution
/// fixed-end-moment method — so this is an exact check, not an
/// approximation, despite using only one element.
#[test]
fn simply_supported_beam_udl_matches_closed_form_end_rotation() {
    let (e, iz, area) = (30000.0, 1000.0, 100.0);
    let length = 100.0;
    let w = -1.0; // local +y transverse load

    let mut domain = Domain::new();
    // Pin at node_i (ux, uy fixed), roller at node_j (uy fixed) — a simply
    // supported beam. Rotations free at both.
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

    let theta_i = analysis.domain().node(node_i).displacement[2];
    let theta_j = analysis.domain().node(node_j).displacement[2];

    let expected = (w.abs() * length.powi(3)) / (24.0 * e * iz);
    assert!(
        (theta_i.abs() - expected).abs() < 1e-9,
        "expected |theta_i|={expected}, got {}",
        theta_i.abs()
    );
    assert!(
        (theta_j.abs() - expected).abs() < 1e-9,
        "expected |theta_j|={expected}, got {}",
        theta_j.abs()
    );
    assert!(
        (theta_i + theta_j).abs() < 1e-9,
        "end rotations should be equal and opposite by symmetry"
    );
}

/// `GeomTransf::PDelta` sanity check: with zero axial force it must reduce
/// exactly to `Linear` (the geometric-stiffness correction is a function of
/// axial force, so it vanishes at N=0), and under a compressive axial force
/// it must soften the bending stiffness relative to `Linear` (the intended
/// P-Delta direction — see beam.rs's `geometric_stiffness` doc comment).
///
/// Each case runs as *two* equal `LoadControl` increments, not one: within
/// a single `Algorithm::Linear` solve from a zero initial state, the axial
/// force used for the geometric-stiffness correction is read from the
/// *pre-step* displacement (always zero on a first step from rest), so a
/// one-shot solve can never actually exercise it (see beam.rs's doc
/// comment on why path-following P-Delta needs M4's Newton iteration).
/// Splitting into two increments lets the first establish a nonzero axial
/// force that the second increment's tangent stiffness can then react to —
/// a legitimate (if approximate) incremental-linear P-Delta technique.
#[test]
fn pdelta_matches_linear_at_zero_axial_force_and_softens_under_compression() {
    let (e, iz, area) = (30000.0, 1000.0, 100.0);
    let length = 100.0;

    let build = |transform: GeomTransf, axial_load: f64| {
        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]).with_load(0, axial_load).with_load(1, -1.0));
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            node_i, node_j, e, area, iz, transform,
        )));
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 0.5 })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);
        analysis.step().expect("cantilever beam should solve (increment 1/2)");
        analysis.step().expect("cantilever beam should solve (increment 2/2)");
        analysis.domain().node(node_j).displacement[1]
    };

    // Zero axial force: PDelta's geometric-stiffness correction vanishes,
    // so the transverse-tip-load response must match Linear exactly.
    let v_linear_no_axial = build(GeomTransf::Linear, 0.0);
    let v_pdelta_no_axial = build(GeomTransf::PDelta, 0.0);
    assert!(
        (v_linear_no_axial - v_pdelta_no_axial).abs() < 1e-9,
        "PDelta with zero axial force should match Linear exactly"
    );

    // Compressive axial force: PDelta should soften the beam, i.e. produce
    // a larger-magnitude transverse tip deflection than Linear for the
    // same transverse load.
    let v_linear_compression = build(GeomTransf::Linear, -1000.0);
    let v_pdelta_compression = build(GeomTransf::PDelta, -1000.0);
    assert!(
        v_pdelta_compression.abs() > v_linear_compression.abs(),
        "compressive PDelta should soften the beam relative to Linear: got {v_pdelta_compression} vs {v_linear_compression}"
    );
}
