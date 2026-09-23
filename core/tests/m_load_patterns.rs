use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator};
use carapace_core::model::{Domain, Element, ElasticBeamColumn, GeomTransf, LoadSeries, Node, NodeId};

/// A cantilever `ElasticBeamColumn` with two independent `LoadPattern`s: a
/// compressive axial ("gravity") load on `default_pattern()`, and a
/// separate transverse ("lateral") load on its own pattern. Axial and
/// bending DOFs are fully decoupled for a straight, axis-aligned prismatic
/// beam (`ElasticBeamColumn::geometric_stiffness` only ever populates the
/// bending rows/columns, and the local-to-global transform here is the
/// identity, since the beam runs along +x) — so the axial DOF's
/// displacement is a clean probe for whether the gravity pattern's
/// contribution actually stayed fixed through a second phase that never
/// touches it directly.
fn build(
    transform: GeomTransf,
    gravity_axial: f64,
    node_j: &mut Option<NodeId>,
) -> (Domain, carapace_core::model::LoadPatternId) {
    let (e, iz, area, length) = (30000.0, 1000.0, 100.0, 100.0);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let j = domain.add_node(Node::new([length, 0.0]));
    *node_j = Some(j);

    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        node_i, j, e, area, iz, transform,
    )));

    let gravity = domain.default_pattern();
    domain.add_nodal_load(gravity, j, 0, gravity_axial);

    // The lateral pattern is created now but deliberately left with no
    // nodal load until phase 2 (see `run_two_phase`) — adding the load
    // itself only after gravity's phase has already converged is what
    // keeps phase 1 genuinely gravity-only, since an *unfrozen* pattern
    // with no reference load contributes exactly zero regardless of how
    // far its series has ramped.
    let lateral_pattern = domain.add_load_pattern(LoadSeries::Linear { slope: 1.0 });

    (domain, lateral_pattern)
}

/// Runs gravity to full scale (two `LoadControl` increments — `PDelta`'s
/// geometric-stiffness correction is read from the *pre-step* axial force,
/// so it needs a first step to establish a nonzero force before a second
/// step's tangent can react to it, same subtlety `m3_beam.rs` documents),
/// freezes gravity's pattern, then runs a second phase (a fresh `Analysis`
/// built from the same `Domain`, pseudo-time restarting at 0) that ramps
/// only the lateral pattern. Returns `(axial_after_gravity,
/// axial_after_lateral, transverse_after_lateral)`.
fn run_two_phase(
    domain: Domain,
    gravity: carapace_core::model::LoadPatternId,
    lateral_pattern: carapace_core::model::LoadPatternId,
    lateral: f64,
    node_j: NodeId,
) -> (f64, f64, f64) {
    let mut gravity_phase = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.5 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    gravity_phase.step().expect("gravity increment 1/2 should solve");
    let gravity_time = gravity_phase.step().expect("gravity increment 2/2 should solve").load_factor;
    let axial_after_gravity = gravity_phase.domain().node(node_j).displacement[0];

    let mut domain = gravity_phase.into_domain();
    domain.hold_pattern_constant(gravity, gravity_time);
    domain.add_nodal_load(lateral_pattern, node_j, 1, lateral);

    let mut lateral_phase = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);
    lateral_phase.step().expect("lateral phase should solve");

    let axial_after_lateral = lateral_phase.domain().node(node_j).displacement[0];
    let transverse_after_lateral = lateral_phase.domain().node(node_j).displacement[1];
    (axial_after_gravity, axial_after_lateral, transverse_after_lateral)
}

#[test]
fn frozen_gravity_pattern_survives_a_second_phase_untouched() {
    let mut node_j = None;
    let (domain, lateral_pattern) = build(GeomTransf::Linear, -1000.0, &mut node_j);
    let gravity = domain.default_pattern();
    let (axial_after_gravity, axial_after_lateral, _) =
        run_two_phase(domain, gravity, lateral_pattern, -1.0, node_j.unwrap());

    assert!(axial_after_gravity.abs() > 1e-9, "gravity phase should have produced nonzero axial displacement");
    assert!(
        (axial_after_gravity - axial_after_lateral).abs() < 1e-9,
        "frozen gravity pattern's axial contribution must be unchanged by the lateral phase: {axial_after_gravity} vs {axial_after_lateral}"
    );
}

/// Same two-phase (frozen gravity → lateral) protocol, run once with
/// `GeomTransf::Linear` and once with `PDelta`: with zero gravity the two
/// must match exactly (no axial force means `PDelta`'s geometric-stiffness
/// correction vanishes), and with compressive gravity `PDelta` must soften
/// the beam relative to `Linear` — the same physical checks `m3_beam.rs`
/// makes within a single `Analysis`, but here driven by two genuinely
/// separate analysis phases sharing one evolving `Domain`.
#[test]
fn pdelta_softening_holds_across_a_frozen_gravity_phase_boundary() {
    let mut nj = None;
    let (d, lp) = build(GeomTransf::Linear, 0.0, &mut nj);
    let g = d.default_pattern();
    let (_, _, v_linear_no_axial) = run_two_phase(d, g, lp, -1.0, nj.unwrap());

    let mut nj = None;
    let (d, lp) = build(GeomTransf::PDelta, 0.0, &mut nj);
    let g = d.default_pattern();
    let (_, _, v_pdelta_no_axial) = run_two_phase(d, g, lp, -1.0, nj.unwrap());

    assert!(
        (v_linear_no_axial - v_pdelta_no_axial).abs() < 1e-9,
        "PDelta with zero gravity axial force should match Linear exactly"
    );

    let mut nj = None;
    let (d, lp) = build(GeomTransf::Linear, -1000.0, &mut nj);
    let g = d.default_pattern();
    let (_, _, v_linear_compression) = run_two_phase(d, g, lp, -1.0, nj.unwrap());

    let mut nj = None;
    let (d, lp) = build(GeomTransf::PDelta, -1000.0, &mut nj);
    let g = d.default_pattern();
    let (_, _, v_pdelta_compression) = run_two_phase(d, g, lp, -1.0, nj.unwrap());

    assert!(
        v_pdelta_compression.abs() > v_linear_compression.abs(),
        "compressive frozen gravity should still soften the PDelta lateral response relative to Linear: got {v_pdelta_compression} vs {v_linear_compression}"
    );
}
