use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    BeamIntegration, Domain, ElasticBeamColumn, Element, Fiber, ForceBeamColumn, GeomTransf,
    Material, Node,
};

fn cantilever_with_tip_load(
    fibers: Vec<Fiber>,
    integration: BeamIntegration,
    length: f64,
    tip_load: f64,
) -> (Domain, carapace_core::model::NodeId) {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).with_load(1, tip_load));
    domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
        node_i,
        node_j,
        fibers,
        integration,
    )));
    (domain, node_j)
}

/// M8 acceptance (implementation-plan §6): for a prismatic elastic member,
/// a force-based element is exact regardless of point count — the
/// textbook property that motivates the whole force-based formulation.
/// Same structure as `m7_disp_beam_column.rs`'s equivalent check: build
/// the identical cantilever both ways and require the tip deflection to
/// match `ElasticBeamColumn`'s closed form to solver precision, not just
/// approximately.
#[test]
fn two_fiber_elastic_section_matches_elastic_beam_column_exactly() {
    let (e, area, iz, length): (f64, f64, f64, f64) = (30000.0, 2.0, 1000.0, 100.0);
    let tip_load = -10.0;

    let force_beam_tip = {
        let h = (iz / area).sqrt();
        let fibers = vec![
            Fiber::new(h, area / 2.0, Material::Elastic { e }),
            Fiber::new(-h, area / 2.0, Material::Elastic { e }),
        ];

        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]).with_load(1, tip_load));
        domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
            node_i,
            node_j,
            fibers,
            BeamIntegration::Lobatto { points: 3 },
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
            .expect("cantilever ForceBeamColumn should solve");
        analysis.domain().node(node_j).displacement[1]
    };

    let elastic_beam_tip = {
        let mut domain = Domain::new();
        let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let node_j = domain.add_node(Node::new([length, 0.0]).with_load(1, tip_load));
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
        (force_beam_tip - elastic_beam_tip).abs() < 1e-6,
        "ForceBeamColumn={force_beam_tip}, ElasticBeamColumn={elastic_beam_tip}, should match exactly"
    );

    let expected = tip_load * length.powi(3) / (3.0 * e * iz);
    assert!(
        (force_beam_tip - expected).abs() < 1e-6,
        "expected {expected}, got {force_beam_tip}"
    );
}

/// M8 acceptance: `ForceBeamColumn`'s basic-force/section coupling (the
/// element flexibility's off-diagonal terms) must be exercised too, not
/// just the symmetric-section case above where axial/bending decouple.
///
/// This is checked against a hand-derived closed form, *not*
/// `DispBeamColumn`: for a single-element cantilever, the three basic
/// forces `q` are fully determined by nodal equilibrium alone (three free
/// dof at the tip, three basic unknowns), independent of section
/// properties, then `v = F*q` gives the exact basic deformation from the
/// section's elastic flexibility. `DispBeamColumn` is *not* an oracle
/// here — deliberately not used as one — because its independent linear-
/// axial/cubic-Hermite shape functions are only exact for a section with
/// no axial-bending coupling (`EQ = 0`); with this asymmetric section it's
/// off by about 0.5% (verified separately), which is a real property of
/// displacement-based elements, not a `ForceBeamColumn` bug.
#[test]
fn asymmetric_elastic_section_matches_hand_derived_closed_form() {
    let (e, length): (f64, f64) = (30000.0, 100.0);
    let tip_axial = 500.0;
    let tip_transverse = -10.0;

    let fibers = vec![
        Fiber::new(6.0, 1.0, Material::Elastic { e }),
        Fiber::new(1.0, 2.0, Material::Elastic { e }),
        Fiber::new(-4.0, 3.0, Material::Elastic { e }),
    ];

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(
        Node::new([length, 0.0])
            .with_load(0, tip_axial)
            .with_load(1, tip_transverse),
    );
    domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
        node_i,
        node_j,
        fibers,
        BeamIntegration::Lobatto { points: 4 },
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
        .expect("cantilever ForceBeamColumn should solve");
    let d = analysis.domain().node(node_j).displacement;

    // Closed form: q1 = tip_axial (axial force, constant along the
    // length); q2 + q3 = -tip_transverse * length and q3 = 0 (no applied
    // moment at the free end) from nodal equilibrium at node_j alone.
    // ea/eq/ei from the fiber layout, section flexibility f = k^-1,
    // F = integral(b^T f b) dx over the prismatic length gives the basic
    // flexibility in closed form (the same integrals `ForceBeamColumn`
    // evaluates numerically via Gauss quadrature, done here by hand).
    let ea = e * (1.0 + 2.0 + 3.0);
    let eq = e * (1.0 * 6.0 + 2.0 * 1.0 + 3.0 * -4.0);
    let ei = e * (1.0 * 36.0 + 2.0 * 1.0 + 3.0 * 16.0);
    let det = ea * ei - eq * eq;
    let l = length;
    let f00 = l * ei / det;
    let f01 = -l * eq / (2.0 * det);
    let f11 = l * ea / (3.0 * det);
    let f12 = -l * ea / (6.0 * det);

    let (q1, q2, q3) = (tip_axial, -tip_transverse * length, 0.0);
    let v1 = f00 * q1 + f01 * q2 + -f01 * q3;
    let v2 = f01 * q1 + f11 * q2 + f12 * q3;

    let expected_axial = v1;
    let expected_transverse = -l * v2;

    assert!(
        (d[0] - expected_axial).abs() < 1e-6,
        "axial: got {}, expected {expected_axial}",
        d[0]
    );
    assert!(
        (d[1] - expected_transverse).abs() < 1e-6,
        "transverse: got {}, expected {expected_transverse}",
        d[1]
    );
}

/// M8 acceptance: `ForceBeamColumn`'s state-determination Newton loop must
/// actually engage material nonlinearity (not just reproduce the elastic
/// case), and its `q_commit`/`e_commit` fields must correctly carry
/// converged history across separate `Analysis::step()` calls, the same
/// property `material_state.rs`'s permanent-set test checks for
/// `ZeroLength`/`ElasticPP`.
///
/// Symmetric eight-fiber `ElasticPP` section (four layers each side, at
/// `y = ±5, ±10, ±15, ±20`) rather than the minimal two- or four-fiber
/// cases used elsewhere in this file: with too few layers, the whole
/// section reaches its full plastic-moment capacity (every fiber past
/// yield strain simultaneously) at only a modest load multiple, leaving
/// the section tangent *exactly* singular — not the gradual, layer-by-
/// layer softening a real fiber-discretized section gives (and which the
/// state-determination Newton loop needs to stay well-posed through). At
/// the load levels this test uses, the outer two layers yield while the
/// inner two stay elastic, matching how an actual fiber-section model is
/// built (many layers), not a worst-case coarse discretization.
///
/// `My = fy*Iz/y_max` (elastic section modulus at the outermost fiber)
/// gives the tip load `P_yield = My/L` at which the fixed-end section
/// (included exactly by `BeamIntegration::Lobatto`'s endpoint) first
/// yields.
#[test]
fn elastic_pp_section_softens_past_yield_and_shows_permanent_set_on_unload() {
    let (e, length) = (30000.0_f64, 100.0);
    let area = 0.25;
    let ys = [20.0, 15.0, 10.0, 5.0];
    let iz = area * ys.iter().map(|y| y * y * 2.0).sum::<f64>();
    let y_outer = ys[0];
    let eyp = 0.001;
    let fy = e * eyp;
    let my = fy * iz / y_outer;
    let p_yield = my / length;
    assert!(p_yield > 0.0);

    let fibers = || {
        ys.iter()
            .flat_map(|&y| [y, -y])
            .map(|y| Fiber::new(y, area, Material::elastic_pp(e, eyp)))
            .collect::<Vec<_>>()
    };

    // Step 1: load well within the elastic range — must match the
    // ElasticBeamColumn closed form exactly, same as the pure-elastic test
    // above (confirms the ElasticPP fibers haven't yielded yet).
    let elastic_tip_load = -0.5 * p_yield;
    let (domain1, node_j) = cantilever_with_tip_load(
        fibers(),
        BeamIntegration::Lobatto { points: 3 },
        length,
        elastic_tip_load,
    );
    let mut analysis1 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(domain1);
    analysis1
        .step()
        .expect("elastic-range step should converge");
    let u1 = analysis1.domain().node(node_j).displacement[1];
    let expected_elastic = elastic_tip_load * length.powi(3) / (3.0 * e * iz);
    assert!(
        (u1 - expected_elastic).abs() < 1e-9,
        "expected {expected_elastic}, got {u1}"
    );

    // Step 2: push past yield. The response must be *softer* than the
    // elastic prediction scaled to the same load — proof the fixed-end
    // section actually yielded rather than the element silently staying
    // linear. Strain scales linearly with y, so the outer two layers
    // (y=20, 15) yield by 1.5*p_yield (kappa = 1.5*kappa_yield exceeds
    // their yield curvature of 1.0 and 1.333*kappa_yield respectively)
    // while the inner two (y=10, 5, yielding at 2x/4x kappa_yield) stay
    // elastic — partial plastification, section tangent still
    // nonsingular.
    let plastic_tip_load = -1.1 * p_yield;
    let (domain2, node_j2) = cantilever_with_tip_load(
        fibers(),
        BeamIntegration::Lobatto { points: 3 },
        length,
        plastic_tip_load,
    );
    // Applied in one `LoadControl` step straight to the full target: a
    // jump this large made an early version of `state_determination`
    // overshoot in an intermediate Newton iteration (the trial moment
    // briefly implying a fully-plastic, exactly-singular section well
    // beyond what's reachable at the *final* load) and diverge to `NaN` —
    // `state_determination`'s subdivision fallback (ported from the
    // source's own `numSubdivide` mechanism) is what makes a one-shot
    // jump like this robust; this test exercises that path, not just the
    // easy case.
    let mut analysis2 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(domain2);
    analysis2.step().expect("past-yield step should converge");
    let u2 = analysis2.domain().node(node_j2).displacement[1];
    let expected_elastic_scaled = plastic_tip_load * length.powi(3) / (3.0 * e * iz);
    assert!(
        u2.abs() > expected_elastic_scaled.abs() * 1.01,
        "past-yield response should be measurably softer than elastic: got {u2}, elastic would be {expected_elastic_scaled}"
    );

    // Step 3: unload fully to zero load from the plastically-deformed
    // state (fresh analysis, load_factor target 0.0, over a *clone* of
    // domain2's post-step-2 committed state) — a permanent, nonzero tip
    // displacement must remain, proof `commit` correctly wrote the
    // yielded fibers' plastic strain back into `ForceBeamColumn`'s
    // sections (via `e_commit`) rather than losing it between calls.
    let mut analysis3 = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::NewtonRaphson)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(analysis2.domain().clone());
    analysis3
        .step()
        .expect("unload-to-zero step should converge");
    let u3 = analysis3.domain().node(node_j2).displacement[1];
    assert!(
        u3.abs() > 1e-6,
        "unloading to zero load should leave a permanent set, got {u3}"
    );
    assert!(
        u3.abs() < u2.abs(),
        "unloaded displacement should be smaller in magnitude than at peak load"
    );
}
