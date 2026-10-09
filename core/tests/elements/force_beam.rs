//! `ForceBeamColumn` / `ForceBeamColumn3`: exact for prismatic elastic members, section coupling,
//! and yielding with permanent set, in 2D and both 3D bending planes.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    BeamIntegration, Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3, Element, Element3,
    Fiber, Fiber3, ForceBeamColumn, ForceBeamColumn3, GeomTransf, GeomTransf3, Material, Node,
    Node3, SpatialDof,
};

fn cantilever_with_tip_load(
    fibers: Vec<Fiber>,
    integration: BeamIntegration,
    length: f64,
    tip_load: f64,
) -> (Domain, carapace_core::model::NodeId) {
    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 1, tip_load);
    domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
        node_i,
        node_j,
        fibers,
        integration,
    )));
    (domain, node_j)
}

/// For a prismatic elastic member,
/// a force-based element is exact regardless of point count — the
/// textbook property that motivates the whole force-based formulation.
/// Same structure as `elements/disp_beam.rs`'s equivalent check: build
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
        let node_j = domain.add_node(Node::new([length, 0.0]));
        domain.load_node(node_j, 1, tip_load);
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
        (force_beam_tip - elastic_beam_tip).abs() < 1e-6,
        "ForceBeamColumn={force_beam_tip}, ElasticBeamColumn={elastic_beam_tip}, should match exactly"
    );

    let expected = tip_load * length.powi(3) / (3.0 * e * iz);
    assert!(
        (force_beam_tip - expected).abs() < 1e-6,
        "expected {expected}, got {force_beam_tip}"
    );
}

/// `ForceBeamColumn`'s basic-force/section coupling (the
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
    let node_j = domain.add_node(Node::new([length, 0.0]));
    domain.load_node(node_j, 0, tip_axial);
    domain.load_node(node_j, 1, tip_transverse);
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

/// `ForceBeamColumn`'s state-determination Newton loop must
/// actually engage material nonlinearity (not just reproduce the elastic
/// case), and its `q_commit`/`e_commit` fields must correctly carry
/// converged history across separate `Analysis::step()` calls, the same
/// property `materials/state.rs`'s permanent-set test checks for
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
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
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
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
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
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
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

/// For a prismatic elastic
/// member, `ForceBeamColumn3` is exact regardless of point count — same
/// property `elements/force_beam.rs` checks planarly, now exercising both
/// bending planes plus axial and (decoupled) torsion. This is the actual
/// verification for `ForceBeamColumn3`'s basic system (`v4`/`v5`'s
/// rotation-coefficient sign flip relative to `v2`/`v3`) — see that type's
/// doc comment for why this test, not inspection, is the source of truth.
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
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            })
            .build(domain);
        analysis
            .step()
            .expect("cantilever ForceBeamColumn3 should solve");
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

/// `ForceBeamColumn3`'s state-determination Newton loop
/// must engage material nonlinearity in *both* bending planes independently
/// (not just reproduce the elastic case, and not just one plane) — the
/// spatial analogue of `elements/force_beam.rs`'s
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

    for (tip_load_dof, disp_dof, ei) in [
        (SpatialDof::Uy as usize, SpatialDof::Uy as usize, iz),
        (SpatialDof::Uz as usize, SpatialDof::Uz as usize, iy),
    ] {
        // Elastic range: must match the closed form exactly.
        let elastic_tip_load = -0.5 * p_yield;
        let (domain1, node_j1) = build(tip_load_dof, elastic_tip_load);
        let mut analysis1 = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 30,
            })
            .build(domain1);
        analysis1
            .step()
            .expect("elastic-range step should converge");
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
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 30,
            })
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
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 30,
            })
            .build(analysis2.domain().clone());
        analysis3
            .step()
            .expect("unload-to-zero step should converge");
        let u3 = analysis3.domain().node(node_j2).displacement[disp_dof];
        assert!(
            u3.abs() > 1e-6,
            "dof {disp_dof}: unloading to zero load should leave a permanent set, got {u3}"
        );
        assert!(u3.abs() < u2.abs(), "dof {disp_dof}: unloaded displacement should be smaller in magnitude than at peak load");
    }
}

/// The two tests above
/// both use doubly-symmetric fiber layouts with zero product of inertia
/// (`Iyz = 0`), so neither exercises `FiberSection3::trial`'s `eiyz` tangent
/// term or the resulting cross-coupling between the two bending planes that
/// a real, non-doubly-symmetric section has. This is the spatial analogue of
/// `elements/force_beam.rs`'s `asymmetric_elastic_section_matches_hand_
/// derived_closed_form`, generalized from a 2x2 to a 3x3 section flexibility.
///
/// Same technique as that planar test: for a single-element cantilever with
/// only tip forces applied (no end moments), all five basic forces
/// `q = [q1..q5]` are fixed by nodal equilibrium alone, independent of
/// section properties (`q1` = axial, `q3 = q5 = 0` since there's no applied
/// end moment in either plane, and `q2`/`q4` follow from the tip shear
/// balance in each plane via `basic_deformation_matrix`'s `A^T` map, exactly
/// as the planar test derives `q2 = -tip_transverse * length`). Then
/// `v = F*q` gives the exact basic deformation from the section's *elastic*
/// flexibility (constant along the member, since an elastic material's
/// tangent doesn't depend on strain), with `F = integral(b^T f b) dx`
/// evaluated in closed form from `b_matrix`'s `(xi-1)`/`xi` shape functions
/// (`integral (xi-1)^2 = integral xi^2 = 1/3`, `integral (xi-1)*xi = -1/6`,
/// `integral (xi-1) = -1/2`, `integral xi = 1/2` over `xi` in `[0,1]`) — the
/// same integrals the planar test evaluates by hand, just applied to a 3x3
/// `f` instead of a 2x2 one. `f = k^-1` (`k` `FiberSection3::trial`'s
/// tangent) is inverted here via the explicit 3x3 cofactor/determinant
/// formula, independently of the `nalgebra` inversion the element itself
/// uses. Only `v1`/`v2`/`v4` are needed (not `v3`/`v5`, which fix the tip
/// *rotations* `Rz`/`Ry`): with node_i fully fixed, `basic_deformation_
/// matrix`'s `v2`/`v4` rows depend only on `Uy_j`/`Uz_j`, giving the tip
/// translations directly.
#[test]
fn asymmetric_biaxial_section_matches_hand_derived_closed_form() {
    let (e, g, j, length): (f64, f64, f64, f64) = (30_000.0, 12_000.0, 50.0, 100.0);
    let tip_axial = 500.0;
    let (fy, fz) = (-10.0, 6.0);

    // An asymmetric layout (no symmetry about either axis and no fiber at
    // the centroid), giving nonzero EQz, EQy, *and* EIyz — the section
    // property this test exists to exercise.
    let fibers = vec![
        Fiber3::new(6.0, 1.0, 1.0, Material::Elastic { e }),
        Fiber3::new(1.0, 2.0, 2.0, Material::Elastic { e }),
        Fiber3::new(-4.0, -1.0, 3.0, Material::Elastic { e }),
    ];

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
    domain.load_node(node_j, SpatialDof::Ux as usize, tip_axial);
    domain.load_node(node_j, SpatialDof::Uy as usize, fy);
    domain.load_node(node_j, SpatialDof::Uz as usize, fz);
    domain.add_element(Element3::ForceBeamColumn3(ForceBeamColumn3::new(
        node_i,
        node_j,
        g,
        j,
        [0.0, 0.0, 1.0],
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
        .expect("cantilever ForceBeamColumn3 should solve");
    let d = analysis.domain().node(node_j).displacement;

    // Section stiffness sums (EA, EQz, EQy, EIzz, EIyy, EIyz), matching
    // `FiberSection3::trial`'s tangent exactly but computed independently
    // here from the fiber layout rather than by calling `trial`.
    let ea = e * (1.0 + 2.0 + 3.0);
    let eqz = e * (6.0 * 1.0 + 1.0 * 2.0 + -4.0 * 3.0);
    let eqy = e * (1.0 * 1.0 + 2.0 * 2.0 + -1.0 * 3.0);
    let eizz = e * (36.0 * 1.0 + 1.0 * 2.0 + 16.0 * 3.0);
    let eiyy = e * (1.0 * 1.0 + 4.0 * 2.0 + 1.0 * 3.0);
    let eiyz = e * (6.0 * 1.0 * 1.0 + 1.0 * 2.0 * 2.0 + -4.0 * -1.0 * 3.0);

    // k = [[a,b,c],[b,d,g],[c,g,f]], matching `FiberSection3::trial`'s
    // `[[ea,-eqz,eqy],[-eqz,eizz,-eiyz],[eqy,-eiyz,eiyy]]`. Inverted here by
    // the explicit symmetric-3x3 cofactor/determinant formula.
    let (ka, kb, kc, kd, ke, kf) = (ea, -eqz, eqy, eizz, -eiyz, eiyy);
    let det = ka * kd * kf - ka * ke * ke - kb * kb * kf + 2.0 * kb * kc * ke - kc * kc * kd;

    let f11 = (kd * kf - ke * ke) / det;
    let f12 = (kc * ke - kb * kf) / det;
    let f13 = (kb * ke - kc * kd) / det;
    let f22 = (ka * kf - kc * kc) / det;
    let f23 = (kb * kc - ka * ke) / det;
    let f33 = (ka * kd - kb * kb) / det;

    // Nodal equilibrium at the free tip alone, no applied end moments:
    // q1 = axial, q3 = q5 = 0, and (q2, q4) from the tip shear balance in
    // each plane, the direct biaxial extension of the planar test's
    // `q2 = -tip_transverse * length`.
    let (q1, q2, q4) = (tip_axial, -fy * length, -fz * length);

    let l = length;
    let v1 = l * f11 * q1 - (l * f12 / 2.0) * q2 - (l * f13 / 2.0) * q4;
    let v2 = -(l * f12 / 2.0) * q1 + (l * f22 / 3.0) * q2 + (l * f23 / 3.0) * q4;
    let v4 = -(l * f13 / 2.0) * q1 + (l * f23 / 3.0) * q2 + (l * f33 / 3.0) * q4;

    let expected_axial = v1;
    let expected_uy = -l * v2;
    let expected_uz = -l * v4;

    assert!(
        (d[0] - expected_axial).abs() < 1e-6,
        "axial: got {}, expected {expected_axial}",
        d[0]
    );
    assert!(
        (d[1] - expected_uy).abs() < 1e-6,
        "uy: got {}, expected {expected_uy}",
        d[1]
    );
    assert!(
        (d[2] - expected_uz).abs() < 1e-6,
        "uz: got {}, expected {expected_uz}",
        d[2]
    );
}
