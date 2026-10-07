//! `docs/arclength.md` §9: arc-length acceptance tests beyond the
//! degrading-strength fixtures — elastic closed forms in both profiles, an
//! exactly singular tangent, smooth snap-back, geometric snap-through,
//! bifurcation diagnostics, stop criteria, auto-scaling, cutbacks and
//! rollback, frozen loads, constraints, scale invariance and invalid
//! configurations.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ArcDirection, ArcLength, ArcPredictor,
    ArcScales, ArcSeed, ConstraintHandler, ConvergenceTest, DisplacementTarget, Integrator,
    LineSearch, LoadFactorTarget, StepResult, StopReason, TangentStrategy,
};
use carapace_core::model::{
    two_node_dofs, DofMask, Domain, Domain3, ElasticBeamColumn, Element, Element3, ElementForce,
    ElementOps, GeomTransf, LoadSeries, Material, Node, Node3, NodeId, NodeList, NodeView,
    TangentSink, VectorSink, ZeroLength, ZeroLength3,
};
use nalgebra::{SMatrix, SVector};
use slotmap::new_key_type;

fn newton() -> Algorithm {
    Algorithm::Newton {
        tangent: TangentStrategy::Current,
        line_search: None,
    }
}

fn combined(max_iter: usize) -> ConvergenceTest {
    ConvergenceTest::Combined {
        force_tol: 1e-9,
        moment_tol: 1e-9,
        relative_tol: 1e-9,
        displacement_tol: None,
        max_iter,
    }
}

fn explicit(displacement: f64, load: f64) -> ArcScales {
    ArcScales::Explicit {
        displacement,
        rotation: None,
        load,
    }
}

fn build<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: Domain<NDIM, NDOF, NId, E>,
    integrator: Integrator<NId>,
    test: ConvergenceTest,
) -> Analysis<NDIM, NDOF, NId, E>
where
    NId: slotmap::Key,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    let handler = if domain.has_mp_constraints() {
        ConstraintHandler::Transformation
    } else {
        ConstraintHandler::Plain
    };
    AnalysisBuilder::new()
        .constraint_handler(handler)
        .integrator(integrator)
        .algorithm(newton())
        .test(test)
        .build(domain)
}

/// Elastic spring of stiffness `k` from a fixed node to a node free in x,
/// unit reference load on the free node.
fn elastic_spring(k: f64) -> (Domain, NodeId, NodeId) {
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, free).with_material(0, Material::Elastic { e: k }),
    ));
    domain.load_node(free, 0, 1.0);
    (domain, fixed, free)
}

fn degrading_spring() -> Material {
    Material::hysteretic(
        10.0, 0.01, 6.0, 0.02, 2.0, 0.03, -10.0, -0.01, -6.0, -0.02, -2.0, -0.03, 1.0, 1.0, 0.0,
        0.0, 0.0,
    )
}

fn softening() -> (Domain, NodeId) {
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, free).with_material(0, degrading_spring()),
    ));
    domain.load_node(free, 0, 1.0);
    (domain, free)
}

fn arc(result: &StepResult) -> carapace_core::analysis::ArcStepInfo {
    result.arc.expect("arc step diagnostics")
}

// ---------------------------------------------------------------------------
// Elastic closed forms
// ---------------------------------------------------------------------------

#[test]
fn elastic_closed_form_planar_needs_no_corrector_and_honours_direction() {
    let (k, u_star, lambda_star, s): (f64, f64, f64, f64) = (250.0, 0.02, 4.0, 0.1);
    // Along lambda = k u, each step moves dlambda = s / sqrt((1/(k u*))^2 + (1/lambda*)^2).
    let dlambda = s / ((1.0 / (k * u_star)).powi(2) + (1.0 / lambda_star).powi(2)).sqrt();
    for (direction, sign) in [
        (ArcDirection::Increasing, 1.0),
        (ArcDirection::Decreasing, -1.0),
    ] {
        let (domain, _fixed, free) = elastic_spring(k);
        let mut config = ArcLength::fixed(s, explicit(u_star, lambda_star));
        config.direction = direction;
        let mut analysis = build(domain, Integrator::ArcLength(config), combined(10));
        for step in 1..=5 {
            let result = analysis.step().unwrap();
            let info = arc(&result);
            assert!((result.load_factor - sign * dlambda * step as f64).abs() < 1e-9);
            let u = analysis.domain().node(free).displacement[0];
            assert!((k * u - result.load_factor).abs() < 1e-9);
            assert_eq!(
                info.corrector_solves, 0,
                "exact elastic predictor needs no corrector"
            );
            assert_eq!(info.det_sign, Some(1));
            assert!(!info.bifurcation_suspected);
        }
    }
}

#[test]
fn elastic_closed_form_spatial() {
    let (k, u_star, lambda_star, s): (f64, f64, f64, f64) = (80.0, 0.05, 2.0, 0.2);
    let mut domain = Domain3::new();
    let fixed = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(0)
            .fix(1)
            .fix(2)
            .fix(3)
            .fix(4)
            .fix(5),
    );
    let free = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(0)
            .fix(1)
            .fix(3)
            .fix(4)
            .fix(5),
    );
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(fixed, free).with_material(2, Material::Elastic { e: k }),
    ));
    domain.load_node(free, 2, 1.0);
    let mut analysis = build(
        domain,
        Integrator::ArcLength(ArcLength::fixed(s, explicit(u_star, lambda_star))),
        combined(10),
    );
    let dlambda = s / ((1.0 / (k * u_star)).powi(2) + (1.0 / lambda_star).powi(2)).sqrt();
    for step in 1..=4 {
        let result = analysis.step().unwrap();
        assert!((result.load_factor - dlambda * step as f64).abs() < 1e-9);
        assert!(
            (k * analysis.domain().node(free).displacement[2] - result.load_factor).abs() < 1e-9
        );
        assert_eq!(arc(&result).corrector_solves, 0);
    }
}

// ---------------------------------------------------------------------------
// Test-only conservative exponential spring: F(u) = k u exp(-u/a)
// ---------------------------------------------------------------------------

new_key_type! {
    struct ExpId;
}

/// Uniaxial (x) spring between two planar nodes with the smooth peaked
/// law `F(u) = k u exp(-u/a)`, `u = ux_j - ux_i`; peak `k a / e` at
/// `u = a`, where the tangent is exactly zero. Path-independent.
#[derive(Debug, Clone)]
enum TestElement {
    Exp { nodes: [NodeId; 2], k: f64, a: f64 },
    Linear { nodes: [NodeId; 2], k: f64 },
}

fn exp_force(k: f64, a: f64, u: f64) -> f64 {
    k * u * (-u / a).exp()
}

fn exp_tangent(k: f64, a: f64, u: f64) -> f64 {
    k * (-u / a).exp() * (1.0 - u / a)
}

/// This test catalog has no element loads.
#[derive(Clone, Copy)]
struct NoLoad;

impl std::ops::Add for NoLoad {
    type Output = NoLoad;
    fn add(self, _: NoLoad) -> NoLoad {
        NoLoad
    }
}

impl carapace_core::model::ElementLoadComponents for NoLoad {
    const COUNT: usize = 0;

    fn component(&self, _: usize) -> f64 {
        0.0
    }
}

impl std::ops::Mul<f64> for NoLoad {
    type Output = NoLoad;
    fn mul(self, _: f64) -> NoLoad {
        NoLoad
    }
}

impl TestElement {
    fn pair(&self) -> [NodeId; 2] {
        match self {
            TestElement::Exp { nodes, .. } | TestElement::Linear { nodes, .. } => *nodes,
        }
    }

    fn tangent_and_resistance(
        &self,
        node_i: &Node,
        node_j: &Node,
    ) -> (SMatrix<f64, 6, 6>, SVector<f64, 6>) {
        let u = node_j.displacement[0] - node_i.displacement[0];
        let (force, tangent) = match *self {
            TestElement::Exp { k, a, .. } => (exp_force(k, a, u), exp_tangent(k, a, u)),
            TestElement::Linear { k, .. } => (k * u, k),
        };
        let mut kmat = SMatrix::<f64, 6, 6>::zeros();
        kmat[(0, 0)] = tangent;
        kmat[(3, 3)] = tangent;
        kmat[(0, 3)] = -tangent;
        kmat[(3, 0)] = -tangent;
        let mut r = SVector::<f64, 6>::zeros();
        r[0] = -force;
        r[3] = force;
        (kmat, r)
    }
}

/// A user-defined catalog: the element interface is not tied to the crate's
/// own elements.
impl ElementOps<2, 3, NodeId> for TestElement {
    type Load = NoLoad;
    type Id = ExpId;

    fn nodes(&self) -> NodeList<NodeId> {
        self.pair().into_iter().collect()
    }

    fn dof_mask(&self) -> DofMask {
        DofMask::all(3)
    }

    fn assemble_tangent<S: TangentSink<NodeId>>(
        &self,
        nodes: &NodeView<'_, 2, 3, NodeId>,
        _: Option<&NoLoad>,
        sink: &mut S,
    ) {
        let [i, j] = self.pair();
        let (k, r) = self.tangent_and_resistance(nodes.get(i), nodes.get(j));
        sink.add(&two_node_dofs::<_, 6>(i, j, 3), &k, &r);
    }

    fn assemble_load<S: VectorSink<NodeId>>(
        &self,
        _: &NodeView<'_, 2, 3, NodeId>,
        _: Option<&NoLoad>,
        _: &mut S,
    ) {
    }

    fn assemble_mass<S: VectorSink<NodeId>>(&self, _: &NodeView<'_, 2, 3, NodeId>, _: &mut S) {}

    fn commit(&mut self, _: &NodeView<'_, 2, 3, NodeId>, _: Option<&NoLoad>) {}

    fn local_force_width(&self) -> usize {
        6
    }

    fn local_force(&self, nodes: &NodeView<'_, 2, 3, NodeId>) -> ElementForce {
        let [i, j] = self.pair();
        self.tangent_and_resistance(nodes.get(i), nodes.get(j))
            .1
            .into()
    }

    fn local_load_force(&self, _: &NodeView<'_, 2, 3, NodeId>, _: Option<&NoLoad>) -> ElementForce {
        SVector::<f64, 6>::zeros().into()
    }

    fn fiber_responses(&self, _: &NodeView<'_, 2, 3, NodeId>) -> Option<Vec<Vec<(f64, f64)>>> {
        None
    }
}

type TestDomain = Domain<2, 3, NodeId, TestElement>;

#[test]
fn exactly_singular_tangent_needs_a_seed_and_then_continues_past_the_peak() {
    let (k, a) = (100.0, 0.5);
    let peak = exp_force(k, a, a);
    assert!((peak - k * a / std::f64::consts::E).abs() < 1e-12);
    assert_eq!(
        exp_tangent(k, a, a),
        0.0,
        "tangent is exactly singular at u = a"
    );

    // Start in equilibrium exactly at the peak: displacement a, a constant
    // (frozen-style) preload carrying the peak force, and a linear
    // reference pattern for continuation.
    let make = || {
        let mut domain = TestDomain::new();
        let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let free = domain.add_node(
            Node::new([0.0, 0.0])
                .fix(1)
                .fix(2)
                .with_initial_displacement(0, a),
        );
        domain.add_element(TestElement::Exp {
            nodes: [fixed, free],
            k,
            a,
        });
        let preload = domain.add_load_pattern(LoadSeries::Constant);
        domain.add_nodal_load(preload, free, 0, peak);
        domain.load_node(free, 0, 1.0);
        (domain, free)
    };
    let scales = explicit(a, peak);

    let (domain, _) = make();
    let mut analysis = build(
        domain,
        Integrator::ArcLength(ArcLength::fixed(0.05, scales)),
        combined(30),
    );
    assert_eq!(
        analysis.step().unwrap_err(),
        AnalysisError::MissingSeedDirection
    );

    let (domain, free) = make();
    let mut config = ArcLength::fixed(0.05, scales);
    config.seed = Some(ArcSeed {
        components: vec![(free, 0, 1.0)],
        load: 0.0,
    });
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
    let mut u_prev = a;
    let mut lambda_prev = 0.0;
    for _ in 0..40 {
        let result = analysis
            .step()
            .expect("bordered system handles the singular peak");
        let u = analysis.domain().node(free).displacement[0];
        let lambda = result.load_factor;
        assert!(u > u_prev, "seed orients travel to increasing u");
        assert!(
            lambda < lambda_prev + 1e-12,
            "past the peak the load only decreases"
        );
        assert!((peak + lambda - exp_force(k, a, u)).abs() < 1e-7);
        (u_prev, lambda_prev) = (u, lambda);
    }
    assert!(u_prev > 1.5 * a);
}

#[test]
fn smooth_snap_back_through_nonsingular_displacement_turning_points() {
    let (k, a) = (100.0, 0.5);
    let k_s = 0.1 * k;
    assert!(k_s < k / std::f64::consts::E.powi(2));
    for predictor in [ArcPredictor::Secant, ArcPredictor::Tangent] {
        let mut domain = TestDomain::new();
        let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let x_node = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
        let y_node = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
        domain.add_element(TestElement::Exp {
            nodes: [fixed, x_node],
            k,
            a,
        });
        domain.add_element(TestElement::Linear {
            nodes: [x_node, y_node],
            k: k_s,
        });
        domain.load_node(y_node, 0, 1.0);

        let mut config = ArcLength::adaptive(0.05, 0.005, 0.1, explicit(a, k * a));
        config.predictor = predictor;
        let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
        let (mut x_prev, mut y_prev) = (0.0, 0.0);
        let mut decreasing_y = 0;
        let mut recovering_y = 0;
        for _ in 0..2000 {
            let result = analysis.step().expect("smooth snap-back traced");
            let info = arc(&result);
            assert_eq!(
                info.retries, 0,
                "{predictor:?}: moderate radius needs no cutbacks"
            );
            assert!(
                !info.bifurcation_suspected,
                "limit points are not bifurcations"
            );
            let x = analysis.domain().node(x_node).displacement[0];
            let y = analysis.domain().node(y_node).displacement[0];
            let lambda = result.load_factor;
            assert!((lambda - exp_force(k, a, x)).abs() < 1e-7);
            assert!((y - (x + exp_force(k, a, x) / k_s)).abs() < 1e-9);
            assert!(x > x_prev);
            if x_prev > 1.45 * a && x < 2.95 * a {
                assert!(
                    y < y_prev,
                    "{predictor:?}: snap-back between turning points"
                );
                decreasing_y += 1;
            }
            if x_prev > 3.05 * a {
                assert!(
                    y > y_prev,
                    "{predictor:?}: recovery after second turning point"
                );
                recovering_y += 1;
            }
            (x_prev, y_prev) = (x, y);
            if x > 4.0 * a {
                break;
            }
        }
        assert!(x_prev > 4.0 * a);
        assert!(decreasing_y >= 5 && recovering_y >= 3);
    }
}

// ---------------------------------------------------------------------------
// Geometric snap-through and bifurcation diagnostics
// ---------------------------------------------------------------------------

/// Shallow two-bar (von Mises) truss of corotational elastic members,
/// pinned at the supports and at the apex (two coincident apex nodes tied
/// in translation by `equal_dof`, so each bar is pin-ended and carries no
/// moment), apex x fixed as the symmetry constraint. Oracle:
/// `P(v) = 2 EA (L0 - L)/L0 * (h - v)/L`, `L = sqrt(b^2 + (h - v)^2)`.
#[test]
fn geometric_snap_through_matches_the_two_bar_oracle() {
    let (b, h, e, area, iz): (f64, f64, f64, f64, f64) = (1.0, 0.1, 1e7, 1e-3, 1e-8);
    let ea = e * area;
    let l0 = (b * b + h * h).sqrt();
    let oracle = |v: f64| {
        let l = (b * b + (h - v).powi(2)).sqrt();
        2.0 * ea * (l0 - l) / l0 * (h - v) / l
    };

    let mut domain = Domain::new();
    let left = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let right = domain.add_node(Node::new([2.0 * b, 0.0]).fix(0).fix(1));
    let apex = domain.add_node(Node::new([b, h]).fix(0));
    let apex_twin = domain.add_node(Node::new([b, h]));
    domain.equal_dof(apex, apex_twin, &[0, 1]);
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        left,
        apex,
        e,
        area,
        iz,
        GeomTransf::Corotational,
    )));
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        apex_twin,
        right,
        e,
        area,
        iz,
        GeomTransf::Corotational,
    )));
    domain.load_node(apex, 1, -1.0);

    let config = ArcLength::adaptive(
        0.05,
        0.005,
        0.1,
        ArcScales::Explicit {
            displacement: h,
            rotation: Some(h),
            load: 5.0,
        },
    );
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
    let mut v_prev = 0.0;
    let (mut saw_negative, mut recovered) = (false, false);
    for _ in 0..1000 {
        let result = analysis.step().expect("snap-through traced");
        let v = -analysis.domain().node(apex).displacement[1];
        let lambda = result.load_factor;
        assert!(
            (lambda - oracle(v)).abs() <= 1e-6 * (1.0 + lambda.abs()),
            "v={v}: lambda {lambda} vs oracle {}",
            oracle(v)
        );
        assert!(v > v_prev, "apex keeps moving down");
        saw_negative |= lambda < -0.1;
        recovered |= saw_negative && v > 2.0 * h && lambda > 0.0;
        v_prev = v;
        if v > 2.5 * h {
            break;
        }
    }
    assert!(saw_negative, "descending branch reaches negative load");
    assert!(recovered, "recovering branch after inversion");
}

/// A perfectly straight corotational cantilever under axial compression:
/// the primary (straight) path keeps loading through the Euler load, where
/// one eigenvalue of `K` crosses zero. `lambda` never reverses, so the
/// determinant sign change must be flagged as a suspected bifurcation.
#[test]
fn determinant_sign_change_without_load_reversal_flags_a_bifurcation() {
    let (length, e, area, iz, elements) = (2.0, 2e5, 0.01, 1e-5, 4);
    let p_cr = std::f64::consts::PI.powi(2) * e * iz / (4.0 * length * length);
    let mut domain = Domain::new();
    let mut nodes = Vec::new();
    for i in 0..=elements {
        let y = length * i as f64 / elements as f64;
        let node = Node::new([0.0, y]);
        nodes.push(domain.add_node(if i == 0 {
            node.fix(0).fix(1).fix(2)
        } else {
            node
        }));
    }
    for pair in nodes.windows(2) {
        domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
            pair[0],
            pair[1],
            e,
            area,
            iz,
            GeomTransf::Corotational,
        )));
    }
    let tip = *nodes.last().unwrap();
    domain.load_node(tip, 1, -1.0);

    let shortening = p_cr * length / (e * area);
    let config = ArcLength::fixed(
        0.1,
        ArcScales::Explicit {
            displacement: shortening,
            rotation: Some(0.01),
            load: p_cr,
        },
    );
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
    let mut flagged_at = None;
    let mut lambda_prev = 0.0;
    for step in 1..=40 {
        let result = analysis.step().expect("primary path traced");
        let info = arc(&result);
        assert!(
            result.load_factor > lambda_prev,
            "load keeps increasing on the primary path"
        );
        assert!(
            analysis.domain().node(tip).displacement[0].abs() < 1e-9,
            "stays straight"
        );
        if info.bifurcation_suspected {
            assert!(flagged_at.is_none(), "only one crossing in range");
            flagged_at = Some((step, lambda_prev, result.load_factor));
        }
        lambda_prev = result.load_factor;
        if lambda_prev > 1.4 * p_cr {
            break;
        }
    }
    let (_, before, after) = flagged_at.expect("crossing of the critical load is flagged");
    // The discrete 4-element critical load is within a few percent of Euler.
    assert!(
        before < 1.05 * p_cr && after > 0.95 * p_cr,
        "{before}..{after} vs {p_cr}"
    );
}

// ---------------------------------------------------------------------------
// Stop criteria
// ---------------------------------------------------------------------------

fn softening_config() -> ArcLength {
    let mut config = ArcLength::fixed(0.05, explicit(0.01, 10.0));
    config.arc_tolerance = 1e-9;
    config
}

#[test]
fn displacement_target_lands_exactly_and_further_steps_are_refused() {
    let (domain, free) = softening();
    let mut config = softening_config();
    config.stop.displacement = Some(DisplacementTarget {
        node: free,
        dof: 0,
        value: 0.02,
        exact: true,
    });
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
    let mut stop = None;
    for _ in 0..200 {
        let result = analysis.step().unwrap();
        if let Some(s) = arc(&result).stop {
            stop = Some((s, result.load_factor));
            break;
        }
    }
    let (stop, lambda) = stop.expect("target reached");
    assert_eq!(stop.reason, StopReason::DisplacementTarget);
    assert!(stop.landed_exactly);
    let u = analysis.domain().node(free).displacement[0];
    assert!((u - 0.02).abs() < 1e-9, "u = {u}");
    assert!((lambda - 6.0).abs() < 1e-6, "lambda = {lambda}");
    assert!(stop.overshoot.abs() < 1e-9);

    let before = analysis.domain().node(free).displacement[0];
    assert_eq!(
        analysis.step().unwrap_err(),
        AnalysisError::ContinuationComplete
    );
    assert_eq!(analysis.domain().node(free).displacement[0], before);
}

#[test]
fn load_factor_target_and_zero_crossing() {
    let (domain, _fixed, free) = elastic_spring(100.0);
    let mut config = ArcLength::fixed(0.07, explicit(0.01, 1.0));
    config.stop.load_factor = Some(LoadFactorTarget {
        value: 1.0,
        exact: true,
    });
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(10));
    let mut last = None;
    for _ in 0..100 {
        let result = analysis.step().unwrap();
        if arc(&result).stop.is_some() {
            last = Some(result);
            break;
        }
    }
    let last = last.expect("load target reached");
    let stop = arc(&last).stop.unwrap();
    assert_eq!(stop.reason, StopReason::LoadFactorTarget);
    assert!(stop.landed_exactly);
    assert!((last.load_factor - 1.0).abs() < 1e-9);

    // A new direction changes the path definition: history restarts at
    // lambda = 1 and the zero crossing ends the phase.
    let mut reverse = ArcLength::fixed(0.07, explicit(0.01, 1.0));
    reverse.direction = ArcDirection::Decreasing;
    reverse.stop.load_factor_zero_crossing = true;
    analysis.set_integrator(Integrator::ArcLength(reverse));
    let mut crossed = None;
    for _ in 0..100 {
        let result = analysis.step().unwrap();
        assert!(analysis.domain().node(free).displacement[0] < 0.0101);
        if let Some(stop) = arc(&result).stop {
            crossed = Some((stop, result.load_factor));
            break;
        }
    }
    let (stop, lambda) = crossed.expect("zero crossing reached");
    assert_eq!(stop.reason, StopReason::LoadFactorZeroCrossing);
    assert!(lambda <= 0.0);
}

#[test]
fn chord_and_step_caps_and_reopening_with_new_criteria() {
    let (domain, _fixed, free) = elastic_spring(100.0);
    let mut config = ArcLength::fixed(0.05, explicit(0.01, 1.0));
    config.stop.max_chord_length = Some(0.48);
    let mut analysis = build(domain, Integrator::ArcLength(config.clone()), combined(10));
    let mut steps = 0;
    loop {
        steps += 1;
        let result = analysis.step().unwrap();
        if let Some(stop) = arc(&result).stop {
            assert_eq!(stop.reason, StopReason::ChordLength);
            assert!(arc(&result).chord_length >= 0.48);
            break;
        }
    }
    assert_eq!(steps, 10);

    // Same path definition, new criteria: history (direction, chord) kept.
    let u_before = analysis.domain().node(free).displacement[0];
    let mut more = config;
    more.stop.max_chord_length = None;
    more.stop.max_steps = Some(3);
    analysis.set_integrator(Integrator::ArcLength(more));
    for step in 1..=3 {
        let result = analysis.step().unwrap();
        let stop = arc(&result).stop;
        assert_eq!(
            stop.map(|s| s.reason),
            (step == 3).then_some(StopReason::StepCount)
        );
        assert!(arc(&result).chord_length > 0.48);
    }
    assert!(analysis.domain().node(free).displacement[0] > u_before);
    assert_eq!(
        analysis.step().unwrap_err(),
        AnalysisError::ContinuationComplete
    );
}

// ---------------------------------------------------------------------------
// Auto-scale
// ---------------------------------------------------------------------------

#[test]
fn auto_scale_derives_the_hand_computed_scale_and_reproduces_the_path() {
    let run = |scales: ArcScales| {
        let (domain, free) = softening();
        let mut config = softening_config();
        config.scales = scales;
        let mut analysis = build(domain, Integrator::ArcLength(config), combined(30));
        let mut path = Vec::new();
        let mut scale = 0.0;
        for _ in 0..40 {
            let result = analysis.step().unwrap();
            scale = arc(&result).displacement_scale;
            path.push((
                analysis.domain().node(free).displacement[0],
                result.load_factor,
            ));
        }
        (path, scale)
    };
    let (explicit_path, explicit_scale) = run(explicit(0.01, 10.0));
    let (auto_path, auto_scale) = run(ArcScales::Auto { load: 10.0 });
    // K v = p with K = 1000: u* = 10 * 1/1000 = 0.01.
    assert_eq!(explicit_scale, 0.01);
    assert!((auto_scale - 0.01).abs() < 1e-15);
    for ((u_e, l_e), (u_a, l_a)) in explicit_path.iter().zip(&auto_path) {
        assert!((u_e - u_a).abs() < 1e-12 && (l_e - l_a).abs() < 1e-9);
    }
}

#[test]
fn auto_scale_with_singular_first_tangent_fails_initialization() {
    let (k, a) = (100.0, 0.5);
    let peak = exp_force(k, a, a);
    let mut domain = TestDomain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(
        Node::new([0.0, 0.0])
            .fix(1)
            .fix(2)
            .with_initial_displacement(0, a),
    );
    domain.add_element(TestElement::Exp {
        nodes: [fixed, free],
        k,
        a,
    });
    let preload = domain.add_load_pattern(LoadSeries::Constant);
    domain.add_nodal_load(preload, free, 0, peak);
    domain.load_node(free, 0, 1.0);
    let mut analysis = build(
        domain,
        Integrator::ArcLength(ArcLength::fixed(0.05, ArcScales::Auto { load: peak })),
        combined(30),
    );
    assert_eq!(analysis.step().unwrap_err(), AnalysisError::SingularSystem);
}

// ---------------------------------------------------------------------------
// Cutbacks, rollback and material-history isolation
// ---------------------------------------------------------------------------

#[test]
fn oversized_radius_cuts_back_and_the_accepted_path_stays_on_the_oracle() {
    // Smooth nonlinearity with a deliberately tiny budget of two corrector
    // solves: a large radius genuinely needs more than that, a smaller one
    // doesn't, so cutbacks must recover every step that starts oversized.
    let (k, a) = (100.0, 0.5);
    let mut domain = TestDomain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(TestElement::Exp {
        nodes: [fixed, free],
        k,
        a,
    });
    domain.load_node(free, 0, 1.0);
    let config = ArcLength::adaptive(1.0, 0.01, 1.0, explicit(a, k * a));
    let mut analysis = build(domain, Integrator::ArcLength(config), combined(2));
    let mut total_retries = 0;
    let mut u_prev = 0.0;
    for _ in 0..500 {
        let result = analysis.step().expect("cutbacks recover");
        let info = arc(&result);
        total_retries += info.retries;
        assert!(info.radius <= 1.0 && info.radius >= 0.01);
        let u = analysis.domain().node(free).displacement[0];
        assert!(u > u_prev);
        assert!((result.load_factor - exp_force(k, a, u)).abs() < 1e-7);
        u_prev = u;
        if u > 2.0 * a {
            break;
        }
    }
    assert!(u_prev > 2.0 * a);
    assert!(
        total_retries > 0,
        "the oversized radius must have provoked retries"
    );
}

#[test]
fn exhausted_cutbacks_restore_everything_and_leave_history_untouched() {
    // Clean reference run.
    let (domain, free) = softening();
    let mut clean = build(
        domain,
        Integrator::ArcLength(softening_config()),
        combined(30),
    );
    let mut clean_path = Vec::new();
    for _ in 0..40 {
        let result = clean.step().unwrap();
        clean_path.push((
            clean.domain().node(free).displacement[0],
            result.load_factor,
        ));
    }

    // Same run, but step 29 (the step that crosses the kink, where the
    // secant predictor is inexact) is first attempted with an unreachable
    // correction tolerance and no retries, so it must fail.
    let (domain, free) = softening();
    let mut faulty = build(
        domain,
        Integrator::ArcLength(softening_config()),
        combined(30),
    );
    for expected in clean_path.iter().take(28) {
        let result = faulty.step().unwrap();
        assert_eq!(result.load_factor, expected.1);
    }
    let mut impossible = softening_config();
    impossible.correction_tolerance = Some(1e-300);
    impossible.max_retries = 0;
    faulty.set_integrator(Integrator::ArcLength(impossible));
    let before = faulty.domain().node(free).displacement[0];
    let error = faulty.step().unwrap_err();
    assert!(
        matches!(error, AnalysisError::CutbacksExhausted { step: 29, .. }),
        "{error:?}"
    );
    assert_eq!(faulty.domain().node(free).displacement[0], before);

    // Back to the original settings (same path definition, so history is
    // kept): the rest of the path is bit-identical to the clean run, which
    // proves the failed attempts left material history, load factor, step
    // count, direction and radius untouched.
    faulty.set_integrator(Integrator::ArcLength(softening_config()));
    for (i, expected) in clean_path.iter().enumerate().skip(28) {
        let result = faulty.step().unwrap();
        assert_eq!(result.step, i + 1);
        assert_eq!(
            (
                faulty.domain().node(free).displacement[0],
                result.load_factor
            ),
            *expected
        );
    }
}

// ---------------------------------------------------------------------------
// Frozen loads, phase transitions, invalid states
// ---------------------------------------------------------------------------

#[test]
fn frozen_gravity_stays_constant_while_lateral_load_is_continued() {
    let (kx, ky, gravity) = (200.0, 500.0, -50.0);
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mass = domain.add_node(Node::new([0.0, 0.0]).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, mass)
            .with_material(0, Material::Elastic { e: kx })
            .with_material(1, Material::Elastic { e: ky }),
    ));
    let gravity_pattern = domain.default_pattern();
    domain.load_node(mass, 1, gravity);
    let lateral = domain.add_load_pattern(LoadSeries::Linear { slope: 1.0 });

    let mut gravity_phase = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(newton())
        .test(combined(10))
        .build(domain);
    let result = gravity_phase.step().unwrap();
    gravity_phase
        .domain_mut()
        .hold_pattern_constant(gravity_pattern, result.load_factor);
    let mut domain = gravity_phase.into_domain();
    domain.add_nodal_load(lateral, mass, 0, 1.0);
    let uy = domain.node(mass).displacement[1];
    assert!((uy - gravity / ky).abs() < 1e-12);

    let mut push = build(
        domain,
        Integrator::ArcLength(ArcLength::fixed(0.1, explicit(0.01, 1.0))),
        combined(10),
    );
    for _ in 0..5 {
        let result = push.step().unwrap();
        let d = push.domain();
        assert_eq!(d.node(mass).displacement[1], uy, "gravity is frozen");
        assert!((kx * d.node(mass).displacement[0] - result.load_factor).abs() < 1e-9);
        assert!((d.reaction(ground, 0, result.load_factor) + result.load_factor).abs() < 1e-9);
        assert!((d.reaction(ground, 1, result.load_factor) + gravity).abs() < 1e-9);
    }
}

#[test]
fn invalid_entry_states_and_load_definitions_are_rejected() {
    let arc = || Integrator::ArcLength(ArcLength::fixed(0.1, explicit(0.01, 1.0)));

    // Not in equilibrium at entry: a constant load with nothing resisting it yet.
    let (mut domain, _fixed, free) = elastic_spring(100.0);
    let constant = domain.add_load_pattern(LoadSeries::Constant);
    domain.add_nodal_load(constant, free, 0, 5.0);
    let mut analysis = build(domain, arc(), combined(10));
    assert!(matches!(
        analysis.step().unwrap_err(),
        AnalysisError::InitialStateNotInEquilibrium { .. }
    ));

    // An active Path series.
    let (mut domain, _fixed, free) = elastic_spring(100.0);
    let path = domain.add_load_pattern(LoadSeries::Path {
        times: vec![0.0, 1.0],
        factors: vec![0.0, 1.0],
    });
    domain.add_nodal_load(path, free, 0, 1.0);
    let mut analysis = build(domain, arc(), combined(10));
    assert_eq!(
        analysis.step().unwrap_err(),
        AnalysisError::UnsupportedLoadSeries
    );

    // No continued load at all.
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, free).with_material(0, Material::Elastic { e: 1.0 }),
    ));
    let mut analysis = build(domain, arc(), combined(10));
    assert_eq!(
        analysis.step().unwrap_err(),
        AnalysisError::ZeroLoadSensitivity
    );
}

#[test]
fn invalid_options_and_unsupported_combinations_are_rejected() {
    let invalid = |field| AnalysisError::InvalidOption { field };
    let run = |config: ArcLength, algorithm: Algorithm, test: ConvergenceTest| {
        let (domain, _fixed, _free) = elastic_spring(100.0);
        AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::ArcLength(config))
            .algorithm(algorithm)
            .test(test)
            .build(domain)
            .step()
            .unwrap_err()
    };
    let base = || ArcLength::fixed(0.1, explicit(0.01, 1.0));

    let mut bad = ArcLength::adaptive(0.1, 0.2, 0.3, explicit(0.01, 1.0));
    assert_eq!(
        run(bad.clone(), newton(), combined(5)),
        invalid("integrator.initialRadius")
    );
    bad = base();
    bad.scales = explicit(f64::NAN, 1.0);
    assert_eq!(
        run(bad, newton(), combined(5)),
        invalid("integrator.displacementScale")
    );
    bad = base();
    bad.arc_tolerance = 0.0;
    assert_eq!(
        run(bad, newton(), combined(5)),
        invalid("integrator.arcTolerance")
    );

    assert_eq!(
        run(base(), Algorithm::Linear, combined(5)),
        invalid("algorithm")
    );
    assert_eq!(
        run(
            base(),
            Algorithm::Newton {
                tangent: TangentStrategy::ReuseAtStepStart,
                line_search: None
            },
            combined(5)
        ),
        invalid("algorithm")
    );
    assert_eq!(
        run(
            base(),
            Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: Some(LineSearch::Bisection {
                    tol: 0.5,
                    max_iter: 5,
                    max_eta: 4.0
                })
            },
            combined(5)
        ),
        invalid("algorithm")
    );
    assert_eq!(
        run(
            base(),
            Algorithm::KrylovNewton {
                tangent: TangentStrategy::Initial,
                max_dimension: 3
            },
            combined(5)
        ),
        invalid("algorithm")
    );
    for legacy in [
        ConvergenceTest::NormDispIncr {
            tol: 1e-9,
            max_iter: 5,
        },
        ConvergenceTest::EnergyIncr {
            tol: 1e-9,
            max_iter: 5,
        },
    ] {
        assert_eq!(run(base(), newton(), legacy), invalid("convergence"));
    }
    // NormUnbalance is accepted (plus the mandatory arc gate).
    let (domain, _fixed, _free) = elastic_spring(100.0);
    let mut ok = build(
        domain,
        Integrator::ArcLength(base()),
        ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        },
    );
    assert!(ok.step().is_ok());

    // A beam model has rotational DOFs, so explicit scales need a rotation scale.
    let mut domain = Domain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let b = domain.add_node(Node::new([1.0, 0.0]));
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        a,
        b,
        1e3,
        1.0,
        1.0,
        GeomTransf::Linear,
    )));
    domain.load_node(b, 1, 1.0);
    let mut beam = build(domain, Integrator::ArcLength(base()), combined(5));
    assert_eq!(
        beam.step().unwrap_err(),
        invalid("integrator.rotationScale")
    );
}

// ---------------------------------------------------------------------------
// Constraints and scale invariance
// ---------------------------------------------------------------------------

#[test]
fn equal_dof_reduced_coordinates_match_the_equivalent_single_node_model() {
    // Two parallel degrading springs to two nodes tied in x is one spring of
    // twice the strength on one node.
    let mut tied = Domain::new();
    let fixed = tied.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let a = tied.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let b = tied.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    tied.equal_dof(a, b, &[0]);
    for node in [a, b] {
        tied.add_element(Element::ZeroLength(
            ZeroLength::new(fixed, node).with_material(0, degrading_spring()),
        ));
    }
    tied.load_node(a, 0, 2.0);
    let (single, free) = softening();

    let mut tied_run = build(
        tied,
        Integrator::ArcLength(softening_config()),
        combined(30),
    );
    let mut single_run = build(
        single,
        Integrator::ArcLength(softening_config()),
        combined(30),
    );
    for _ in 0..50 {
        let t = tied_run.step().unwrap();
        let s = single_run.step().unwrap();
        assert!((t.load_factor - s.load_factor).abs() < 1e-9);
        let (ua, ub) = (
            tied_run.domain().node(a).displacement[0],
            tied_run.domain().node(b).displacement[0],
        );
        assert_eq!(ua, ub);
        assert!((ua - single_run.domain().node(free).displacement[0]).abs() < 1e-12);
    }
}

#[test]
fn unit_and_reference_load_changes_with_matching_scales_recover_the_same_path() {
    // Baseline: metres, unit reference load.
    let (domain, free) = softening();
    let mut base = build(
        domain,
        Integrator::ArcLength(softening_config()),
        combined(30),
    );

    // Millimetres (backbone rotations x1000) with u* = 10 mm, and a
    // reference load of 4 with lambda* = 10/4.
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mm = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, mm).with_material(
            0,
            Material::hysteretic(
                10.0, 10.0, 6.0, 20.0, 2.0, 30.0, -10.0, -10.0, -6.0, -20.0, -2.0, -30.0, 1.0, 1.0,
                0.0, 0.0, 0.0,
            ),
        ),
    ));
    domain.load_node(mm, 0, 4.0);
    let mut config = softening_config();
    config.scales = explicit(10.0, 2.5);
    let mut scaled = build(
        domain,
        Integrator::ArcLength(config),
        ConvergenceTest::Combined {
            force_tol: 1e-9,
            moment_tol: 1e-9,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 30,
        },
    );

    for _ in 0..50 {
        let b = base.step().unwrap();
        let s = scaled.step().unwrap();
        assert!((b.load_factor - 4.0 * s.load_factor).abs() < 1e-8);
        let (u_m, u_mm) = (
            base.domain().node(free).displacement[0],
            scaled.domain().node(mm).displacement[0],
        );
        assert!((u_m - u_mm / 1000.0).abs() < 1e-12, "{u_m} vs {u_mm} mm");
    }
}

#[test]
fn domain_mut_discards_history_and_revalidates() {
    let (domain, _fixed, free) = elastic_spring(100.0);
    let mut analysis = build(
        domain,
        Integrator::ArcLength(ArcLength::fixed(0.1, explicit(0.01, 1.0))),
        combined(10),
    );
    analysis.step().unwrap();
    // A load added mid-phase breaks equilibrium; continuation must notice.
    let extra = analysis.domain_mut().add_load_pattern(LoadSeries::Constant);
    analysis.domain_mut().add_nodal_load(extra, free, 0, 3.0);
    assert!(matches!(
        analysis.step().unwrap_err(),
        AnalysisError::InitialStateNotInEquilibrium { .. }
    ));
}
