//! `docs/arclength.md` §8: arc-length continuation through a degrading-
//! strength peak (one `Hysteretic` spring) and through an actual snap-back
//! (the same spring in series with an elastic one), against closed-form
//! backbone oracles at every accepted step.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ArcLength, ArcScales, ArcStepInfo,
    ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, NodeId, ZeroLength};

const U_STAR: f64 = 0.01;
const LAMBDA_STAR: f64 = 10.0;
const RADIUS: f64 = 0.05;

/// Monotonic positive backbone of the degrading spring.
fn backbone(u: f64) -> f64 {
    if u <= 0.01 {
        1000.0 * u
    } else {
        10.0 - 400.0 * (u - 0.01)
    }
}

fn degrading_spring() -> Material {
    Material::hysteretic(
        10.0, 0.01, 6.0, 0.02, 2.0, 0.03, -10.0, -0.01, -6.0, -0.02, -2.0, -0.03, 1.0, 1.0, 0.0,
        0.0, 0.0,
    )
}

fn combined() -> ConvergenceTest {
    ConvergenceTest::Combined {
        force_tol: 1e-9,
        moment_tol: 1e-9,
        relative_tol: 1e-9,
        displacement_tol: None,
        max_iter: 30,
    }
}

fn fixture_arc() -> ArcLength {
    let mut arc = ArcLength::fixed(
        RADIUS,
        ArcScales::Explicit {
            displacement: U_STAR,
            rotation: None,
            load: LAMBDA_STAR,
        },
    );
    arc.arc_tolerance = 1e-9;
    arc
}

fn build(domain: Domain, integrator: Integrator) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(integrator)
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(combined())
        .build(domain)
}

/// Two coincident nodes, one fully fixed, the other free only in x, a unit
/// reference force on the free x DOF.
fn single_spring() -> (Domain, NodeId, NodeId) {
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, free).with_material(0, degrading_spring()),
    ));
    domain.load_node(free, 0, 1.0);
    (domain, fixed, free)
}

fn info(result: &carapace_core::analysis::StepResult) -> ArcStepInfo {
    result
        .arc
        .expect("arc steps report continuation diagnostics")
}

/// Arc-length continuation traces a degrading-strength spring through its peak and down the
/// softening branch. Every accepted point lies on the independent backbone curve (lambda =
/// backbone(u)), the support reaction balances the applied load, and the load decreases after the
/// peak.
#[test]
fn degrading_spring_is_traced_past_peak_strength() {
    let (domain, fixed, free) = single_spring();
    let mut analysis = build(domain, Integrator::ArcLength(fixture_arc()));

    let (mut u_prev, mut lambda_prev) = (0.0, 0.0);
    let mut lambdas = Vec::new();
    let mut us = Vec::new();
    let mut post_peak_decreases = 0;
    for _ in 0..100 {
        let result = analysis
            .step()
            .expect("arc length should trace the softening branch");
        let arc = info(&result);
        let u = analysis.domain().node(free).displacement[0];
        let lambda = result.load_factor;

        assert!(u.is_finite() && lambda.is_finite());
        assert!(arc.force_measure <= 1.0 && arc.arc_error <= 1e-9);
        // Equilibrium against the independent backbone oracle, and the
        // support reaction balances the applied load.
        assert!(
            (lambda - backbone(u)).abs() <= 1e-6,
            "lambda {lambda} vs backbone({u}) = {}",
            backbone(u)
        );
        let reaction = analysis.domain().reaction(fixed, 0, lambda);
        assert!(
            (reaction + lambda).abs() <= 1e-6,
            "reaction {reaction} vs {lambda}"
        );
        // The measured chord (total step increments) is on the sphere.
        let chord =
            ((u - u_prev) / U_STAR).powi(2) + ((lambda - lambda_prev) / LAMBDA_STAR).powi(2);
        assert!(
            (chord - RADIUS * RADIUS).abs() / (RADIUS * RADIUS) <= 1e-8,
            "chord {chord}"
        );
        assert!(
            u > u_prev,
            "u must increase through the peak: {u_prev} -> {u}"
        );
        if u_prev > 0.01 && lambda < lambda_prev {
            post_peak_decreases += 1;
        }

        us.push(u);
        lambdas.push(lambda);
        (u_prev, lambda_prev) = (u, lambda);
        if u >= 0.028 {
            break;
        }
    }

    let final_u = *us.last().unwrap();
    let final_lambda = *lambdas.last().unwrap();
    assert!((0.028..0.029).contains(&final_u), "final u {final_u}");
    assert!(final_lambda < 3.0, "final force {final_lambda}");
    assert!(post_peak_decreases >= 10);
    // An interior maximum with samples on both sides of the peak, within
    // the radius-bounded undersampling of 1000 * u* * s = 0.5.
    let peak = lambdas.iter().cloned().fold(f64::MIN, f64::max);
    assert!(us.iter().any(|&u| u < 0.01) && us.iter().any(|&u| u > 0.01));
    assert!(peak > 10.0 - 0.5 && peak <= 10.0 + 1e-9);

    // Landmarks, by bracketing interpolation (the path needn't sample them).
    for (u_mark, f_mark) in [(0.01, 10.0), (0.02, 6.0), (0.028, 2.8)] {
        let i = us
            .iter()
            .position(|&u| u >= u_mark)
            .expect("path passes landmark");
        if i == 0 {
            continue;
        }
        let (u0, u1, l0, l1) = (us[i - 1], us[i], lambdas[i - 1], lambdas[i]);
        let interpolated = l0 + (l1 - l0) * (u_mark - u0) / (u1 - u0);
        // Linear interpolation across the kink at 0.01 can undershoot.
        let tolerance = if u_mark == 0.01 { 0.5 } else { 1e-6 };
        assert!(
            (interpolated - f_mark).abs() <= tolerance,
            "landmark u={u_mark}: {interpolated} vs {f_mark}"
        );
    }
}

/// Baseline for the arc-length tests: load control converges before the peak but fails to converge
/// past it, and the failed step rolls the node displacement back exactly.
#[test]
fn load_control_fails_past_the_peak_and_rolls_back_exactly() {
    let (domain, _fixed, free) = single_spring();
    let mut analysis = build(domain, Integrator::LoadControl { increment: 2.0 });
    for _ in 0..5 {
        analysis.step().expect("pre-peak load control converges");
    }
    let before = analysis.domain().node(free).displacement[0];
    let result = analysis.step();
    assert!(matches!(
        result,
        Err(AnalysisError::FailedToConverge { .. })
    ));
    assert_eq!(analysis.domain().node(free).displacement[0], before);
}

/// The degrading spring in series with a stiff elastic spring produces a true snap-back (both load
/// and displacement decrease after the peak). Arc-length continuation traces it, including the step
/// that crosses the backbone kink.
#[test]
fn series_spring_snap_back_is_traced_with_decreasing_load_and_displacement() {
    const K_S: f64 = 380.0;
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let x_node = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let y_node = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, x_node).with_material(0, degrading_spring()),
    ));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(x_node, y_node).with_material(0, Material::Elastic { e: K_S }),
    ));
    domain.load_node(y_node, 0, 1.0);
    let mut analysis = build(domain, Integrator::ArcLength(fixture_arc()));

    let (mut x_prev, mut y_prev, mut lambda_prev) = (0.0, 0.0, 0.0);
    let mut post_peak_steps = 0;
    let mut saw_kink_step = false;
    for _ in 0..250 {
        let result = analysis
            .step()
            .expect("arc length should trace the snap-back");
        let arc = info(&result);
        let x = analysis.domain().node(x_node).displacement[0];
        let y = analysis.domain().node(y_node).displacement[0];
        let lambda = result.load_factor;

        // Series oracle: both springs carry lambda.
        assert!(
            (lambda - backbone(x)).abs() <= 1e-6,
            "spring force at x={x}"
        );
        let y_oracle = x + backbone(x) / K_S;
        assert!((y - y_oracle).abs() <= 1e-9, "y {y} vs {y_oracle}");
        assert!((analysis.domain().reaction(fixed, 0, lambda) + lambda).abs() <= 1e-6);
        let chord = ((x - x_prev) / U_STAR).powi(2)
            + ((y - y_prev) / U_STAR).powi(2)
            + ((lambda - lambda_prev) / LAMBDA_STAR).powi(2);
        assert!((chord - RADIUS * RADIUS).abs() / (RADIUS * RADIUS) <= 1e-8);
        assert!(x > x_prev, "x keeps increasing through the snap-back");
        // The step straddling the 84.4-degree kink must not need cutbacks.
        if x_prev < 0.01 && x > 0.01 {
            saw_kink_step = true;
            assert_eq!(arc.retries, 0);
        }
        if x_prev > 0.01 {
            post_peak_steps += 1;
            assert!(lambda < lambda_prev, "load decreases after the peak");
            assert!(
                y < y_prev,
                "loaded displacement decreases: actual snap-back"
            );
        }
        (x_prev, y_prev, lambda_prev) = (x, y, lambda);
        if x >= 0.028 {
            break;
        }
    }
    assert!(saw_kink_step);
    assert!(post_peak_steps >= 10);
    assert!(x_prev >= 0.028);
    // Landmarks from the exact series formula.
    for (x_mark, y_mark, l_mark) in [
        (0.01, 0.01 + 10.0 / K_S, 10.0),
        (0.02, 0.02 + 6.0 / K_S, 6.0),
        (0.028, 0.028 + 2.8 / K_S, 2.8),
    ] {
        assert!((y_mark - (x_mark + backbone(x_mark) / K_S)).abs() < 1e-15);
        assert!((backbone(x_mark) - l_mark).abs() < 1e-12);
    }
}

/// Starting arc-length in the decreasing direction follows the mirror-image negative backbone:
/// displacement decreases monotonically and every point satisfies lambda = -backbone(-u).
#[test]
fn reverse_initial_direction_follows_the_negative_backbone() {
    let (domain, _fixed, free) = single_spring();
    let mut arc = fixture_arc();
    arc.direction = carapace_core::analysis::ArcDirection::Decreasing;
    let mut analysis = build(domain, Integrator::ArcLength(arc));
    let mut u_prev = 0.0;
    for _ in 0..100 {
        let result = analysis.step().expect("negative branch traced");
        let u = analysis.domain().node(free).displacement[0];
        assert!(u < u_prev);
        assert!((result.load_factor + backbone(-u)).abs() <= 1e-6);
        u_prev = u;
        if u <= -0.028 {
            break;
        }
    }
    assert!(u_prev <= -0.028);
}
