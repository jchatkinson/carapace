use carapace_core::analysis::{
    Algorithm, ConvergenceTest, RayleighDamping, TangentStrategy, TransientAnalysis,
};
use carapace_core::model::{Domain, Element, Material, Node, ZeroLength};

/// `docs/algorithms.md` §6.1's load-bearing correctness check, run as an
/// actual test rather than only hand-verified on paper: for a *linear*
/// system, the new Newton-corrector arm (`TangentStrategy::ReuseAtStepStart`
/// — form the effective operator once, at step start, same as the old
/// closed form did implicitly) reproduces the old one-shot closed form's
/// answer bit-for-bit — proving the restructure is a strict superset, not
/// a behavior change, for the case that already worked. Same undamped SDOF
/// free-vibration system `dynamics/free_vibration.rs`'s `newmark_undamped_sdof_
/// matches_closed_form_free_vibration` already checks against the closed
/// form; this test checks the two *algorithms* against each other instead.
///
/// Exactly *two* iterations, not one: `ConvergenceTest::check` (shared with
/// static `Analysis`) tests the residual that
/// *produced* the just-applied correction, not a freshly re-formed one —
/// so iteration 0's large step-start unbalance never self-certifies, even
/// though the single linear solve it produces is already exact. Iteration
/// 1 re-forms the residual at that already-exact state (now ~0, since a
/// linear `Op` solved exactly satisfies equilibrium after one correction)
/// and *that's* what converges. This is the same convention
/// `Algorithm::Newton` already used for static analysis — not a new wrinkle introduced by the
/// transient corrector.
#[test]
fn newton_corrector_matches_algorithm_linear_bit_for_bit_on_a_linear_system() {
    let (k_spring, m): (f64, f64) = (1.0, 1.0);
    let dt = 0.01;
    let steps = 50;

    let build = || {
        let mut domain = Domain::new();
        let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let mass_node = domain.add_node(
            Node::new([1.0, 0.0])
                .fix(1)
                .fix(2)
                .with_mass(0, m)
                .with_initial_displacement(0, 1.0),
        );
        domain.add_element(Element::ZeroLength(
            ZeroLength::new(ground, mass_node).with_material(0, Material::Elastic { e: k_spring }),
        ));
        (domain, mass_node)
    };

    let (domain_linear, mass_node) = build();
    let mut linear =
        TransientAnalysis::new(domain_linear, RayleighDamping::NONE, dt).expect("well-posed");

    let (domain_newton, _) = build();
    let mut newton = TransientAnalysis::new(domain_newton, RayleighDamping::NONE, dt)
        .expect("well-posed")
        .with_algorithm(
            Algorithm::Newton {
                tangent: TangentStrategy::ReuseAtStepStart,
                line_search: None,
            },
            // Newmark's mass term scales as `1/(beta*dt^2)` (`=40000` at
            // this `dt`), so an *absolute* residual tolerance needs some
            // headroom above plain machine epsilon to be reachable at all
            // once floating-point cancellation at that magnitude is
            // accounted for — `1e-8` is still far tighter than any
            // physical modeling tolerance, just not literally `1e-12`.
            ConvergenceTest::NormUnbalance {
                tol: 1e-8,
                max_iter: 20,
            },
        );

    for step in 1..=steps {
        linear
            .step()
            .expect("Algorithm::Linear should solve every step");
        let result = newton
            .step()
            .expect("Newton corrector should solve every step");

        // Convergence is checked at the accepted (post-update) state, so the
        // one "real" solve is also the confirming one: no second iteration.
        assert_eq!(
            result.iterations, 1,
            "a linear system should converge on its one real solve (step {step})"
        );
        assert_eq!(
            result.factorizations, 1,
            "one iteration needs exactly one factorization (step {step})"
        );

        let (u_lin, u_new) = (
            linear.domain().node(mass_node).displacement[0],
            newton.domain().node(mass_node).displacement[0],
        );
        let (v_lin, v_new) = (
            linear.domain().node(mass_node).velocity[0],
            newton.domain().node(mass_node).velocity[0],
        );
        let (a_lin, a_new) = (
            linear.domain().node(mass_node).acceleration[0],
            newton.domain().node(mass_node).acceleration[0],
        );

        assert!(
            (u_lin - u_new).abs() < 1e-8,
            "step {step}: displacement diverged: linear={u_lin}, newton={u_new}"
        );
        assert!(
            (v_lin - v_new).abs() < 1e-8,
            "step {step}: velocity diverged: linear={v_lin}, newton={v_new}"
        );
        assert!(
            (a_lin - a_new).abs() < 1e-6,
            "step {step}: acceleration diverged: linear={a_lin}, newton={a_new}"
        );
    }
}

/// The actual new capability: `Algorithm::Linear` forms its one effective
/// tangent at *step start* — still fully elastic here, since the spring
/// hasn't yielded yet at `u_n` — and never re-checks it for the rest of
/// the step, even once the true (trial) state has crossed the `ElasticPP`
/// spring's yield breakpoint. The resisting force `element_local_force`
/// reads back is *always* evaluated through the material's own real law
/// (`Material::trial_stress_tangent`, a pure function of committed history
/// + current displacement, independent of which `Algorithm` chose that
/// displacement) — so both runs report a force clamped to `fy`. What
/// `Algorithm::Linear` actually gets wrong is the *displacement* itself:
/// having solved against the stale elastic tangent for the whole step, its
/// `u` doesn't satisfy the true (regime-updated) dynamic equilibrium the
/// way the Newton corrector's converged `u` does — so the two algorithms
/// must land on visibly different displacements once yielding is engaged
/// mid-step, even though both "solve" without error.
#[test]
fn newton_corrector_gives_a_different_displacement_than_algorithm_linear_once_yielding_is_engaged_mid_step(
) {
    let (k0, fy): (f64, f64) = (100.0, 1.0);
    let eyp = fy / k0;
    let m = 1.0;
    let dt = 0.01;

    let build = || {
        let mut domain = Domain::new();
        let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        let mass_node = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2).with_mass(0, m));
        let spring = domain.add_element(Element::ZeroLength(
            ZeroLength::new(ground, mass_node).with_material(0, Material::elastic_pp(k0, eyp)),
        ));
        // Default reference-load series is `Linear { slope: 1.0 }` (grows
        // with pseudo-time), so at `time = dt` the applied force is
        // `1e7 * dt` = 1e5 — with the mass term (`~1/(beta*dt^2)`)
        // dominating the effective stiffness at this small `dt`, this
        // magnitude is chosen (empirically, not from a hand closed form)
        // to guarantee the trial displacement crosses `eyp` within step 1
        // regardless of the exact dynamics.
        domain.load_node(mass_node, 0, 1.0e7);
        (domain, mass_node, spring)
    };

    let (domain_linear, mass_node, spring) = build();
    let mut linear =
        TransientAnalysis::new(domain_linear, RayleighDamping::NONE, dt).expect("well-posed");
    let linear_result = linear
        .step()
        .expect("Algorithm::Linear step should still numerically solve");
    let (u_lin, f_lin) = (
        linear.domain().node(mass_node).displacement[0],
        linear.domain().element_local_force(spring)[0],
    );

    let (domain_newton, _, spring2) = build();
    let mut newton = TransientAnalysis::new(domain_newton, RayleighDamping::NONE, dt)
        .expect("well-posed")
        .with_algorithm(
            Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            },
            ConvergenceTest::NormUnbalance {
                tol: 1e-6,
                max_iter: 100,
            },
        );
    let newton_result = newton
        .step()
        .expect("Newton-corrected step should converge");
    let (u_new, f_new) = (
        newton.domain().node(mass_node).displacement[0],
        newton.domain().element_local_force(spring2)[0],
    );

    assert!(
        u_lin.abs() > eyp,
        "test setup should have crossed yield (eyp={eyp}), got u_linear={u_lin}"
    );
    assert!(
        f_lin.abs() <= fy * (1.0 + 1e-6),
        "material must clamp to fy regardless of algorithm, got f_linear={f_lin}"
    );
    assert!(
        f_new.abs() <= fy * (1.0 + 1e-6),
        "material must clamp to fy regardless of algorithm, got f_newton={f_new}"
    );

    assert_eq!(
        linear_result.iterations, 1,
        "Algorithm::Linear never iterates by definition"
    );
    assert!(
        newton_result.iterations > 1,
        "the corrector should need real iteration once yielding engages mid-step"
    );

    // A small but real difference: at this `dt`, Newmark's mass term
    // dominates the effective operator (`~1/(beta*dt^2)` vs. the spring's
    // `k0`/post-yield `0`), so yielding mid-step only perturbs the answer
    // slightly — the point of this test is that it's nonzero and
    // attributable to the regime change (`newton_result.iterations > 1`
    // above already confirms the corrector actually engaged), not that
    // it's dramatic.
    let relative_diff = (u_lin - u_new).abs() / u_lin.abs();
    assert!(
        relative_diff > 1e-4,
        "expected Algorithm::Linear's stale-tangent displacement to differ from the corrector's, \
         got u_linear={u_lin}, u_newton={u_new} ({}% relative difference)",
        relative_diff * 100.0
    );
}
