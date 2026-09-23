use carapace_core::analysis::{GroundMotion, RayleighDamping, TransientAnalysis};
use carapace_core::model::{Domain, Element, LoadSeries, Material, Node, NodeId, ZeroLength};

fn build_sdof(m: f64, k: f64) -> (Domain, NodeId) {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mass_node = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2).with_mass(0, m));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, mass_node).with_material(0, Material::Elastic { e: k }),
    ));
    (domain, mass_node)
}

/// Undamped SDOF (`m=1, k=1, omega=1`), at rest, under a *constant* ground
/// acceleration `ag0` starting at `t=0` (`LoadSeries::Constant` — `factor(t)
/// = 1.0` for every `t`, so `GroundMotion::acceleration` is just `ag0`
/// throughout). Standard step-response closed form (e.g. Chopra):
/// `u(t) = -(ag0/omega^2) * (1 - cos(omega*t))`. Same closed-form-vs-Newmark
/// comparison style as `m6_dynamics.rs`'s `newmark_undamped_sdof_matches_
/// closed_form_free_vibration`, but forced (nonzero right-hand side) rather
/// than free vibration — isolates that `GroundMotion`'s effective force
/// term (`-M*ι*ag(t)`) is wired correctly into both the initial-
/// acceleration equilibrium and every Newmark step's `f_eff`.
#[test]
fn constant_ground_acceleration_matches_undamped_step_response_closed_form() {
    let (m, k): (f64, f64) = (1.0, 1.0);
    let omega: f64 = (k / m).sqrt();
    let ag0 = 2.0;

    let (domain, mass_node) = build_sdof(m, k);
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, 0.01)
        .expect("undamped SDOF should be well-posed")
        .with_ground_motion(GroundMotion::new(0, LoadSeries::Constant).with_scale_factor(ag0));

    let steps = 100;
    for _ in 0..steps {
        analysis.step().expect("forced undamped SDOF should solve every step");
    }

    let t = 0.01 * steps as f64;
    let expected = -(ag0 / (omega * omega)) * (1.0 - (omega * t).cos());
    let u = analysis.domain().node(mass_node).displacement[0];
    assert!((u - expected).abs() < 1e-3, "expected u({t})={expected}, got {u}");
}

/// `GroundMotion`'s effective force (`-M*ι*ag(t)`, mass-proportional) must
/// be mathematically identical to applying an ordinary nodal `LoadPattern`
/// whose reference load at the mass DOF is `-m` and whose series *is*
/// `ag(t)` directly (`load(t) = series.factor(t) * (-m) = -m*ag(t)`) — the
/// same physics, two different code paths (`GroundMotion`'s `direction_incidence`
/// vs. an ordinary pattern). Checked here with a `LoadSeries::Path`
/// (piecewise-linear — an accelerogram shape), rather than re-deriving a
/// closed form for a forced-vibration response to a nontrivial pulse: if
/// both paths agree to near machine precision at every step, `GroundMotion`
/// is correctly wired.
#[test]
fn ground_motion_matches_an_equivalent_mass_proportional_load_pattern() {
    let (m, k) = (2.0, 5.0);
    let ag_path = LoadSeries::Path {
        times: vec![0.0, 0.1, 0.2, 0.3, 0.5],
        factors: vec![0.0, 3.0, -1.5, 0.8, 0.8],
    };

    let (domain_gm, node_gm) = build_sdof(m, k);
    let mut analysis_gm = TransientAnalysis::new(domain_gm, RayleighDamping::NONE, 0.01)
        .expect("SDOF should be well-posed")
        .with_ground_motion(GroundMotion::new(0, ag_path.clone()));

    let (mut domain_pattern, node_pattern) = build_sdof(m, k);
    let pattern = domain_pattern.add_load_pattern(ag_path);
    domain_pattern.add_nodal_load(pattern, node_pattern, 0, -m);
    let mut analysis_pattern =
        TransientAnalysis::new(domain_pattern, RayleighDamping::NONE, 0.01).expect("SDOF should be well-posed");

    for _ in 0..60 {
        let r_gm = analysis_gm.step().expect("ground-motion-driven SDOF should solve every step");
        let r_pattern = analysis_pattern
            .step()
            .expect("equivalent-load-pattern SDOF should solve every step");
        assert_eq!(r_gm.step, r_pattern.step);

        let u_gm = analysis_gm.domain().node(node_gm).displacement[0];
        let u_pattern = analysis_pattern.domain().node(node_pattern).displacement[0];
        assert!(
            (u_gm - u_pattern).abs() < 1e-9,
            "step {}: GroundMotion (u={u_gm}) should exactly match the equivalent load pattern (u={u_pattern})",
            r_gm.step
        );
    }
}
