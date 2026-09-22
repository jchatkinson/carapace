use carapace_core::analysis::{modal_analysis, RayleighDamping, TransientAnalysis};
use carapace_core::model::{Domain, Element, GeomTransf, Material, Node, Truss, ZeroLength};

/// M6 acceptance, part 1 (implementation-plan §6): element-consistent
/// lumped mass. A single `Truss` with a nonzero `density`, fixed at one
/// end, has its total mass (`density*area*length`) split half to each
/// node's translational DOFs (`Element::form_mass`) — since the fixed
/// node's half never enters the free-DOF mass diagonal, the free node
/// carries exactly half the truss's total mass. This is a genuine SDOF
/// (stiffness `k = E*A/L`, mass `m = density*A*L/2`), so its fundamental
/// frequency has an exact closed form: `omega = sqrt(k/m)`.
#[test]
fn truss_element_mass_matches_sdof_closed_form_frequency() {
    let (e, area, length, density) = (30000.0, 2.0, 100.0, 0.5);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.add_element(Element::Truss(
        Truss::new(node_i, node_j, area, Material::Elastic { e }).with_density(density),
    ));

    let modes = modal_analysis(&mut domain, 1).expect("SDOF truss should have a well-posed eigenproblem");

    let k = e * area / length;
    let m = density * area * length / 2.0;
    let expected = (k / m).sqrt();
    assert!(
        (modes[0].frequency - expected).abs() < 1e-9,
        "expected omega={expected}, got {}",
        modes[0].frequency
    );
}

/// Same check for `ElasticBeamColumn`'s axial mode (transverse bending DOFs
/// are fixed here so this exercises the same `k = E*A/L` axial stiffness
/// against the beam's `form_mass`, independent of the truss test above).
#[test]
fn beam_element_mass_matches_sdof_closed_form_frequency() {
    let (e, area, iz, length, density) = (30000.0, 2.0, 1000.0, 100.0, 0.5);

    let mut domain = Domain::new();
    let node_i = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node_j = domain.add_node(Node::new([length, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ElasticBeamColumn(
        carapace_core::model::ElasticBeamColumn::new(node_i, node_j, e, area, iz, GeomTransf::Linear)
            .with_density(density),
    ));

    let modes = modal_analysis(&mut domain, 1).expect("SDOF beam should have a well-posed eigenproblem");

    let k = e * area / length;
    let m = density * area * length / 2.0;
    let expected = (k / m).sqrt();
    assert!(
        (modes[0].frequency - expected).abs() < 1e-9,
        "expected omega={expected}, got {}",
        modes[0].frequency
    );
}

/// M6 acceptance, part 2: Newmark + Rayleigh damping against the classical
/// damped-SDOF free-vibration closed form. `ZeroLength`+`Elastic` (k=1)
/// between a fixed ground and a node carrying `Node::mass = 1` (m=1, so
/// omega=1 rad/s), mass-proportional-only Rayleigh damping
/// (`alpha_m = 2*xi*omega`, `beta_k = 0`) for an exact damping ratio
/// `xi = 5%`, released from `u0 = 1`, `v0 = 0`.
///
/// Closed form: `u(t) = exp(-xi*omega*t) * u0 * [cos(omega_d*t) +
/// (xi*omega/omega_d)*sin(omega_d*t)]`, `omega_d = omega*sqrt(1-xi^2)`.
///
/// Newmark's average-acceleration scheme is second-order accurate, not
/// exact, so the tolerance here (1e-4) reflects real (small, at this
/// dt/period ratio) discretization error — unlike the milestones through
/// M5, which check exact closed forms to solver precision.
#[test]
fn newmark_rayleigh_damped_sdof_matches_closed_form_free_vibration() {
    let (k_spring, m, xi): (f64, f64, f64) = (1.0, 1.0, 0.05);
    let omega = (k_spring / m).sqrt();
    let alpha_m = 2.0 * xi * omega;

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

    let dt = 0.01;
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::new(alpha_m, 0.0), dt)
        .expect("SDOF with assigned mass should have a well-posed transient system");

    let steps = 100;
    for _ in 0..steps {
        analysis.step().expect("linear damped SDOF should solve every step");
    }

    let t = dt * steps as f64;
    let omega_d = omega * (1.0 - xi * xi).sqrt();
    let expected = (-xi * omega * t).exp() * (omega_d * t).cos()
        + (-xi * omega * t).exp() * (xi * omega / omega_d) * (omega_d * t).sin();

    let u = analysis.domain().node(mass_node).displacement[0];
    assert!(
        (u - expected).abs() < 1e-4,
        "expected u({t})={expected}, got {u}"
    );
}

/// Zero damping (`RayleighDamping::NONE`) should match undamped SDOF free
/// vibration `u(t) = u0*cos(omega*t)` exactly (to Newmark discretization
/// error) — a sanity check that `RayleighDamping::NONE` really contributes
/// nothing, isolated from the damped case above.
#[test]
fn newmark_undamped_sdof_matches_closed_form_free_vibration() {
    let (k_spring, m): (f64, f64) = (1.0, 1.0);
    let omega = (k_spring / m).sqrt();

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

    let dt = 0.01;
    let mut analysis =
        TransientAnalysis::new(domain, RayleighDamping::NONE, dt).expect("undamped SDOF should be well-posed");

    let steps = 100;
    for _ in 0..steps {
        analysis.step().expect("undamped SDOF should solve every step");
    }

    let t = dt * steps as f64;
    let expected = (omega * t).cos();
    let u = analysis.domain().node(mass_node).displacement[0];
    assert!((u - expected).abs() < 1e-4, "expected u({t})={expected}, got {u}");
}

