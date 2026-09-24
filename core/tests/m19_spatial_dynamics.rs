use carapace_core::analysis::{modal_analysis, GroundMotion, RayleighDamping, TransientAnalysis};
use carapace_core::model::{Domain3, Element3, LoadSeries, Material, Node3, SpatialDof, Truss3, ZeroLength3};

/// Spatial counterpart to `m6_dynamics.rs`'s `truss_element_mass_matches_
/// sdof_closed_form_frequency` — `modal_analysis`/`TransientAnalysis`
/// generalized to the same `<NDIM, NDOF, ELEMENT_DOF, NId, E>` profile as
/// `Domain`/`Analysis` (see `spatial-architecture.md`'s "Spatial dynamics"
/// milestone), run through the real `Domain3`/`Analysis3` stack rather than
/// asserted by analogy. Every non-axial DOF at the free node is fixed, so
/// this is a genuine SDOF (`k = E*A/L`, `m = density*A*L/2`, half the
/// truss's total mass — the fixed node's half never enters the free-DOF
/// mass diagonal) with the exact same closed form as the planar test.
#[test]
fn spatial_truss_mass_matches_sdof_closed_form_frequency_via_modal_analysis3() {
    let (e, area, length, density) = (30_000.0, 2.0, 100.0, 0.5);

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
    let node_j = domain.add_node(
        Node3::new([length, 0.0, 0.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    domain.add_element(Element3::Truss3(
        Truss3::new(node_i, node_j, area, Material::Elastic { e }).with_density(density),
    ));

    let modes = modal_analysis(&mut domain, 1).expect("spatial SDOF truss should have a well-posed eigenproblem");

    let k = e * area / length;
    let m = density * area * length / 2.0;
    let expected = (k / m).sqrt();
    assert!(
        (modes[0].frequency - expected).abs() < 1e-9,
        "expected omega={expected}, got {}",
        modes[0].frequency
    );
}

/// Spatial counterpart to `m6_dynamics.rs`'s `newmark_undamped_sdof_
/// matches_closed_form_free_vibration` — `TransientAnalysis3` (Newmark
/// stepping through `Domain3`) against the same undamped free-vibration
/// closed form, `u(t) = u0*cos(omega*t)`, released from rest at `u0=1`.
/// Every DOF but `ux` at the mass node is fixed, isolating one spatial
/// translational DOF exactly like the planar SDOF case.
#[test]
fn spatial_zero_length_sdof_matches_closed_form_free_vibration_via_transient_analysis3() {
    let (k_spring, m): (f64, f64) = (1.0, 1.0);
    let omega = (k_spring / m).sqrt();

    let mut domain = Domain3::new();
    let ground = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    let mass_node = domain.add_node(
        Node3::new([1.0, 0.0, 0.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Ux as usize, m)
            .with_initial_displacement(SpatialDof::Ux as usize, 1.0),
    );
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, mass_node).with_material(SpatialDof::Ux as usize, Material::Elastic { e: k_spring }),
    ));

    let dt = 0.01;
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, dt)
        .expect("spatial undamped SDOF should be well-posed");

    let steps = 100;
    for _ in 0..steps {
        analysis.step().expect("spatial undamped SDOF should solve every step");
    }

    let t = dt * steps as f64;
    let expected = (omega * t).cos();
    let u = analysis.domain().node(mass_node).displacement[SpatialDof::Ux as usize];
    assert!((u - expected).abs() < 1e-4, "expected u({t})={expected}, got {u}");
}

/// Spatial Rayleigh-damped counterpart — same closed form and tolerance as
/// `m6_dynamics.rs`'s `newmark_rayleigh_damped_sdof_matches_closed_form_
/// free_vibration`, run on a single spatial translational DOF (`uy` this
/// time, not `ux`, so the two Newmark tests above don't both happen to
/// exercise the same DOF index).
#[test]
fn spatial_rayleigh_damped_sdof_matches_closed_form_free_vibration_via_transient_analysis3() {
    let (k_spring, m, xi): (f64, f64, f64) = (1.0, 1.0, 0.05);
    let omega = (k_spring / m).sqrt();
    let alpha_m = 2.0 * xi * omega;

    let mut domain = Domain3::new();
    let ground = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    let mass_node = domain.add_node(
        Node3::new([0.0, 1.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Uy as usize, m)
            .with_initial_displacement(SpatialDof::Uy as usize, 1.0),
    );
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, mass_node).with_material(SpatialDof::Uy as usize, Material::Elastic { e: k_spring }),
    ));

    let dt = 0.01;
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::new(alpha_m, 0.0), dt)
        .expect("spatial damped SDOF should be well-posed");

    let steps = 100;
    for _ in 0..steps {
        analysis.step().expect("spatial damped SDOF should solve every step");
    }

    let t = dt * steps as f64;
    let omega_d = omega * (1.0 - xi * xi).sqrt();
    let expected = (-xi * omega * t).exp()
        * ((omega_d * t).cos() + (xi * omega / omega_d) * (omega_d * t).sin());

    let u = analysis.domain().node(mass_node).displacement[SpatialDof::Uy as usize];
    assert!((u - expected).abs() < 1e-4, "expected u({t})={expected}, got {u}");
}

/// "3D ground motion": two independent `GroundMotion`s active at once, one
/// per translational direction (`ux`, `uz`), driving one spatial mass node
/// whose spring stiffnesses differ per direction and which has no coupling
/// between them (independent `ZeroLength3` materials per DOF — see
/// `ZeroLength3`'s doc comment). Because the assembled system is diagonal,
/// each direction's response must match the *same* single-direction
/// step-response closed form `m_ground_motion.rs` already verifies
/// (`u(t) = -(ag0/omega^2)*(1-cos(omega*t))`), independently and
/// simultaneously — the spatial generalization of "ground motion in more
/// than one direction doesn't cross-talk between DOFs".
#[test]
fn orthogonal_ground_motions_drive_independent_spatial_sdof_responses() {
    let (m_x, k_x, ag_x): (f64, f64, f64) = (1.0, 1.0, 2.0);
    let (m_z, k_z, ag_z): (f64, f64, f64) = (2.0, 8.0, -1.5);
    let omega_x = (k_x / m_x).sqrt();
    let omega_z = (k_z / m_z).sqrt();

    let mut domain = Domain3::new();
    let ground = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    let mass_node = domain.add_node(
        Node3::new([1.0, 0.0, 1.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Ux as usize, m_x)
            .with_mass(SpatialDof::Uz as usize, m_z),
    );
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, mass_node)
            .with_material(SpatialDof::Ux as usize, Material::Elastic { e: k_x })
            .with_material(SpatialDof::Uz as usize, Material::Elastic { e: k_z }),
    ));

    let dt = 0.01;
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, dt)
        .expect("orthogonal spatial SDOF pair should be well-posed")
        .with_ground_motion(GroundMotion::new(SpatialDof::Ux as usize, LoadSeries::Constant).with_scale_factor(ag_x))
        .with_ground_motion(GroundMotion::new(SpatialDof::Uz as usize, LoadSeries::Constant).with_scale_factor(ag_z));

    let steps = 100;
    for _ in 0..steps {
        analysis
            .step()
            .expect("orthogonal ground-motion-driven spatial SDOF pair should solve every step");
    }

    let t = dt * steps as f64;
    let expected_x = -(ag_x / (omega_x * omega_x)) * (1.0 - (omega_x * t).cos());
    let expected_z = -(ag_z / (omega_z * omega_z)) * (1.0 - (omega_z * t).cos());

    let node = analysis.domain().node(mass_node);
    let u_x = node.displacement[SpatialDof::Ux as usize];
    let u_z = node.displacement[SpatialDof::Uz as usize];
    assert!(
        (u_x - expected_x).abs() < 1e-3,
        "ux: expected {expected_x}, got {u_x}"
    );
    assert!(
        (u_z - expected_z).abs() < 1e-3,
        "uz: expected {expected_z}, got {u_z}"
    );
}
