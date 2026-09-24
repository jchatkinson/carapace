use carapace_core::analysis::{Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, RayleighDamping, TransientAnalysis};
use carapace_core::model::{Axis3, Domain3, Element3, Material, Node3, SpatialDof, ZeroLength3};

/// The classical "asymmetric rigid diaphragm" acceptance case: two lateral
/// springs ("shear walls") at different offsets from the diaphragm's master
/// (retained) node, tied to it via `Domain3::rigid_diaphragm` (default
/// normal `Axis3::Y`, so the diaphragm plane is `x`-`z` — see
/// `Domain3::rigid_diaphragm`'s doc comment for the kinematics). A lateral
/// load applied at the master node (the diaphragm's center) still induces
/// rotation whenever the walls' combined stiffness isn't itself centered
/// there — exactly the mechanism the lever-arm term exists for; an identity
/// alias (what planar `rigid_diaphragm`/`equal_dof` do) could never produce
/// this coupling.
///
/// Each wall is a `ZeroLength3` spring (ground → wall node) stiff only in
/// local `x` (`SpatialDof::Ux`), attached to a diaphragm-constrained wall
/// node offset from the master node in `z`. Assembling this through
/// `Domain3`'s general affine-constraint machinery (`dof_terms`'s
/// coefficient-product distribution in `assemble_stiffness_triplets`)
/// produces exactly the standard rigid-diaphragm torsional stiffness
/// matrix used throughout seismic engineering practice:
///
/// ```text
/// [ sum(k_i)         sum(k_i*z_i)   ] [U]     [P]
/// [ sum(k_i*z_i)      sum(k_i*z_i^2)] [theta] = [M]
/// ```
///
/// (`U` = master's `Ux`, `theta` = master's `Ry`, `z_i` each wall's `z`
/// offset from the master) — derived independently here from first
/// principles (not copied from the implementation) and solved by hand
/// (`M = 0`: the load is applied at the master node itself, not offset).
#[test]
fn asymmetric_shear_walls_reproduce_the_standard_torsional_diaphragm_stiffness_matrix() {
    let (k_a, z_a) = (10.0, 5.0);
    let (k_b, z_b) = (6.0, -3.0);
    let load = 100.0;

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

    // Master (retained) node at the diaphragm's geometric center. Only its
    // in-plane x-translation and rotation-about-y are free — the two DOFs
    // this rigid-diaphragm problem actually has; everything else (vertical
    // translation, both other rotations) has no source of stiffness in
    // this model and must be fixed directly, same as any other DOF a
    // model doesn't otherwise restrain.
    let master = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Rz as usize),
    );
    domain.load_node(master, SpatialDof::Ux as usize, load);

    let wall = |domain: &mut Domain3, z: f64, k: f64| {
        let node = domain.add_node(
            Node3::new([0.0, 0.0, z])
                .fix(SpatialDof::Uy as usize)
                .fix(SpatialDof::Rx as usize)
                .fix(SpatialDof::Ry as usize)
                .fix(SpatialDof::Rz as usize),
        );
        domain.add_element(Element3::ZeroLength3(
            ZeroLength3::new(ground, node).with_material(SpatialDof::Ux as usize, Material::Elastic { e: k }),
        ));
        node
    };
    let wall_a = wall(&mut domain, z_a, k_a);
    let wall_b = wall(&mut domain, z_b, k_b);

    domain.rigid_diaphragm(master, &[wall_a, wall_b]);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol: 1e-9, max_iter: 10 })
        .build(domain);

    analysis.step().expect("linear elastic torsional rigid diaphragm should solve");

    let k_uu = k_a + k_b;
    let k_ut = k_a * z_a + k_b * z_b;
    let k_tt = k_a * z_a * z_a + k_b * z_b * z_b;
    let det = k_uu * k_tt - k_ut * k_ut;
    let expected_u = k_tt * load / det;
    let expected_theta = -k_ut * load / det;

    let u = analysis.domain().node(master).displacement[SpatialDof::Ux as usize];
    let theta = analysis.domain().node(master).displacement[SpatialDof::Ry as usize];
    assert!((u - expected_u).abs() < 1e-9, "expected U={expected_u}, got {u}");
    assert!(
        (theta - expected_theta).abs() < 1e-9,
        "expected theta={expected_theta}, got {theta}"
    );

    // Each wall's own x-displacement must match the rigid-body kinematic
    // relation the diaphragm enforces: `u_c[x] = U + theta*(z_c - z_master)`
    // (this problem's `a = z`, `b = x` per `Domain3::rigid_diaphragm`'s
    // general formula, with `z_master = 0`) — checked independently of the
    // stiffness-matrix closed form above, against the raw kinematic
    // definition.
    let ux_a = analysis.domain().node(wall_a).displacement[SpatialDof::Ux as usize];
    let ux_b = analysis.domain().node(wall_b).displacement[SpatialDof::Ux as usize];
    assert!(
        (ux_a - (expected_u + expected_theta * z_a)).abs() < 1e-9,
        "wall A kinematics: expected {}, got {ux_a}",
        expected_u + expected_theta * z_a
    );
    assert!(
        (ux_b - (expected_u + expected_theta * z_b)).abs() < 1e-9,
        "wall B kinematics: expected {}, got {ux_b}",
        expected_u + expected_theta * z_b
    );
}

/// Rigid-rotation invariance: a diaphragm-constrained node's *out-of-plane*
/// translation (`uy`, never tied by `rigid_diaphragm` — see its doc
/// comment: "never vertical") and every one of its rotational DOFs
/// (likewise never tied — "never rotational") must be completely
/// unaffected by the master node's rotation, however large. Only the two
/// in-plane translations respond to the lever arm. This is the direct
/// acceptance check for the "translational constraint only, never
/// rotational, never vertical" scope the constraint is deliberately
/// limited to.
///
/// Driven through `TransientAnalysis::new` (public API) rather than the
/// crate-private `scatter_state` directly: it computes a consistent
/// initial acceleration from equilibrium, which exercises `Domain::
/// scatter_state`'s affine-derived-dof pass (`apply_displacement_
/// increment`'s twin — see its doc comment) with zero stiffness coupling
/// in play (`ZeroLength3` with no material assigned on any of the DOFs
/// under test has zero stiffness there), isolating pure kinematics. The
/// constrained node's untied DOFs (`uy`, `rx`, `ry`, `rz`) start at
/// distinct nonzero values via `Node3::with_initial_displacement`; if the
/// diaphragm's affine transform ever touched them, they'd be overwritten
/// to something derived from the master node's state instead of surviving
/// unchanged.
#[test]
fn rigid_diaphragm_never_couples_vertical_translation_or_any_rotation() {
    let mut domain = Domain3::new();

    let (u_r, w_r, theta_r) = (2.0, -1.5, 0.3);
    let master = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .with_initial_displacement(SpatialDof::Ux as usize, u_r)
            .with_initial_displacement(SpatialDof::Uz as usize, w_r)
            .with_initial_displacement(SpatialDof::Ry as usize, theta_r)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Ux as usize, 1.0)
            .with_mass(SpatialDof::Uz as usize, 1.0)
            .with_mass(SpatialDof::Ry as usize, 1.0),
    );
    let (x_c, z_c) = (4.0, -6.0);
    let constrained = domain.add_node(
        Node3::new([x_c, 0.0, z_c])
            .with_initial_displacement(SpatialDof::Uy as usize, 7.0)
            .with_initial_displacement(SpatialDof::Rx as usize, 1.0)
            .with_initial_displacement(SpatialDof::Ry as usize, 2.0)
            .with_initial_displacement(SpatialDof::Rz as usize, 3.0)
            .with_mass(SpatialDof::Uy as usize, 1.0)
            .with_mass(SpatialDof::Rx as usize, 1.0)
            .with_mass(SpatialDof::Ry as usize, 1.0)
            .with_mass(SpatialDof::Rz as usize, 1.0),
    );
    domain.rigid_diaphragm(master, &[constrained]);

    let analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, 0.01)
        .expect("diaphragm-constrained system with mass on every free dof should be well-posed");

    let c = analysis.domain().node(constrained);
    assert_eq!(c.displacement[SpatialDof::Uy as usize], 7.0, "out-of-plane translation must stay untouched");
    assert_eq!(c.displacement[SpatialDof::Rx as usize], 1.0, "no rotational dof is ever tied");
    assert_eq!(c.displacement[SpatialDof::Ry as usize], 2.0, "not even rotation about the diaphragm normal itself");
    assert_eq!(c.displacement[SpatialDof::Rz as usize], 3.0, "no rotational dof is ever tied");

    // The two in-plane translations *do* respond to the lever arm:
    // `a = z`, `b = x` for `normal = Y` (see `Domain3::rigid_diaphragm`'s
    // doc comment), with `lever_a = z_c - z_r`, `lever_b = x_c - x_r`, and
    // `(x_r, z_r) = (0, 0)` here.
    let expected_z = w_r - theta_r * x_c;
    let expected_x = u_r + theta_r * z_c;
    assert!((c.displacement[SpatialDof::Uz as usize] - expected_z).abs() < 1e-12);
    assert!((c.displacement[SpatialDof::Ux as usize] - expected_x).abs() < 1e-12);
}

/// `rigid_diaphragm_about`'s explicit `normal` parameter, checked against
/// Xara/OpenSees's own documented `RigidDiaphragm` constraint equations for
/// their default `perpDirn = 3` (z-normal, diaphragm plane `x`-`y`): `Ux_c
/// = Ux_r - Rz_r*(Yc-Yr)`, `Uy_c = Uy_r + Rz_r*(Xc-Xr)` — literally, not
/// via this crate's general `a`/`b`/`normal` formula (`Domain3::
/// rigid_diaphragm`'s doc comment shows the two agree, but this test
/// doesn't rely on that derivation being right). Same kinematics-only
/// approach as the vertical/rotation test above.
#[test]
fn rigid_diaphragm_about_z_normal_matches_opensees_rigid_diaphragm_equations() {
    let mut domain = Domain3::new();

    let (u_r, v_r, theta_r) = (1.2, -0.8, 0.15); // Ux_r, Uy_r, Rz_r
    let master = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .with_initial_displacement(SpatialDof::Ux as usize, u_r)
            .with_initial_displacement(SpatialDof::Uy as usize, v_r)
            .with_initial_displacement(SpatialDof::Rz as usize, theta_r)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .with_mass(SpatialDof::Ux as usize, 1.0)
            .with_mass(SpatialDof::Uy as usize, 1.0)
            .with_mass(SpatialDof::Rz as usize, 1.0),
    );
    let (x_c, y_c) = (3.0, 7.0);
    let constrained = domain.add_node(
        Node3::new([x_c, y_c, 0.0])
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    domain.rigid_diaphragm_about(master, &[constrained], Axis3::Z);

    let analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, 0.01)
        .expect("z-normal diaphragm-constrained system should be well-posed");

    let c = analysis.domain().node(constrained);
    let expected_ux = u_r - theta_r * y_c;
    let expected_uy = v_r + theta_r * x_c;
    assert!(
        (c.displacement[SpatialDof::Ux as usize] - expected_ux).abs() < 1e-12,
        "expected Ux_c={expected_ux}, got {}",
        c.displacement[SpatialDof::Ux as usize]
    );
    assert!(
        (c.displacement[SpatialDof::Uy as usize] - expected_uy).abs() < 1e-12,
        "expected Uy_c={expected_uy}, got {}",
        c.displacement[SpatialDof::Uy as usize]
    );
}

/// Documented limitation: nodal mass directly on a rigid-diaphragm-tied
/// translation can't be represented by the lumped mass *diagonal*
/// `assemble_mass_diagonal` returns — a real rigid diaphragm's tributary
/// mass there would need an off-diagonal (rotational-inertia-at-the-
/// master) term, which this crate's dynamics don't support yet (see
/// `spatial-architecture.md`). A clear panic, not a silently-wrong lumped
/// mass, is the intended failure mode — same "loud, not silent" reasoning
/// `GeomTransf3::local_axes`'s degenerate-input asserts already use.
#[test]
#[should_panic(expected = "rigid-diaphragm-affine-constrained dof")]
fn nodal_mass_on_a_diaphragm_tied_translation_panics_rather_than_silently_dropping_inertia() {
    let mut domain = Domain3::new();
    let master = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Ux as usize, 1.0)
            .with_mass(SpatialDof::Uz as usize, 1.0)
            .with_mass(SpatialDof::Ry as usize, 1.0),
    );
    let constrained = domain.add_node(
        Node3::new([4.0, 0.0, -6.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize)
            .with_mass(SpatialDof::Ux as usize, 2.0), // on a diaphragm-tied dof — not representable
    );
    domain.rigid_diaphragm(master, &[constrained]);

    let _ = domain.assemble_mass_diagonal();
}
