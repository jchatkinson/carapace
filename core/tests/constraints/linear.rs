//! General linear multi-point constraints: `u_slave = sum(coeff * u_master)`
//! with masters on any nodes. Covers the generated kinds (`equal_dof`, rigid
//! diaphragms, rigid links), arbitrary coefficients, chains of any depth,
//! normalization, every structural error, the refusal to prescribe a
//! constraint master, and reading slave motion back out of free-DOF vectors.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, AnalysisError, ConstraintHandler,
    ConvergenceTest, Integrator, RayleighDamping, TransientAnalysis,
};
use carapace_core::model::{
    Domain, Domain3, Element, Element3, Material, ModelError, Node, Node3, NodeId, SpatialDof,
    ZeroLength, ZeroLength3,
};
use nalgebra::{DMatrix, DVector};

fn spring(e: f64) -> Material {
    Material::Elastic { e }
}

fn analysis<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: Domain<NDIM, NDOF, NId, E>,
) -> Analysis<NDIM, NDOF, NId, E>
where
    NId: slotmap::Key,
    E: carapace_core::model::ElementOps<NDIM, NDOF, NId>,
{
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain)
}

/// A ground spring on one DOF of a node: `ground` is fully fixed.
fn grounded(domain: &mut Domain, node: NodeId, dof: usize, k: f64) {
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, node).with_material(dof, spring(k)),
    ));
}

// ---------------------------------------------------------------------------
// Rigid links against hand-solved rigid-body kinematics.
// ---------------------------------------------------------------------------

/// 2D: springs `kx`, `ky` on a slave offset `(dx, dy)` from a loaded master.
/// Master stiffness is `T^T K T` with `u_s = [ux - theta*dy, uy + theta*dx]`.
#[test]
fn planar_rigid_link_matches_the_hand_assembled_master_stiffness() {
    let (kx, ky, kr, dx, dy) = (400.0, 250.0, 1000.0, 3.0, 4.0);
    let load = [10.0, -6.0, 20.0];

    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]));
    let slave = domain.add_node(Node::new([dx, dy]));
    grounded(&mut domain, slave, 0, kx);
    grounded(&mut domain, slave, 1, ky);
    // The slave's rotation is the master's, so a rotational spring there
    // restrains the master's rotation (without it the system has a mechanism).
    grounded(&mut domain, slave, 2, kr);
    domain.rigid_link(master, slave);
    for (dof, value) in load.into_iter().enumerate() {
        domain.load_node(master, dof, value);
    }

    let mut analysis = analysis(domain);
    analysis.step().expect("rigid link solves");

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(3, 3, &[
        kx,          0.0,        -kx * dy,
        0.0,         ky,          ky * dx,
        -kx * dy,    ky * dx,     kx * dy * dy + ky * dx * dx + kr,
    ]);
    let expected = k
        .lu()
        .solve(&DVector::from_row_slice(&load))
        .expect("nonsingular");
    let m = analysis.domain().node(master).displacement;
    for i in 0..3 {
        assert!(
            (m[i] - expected[i]).abs() < 1e-12,
            "master dof {i}: {} vs {}",
            m[i],
            expected[i]
        );
    }
    // The slave follows rigidly: translation from the lever arm, same rotation.
    let s = analysis.domain().node(slave).displacement;
    assert!((s[0] - (m[0] - m[2] * dy)).abs() < 1e-12);
    assert!((s[1] - (m[1] + m[2] * dx)).abs() < 1e-12);
    assert!((s[2] - m[2]).abs() < 1e-15);
}

/// 3D: translational springs on an offset slave, `u_s = u_m + theta x d`, so
/// the master's 6x6 stiffness is `T^T K T` with `T = [I | -[d]x]`.
#[test]
fn spatial_rigid_link_matches_the_hand_assembled_master_stiffness() {
    let k = [300.0, 200.0, 100.0];
    let kr = [500.0, 600.0, 700.0];
    let d = [2.0, -1.0, 3.0];
    let load = [5.0, -3.0, 8.0, 1.0, 2.0, -4.0];

    let mut domain = Domain3::new();
    let master = domain.add_node(Node3::new([0.0, 0.0, 0.0]));
    let slave = domain.add_node(Node3::new(d));
    let ground = domain.add_node(Node3::new(d).fix(0).fix(1).fix(2).fix(3).fix(4).fix(5));
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, slave)
            .with_material(0, spring(k[0]))
            .with_material(1, spring(k[1]))
            .with_material(2, spring(k[2]))
            // Rotational springs: the slave's rotations are the master's, and
            // without them the master's rotations are a mechanism.
            .with_material(3, spring(kr[0]))
            .with_material(4, spring(kr[1]))
            .with_material(5, spring(kr[2])),
    ));
    domain.rigid_link(master, slave);
    for (dof, value) in load.into_iter().enumerate() {
        domain.load_node(master, dof, value);
    }
    let mut analysis = analysis(domain);
    analysis.step().expect("3D rigid link solves");

    // T: 3 x 6, u_s = [I | -[d]x] [u_m; theta_m] with [d]x the skew matrix.
    #[rustfmt::skip]
    let t = DMatrix::from_row_slice(3, 6, &[
        1.0, 0.0, 0.0,   0.0,  d[2], -d[1],
        0.0, 1.0, 0.0, -d[2],  0.0,   d[0],
        0.0, 0.0, 1.0,  d[1], -d[0],  0.0,
    ]);
    let kdiag = DMatrix::from_diagonal(&DVector::from_row_slice(&k));
    let mut kmm = t.transpose() * kdiag * &t;
    for i in 0..3 {
        kmm[(3 + i, 3 + i)] += kr[i];
    }
    let expected = kmm
        .lu()
        .solve(&DVector::from_row_slice(&load))
        .expect("nonsingular");
    let m = analysis.domain().node(master).displacement;
    for i in 0..6 {
        assert!(
            (m[i] - expected[i]).abs() < 1e-11,
            "master dof {i}: {} vs {}",
            m[i],
            expected[i]
        );
    }
    let s = analysis.domain().node(slave).displacement;
    let slave_translation = &t * DVector::from_row_slice(&m);
    for i in 0..3 {
        assert!(
            (s[i] - slave_translation[i]).abs() < 1e-11,
            "slave translation {i}"
        );
    }
    for i in 3..6 {
        assert!(
            (s[i] - m[i]).abs() < 1e-15,
            "slave rotation {i} follows the master"
        );
    }
}

/// A rigid translation of the whole model satisfies every constraint, so the
/// ground-motion influence vector must read back as exactly 1 on every
/// translation of every node — including slaves — and 0 on rotations.
#[test]
fn rigid_translation_field_satisfies_links_diaphragms_and_chains() {
    let mut domain = Domain3::new();
    let master = domain.add_node(Node3::new([0.0, 0.0, 0.0]));
    let link = domain.add_node(Node3::new([1.0, 2.0, 3.0]));
    let chained = domain.add_node(Node3::new([4.0, -1.0, 2.0]));
    let floor = domain.add_node(Node3::new([5.0, 0.0, -2.0]));
    domain.rigid_link(master, link);
    domain.rigid_link(link, chained); // a chain of two links
    domain.rigid_diaphragm(master, &[floor]);
    domain.validate().expect("valid model");

    for direction in 0..3 {
        let iota = domain.direction_incidence(direction);
        for node in [master, link, chained] {
            for dof in 0..6 {
                let expected = if dof == direction { 1.0 } else { 0.0 };
                assert!(
                    (domain.value_at(&iota, node, dof) - expected).abs() < 1e-15,
                    "direction {direction}, dof {dof}"
                );
            }
        }
        // The diaphragm only ties its two in-plane translations.
        for dof in [0, 2] {
            let expected = if dof == direction { 1.0 } else { 0.0 };
            assert!((domain.value_at(&iota, floor, dof) - expected).abs() < 1e-15);
        }
    }
}

// ---------------------------------------------------------------------------
// Arbitrary coefficients, normalization and chains.
// ---------------------------------------------------------------------------

/// `u_slave = 2 * u_master`: a ground spring `k` on the slave and a load `P`
/// on the master give a master stiffness of `4k`, so `u_m = P / (4k)`.
#[test]
fn a_coefficient_of_two_is_honored_not_treated_as_an_alias() {
    let (k, p) = (50.0, 8.0);
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let slave = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    grounded(&mut domain, slave, 0, k);
    domain.add_constraint((slave, 0), &[(master, 0, 2.0)]);
    domain.load_node(master, 0, p);

    let mut analysis = analysis(domain);
    analysis.step().expect("scaled tie solves");
    let (um, us) = (
        analysis.domain().node(master).displacement[0],
        analysis.domain().node(slave).displacement[0],
    );
    assert!((um - p / (4.0 * k)).abs() < 1e-14);
    assert!((us - 2.0 * um).abs() < 1e-14);
    // The slave is not an unknown of its own and has no single equation.
    assert!(analysis.domain().equation_of(slave, 0).is_none());
    assert!(analysis.domain().equation_of(master, 0).is_some());
}

/// Repeated masters merge and zero coefficients vanish, so a term with a zero
/// coefficient does not even activate its DOF.
#[test]
fn repeated_masters_merge_and_zero_terms_do_not_activate_their_dof() {
    let (k, p) = (50.0, 8.0);
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]).fix(1));
    let slave = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    grounded(&mut domain, slave, 0, k);
    // 1 + 1 on the same master dof is a coefficient of 2; the zero term on rz drops out.
    domain.add_constraint(
        (slave, 0),
        &[(master, 0, 1.0), (master, 0, 1.0), (master, 2, 0.0)],
    );
    domain.load_node(master, 0, p);

    let mut analysis = analysis(domain);
    analysis.step().expect("merged tie solves");
    assert!((analysis.domain().node(master).displacement[0] - p / (4.0 * k)).abs() < 1e-14);
    assert!(
        !analysis.domain().is_active(master, 2),
        "a zero-coefficient master is not activated"
    );
}

/// The same tie written as `equal_dof` and as a general constraint solves
/// identically, and chains resolve whatever order the constraints were added in.
#[test]
fn equal_dof_chains_and_general_constraints_agree() {
    let solve = |chain_in_reverse: bool, general: bool| {
        let mut domain: Domain = Domain::new();
        let a = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
        let b = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
        let c = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
        grounded(&mut domain, c, 0, 70.0);
        domain.load_node(a, 0, 14.0);
        let tie = |domain: &mut Domain, retained: NodeId, constrained: NodeId| {
            if general {
                domain.add_constraint((constrained, 0), &[(retained, 0, 1.0)]);
            } else {
                domain.equal_dof(retained, constrained, &[0]);
            }
        };
        if chain_in_reverse {
            tie(&mut domain, b, c);
            tie(&mut domain, a, b);
        } else {
            tie(&mut domain, a, b);
            tie(&mut domain, b, c);
        }
        let mut analysis = analysis(domain);
        analysis.step().expect("chain solves");
        [a, b, c].map(|n| analysis.domain().node(n).displacement[0])
    };
    let reference = solve(false, false);
    assert!((reference[0] - 14.0 / 70.0).abs() < 1e-14);
    assert_eq!(reference[0], reference[1]);
    assert_eq!(reference[1], reference[2]);
    for (reverse, general) in [(true, false), (false, true), (true, true)] {
        assert_eq!(
            solve(reverse, general),
            reference,
            "reverse={reverse} general={general}"
        );
    }
}

/// A rigid diaphragm whose retained node is itself a slave of another rigid
/// link (a chained diaphragm) used to panic.
#[test]
fn a_diaphragm_can_hang_off_a_slave() {
    let mut domain = Domain3::new();
    let top = domain.add_node(Node3::new([0.0, 0.0, 0.0]));
    let retained = domain.add_node(Node3::new([1.0, 0.0, 1.0]));
    let tied = domain.add_node(Node3::new([3.0, 0.0, -1.0]));
    domain.rigid_link(top, retained);
    domain.rigid_diaphragm(retained, &[tied]);
    domain
        .validate()
        .expect("a diaphragm retained node may be a slave");
    let mut q = DVector::zeros(domain.num_free_dofs());
    // theta_y at the top rotates retained and, through it, tied.
    let ry = domain
        .equation_of(top, SpatialDof::Ry as usize)
        .expect("ry is a free unknown");
    q[ry] = 0.01;
    let ux_tied = domain.value_at(&q, tied, SpatialDof::Ux as usize);
    // The link gives `retained` ux = theta_y * dz = 0.01 * 1 and ry = 0.01. The
    // diaphragm (normal Y) then adds ry * (z_tied - z_retained) = 0.01 * (-1 - 1)
    // to the tied node's ux.
    let expected = 0.01 * 1.0 + 0.01 * (-1.0 - 1.0);
    assert!(
        (ux_tied - expected).abs() < 1e-15,
        "{ux_tied} vs {expected}"
    );
}

// ---------------------------------------------------------------------------
// Structural errors.
// ---------------------------------------------------------------------------

fn pair() -> (Domain, NodeId, NodeId) {
    let mut domain: Domain = Domain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let b = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    grounded(&mut domain, a, 0, 10.0);
    (domain, a, b)
}

/// Constraining the same slave DOF twice is rejected with `DuplicateSlave`.
#[test]
fn a_slave_defined_twice_is_rejected() {
    let (mut domain, a, b) = pair();
    domain.add_constraint((b, 0), &[(a, 0, 1.0)]);
    domain.add_constraint((b, 0), &[(a, 0, 2.0)]);
    assert_eq!(
        domain.validate(),
        Err(ModelError::DuplicateSlave { node: 1, dof: 0 })
    );
}

/// A constraint whose slave DOF is also fixed is rejected by `validate` with `SlaveIsFixed`.
#[test]
fn a_fixed_slave_is_rejected() {
    let (mut domain, a, _) = pair();
    let fixed = domain.add_node(Node::new([2.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_constraint((fixed, 0), &[(a, 0, 1.0)]);
    // Nodes in insertion order: a, b, the ground `pair` adds, then `fixed`.
    assert_eq!(
        domain.validate(),
        Err(ModelError::SlaveIsFixed { node: 3, dof: 0 })
    );
}

/// Constraint cycles (A depends on B and B on A) and a slave that references itself are rejected
/// with `ConstraintCycle`.
#[test]
fn cycles_and_self_references_are_rejected() {
    let (mut domain, a, b) = pair();
    domain.add_constraint((a, 0), &[(b, 0, 1.0)]);
    domain.add_constraint((b, 0), &[(a, 0, 1.0)]);
    assert!(matches!(
        domain.validate(),
        Err(ModelError::ConstraintCycle { dof: 0, .. })
    ));

    let (mut domain, a, _) = pair();
    domain.add_constraint((a, 0), &[(a, 0, 1.0)]);
    assert!(matches!(
        domain.validate(),
        Err(ModelError::ConstraintCycle { dof: 0, .. })
    ));
}

/// A master DOF fixed at a nonzero initial displacement is rejected with
/// `ConstraintOnPrescribedDof`, since constrained slaves cannot follow a prescribed value.
#[test]
fn a_master_fixed_at_a_nonzero_value_is_rejected() {
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(
        Node::new([0.0, 0.0])
            .fix(0)
            .fix(1)
            .fix(2)
            .with_initial_displacement(0, 0.5),
    );
    let slave = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    domain.equal_dof(master, slave, &[0]);
    assert_eq!(
        domain.validate(),
        Err(ModelError::ConstraintOnPrescribedDof { node: 0, dof: 0 })
    );
}

/// Nodal mass on a slave DOF that is a combination of several unknowns (a rigid-diaphragm slave
/// with a lever arm) cannot be lumped. `validate` and `TransientAnalysis::new` report
/// `MassOnConstrainedDof` as an error instead of panicking.
#[test]
fn mass_on_a_slave_that_combines_several_unknowns_is_an_error_not_a_panic() {
    let mut domain = Domain3::new();
    let master = domain.add_node(Node3::new([0.0, 0.0, 0.0]).fix(1).fix(3).fix(5));
    let slave = domain.add_node(Node3::new([4.0, 0.0, -6.0]).fix(1).with_mass(0, 2.0));
    domain.rigid_diaphragm(master, &[slave]);
    assert_eq!(
        domain.validate(),
        Err(ModelError::MassOnConstrainedDof { node: 1, dof: 0 })
    );
    let transient = TransientAnalysis::new(domain, RayleighDamping::NONE, 0.01);
    assert!(matches!(
        transient,
        Err(AnalysisError::InvalidModel(
            ModelError::MassOnConstrainedDof { .. }
        ))
    ));
}

/// A single scaled unknown (`u = c * q`) can carry mass: `c^2 m` lands on `q`.
#[test]
fn mass_on_a_scaled_single_term_slave_contributes_coefficient_squared() {
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2).with_mass(0, 1.0));
    let slave = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2).with_mass(0, 3.0));
    domain.add_constraint((slave, 0), &[(master, 0, 2.0)]);
    domain.validate().expect("single-term slave may carry mass");
    let mass = domain.assemble_mass_diagonal();
    assert_eq!(mass.len(), 1);
    assert!((mass[0] - (1.0 + 4.0 * 3.0)).abs() < 1e-15);
}

// ---------------------------------------------------------------------------
// Initial state.
// ---------------------------------------------------------------------------

/// A slave whose initial displacement disagrees with its constraint equation is rejected with
/// `InconsistentInitialState`.
#[test]
fn an_inconsistent_nonzero_initial_slave_state_is_rejected() {
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let slave = domain.add_node(
        Node::new([1.0, 0.0])
            .fix(1)
            .fix(2)
            .with_initial_displacement(0, 1.0),
    );
    domain.equal_dof(master, slave, &[0]);
    assert_eq!(
        domain.validate(),
        Err(ModelError::InconsistentInitialState { node: 1, dof: 0 })
    );
}

/// A slave left at the default zero takes the value its masters imply (in
/// dependency order along a chain); a consistent nonzero value is accepted.
#[test]
fn zero_slave_state_is_derived_from_the_masters_and_consistent_state_passes() {
    let mut domain: Domain = Domain::new();
    let a = domain.add_node(
        Node::new([0.0, 0.0])
            .fix(1)
            .fix(2)
            .with_initial_displacement(0, 0.25)
            .with_initial_velocity(0, -2.0),
    );
    let b = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    let c = domain.add_node(
        Node::new([2.0, 0.0])
            .fix(1)
            .fix(2)
            .with_initial_displacement(0, 0.5), // 2 * 0.25: consistent
    );
    // c's constraint is registered before b's, so the chain must resolve in order.
    domain.add_constraint((c, 0), &[(b, 0, 2.0)]);
    domain.equal_dof(a, b, &[0]);
    // `a` is the only unknown, so give it something to push against.
    grounded(&mut domain, a, 0, 5.0);
    domain.validate().expect("consistent state");
    assert_eq!(domain.node(b).displacement[0], 0.25);
    assert_eq!(domain.node(b).velocity[0], -2.0);
    assert_eq!(
        domain.node(c).velocity[0],
        -4.0,
        "velocity derived through the chain"
    );
}

// ---------------------------------------------------------------------------
// Prescribed displacements.
// ---------------------------------------------------------------------------

/// A fixed DOF that is a constraint master cannot be prescribed (its slaves
/// would not follow), including the first link of a chain, and a refused step
/// leaves the domain untouched.
#[test]
fn prescribing_a_constraint_master_is_refused_and_leaves_state_unchanged() {
    let mut domain: Domain = Domain::new();
    let support = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let s1 = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2));
    let s2 = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2));
    let free = domain.add_node(Node::new([3.0, 0.0]).fix(1).fix(2));
    let other = domain.add_node(Node::new([4.0, 0.0]).fix(0).fix(1).fix(2));
    // support.ux (fixed) -> s1.ux -> s2.ux: a chain on one DOF.
    domain.equal_dof(support, s1, &[0]);
    domain.equal_dof(s1, s2, &[0]);
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(s2, free).with_material(0, spring(10.0)),
    ));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(other, free).with_material(0, spring(10.0)),
    ));

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);

    assert!(
        !analysis.domain().can_prescribe(support, 0, 1.0),
        "first link of a chain"
    );
    assert!(
        analysis.domain().can_prescribe(other, 0, 1.0),
        "an ordinary support"
    );

    let nodes = [support, s1, s2, free, other];
    let snapshot = |a: &Analysis| nodes.map(|n| a.domain().node(n).displacement);
    let before = snapshot(&analysis);
    let refused = analysis.step_prescribed(&[(other, 0, 0.5), (support, 0, 1.0)]);
    assert_eq!(refused.unwrap_err(), AnalysisError::InvalidConstraint);
    assert_eq!(
        before,
        snapshot(&analysis),
        "no target is applied when any target is refused"
    );

    // The legal target alone is accepted: `free` sits halfway between 0 and 0.5.
    analysis
        .step_prescribed(&[(other, 0, 0.5)])
        .expect("ordinary settlement");
    assert_eq!(analysis.domain().node(other).displacement[0], 0.5);
    assert!((analysis.domain().node(free).displacement[0] - 0.25).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// Reading slave motion out of free-DOF vectors.
// ---------------------------------------------------------------------------

/// A constrained node's motion in a mode shape comes from `value_at`, not
/// from the (absent) equation of the slave: a mass on a rigid-link master with
/// a spring on the slave moves the slave by the lever-arm kinematics.
#[test]
fn mode_shapes_report_slave_motion_through_the_constraint() {
    let (dx, dy) = (3.0, 4.0);
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(
        Node::new([0.0, 0.0])
            .with_mass(0, 2.0)
            .with_mass(1, 2.0)
            .with_mass(2, 5.0),
    );
    let slave = domain.add_node(Node::new([dx, dy]));
    grounded(&mut domain, slave, 0, 400.0);
    grounded(&mut domain, slave, 1, 250.0);
    grounded(&mut domain, slave, 2, 1000.0);
    domain.rigid_link(master, slave);

    let modes = modal_analysis(&mut domain, 3).expect("three modes");
    for mode in &modes {
        let (ux, uy, rz) = (
            domain.value_at(&mode.shape, master, 0),
            domain.value_at(&mode.shape, master, 1),
            domain.value_at(&mode.shape, master, 2),
        );
        let sx = domain.value_at(&mode.shape, slave, 0);
        let sy = domain.value_at(&mode.shape, slave, 1);
        let sr = domain.value_at(&mode.shape, slave, 2);
        assert!(
            (sx - (ux - rz * dy)).abs() < 1e-12,
            "slave ux in mode {}",
            mode.frequency
        );
        assert!(
            (sy - (uy + rz * dx)).abs() < 1e-12,
            "slave uy in mode {}",
            mode.frequency
        );
        assert!((sr - rz).abs() < 1e-12);
        // And the slave genuinely moves (it is not read as zero).
        assert!(sx.abs() + sy.abs() > 1e-9);
    }
}
