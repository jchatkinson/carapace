use carapace_core::analysis::modal_analysis;
use carapace_core::model::{Domain, Element, Material, Node, ZeroLength};

/// M5 acceptance (implementation-plan §6): eigenvalues of a small known
/// model match closed form. This is the classic 2-DOF "1-1-1-1" mass-spring
/// chain (unit masses, unit spring stiffnesses in series): ground -k1- m1
/// -k2- m2, giving `K = [[2,-1],[-1,1]]`, `M = I`. Its natural frequencies
/// have an exact closed form in the golden ratio: `omega = 1/phi` and
/// `omega = phi`, where `phi = (1+sqrt(5))/2` — from the characteristic
/// polynomial `lambda^2 - 3*lambda + 1 = 0`.
#[test]
fn mass_spring_chain_matches_golden_ratio_closed_form() {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let m1 = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2).with_mass(0, 1.0));
    let m2 = domain.add_node(Node::new([2.0, 0.0]).fix(1).fix(2).with_mass(0, 1.0));

    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, m1).with_material(0, Material::Elastic { e: 1.0 }),
    ));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(m1, m2).with_material(0, Material::Elastic { e: 1.0 }),
    ));

    let modes = modal_analysis(&mut domain, 2).expect("2-DOF chain should have a well-posed eigenproblem");
    assert_eq!(modes.len(), 2);

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let expected = [1.0 / phi, phi];

    for (mode, &expected_omega) in modes.iter().zip(expected.iter()) {
        assert!(
            (mode.frequency - expected_omega).abs() < 1e-9,
            "expected omega={expected_omega}, got {}",
            mode.frequency
        );
    }
}

/// Asking for fewer modes than the model has free DOFs is the whole point
/// of the Lanczos switch (implementation-plan §5 decision #2, revisited
/// post-`SparseSolver` sparse switch) — a dense full-spectrum eigensolve
/// can't do partial-spectrum extraction at all. A 5-mass chain (5 free
/// DOFs) asking for just the lowest 2 modes should return exactly those 2,
/// still matching this chain's true lowest frequencies (checked against a
/// `num_modes = 5` full solve of the same system).
#[test]
fn requesting_fewer_modes_than_free_dofs_returns_the_lowest_ones() {
    let build = || {
        let mut domain = Domain::new();
        let mut prev = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
        for i in 1..=5 {
            let node = domain.add_node(
                Node::new([i as f64, 0.0]).fix(1).fix(2).with_mass(0, 1.0),
            );
            domain.add_element(Element::ZeroLength(
                ZeroLength::new(prev, node).with_material(0, Material::Elastic { e: 1.0 }),
            ));
            prev = node;
        }
        domain
    };

    let mut full_domain = build();
    let full = modal_analysis(&mut full_domain, 5).expect("5-DOF chain should solve");

    let mut partial_domain = build();
    let partial = modal_analysis(&mut partial_domain, 2).expect("partial request should solve");

    assert_eq!(partial.len(), 2);
    for (p, f) in partial.iter().zip(full.iter()) {
        assert!(
            (p.frequency - f.frequency).abs() < 1e-7,
            "partial-spectrum frequency {} should match full-spectrum {}",
            p.frequency,
            f.frequency
        );
    }
}

/// A node with a free DOF but no assigned mass makes the generalized
/// eigenproblem singular (physically: a massless DOF has infinite natural
/// frequency) — `modal_analysis` should report that, not silently divide
/// by zero.
#[test]
fn massless_free_dof_is_reported_as_a_singular_system() {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let m1 = domain.add_node(Node::new([1.0, 0.0]).fix(1).fix(2)); // no mass assigned

    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, m1).with_material(0, Material::Elastic { e: 1.0 }),
    ));

    let result = modal_analysis(&mut domain, 1);
    assert!(result.is_err(), "a massless free DOF should be reported as an error");
}
