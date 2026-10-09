//! DOF activation: a node DOF is an equation only if an element stiffens it,
//! a constraint uses it as a master, it carries a nodal mass, or the user
//! fixed it. Models therefore no longer need to fix the DOFs an element type
//! never uses (a truss's rotations), and a load on a DOF nothing resists is
//! reported instead of silently dropped.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, ArcLength, ArcScales, ConstraintHandler,
    ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    Domain, Domain3, ElasticBeamColumn, Element, Element3, GeomTransf, Material, ModelError, Node,
    Node3, NodeId, SpatialDof, Truss, Truss3, ZeroLength,
};

fn linear(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(if domain.has_mp_constraints() {
            ConstraintHandler::Transformation
        } else {
            ConstraintHandler::Plain
        })
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain)
}

/// A two-bar truss, pinned at one support and on a roller at the other, with
/// a vertical load at the apex. Rotations are fixed only when asked.
fn triangle_truss(fix_rotations: bool) -> (Domain, [NodeId; 3]) {
    let mut domain = Domain::new();
    let pin = |n: Node| if fix_rotations { n.fix(2) } else { n };
    let a = domain.add_node(pin(Node::new([0.0, 0.0]).fix(0).fix(1)));
    let b = domain.add_node(pin(Node::new([4000.0, 0.0]).fix(1)));
    let c = domain.add_node(pin(Node::new([2000.0, 1500.0])));
    let material = || Material::Elastic { e: 200_000.0 };
    for (i, j) in [(a, b), (a, c), (b, c)] {
        domain.add_element(Element::Truss(Truss::new(i, j, 500.0, material())));
    }
    domain.load_node(c, 1, -1.0e4);
    (domain, [a, b, c])
}

/// A truss-only 2D model solves without fixing rotations: only the 3 translational DOFs become
/// equations, and the displacements are identical to the same model with rotations fixed by hand.
#[test]
fn truss_only_model_needs_no_fixed_rotations_and_matches_the_fixed_model() {
    let (domain, nodes) = triangle_truss(false);
    let mut free = linear(domain);
    free.step()
        .expect("truss-only model solves with unfixed rotations");

    let (domain, _) = triangle_truss(true);
    let mut fixed = linear(domain);
    fixed
        .step()
        .expect("same model with rotations fixed by hand");

    // Only the translations are unknowns: 2 + 1 + 2 DOFs, not 3 per node.
    assert_eq!(free.domain().num_free_dofs(), 3);
    for node in nodes {
        for dof in 0..2 {
            let (u1, u2) = (
                free.domain().node(node).displacement[dof],
                fixed.domain().node(node).displacement[dof],
            );
            assert_eq!(u1, u2, "node displacement must be identical");
        }
        assert_eq!(free.domain().node(node).displacement[2], 0.0);
    }
}

/// A 3D skew truss with no rotations fixed solves, with only the translational DOFs as unknowns.
#[test]
fn spatial_truss_needs_no_fixed_rotations() {
    // The skew cantilever truss from the domain tests, with no rotation fixed.
    let mut domain = Domain3::new();
    let fixed = domain.add_node(
        Node3::new([0.0, 0.0, 0.0])
            .fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize),
    );
    let free = domain.add_node(Node3::new([3.0, 4.0, 0.0]).fix(SpatialDof::Uz as usize));
    let (area, e, force) = (2.0, 1000.0, 100.0);
    domain.add_element(Element3::Truss3(Truss3::new(
        fixed,
        free,
        area,
        Material::Elastic { e },
    )));
    domain.load_node(free, SpatialDof::Ux as usize, force * 0.6);
    domain.load_node(free, SpatialDof::Uy as usize, force * 0.8);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis
        .step()
        .expect("3D truss solves with unfixed rotations");
    let node = analysis.domain().node(free);
    let axial = node.displacement[0] * 0.6 + node.displacement[1] * 0.8;
    assert!((axial - force * 5.0 / (area * e)).abs() < 1e-12);
    assert_eq!(analysis.domain().num_free_dofs(), 2);
}

/// A beam shares a node with a brace truss: the joint's rotation is an
/// unknown (the beam stiffens it) while a truss-only node has none.
#[test]
fn joint_of_beam_and_truss_keeps_rotation_only_where_a_beam_needs_it() {
    let mut domain = Domain::new();
    let wall = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let tip = domain.add_node(Node::new([3000.0, 0.0]));
    let anchor = domain.add_node(Node::new([3000.0, 2000.0]).fix(0).fix(1));
    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
        wall,
        tip,
        200_000.0,
        8000.0,
        1.0e8,
        GeomTransf::Linear,
    )));
    domain.add_element(Element::Truss(Truss::new(
        tip,
        anchor,
        400.0,
        Material::Elastic { e: 200_000.0 },
    )));
    domain.load_node(tip, 0, 5.0e4);

    let mut analysis = linear(domain);
    analysis.step().expect("beam plus brace solves");
    let d = analysis.domain();
    assert!(d.is_active(tip, 2), "the beam stiffens the tip rotation");
    assert!(d.equation_of(tip, 2).is_some());
    assert!(!d.is_active(anchor, 2), "a truss-only node has no rotation");
    assert_eq!(d.equation_of(anchor, 2), None);
    // Axial load on a beam in series with an orthogonal brace: the beam alone
    // carries it axially, u = P L / (E A).
    let expected = 5.0e4 * 3000.0 / (200_000.0 * 8000.0);
    assert!((d.node(tip).displacement[0] - expected).abs() < 1e-9);
}

/// A nodal load on a DOF that no element stiffens (a moment on a truss-only node) is a model error,
/// `LoadOnInactiveDof`, reported by `validate` and returned as a value by `try_build`.
#[test]
fn a_nodal_load_on_a_dof_nothing_stiffens_is_a_model_error() {
    let (mut domain, [_, _, apex]) = triangle_truss(false);
    // Rotations are inactive in a truss-only model; a moment has nothing to push on.
    domain.load_node(apex, 2, 1.0);
    let error = domain.validate().expect_err("moment on a truss-only node");
    assert_eq!(error, ModelError::LoadOnInactiveDof { node: 2, dof: 2 });

    // `try_build` reports it as a value (what the wasm session needs).
    let built = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .try_build(domain);
    assert!(matches!(
        built,
        Err(ModelError::LoadOnInactiveDof { node: 2, dof: 2 })
    ));
}

/// `build` panics on a model that fails validation (a moment applied to a truss-only node), since
/// that is a model-construction error rather than a runtime failure.
#[test]
#[should_panic(expected = "invalid model")]
fn build_panics_on_an_invalid_model() {
    let (mut domain, [_, _, apex]) = triangle_truss(false);
    domain.load_node(apex, 2, 1.0);
    linear(domain);
}

/// Explicit fixity wins: fixing a DOF that no element uses (a truss node's rotation) is legal,
/// keeps the DOF active, and gives it no equation number.
#[test]
fn fixed_dofs_are_active_even_when_nothing_stiffens_them() {
    // Explicit fixity wins: boundary conditions (and their reactions) on a
    // DOF no element uses are legal and unchanged.
    let (mut domain, [a, ..]) = triangle_truss(true);
    assert!(domain.validate().is_ok());
    assert!(domain.is_active(a, 2));
    assert_eq!(domain.equation_of(a, 2), None);
    assert_eq!(domain.num_free_dofs(), 3);
}

/// A diaphragm's retained node typically has no element at all; the DOFs the
/// constraint ties to must still be equations.
#[test]
fn constraint_masters_are_active_without_any_element() {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let wall = domain.add_node(Node::new([0.0, 3000.0]));
    let master = domain.add_node(Node::new([0.0, 3000.0]));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, wall).with_material(0, Material::Elastic { e: 50.0 }),
    ));
    domain.equal_dof(master, wall, &[0]);
    domain.load_node(master, 0, 25.0);

    let mut analysis = linear(domain);
    analysis
        .step()
        .expect("diaphragm-style tie to an element-free master");
    let d = analysis.domain();
    assert!(d.is_active(master, 0));
    assert!(!d.is_active(master, 1));
    assert!(!d.is_active(master, 2));
    assert!((d.node(master).displacement[0] - 0.5).abs() < 1e-12);
    assert_eq!(d.node(master).displacement[0], d.node(wall).displacement[0]);
}

/// A nodal mass is inertia even with no stiffness, so it activates its DOF
/// (a free mass is legal in a transient analysis); a mass on a stiffened DOF
/// is unchanged.
#[test]
fn nodal_mass_activates_its_dof() {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mass = domain.add_node(Node::new([1.0, 0.0]).with_mass(0, 2.0).with_mass(1, 3.0));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, mass).with_material(0, Material::Elastic { e: 8.0 }),
    ));
    assert!(domain.validate().is_ok());
    assert!(domain.is_active(mass, 1), "mass on uy gives it inertia");
    assert!(!domain.is_active(mass, 2), "no mass, no stiffness on rz");
    assert_eq!(domain.num_free_dofs(), 2);
}

/// A modal analysis of a 2D truss SDOF with its rotations left free works and gives the closed-form
/// frequency sqrt(k/m), with k = EA/L and m = rho*A*L/2.
#[test]
fn modal_analysis_of_a_truss_needs_no_fixed_rotations() {
    let (e, area, length, density) = (30_000.0, 2.0, 100.0, 0.5);
    let mut domain = Domain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let b = domain.add_node(Node::new([length, 0.0]).fix(1));
    domain.add_element(Element::Truss(
        Truss::new(a, b, area, Material::Elastic { e }).with_density(density),
    ));
    let modes = modal_analysis(&mut domain, 1).expect("SDOF truss, rotations left free");
    let (k, m) = (e * area / length, density * area * length / 2.0);
    assert!((modes[0].frequency - (k / m).sqrt()).abs() < 1e-9);
}

/// Arc length needs a rotation scale only when the model has rotational
/// unknowns; with the rotations inactive it must not demand one.
#[test]
fn arc_length_needs_no_rotation_scale_when_rotations_are_inactive() {
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, free).with_material(0, Material::Elastic { e: 100.0 }),
    ));
    domain.load_node(free, 0, 1.0);

    let arc = ArcLength::fixed(
        0.1,
        ArcScales::Explicit {
            displacement: 0.01,
            rotation: None,
            load: 1.0,
        },
    );
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(arc))
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: 1e-9,
            moment_tol: 1e-9,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 20,
        })
        .build(domain);
    let step = analysis.step().expect("no rotation scale is needed");
    assert!(step.load_factor > 0.0);
}
