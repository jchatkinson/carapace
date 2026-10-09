//! Continuum elements together with frame elements, constraints, and the dynamic
//! and nonlinear analysis paths.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, ArcLength, ArcScales, ConstraintHandler,
    ConvergenceTest, Integrator, LineSearch, RayleighDamping, TangentStrategy, TransientAnalysis,
};
use carapace_core::model::{
    Domain, ElasticBeamColumn, Element, ElementId, GeomTransf, Node, NodeId, PlaneMaterial, Quad4,
    Quad4Formulation,
};

/// Adds the nodes of an `nx` by `ny` grid over `[0, l] x [0, h]` (row-major) with `setup`
/// customizing each node, then the Quad4 elements.
fn quad_grid(
    domain: &mut Domain,
    (nx, ny): (usize, usize),
    (l, h): (f64, f64),
    thickness: f64,
    density: f64,
    material: &PlaneMaterial,
    formulation: Quad4Formulation,
    setup: impl Fn(Node, usize, usize) -> Node,
) -> (Vec<NodeId>, Vec<ElementId>) {
    let mut ids = Vec::new();
    let mut elements = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let node = Node::new([l * i as f64 / nx as f64, h * j as f64 / ny as f64]);
            ids.push(domain.add_node(setup(node, i, j)));
        }
    }
    let at = |i: usize, j: usize| ids[j * (nx + 1) + i];
    for j in 0..ny {
        for i in 0..nx {
            let nodes = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            elements.push(
                domain.add_element(Element::Quad4(
                    Quad4::new(nodes, thickness, material.clone())
                        .with_density(density)
                        .with_formulation(formulation),
                )),
            );
        }
    }
    (ids, elements)
}

fn linear_analysis(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain)
}

/// A Q4 panel with a frame element along its bottom edge sharing the panel's nodes:
/// under a uniform axial stretch (nu = 0) the stiffness is `E A_panel / L + E A_beam / L`
/// exactly. Proves activation and assembly across element kinds (the beam's rotation
/// is an equation, the panel's is not).
#[test]
fn panel_and_beam_sharing_nodes_stiffen_in_parallel() {
    let (l, h, t, e, area_beam) = (8.0, 2.0, 0.5, 1000.0, 0.3);
    let stretch = 1e-3;
    let axial_reaction = |with_beam: bool| {
        let mut domain = Domain::new();
        let (ids, _) = quad_grid(
            &mut domain,
            (4, 2),
            (l, h),
            t,
            0.0,
            &PlaneMaterial::plane_stress(e, 0.0).unwrap(),
            Quad4Formulation::Full,
            |node, i, j| {
                let mut node = node;
                if i == 0 {
                    node = node.fix(0);
                    if j == 0 {
                        node = node.fix(1);
                    }
                }
                if i == 4 {
                    node = node.fix(0);
                }
                node
            },
        );
        if with_beam {
            for i in 0..4 {
                domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                    ids[i],
                    ids[i + 1],
                    e,
                    area_beam,
                    0.05,
                    GeomTransf::Linear,
                )));
            }
        }
        let mut analysis = linear_analysis(domain);
        let targets: Vec<_> = (0..=2).map(|j| (ids[j * 5 + 4], 0, stretch * l)).collect();
        analysis.step_prescribed(&targets).unwrap();
        (0..=2)
            .map(|j| analysis.domain().reaction(ids[j * 5], 0, 0.0))
            .sum::<f64>()
    };
    let panel = e * t * h / l * stretch * l;
    let beam = e * area_beam / l * stretch * l;
    assert!((axial_reaction(false).abs() - panel).abs() < 1e-9 * panel);
    assert!((axial_reaction(true).abs() - (panel + beam)).abs() < 1e-9 * (panel + beam));
}

/// A rigid link between the end of a beam (master) and a Q4 corner (slave): every
/// constrained mode shape reports the slave's motion through the constraint, equal to the
/// master's plus the lever-arm term.
#[test]
fn modal_shapes_of_a_beam_rigidly_linked_to_a_quad_edge_include_slave_motion() {
    let mut domain = Domain::new();
    let (block, _) = quad_grid(
        &mut domain,
        (2, 1),
        (2.0, 1.0),
        0.5,
        0.0, // massless: element mass on a rigid-link slave is a structured model error
        &PlaneMaterial::plane_stress(1000.0, 0.2).unwrap(),
        Quad4Formulation::Enhanced,
        |node, i, _| if i == 0 { node.fix(0).fix(1) } else { node },
    );
    let slave = block[block.len() - 1]; // top right corner (2, 1)
    let master = domain.add_node(Node::new([3.0, 1.0]));
    let tip = domain.add_node(Node::new([6.0, 1.0]).with_mass(0, 2.0).with_mass(1, 2.0));
    domain.add_element(Element::ElasticBeamColumn(
        ElasticBeamColumn::new(master, tip, 1000.0, 0.5, 0.04, GeomTransf::Linear)
            .with_density(1.0),
    ));
    domain.rigid_link(master, slave);

    let modes = modal_analysis(&mut domain, 3).expect("modal analysis");
    let (dx, dy) = (-1.0, 0.0); // slave - master
    let mut moved = false;
    for mode in &modes {
        let q = &mode.shape;
        let m = |dof| domain.value_at(q, master, dof);
        let s = |dof| domain.value_at(q, slave, dof);
        assert!((s(0) - (m(0) - m(2) * dy)).abs() < 1e-12);
        assert!((s(1) - (m(1) + m(2) * dx)).abs() < 1e-12);
        assert!((s(2) - m(2)).abs() < 1e-12);
        moved |= s(0).abs() + s(1).abs() > 1e-6;
    }
    assert!(moved, "the slave must actually move in some mode");
}

const L: f64 = 10.0;
const H: f64 = 1.0;

fn strip(nx: usize, ny: usize, clamped: bool) -> Domain {
    let mut domain = Domain::new();
    quad_grid(
        &mut domain,
        (nx, ny),
        (L, H),
        1.0,
        1.0,
        &PlaneMaterial::plane_stress(1e4, 0.0).unwrap(),
        Quad4Formulation::Enhanced,
        |node, i, _| {
            if clamped && i == 0 {
                node.fix(0).fix(1)
            } else {
                node
            }
        },
    );
    domain
}

/// Cantilever strip: the fundamental frequency approaches Euler-Bernoulli theory from
/// below (shear and rotary inertia lower it by about one percent at this slenderness),
/// and the total mass is `rho t A` in each direction.
#[test]
fn cantilever_strip_frequency_approaches_beam_theory_under_refinement() {
    let euler =
        1.875_104_068_7_f64.powi(2) * (1e4 * (H.powi(3) / 12.0) / (1.0 * H * L.powi(4))).sqrt();
    let mut errors = Vec::new();
    for (nx, ny) in [(5, 1), (10, 2), (20, 4), (40, 8)] {
        let mut domain = strip(nx, ny, true);
        let omega = modal_analysis(&mut domain, 1).unwrap()[0].frequency;
        errors.push((omega - euler) / euler);
    }
    assert!(errors.iter().all(|e| *e < 0.0), "{errors:?}");
    // The error is dominated by Timoshenko effects and settles; refinement never makes it worse.
    for pair in errors.windows(2) {
        assert!(pair[1].abs() <= pair[0].abs() + 1e-4, "{errors:?}");
    }
    assert!(errors[3].abs() < 0.02, "{errors:?}");

    // Total mass: a free strip has no supports, so the diagonal covers every DOF.
    let mut free = strip(10, 2, false);
    let mass = free.assemble_mass_diagonal();
    assert!(
        (mass.sum() - 2.0 * 1.0 * (L * H)).abs() < 1e-12,
        "{}",
        mass.sum()
    );
}

/// Undamped free vibration of a Q4 block under Newmark average acceleration conserves
/// `1/2 v.M.v + 1/2 u.K.u`.
#[test]
fn undamped_block_conserves_energy_under_newmark() {
    let mut domain = Domain::new();
    let (ids, elements) = quad_grid(
        &mut domain,
        (4, 4),
        (4.0, 4.0),
        1.0,
        1.0,
        &PlaneMaterial::plane_stress(100.0, 0.3).unwrap(),
        Quad4Formulation::Enhanced,
        |node, i, j| {
            let mut node = if i == 0 { node.fix(0).fix(1) } else { node };
            if i == 4 {
                node = node
                    .with_initial_velocity(0, 0.3 * (j as f64 + 1.0))
                    .with_initial_velocity(1, -0.2);
            }
            node
        },
    );
    let mass = domain.assemble_mass_diagonal();
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::NONE, 0.05).unwrap();
    let energy = |analysis: &TransientAnalysis<2, 3, NodeId, Element>| {
        let d = analysis.domain();
        let mut kinetic = 0.0;
        for &id in &ids {
            for dof in 0..2 {
                if let Some(eq) = d.equation_of(id, dof) {
                    kinetic += 0.5 * mass[eq] * d.node(id).velocity[dof].powi(2);
                }
            }
        }
        // Strain energy: 1/2 u . f_int, element by element (global DOF order, 4 nodes x 2).
        let mut strain = 0.0;
        let nodes_of = |index: usize| {
            let (i, j) = (index % 4, index / 4);
            [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)].map(|(a, b)| ids[b * 5 + a])
        };
        for (index, id) in elements.iter().copied().enumerate() {
            let force = d.element_local_force(id);
            for (n, node) in nodes_of(index).into_iter().enumerate() {
                for dof in 0..2 {
                    strain += 0.5 * d.node(node).displacement[dof] * force[2 * n + dof];
                }
            }
        }
        kinetic + strain
    };
    let initial = energy(&analysis);
    assert!(initial > 0.0);
    let mut peak_displacement = 0.0_f64;
    for _ in 0..200 {
        analysis.step().unwrap();
        peak_displacement =
            peak_displacement.max(analysis.domain().node(ids[24]).displacement[0].abs());
        let e = energy(&analysis);
        assert!((e - initial).abs() < 1e-9 * initial, "{e} vs {initial}");
    }
    assert!(
        peak_displacement > 1e-3,
        "the block should actually vibrate"
    );
}

fn loaded_block() -> (Domain, NodeId) {
    let mut domain = Domain::new();
    let (ids, _) = quad_grid(
        &mut domain,
        (3, 2),
        (3.0, 2.0),
        1.0,
        0.0,
        &PlaneMaterial::plane_stress(100.0, 0.3).unwrap(),
        Quad4Formulation::Full,
        |node, i, _| if i == 0 { node.fix(0).fix(1) } else { node },
    );
    for j in 0..=2 {
        domain.load_node(ids[j * 4 + 3], 1, -0.5);
        domain.load_node(ids[j * 4 + 3], 0, 0.2);
    }
    (domain, ids[11])
}

/// A linear continuum model through every nonlinear path needs one corrector iteration (or
/// none for the arc-length predictor): Newton in each tangent mode, with line search, and
/// the arc-length integrator.
#[test]
fn a_linear_continuum_model_converges_immediately_on_every_solver_path() {
    let test = ConvergenceTest::NormUnbalance {
        tol: 1e-9,
        max_iter: 10,
    };
    let (reference_domain, tip) = loaded_block();
    let mut reference = linear_analysis_loaded(reference_domain);
    reference.step().unwrap();
    let expected = reference.domain().node(tip).displacement;

    let line_search = Some(LineSearch::Bisection {
        tol: 1e-10,
        max_iter: 20,
        max_eta: 16.0,
    });
    for algorithm in [
        Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        },
        Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search,
        },
        Algorithm::Newton {
            tangent: TangentStrategy::Initial,
            line_search: None,
        },
    ] {
        let (domain, tip) = loaded_block();
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1.0 })
            .algorithm(algorithm)
            .test(test)
            .build(domain);
        let result = analysis.step().unwrap();
        assert!(
            result.iterations <= 2,
            "{algorithm:?}: {}",
            result.iterations
        );
        let got = analysis.domain().node(tip).displacement;
        assert!((got[0] - expected[0]).abs() < 1e-9 && (got[1] - expected[1]).abs() < 1e-9);
    }

    let (domain, tip) = loaded_block();
    let u_star = expected[1].abs();
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(ArcLength::fixed(
            0.5 * u_star,
            ArcScales::Explicit {
                displacement: u_star,
                rotation: None,
                load: 1.0,
            },
        )))
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(test)
        .build(domain);
    for _ in 0..3 {
        let result = analysis.step().unwrap();
        assert_eq!(result.arc.expect("arc diagnostics").corrector_solves, 0);
        // Linear response: displacement is proportional to the load factor.
        let u = analysis.domain().node(tip).displacement;
        assert!((u[1] - expected[1] * result.load_factor).abs() < 1e-9 * u_star);
    }
}

fn linear_analysis_loaded(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain)
}

/// A mass-bearing element on a rigid-link slave cannot be lumped onto a diagonal mass
/// vector; that is a `ModelError`, not a panic inside the mass assembly.
#[test]
fn element_mass_on_a_rigid_link_slave_is_a_structured_error() {
    let mut domain = Domain::new();
    let (block, _) = quad_grid(
        &mut domain,
        (1, 1),
        (1.0, 1.0),
        1.0,
        1.0,
        &PlaneMaterial::plane_stress(100.0, 0.2).unwrap(),
        Quad4Formulation::Full,
        |node, i, _| if i == 0 { node.fix(0).fix(1) } else { node },
    );
    let master = domain.add_node(Node::new([2.0, 1.0]).with_mass(0, 1.0).with_mass(1, 1.0));
    domain.rigid_link(master, block[3]);
    assert!(matches!(
        domain.validate(),
        Err(carapace_core::model::ModelError::MassOnConstrainedDof { .. })
    ));
    assert!(matches!(
        modal_analysis(&mut domain, 1),
        Err(carapace_core::analysis::AnalysisError::InvalidModel(_))
    ));
}
