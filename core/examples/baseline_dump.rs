//! Result dump for before/after comparison of the domain refactor
//! (`docs/domain-refactor-plan.md`, Phase 0).
//!
//! Builds a fixed set of representative models using only API that is meant
//! to survive the refactor (the `Domain`/`Analysis` aliases, element
//! constructors, `reaction`, `element_local_force`, `modal_analysis`,
//! `TransientAnalysis`), runs them, and prints every result as JSON on
//! stdout. Run it on the `pre-domain-refactor` tag and on the working branch
//! and diff the outputs with `comparison/compare_dumps.py`.
//!
//! Usage: `cargo run --release -p carapace-core --example baseline_dump > dump.json`
//!
//! Every model here fixes its unused DOFs explicitly, so it is valid before
//! and after DOF activation lands. Do not add models that depend on a
//! behavior the refactor intentionally changes.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, ArcLength, ArcScales, ConstraintHandler,
    ConvergenceTest, GroundMotion, Integrator, RayleighDamping, TangentStrategy, TransientAnalysis,
};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, Domain3, ElasticBeamColumn, ElasticBeamColumn3,
    Element, Element3, Fiber, ForceBeamColumn, Friction, GeomTransf, GeomTransf3, LoadSeries,
    Material, Node, Node3, NodeId, Orientation, SpatialDof, Truss, ZeroLength, ZeroLength3,
};

/// Minimal JSON writer: a flat map from a dotted key to an array of numbers.
struct Dump {
    entries: Vec<(String, Vec<f64>)>,
}

impl Dump {
    fn new() -> Self {
        Dump {
            entries: Vec::new(),
        }
    }
    fn put(&mut self, key: impl Into<String>, values: impl IntoIterator<Item = f64>) {
        self.entries
            .push((key.into(), values.into_iter().collect()));
    }
    fn print(&self) {
        println!("{{");
        for (i, (key, values)) in self.entries.iter().enumerate() {
            let body: Vec<String> = values.iter().map(|v| format!("{v:.17e}")).collect();
            let comma = if i + 1 == self.entries.len() { "" } else { "," };
            println!("  \"{key}\": [{}]{comma}", body.join(", "));
        }
        println!("}}");
    }
}

fn newton() -> Algorithm {
    Algorithm::Newton {
        tangent: TangentStrategy::Current,
        line_search: None,
    }
}

fn norm_test(tol: f64) -> ConvergenceTest {
    ConvergenceTest::NormUnbalance { tol, max_iter: 40 }
}

// ---------------------------------------------------------------------------
// 1. Two-story 2D frame, P-Delta, rigid floor diaphragms, two-phase loading.
// ---------------------------------------------------------------------------

fn frame2d(dump: &mut Dump) {
    let (e, a, iz) = (200_000.0, 8_000.0, 1.0e8);
    let mut domain = Domain::new();
    let mut nodes: Vec<NodeId> = Vec::new();
    for (i, (x, y)) in [
        (0.0, 0.0),
        (6000.0, 0.0),
        (0.0, 4000.0),
        (6000.0, 4000.0),
        (0.0, 8000.0),
        (6000.0, 8000.0),
    ]
    .into_iter()
    .enumerate()
    {
        let node = if i < 2 {
            Node::new([x, y]).fix(0).fix(1).fix(2)
        } else {
            Node::new([x, y])
        };
        nodes.push(domain.add_node(node));
    }
    let mut elements = Vec::new();
    let connect = |domain: &mut Domain, i: usize, j: usize, elements: &mut Vec<_>| {
        elements.push(
            domain.add_element(Element::ElasticBeamColumn(
                ElasticBeamColumn::new(nodes[i], nodes[j], e, a, iz, GeomTransf::PDelta)
                    .with_density(7.85e-9),
            )),
        );
    };
    for (i, j) in [(0, 2), (1, 3), (2, 4), (3, 5), (2, 3), (4, 5)] {
        connect(&mut domain, i, j, &mut elements);
    }
    domain.rigid_diaphragm(nodes[2], &[nodes[3]]);
    domain.rigid_diaphragm(nodes[4], &[nodes[5]]);

    // Phase 1: gravity on the default pattern.
    for &n in &nodes[2..] {
        domain.load_node(n, 1, -400_000.0);
    }
    let gravity_pattern = domain.default_pattern();

    let mut analysis: Analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.25 })
        .algorithm(newton())
        .test(norm_test(1e-6))
        .build(domain);

    let mut last = 0.0;
    for step in 0..4 {
        let r = analysis.step().expect("gravity phase converges");
        dump.put(
            format!("frame2d.gravity.step{step}"),
            [r.load_factor, r.iterations as f64, r.factorizations as f64],
        );
        last = r.load_factor;
    }
    analysis
        .domain_mut()
        .hold_pattern_constant(gravity_pattern, last);

    // Phase 2: lateral pattern ramps while gravity stays frozen.
    let lateral = analysis
        .domain_mut()
        .add_load_pattern(LoadSeries::Linear { slope: 1.0 });
    analysis
        .domain_mut()
        .add_nodal_load(lateral, nodes[2], 0, 50_000.0);
    analysis
        .domain_mut()
        .add_nodal_load(lateral, nodes[4], 0, 100_000.0);
    analysis.set_integrator(Integrator::LoadControl { increment: 0.25 });
    for step in 0..6 {
        let r = analysis.step().expect("lateral phase converges");
        dump.put(
            format!("frame2d.lateral.step{step}"),
            [r.load_factor, r.iterations as f64, r.factorizations as f64],
        );
        last = r.load_factor;
    }

    for (i, &n) in nodes.iter().enumerate() {
        dump.put(
            format!("frame2d.disp.node{i}"),
            analysis.domain().node(n).displacement,
        );
        for dof in 0..3 {
            dump.put(
                format!("frame2d.reaction.node{i}.dof{dof}"),
                [analysis.domain().reaction(n, dof, last)],
            );
        }
    }
    for (i, &el) in elements.iter().enumerate() {
        let f = analysis.domain().element_local_force(el);
        dump.put(format!("frame2d.force.el{i}"), (0..6).map(|c| f[c]));
    }

    // Modal analysis of the same frame (mass from element density).
    let mut domain = analysis.into_domain();
    let modes = modal_analysis(&mut domain, 3).expect("frame modes");
    for (m, mode) in modes.iter().enumerate() {
        dump.put(format!("frame2d.modal.freq{m}"), [mode.frequency]);
        let sign = mode
            .shape
            .iter()
            .copied()
            .max_by(|x, y| x.abs().partial_cmp(&y.abs()).unwrap())
            .map(|v| v.signum())
            .unwrap_or(1.0);
        dump.put(
            format!("frame2d.modal.shape{m}"),
            mode.shape.iter().map(|v| v * sign),
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Fiber cantilevers (displacement- and force-based), displacement control.
// ---------------------------------------------------------------------------

fn fibers() -> Vec<Fiber> {
    let steel = || Material::steel01(355.0, 200_000.0, 0.02, 0.0, 1.0, 0.0, 1.0);
    (0..6)
        .map(|i| {
            let y = -150.0 + 60.0 * i as f64;
            Fiber::new(y, 1000.0, steel())
        })
        .collect()
}

fn fiber_cantilever(dump: &mut Dump, force_based: bool) {
    let tag = if force_based { "forcebeam" } else { "dispbeam" };
    let mut domain = Domain::new();
    let base = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let tip = domain.add_node(Node::new([0.0, 2000.0]));
    let integration = BeamIntegration::Lobatto { points: 5 };
    let element = if force_based {
        domain.add_element(Element::ForceBeamColumn(ForceBeamColumn::new(
            base,
            tip,
            fibers(),
            integration,
        )))
    } else {
        domain.add_element(Element::DispBeamColumn(DispBeamColumn::new(
            base,
            tip,
            fibers(),
            integration,
        )))
    };
    domain.load_node(tip, 0, 1.0);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::DisplacementControl {
            node: tip,
            dof: 0,
            increment: 2.0,
        })
        .algorithm(newton())
        .test(norm_test(1e-6))
        .build(domain);
    for step in 0..30 {
        if step == 15 {
            analysis.set_integrator(Integrator::DisplacementControl {
                node: tip,
                dof: 0,
                increment: -2.0,
            });
        }
        let r = analysis.step().expect("fiber cantilever converges");
        dump.put(
            format!("{tag}.step{step}"),
            [
                r.load_factor,
                r.iterations as f64,
                r.factorizations as f64,
                analysis.domain().node(tip).displacement[0],
            ],
        );
    }
    let f = analysis.domain().element_local_force(element);
    dump.put(format!("{tag}.force"), (0..6).map(|c| f[c]));
    dump.put(
        format!("{tag}.reaction"),
        (0..3).map(|d| analysis.domain().reaction(base, d, 0.0)),
    );
    if let Some(points) = analysis.domain().element_fiber_responses(element) {
        for (p, fibers) in points.iter().enumerate() {
            dump.put(
                format!("{tag}.fibers.point{p}"),
                fibers.iter().flat_map(|&(s, t)| [s, t]),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 3. Zero-length springs: oriented local axes, and Coulomb friction.
// ---------------------------------------------------------------------------

fn zero_length(dump: &mut Dump) {
    // Oriented elastic springs under a global load.
    let (c, s) = (30.0_f64.to_radians().cos(), 30.0_f64.to_radians().sin());
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(2));
    let orientation = Orientation::new([c, s, 0.0], [-s, c, 0.0]).expect("valid orientation");
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, free)
            .with_material(0, Material::Elastic { e: 400.0 })
            .with_material(1, Material::Elastic { e: 100.0 })
            .with_orientation(orientation)
            .expect("2D orientation"),
    ));
    domain.load_node(free, 0, 30.0);
    domain.load_node(free, 1, -10.0);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.5 })
        .algorithm(Algorithm::Linear)
        .test(norm_test(1e-9))
        .build(domain);
    for _ in 0..2 {
        analysis.step().expect("oriented spring solves");
    }
    dump.put(
        "zerolength.oriented.disp",
        analysis.domain().node(free).displacement,
    );

    // Friction: normal spring + sliding shear, cyclic shear via displacement control.
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let slider = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, slider)
            .with_material(0, Material::Elastic { e: 1000.0 })
            .with_friction(Friction::new(0, 1, 0.3, 500.0, 0.01)),
    ));
    domain.load_node(slider, 1, 1.0);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::DisplacementControl {
            node: slider,
            dof: 1,
            increment: 0.1,
        })
        .algorithm(newton())
        .test(norm_test(1e-8))
        .build(domain);
    // Normal force: prescribed compression through the fixed normal DOF is
    // not available here, so the friction element sees zero normal force and
    // stays on its (degenerate) yield surface; still a deterministic path.
    let mut ok = Vec::new();
    for step in 0..12 {
        if step == 6 {
            analysis.set_integrator(Integrator::DisplacementControl {
                node: slider,
                dof: 1,
                increment: -0.1,
            });
        }
        match analysis.step() {
            Ok(r) => ok.push([r.load_factor, r.iterations as f64]),
            Err(_) => break,
        }
    }
    dump.put("zerolength.friction.steps", ok.iter().flatten().copied());
}

// ---------------------------------------------------------------------------
// 4. Transient: linear MDOF with ground motion; nonlinear with Newton corrector.
// ---------------------------------------------------------------------------

fn shear_building(springs: [Material; 3]) -> (Domain, [NodeId; 3]) {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let mut prev = ground;
    let mut masses = [ground; 3];
    for (i, spring) in springs.into_iter().enumerate() {
        let node = domain.add_node(
            Node::new([0.0, 3.0 * (i as f64 + 1.0)])
                .fix(1)
                .fix(2)
                .with_mass(0, 2.0 - 0.3 * i as f64),
        );
        domain.add_element(Element::ZeroLength(
            ZeroLength::new(prev, node).with_material(0, spring),
        ));
        masses[i] = node;
        prev = node;
    }
    (domain, masses)
}

fn accelerogram() -> LoadSeries {
    LoadSeries::Path {
        times: vec![0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2],
        factors: vec![0.0, 2.0, -1.5, 1.0, -2.5, 0.5, 0.0],
    }
}

fn transient(dump: &mut Dump) {
    let elastic = |e| Material::Elastic { e };
    let (domain, masses) = shear_building([elastic(900.0), elastic(700.0), elastic(500.0)]);
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::new(0.4, 0.002), 0.01)
        .expect("building is well-posed")
        .with_ground_motion(GroundMotion::new(0, accelerogram()));
    for step in 0..150 {
        let r = analysis.step().expect("linear transient step");
        if step % 10 == 9 {
            dump.put(
                format!("transient.linear.step{step}"),
                masses
                    .iter()
                    .flat_map(|&n| {
                        let node = analysis.domain().node(n);
                        [node.displacement[0], node.velocity[0], node.acceleration[0]]
                    })
                    .chain([r.iterations as f64]),
            );
        }
    }

    // Nonlinear: elastic-perfectly-plastic springs, Newton corrector.
    let epp = |e: f64, fy: f64| Material::elastic_pp(e, fy / e);
    let (domain, masses) = shear_building([epp(900.0, 12.0), epp(700.0, 9.0), epp(500.0, 6.0)]);
    let mut analysis = TransientAnalysis::new(domain, RayleighDamping::new(0.2, 0.0), 0.01)
        .expect("nonlinear building is well-posed")
        .with_algorithm(newton(), norm_test(1e-9))
        .with_ground_motion(GroundMotion::new(0, accelerogram()).with_scale_factor(1.5));
    for step in 0..150 {
        let r = analysis.step().expect("nonlinear transient step");
        if step % 10 == 9 {
            dump.put(
                format!("transient.nonlinear.step{step}"),
                masses
                    .iter()
                    .map(|&n| analysis.domain().node(n).displacement[0])
                    .chain([r.iterations as f64, r.factorizations as f64]),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 5. Arc length through a degrading peak and snap-back (series springs).
// ---------------------------------------------------------------------------

fn arc_length(dump: &mut Dump) {
    let degrading = Material::hysteretic(
        10.0, 0.01, 6.0, 0.02, 2.0, 0.03, -10.0, -0.01, -6.0, -0.02, -2.0, -0.03, 1.0, 1.0, 0.0,
        0.0, 0.0,
    );
    let mut domain = Domain::new();
    let fixed = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let middle = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    let free = domain.add_node(Node::new([0.0, 0.0]).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(fixed, middle).with_material(0, degrading),
    ));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(middle, free).with_material(0, Material::Elastic { e: 500.0 }),
    ));
    domain.load_node(free, 0, 1.0);
    let mut arc = ArcLength::fixed(
        0.05,
        ArcScales::Explicit {
            displacement: 0.01,
            rotation: None,
            load: 10.0,
        },
    );
    arc.arc_tolerance = 1e-9;
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(arc))
        .algorithm(newton())
        .test(ConvergenceTest::Combined {
            force_tol: 1e-9,
            moment_tol: 1e-9,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 30,
        })
        .build(domain);
    for step in 0..60 {
        let Ok(r) = analysis.step() else { break };
        let info = r.arc.expect("arc diagnostics");
        dump.put(
            format!("arc.step{step}"),
            [
                r.load_factor,
                analysis.domain().node(free).displacement[0],
                analysis.domain().node(middle).displacement[0],
                r.iterations as f64,
                info.factorizations as f64,
                info.radius,
                info.retries as f64,
            ],
        );
    }
}

// ---------------------------------------------------------------------------
// 6. 3D frame with rigid diaphragm (static and modal), plus a 3D spring.
// ---------------------------------------------------------------------------

fn frame3d(dump: &mut Dump) {
    let all = |n: Node3| {
        n.fix(SpatialDof::Ux as usize)
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Uz as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize)
    };
    let (e, g, a, j, iy, iz) = (200_000.0, 80_000.0, 8_000.0, 2.0e7, 5.0e7, 1.0e8);
    let mut domain = Domain3::new();
    let corners = [(0.0, 0.0), (6000.0, 0.0), (0.0, 6000.0), (6000.0, 6000.0)];
    let mut top = Vec::new();
    let mut elements = Vec::new();
    let mut bases = Vec::new();
    for &(x, z) in &corners {
        let base = domain.add_node(all(Node3::new([x, 0.0, z])));
        let t = domain.add_node(Node3::new([x, 4000.0, z]));
        elements.push(
            domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
                base,
                t,
                e,
                g,
                a,
                j,
                iy,
                iz,
                GeomTransf3::Linear3 {
                    vec_xz: [1.0, 0.0, 0.0],
                },
            ))),
        );
        bases.push(base);
        top.push(t);
    }
    let master = domain.add_node(
        Node3::new([3000.0, 4000.0, 3000.0])
            .fix(SpatialDof::Uy as usize)
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Rz as usize)
            // Mass lives on the diaphragm master: lumped mass on a
            // diaphragm-constrained slave DOF is unsupported.
            .with_mass(SpatialDof::Ux as usize, 1.0e3)
            .with_mass(SpatialDof::Uz as usize, 1.0e3)
            .with_mass(SpatialDof::Ry as usize, 1.0e9),
    );
    domain.rigid_diaphragm(master, &top);
    domain.load_node(master, SpatialDof::Ux as usize, 100_000.0);
    domain.load_node(master, SpatialDof::Ry as usize, 5.0e7);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(norm_test(1e-9))
        .build(domain);
    analysis.step().expect("3D frame solves");
    dump.put(
        "frame3d.master",
        analysis.domain().node(master).displacement,
    );
    for (i, &t) in top.iter().enumerate() {
        dump.put(
            format!("frame3d.top{i}"),
            analysis.domain().node(t).displacement,
        );
    }
    for (i, &el) in elements.iter().enumerate() {
        let f = analysis.domain().element_local_force(el);
        dump.put(format!("frame3d.force.el{i}"), (0..12).map(|c| f[c]));
    }
    for (i, &b) in bases.iter().enumerate() {
        dump.put(
            format!("frame3d.reaction.base{i}"),
            (0..6).map(|d| analysis.domain().reaction(b, d, 1.0)),
        );
    }

    let mut domain = analysis.into_domain();
    let modes = modal_analysis(&mut domain, 3).expect("3D frame modes");
    for (m, mode) in modes.iter().enumerate() {
        dump.put(format!("frame3d.modal.freq{m}"), [mode.frequency]);
    }

    // A 3D zero-length spring and truss chain, for the 3D catalog.
    let mut domain = Domain3::new();
    let ground = domain.add_node(all(Node3::new([0.0, 0.0, 0.0])));
    let free = domain.add_node(
        Node3::new([1.0, 0.0, 0.0])
            .fix(SpatialDof::Rx as usize)
            .fix(SpatialDof::Ry as usize)
            .fix(SpatialDof::Rz as usize),
    );
    domain.add_element(Element3::ZeroLength3(
        ZeroLength3::new(ground, free)
            .with_material(0, Material::Elastic { e: 300.0 })
            .with_material(1, Material::Elastic { e: 200.0 })
            .with_material(2, Material::Elastic { e: 100.0 }),
    ));
    domain.load_node(free, 0, 30.0);
    domain.load_node(free, 1, 20.0);
    domain.load_node(free, 2, 10.0);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(norm_test(1e-9))
        .build(domain);
    analysis.step().expect("3D spring solves");
    dump.put("spring3d.disp", analysis.domain().node(free).displacement);
}

// ---------------------------------------------------------------------------
// 7. Truss with equal_dof and a prescribed support displacement.
// ---------------------------------------------------------------------------

fn truss_and_prescribed(dump: &mut Dump) {
    let mut domain = Domain::new();
    let a = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let b = domain.add_node(Node::new([3000.0, 0.0]).fix(1).fix(2));
    let c = domain.add_node(Node::new([1500.0, 2000.0]).fix(2));
    let d = domain.add_node(Node::new([1500.0, 2000.0]).fix(2));
    for (i, j) in [(a, b), (b, c), (a, c)] {
        domain.add_element(Element::Truss(Truss::new(
            i,
            j,
            500.0,
            Material::Elastic { e: 200_000.0 },
        )));
    }
    domain.equal_dof(c, d, &[0, 1]);
    domain.load_node(d, 1, -1.0e4);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.5 })
        .algorithm(newton())
        .test(norm_test(1e-9))
        .build(domain);
    for _ in 0..2 {
        analysis.step().expect("truss solves");
    }
    for (i, n) in [a, b, c, d].into_iter().enumerate() {
        dump.put(
            format!("truss.disp.node{i}"),
            analysis.domain().node(n).displacement,
        );
        dump.put(
            format!("truss.reaction.node{i}"),
            (0..2).map(|dof| analysis.domain().reaction(n, dof, 1.0)),
        );
    }

    // Prescribed settlement of a fixed support.
    let mut domain = Domain::new();
    let left = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let right = domain.add_node(Node::new([4000.0, 0.0]).fix(1).fix(2));
    let support = domain.add_node(Node::new([8000.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::Truss(Truss::new(
        left,
        right,
        500.0,
        Material::Elastic { e: 200_000.0 },
    )));
    domain.add_element(Element::Truss(Truss::new(
        right,
        support,
        300.0,
        Material::Elastic { e: 200_000.0 },
    )));
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(newton())
        .test(norm_test(1e-9))
        .build(domain);
    analysis
        .step_prescribed(&[(support, 0, 5.0)])
        .expect("settlement step");
    dump.put(
        "truss.settlement.disp",
        [analysis.domain().node(right).displacement[0]],
    );
    dump.put(
        "truss.settlement.reaction",
        [analysis.domain().reaction(support, 0, 0.0)],
    );
}

fn main() {
    let mut dump = Dump::new();
    frame2d(&mut dump);
    fiber_cantilever(&mut dump, false);
    fiber_cantilever(&mut dump, true);
    zero_length(&mut dump);
    transient(&mut dump);
    arc_length(&mut dump);
    frame3d(&mut dump);
    truss_and_prescribed(&mut dump);
    dump.print();
}
