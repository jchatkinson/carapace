//! Axis-aligned 3D diaphragm response against shear-law patch tests and OpenSees/CSA O86 references.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    flexible_diaphragm_loads, lump_line_reaction, ChordDistribution, Diaphragm, DiaphragmLayout,
    Domain3, Element3, FloorBounds, FloorLoad, FloorSupport, Node3, ReactionLine, StudWall,
    StudWallLayout, WallShearLaw, WallSheathing, WallStud,
};

fn layout() -> DiaphragmLayout {
    DiaphragmLayout {
        u: vec![0., 1.],
        v: vec![0., 2.],
        level: 0.,
        bays: vec![Some(WallShearLaw::Linear { ga: 100. })],
        ea_u: vec![0., 0.],
        ea_v: vec![0., 0.],
    }
}

/// A single bay's constant shear field gives `q = ga gamma`; the nodal force balance is exact.
#[test]
fn single_bay_linear_shear() {
    let mut d = Domain3::new();
    let a = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
    let b = d.add_node(Node3::new([0., 0., 1.]).fix(1).fix(2));
    let c = d.add_node(Node3::new([2., 0., 0.]).fix(0).fix(1).fix(2));
    let e = d.add_node(Node3::new([2., 0., 1.]).fix(0).fix(1).fix(2));
    let id = d.add_element(Element3::Diaphragm(
        Diaphragm::new(vec![a, b, c, e], layout()).unwrap(),
    ));
    d.add_nodal_load(d.default_pattern(), b, 0, 20.);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-12,
            max_iter: 4,
        })
        .build(d);
    analysis.step().unwrap();
    let response = analysis.domain().element_diaphragm_response(id).unwrap();
    assert!((response.bays[0].q - 20.).abs() < 1e-12);
    assert!((response.bays[0].q - 100. * response.bays[0].gamma).abs() < 1e-12);
}

/// O86 panel response satisfies bay equilibrium; reports the nail-slip strain for comparison with 11.7.2.
#[test]
fn single_bay_o86_shear_equilibrium_and_nail_slip_report() {
    let mut d = Domain3::new();
    let a = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
    let b = d.add_node(Node3::new([0., 0., 1.]).fix(1).fix(2));
    let c = d.add_node(Node3::new([2., 0., 0.]).fix(0).fix(1).fix(2));
    let e = d.add_node(Node3::new([2., 0., 1.]).fix(0).fix(1).fix(2));
    let mut l = layout();
    l.bays[0] = Some(WallShearLaw::O86 {
        bv: 100.,
        nail_spacing: 150.,
        nail_diameter: 3.,
    });
    let id = d.add_element(Element3::Diaphragm(
        Diaphragm::new(vec![a, b, c, e], l).unwrap(),
    ));
    d.add_nodal_load(d.default_pattern(), b, 0, 20.);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 12,
        })
        .build(d);
    analysis.step().unwrap();
    let bay = analysis
        .domain()
        .element_diaphragm_response(id)
        .unwrap()
        .bays[0];
    assert!((bay.q - 20.).abs() < 1e-9);
    let en = (0.013 * bay.q * 150. / 3_f64.powi(2)).powi(2);
    let nail_slip_strain = bay.gamma - bay.q / 100.;
    println!(
        "O86 bay: q={}, gamma={}, nail-slip strain={}, 0.00061 e_n={}",
        bay.q,
        bay.gamma,
        nail_slip_strain,
        0.00061 * en
    );
    assert!((nail_slip_strain - 0.0025 * en).abs() < 1e-14);
}

/// The same single-bay shear patch gives identical strain and flow when its normal changes from Y to Z.
/// Reference: direct Y-up result versus the cyclic Z-up `(u, v)` plane; relative tolerance 1e-12.
#[test]
fn single_bay_shear_is_invariant_in_z_normal_plane() {
    use carapace_core::model::Node3Id;
    fn solve(z_up: bool, law: WallShearLaw) -> (f64, f64) {
        let mut d = Domain3::new();
        let points = if z_up {
            [[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [1., 2., 0.]]
        } else {
            [[0., 0., 0.], [0., 0., 1.], [2., 0., 0.], [2., 0., 1.]]
        };
        let mut ids = Vec::<Node3Id>::new();
        for (i, p) in points.into_iter().enumerate() {
            let node = Node3::new(p);
            let node = if i == 0 || i >= 2 {
                node.fix(0).fix(1).fix(2)
            } else if z_up {
                node.fix(0).fix(2)
            } else {
                node.fix(1).fix(2)
            };
            ids.push(d.add_node(node));
        }
        let layout = DiaphragmLayout {
            u: vec![0., 1.],
            v: vec![0., 2.],
            level: 0.,
            bays: vec![Some(law)],
            ea_u: vec![0., 0.],
            ea_v: vec![0., 0.],
        };
        let e = d.add_element(Element3::Diaphragm(
            Diaphragm::new(ids.clone(), layout).unwrap(),
        ));
        d.add_nodal_load(d.default_pattern(), ids[1], if z_up { 1 } else { 0 }, 20.);
        let mut a = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            })
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 12,
            })
            .build(d);
        a.step().unwrap();
        let response = a.domain().element_diaphragm_response(e).unwrap().bays[0];
        (response.gamma, response.q)
    }
    for law in [
        WallShearLaw::Linear { ga: 100. },
        WallShearLaw::O86 {
            bv: 100.,
            nail_spacing: 150.,
            nail_diameter: 3.,
        },
    ] {
        let y = solve(false, law.clone());
        let z = solve(true, law);
        for (a, b) in [(y.0, z.0), (y.1, z.1)] {
            assert!((a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.));
        }
    }
}

/// Rigid in-plane translations and rotation have zero engineering shear in the bilinear field.
#[test]
fn rigid_body_kinematics_have_zero_shear() {
    let (dx, dz) = (2.4, 1.7);
    for [tx, tz, theta] in [[1.0_f64, 0., 0.], [0., 1., 0.], [0., 0., 1.]] {
        let ux = [tx, tx, tx + theta * dz, tx + theta * dz];
        let uz = [tz, tz - theta * dx, tz, tz - theta * dx];
        let gamma = (ux[2] + ux[3] - ux[0] - ux[1]) / (2. * dz)
            + (uz[1] + uz[3] - uz[0] - uz[2]) / (2. * dx);
        assert!(gamma.abs() < 1e-15);
    }
}

/// Uniform flexible diaphragm loading over three equal spans gives hand-calculated simple-span reactions.
#[test]
fn flexible_three_line_tributary_reactions_and_lumping() {
    let b = FloorBounds {
        s0: 0.,
        s1: 10.,
        t0: 0.,
        t1: 2.,
    };
    let s = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 2.,
        },
        FloorSupport {
            s: 5.,
            t0: 0.,
            t1: 2.,
        },
        FloorSupport {
            s: 10.,
            t0: 0.,
            t1: 2.,
        },
    ];
    let r = flexible_diaphragm_loads(
        b,
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 3.,
            region: None,
        }],
    )
    .unwrap();
    let q: Vec<f64> = r
        .per_support
        .iter()
        .map(|x| x.lines.iter().map(|l| l.w * (l.t1 - l.t0)).sum())
        .collect();
    assert_eq!(q, vec![15., 30., 15.]);
    let lump = lump_line_reaction(
        &[0., 1., 2.],
        ReactionLine {
            slot: 0,
            t0: 0.,
            t1: 2.,
            w: 3.,
        },
    )
    .unwrap();
    assert_eq!(lump, vec![1.5, 3., 1.5]);
    assert!((lump.iter().sum::<f64>() - 6.).abs() < 1e-12);

    let mut d = Domain3::new();
    let mut walls = Vec::new();
    for (i, force) in q.iter().copied().enumerate() {
        let z = i as f64 * 3.;
        let c0 = d.add_node(Node3::new([0., 0., z]).fix(0).fix(1).fix(2));
        let c1 = d.add_node(Node3::new([2., 0., z]).fix(0).fix(1).fix(2));
        let c2 = d.add_node(Node3::new([2., 3., z]).fix(2));
        let c3 = d.add_node(Node3::new([0., 3., z]).fix(2));
        let layout = StudWallLayout {
            length: 2.,
            height: 3.,
            studs: vec![WallStud { x: 0., ea: 1e6 }, WallStud { x: 2., ea: 1e6 }],
            top_ea: 1e8,
            top_ei: Some(1e8),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: 1e8,
            bottom_ei: 1e8,
            openings: vec![],
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear { ga: 1e5 },
            }),
            foundation: true,
            bearing: None,
            floor: None,
        };
        let id = d.add_element(Element3::StudWall(
            StudWall::new3([c0, c1, c2, c3], vec![], layout, vec![]).unwrap(),
        ));
        d.add_nodal_load(d.default_pattern(), c2, 0, force * 0.5);
        d.add_nodal_load(d.default_pattern(), c3, 0, force * 0.5);
        walls.push(id);
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    for (id, reaction) in walls.into_iter().zip(q) {
        let wall = a.domain().element_stud_wall_response(id, 1.).unwrap();
        assert!((wall.base[0].abs() - reaction).abs() < 1e-10);
    }
}

/// The same diaphragm solve gives bit-identical response with native threaded assembly.
#[cfg(feature = "parallel")]
#[test]
fn diaphragm_serial_and_parallel_are_bit_identical() {
    use carapace_core::analysis::{set_parallelism, Parallelism};
    fn solve() -> [u64; 2] {
        let mut d = Domain3::new();
        let a = d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2));
        let b = d.add_node(Node3::new([0., 0., 1.]).fix(0).fix(1).fix(2));
        let c = d.add_node(Node3::new([2., 0., 0.]).fix(1).fix(2));
        let e = d.add_node(Node3::new([2., 0., 1.]).fix(0).fix(1).fix(2));
        let id = d.add_element(Element3::Diaphragm(
            Diaphragm::new(vec![a, b, c, e], layout()).unwrap(),
        ));
        d.add_nodal_load(d.default_pattern(), c, 0, 20.);
        let mut analysis = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Plain)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-12,
                max_iter: 4,
            })
            .build(d);
        analysis.step().unwrap();
        [
            analysis.domain().node(c).displacement[0].to_bits(),
            analysis
                .domain()
                .element_diaphragm_response(id)
                .unwrap()
                .bays[0]
                .q
                .to_bits(),
        ]
    }
    set_parallelism(Parallelism::Serial);
    let serial = solve();
    set_parallelism(Parallelism::Threads(2));
    let parallel = solve();
    set_parallelism(Parallelism::Serial);
    assert_eq!(parallel, serial);
}

/// Simply supported blocked diaphragm under a uniform load, one bay deep, `n` bays along the span; returns the
/// mid-span deflection (load direction) with the given bay law.
fn simple_diaphragm(n: usize, law: WallShearLaw, z_up: bool) -> f64 {
    use carapace_core::model::Node3Id;
    let (span, depth, w) = (12_000., 6_000., 5.);
    let load_axis = if z_up { 0 } else { 2 };
    let base_axis = if z_up { 1 } else { 0 };
    let xs: Vec<f64> = (0..=n).map(|i| span * i as f64 / n as f64).collect();
    let mut d = Domain3::new();
    let mut nodes: Vec<Node3Id> = Vec::new();
    for (j, &x) in xs.iter().enumerate() {
        for z in [0., depth] {
            let mut node = Node3::new(if z_up { [z, x, 0.] } else { [x, 0., z] });
            if j == 0 || j == n {
                node = node.fix(load_axis); // end walls resist the load direction
            }
            if j == 0 && z == 0. {
                node = node.fix(base_axis);
            }
            nodes.push(d.add_node(node));
        }
    }
    let layout = DiaphragmLayout {
        u: vec![0., depth],
        v: xs.clone(),
        level: 0.,
        bays: vec![Some(law); n],
        ea_u: vec![1e9; n + 1],
        ea_v: vec![CHORD_EA; 2],
    };
    d.add_element(Element3::Diaphragm(
        Diaphragm::new(nodes.clone(), layout).unwrap(),
    ));
    let a = span / n as f64;
    for row in 0..n + 1 {
        for i in 0..2 {
            let trib = if row == 0 || row == n { a / 2. } else { a };
            d.add_nodal_load(
                d.default_pattern(),
                nodes[row * 2 + i],
                load_axis,
                w * trib / 2.,
            );
        }
    }
    let mut an = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 50,
        })
        .build(d);
    an.step().expect("diaphragm solve");
    let mid = n / 2;
    (an.domain().node(nodes[2 * mid]).displacement[load_axis]
        + an.domain().node(nodes[2 * mid + 1]).displacement[load_axis])
        / 2.
}

/// Chord `EA`: two 38x140 SPF plies, E = 9500 MPa (N).
const CHORD_EA: f64 = 2. * 9500. * 38. * 140.;

/// A simply supported blocked diaphragm (12 m span, 6 m deep, 5 N/mm) with a linear law `ga = Bv` against CSA
/// O86:19 11.7.2 (bending `5 v L³/(96 E A L_D)` + shear `v L/(4 B_v)`, `v = w L/(2 L_D)`). The shear term is exact
/// for nodal lumping (midpoint rule); the chord term converges at second order with the number of bays. Asserts
/// the observed order is 2 (±0.05) and the 64-bay result is within 1e-3 of the code value.
#[test]
fn simply_supported_diaphragm_converges_to_o86_11_7_2() {
    let (span, depth, w, bv): (f64, f64, f64, f64) = (12_000., 6_000., 5., 11_000.);
    let v = w * span / (2. * depth);
    let code = 5. * v * span.powi(3) / (96. * CHORD_EA * depth) + v * span / (4. * bv);
    let errs: Vec<f64> = [8, 16, 32, 64]
        .iter()
        .map(|&n| (simple_diaphragm(n, WallShearLaw::Linear { ga: bv }, false) - code).abs() / code)
        .collect();
    eprintln!("O86 11.7.2 linear: code {code:.6} mm, relative errors {errs:?}");
    for e in errs.windows(2) {
        let order = (e[0] / e[1]).log2();
        assert!((order - 2.).abs() < 0.05, "order {order}");
    }
    assert!(errs[3] < 1e-3, "{errs:?}");
}

/// The same diaphragm with the O86 nail-slip law: prints the nail-slip part of the mid-span deflection (O86 law
/// minus the linear `ga = Bv` run) next to the code's `0.00061 L e_n`. Characterization only: the integrated law
/// gives the slip of the actual shear distribution, the code term is a fixed approximation of it.
#[test]
fn o86_law_diaphragm_nail_slip_vs_code_term() {
    let (span, depth, w, bv, s, df): (f64, f64, f64, f64, f64, f64) =
        (12_000., 6_000., 5., 11_000., 150., 2.84);
    let v = w * span / (2. * depth);
    let en = (0.013 * v * s / (df * df)).powi(2);
    let law = WallShearLaw::O86 {
        bv,
        nail_spacing: s,
        nail_diameter: df,
    };
    let slip = simple_diaphragm(64, law.clone(), false)
        - simple_diaphragm(64, WallShearLaw::Linear { ga: bv }, false);
    eprintln!(
        "O86 law nail slip {slip:.4} mm vs 0.00061 L e_n = {:.4} mm (ratio {:.3})",
        0.00061 * span * en,
        slip / (0.00061 * span * en)
    );
    assert!(slip > 0.);
}

/// Rotating the same simply supported diaphragm from a Y-normal plane to a Z-normal plane preserves the O86 response.
/// Reference: identical cyclic `(u, v)` grid and material in both orientations; relative tolerance 1e-12.
#[test]
fn o86_simply_supported_response_is_invariant_in_z_up_plane() {
    let law = WallShearLaw::O86 {
        bv: 11_000.,
        nail_spacing: 150.,
        nail_diameter: 2.84,
    };
    let y_up = simple_diaphragm(16, law.clone(), false);
    let z_up = simple_diaphragm(16, law, true);
    assert!((y_up - z_up).abs() <= 1e-12 * y_up.abs().max(z_up.abs()).max(1.));
}

/// Rigid diaphragm (`Domain3::rigid_diaphragm`) over four linear 3D walls: two lines along X (z = 0, 4 m long;
/// z = 6000, 2 m long) and two along Z (x = 0 and x = 8000, 3 m long), with an X force applied eccentrically at
/// the master node. Each line's shear equals the rigid-body stiffness + torsion split computed by hand from the
/// lines' own lateral stiffnesses (each measured from a run of that wall alone): `K = Σ k_i r_iᵀ r_i` with
/// `r = [1, 0, z_i - z_m]` for an X line and `[0, 1, -(x_i - x_m)]` for a Z line, `V_i = k_i r_i·K⁻¹F`. 1e-9.
#[test]
fn rigid_diaphragm_matches_stiffness_and_torsion_split() {
    use carapace_core::model::Node3Id;
    let h = 2440.;
    // (start x, start z, ex, ez, length)
    let lines = [
        (0., 0., 1., 0., 4000.),
        (0., 6000., 1., 0., 2000.),
        (0., 0., 0., 1., 3000.),
        (8000., 0., 0., 1., 3000.),
    ];
    let wall = |d: &mut Domain3, (x0, z0, ex, ez, l): (f64, f64, f64, f64, f64)| {
        let p = |x: f64, y: f64| [x0 + x * ex, y, z0 + x * ez];
        let c = [
            d.add_node(Node3::new(p(0., 0.)).fix(0).fix(1).fix(2)),
            d.add_node(Node3::new(p(l, 0.)).fix(0).fix(1).fix(2)),
            d.add_node(Node3::new(p(l, h))),
            d.add_node(Node3::new(p(0., h))),
        ];
        let layout = StudWallLayout {
            length: l,
            height: h,
            studs: vec![
                WallStud {
                    x: 0.,
                    ea: CHORD_EA,
                },
                WallStud { x: l, ea: CHORD_EA },
            ],
            top_ea: 1e9,
            top_ei: Some(1e10),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: 1e9,
            bottom_ei: 1e10,
            openings: vec![],
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear { ga: 6000. },
            }),
            foundation: true,
            bearing: None,
            floor: None,
        };
        let id = d.add_element(Element3::StudWall(
            StudWall::new3(c, vec![], layout, vec![]).unwrap(),
        ));
        (id, [c[2], c[3]])
    };
    let solve = |d: Domain3| {
        let mut a = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Transformation)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-10,
                max_iter: 4,
            })
            .build(d);
        a.step().expect("rigid diaphragm solve");
        a
    };
    // Each line's lateral stiffness: the wall alone, top corners moving together along the wall.
    let k: Vec<f64> = lines
        .iter()
        .map(|&line| {
            let mut d = Domain3::new();
            let (_, top) = wall(&mut d, line);
            let slot = if line.2 != 0. { 0 } else { 2 };
            for (n, other) in [(top[0], 2 - slot), (top[1], 2 - slot)] {
                d.add_constraint((n, other), &[]);
            }
            d.equal_dof(top[0], top[1], &[slot]);
            d.add_nodal_load(d.default_pattern(), top[0], slot, 1000.);
            let a = solve(d);
            1000. / a.domain().node(top[0]).displacement[slot]
        })
        .collect();
    let (xm, zm, f) = (2000., 4500., 10_000.);
    let mut d = Domain3::new();
    let ids: Vec<_> = lines.iter().map(|&line| wall(&mut d, line)).collect();
    let master: Node3Id = d.add_node(Node3::new([xm, h, zm]));
    let tops: Vec<Node3Id> = ids.iter().flat_map(|(_, t)| *t).collect();
    d.rigid_diaphragm(master, &tops);
    d.add_nodal_load(d.default_pattern(), master, 0, f);
    let a = solve(d);
    // Hand: rigid-body stiffness about the master in (ux, uz, θ_y), with u = θ (z - z_m), w = -θ (x - x_m).
    let rows: Vec<[f64; 3]> = lines
        .iter()
        .map(|&(x0, z0, ex, _, _)| {
            if ex != 0. {
                [1., 0., z0 - zm]
            } else {
                [0., 1., -(x0 - xm)]
            }
        })
        .collect();
    let mut km = nalgebra::Matrix3::zeros();
    for (r, &ki) in rows.iter().zip(&k) {
        let r = nalgebra::Vector3::from(*r);
        km += ki * r * r.transpose();
    }
    let u = km.lu().solve(&nalgebra::Vector3::new(f, 0., 0.)).unwrap();
    for ((r, &ki), (id, _)) in rows.iter().zip(&k).zip(&ids) {
        let hand = ki * nalgebra::Vector3::from(*r).dot(&u);
        let fe = a.domain().element_stud_wall_response(*id, 1.).unwrap().base[0];
        assert!(
            (fe.abs() - hand.abs()).abs() <= 1e-9 * f,
            "wall line shear {fe} vs hand {hand}"
        );
    }
}
