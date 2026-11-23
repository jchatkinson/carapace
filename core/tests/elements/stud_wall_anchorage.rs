//! A single-storey shear wall on `Anchorage` pier-end links against the CSA O86:19 11.7.1.2 deflection with the
//! anchorage term built the way WoodWorks Shearwalls builds `d_a`: hold-down elongation `T/k + slack` at the
//! tension end plus bottom-plate crushing `crush(C)` at the compression end, with `T` and `C` from statics
//! (overturning `V H / L` and the dead load on each end post). The wall is code-idealized: end posts only, plates
//! axially rigid stand-ins (`EA × 1e6`), a linear sheathing law, so the bending and shear terms are the code's
//! `2 v H³/(3 E A L)` and `v H / ga`. Tolerance 1e-5 relative (the stand-in plates, as in the multi-storey test).

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    AnchorageBearing, ChordDistribution, Domain, Element, ElementLoad, LoadSeries, Material, Node,
    NodeId, StudWall, StudWallLayout, WallCaseLoad, WallLoad, WallShearLaw, WallSheathing,
    WallStud, ZeroLength,
};

type D2 = Domain<2, 3, NodeId, Element>;
const L: f64 = 2400.;
const H: f64 = 2440.;
const POST: f64 = 2. * 9500. * 38. * 140.;
const GA: f64 = 6000.;
const K_T: f64 = 15_000.;
const QR: f64 = 60_000.;
const W_DEAD: f64 = 3.; // N/mm on the top plate: 3600 N on each end post

/// WoodWorks crushing curve: crush (mm) at compression `c` (N).
fn crush(c: f64) -> f64 {
    let r = c / QR;
    0.5 * if r <= 0.73 {
        r / 0.73
    } else if r <= 1. {
        1. + (r - 0.73) / 0.27
    } else {
        2. * r.powi(3)
    }
}

/// O86 11.7.1.2 with the anchorage term from statics. Positive `v_total` pushes the top towards +x, lifting
/// the left end.
fn hand_drift(v_total: f64, slack: f64) -> f64 {
    let ei = POST * L * L / 2.;
    let couple = v_total.abs() * H / L;
    let p = W_DEAD * L / 2.;
    let (up, down) = (couple - p, couple + p); // tension-end net uplift force, compression-end force
    let tension_end = if up > 0. {
        slack + up / K_T
    } else {
        -crush(-up)
    };
    let rotation = (tension_end + crush(down)) / L;
    let d = v_total.abs() * H.powi(3) / (3. * ei) + v_total.abs() / L * H / GA + H * rotation;
    d * v_total.signum()
}

struct Out {
    drift: f64,
    moment: f64,
}

fn solve(v_total: f64, slack: f64) -> Out {
    let mut d = D2::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let anchor = Material::anchorage(K_T, slack, AnchorageBearing::O86Crush { qr: QR }).unwrap();
    let base = [0., L].map(|x| {
        let n = d.add_node(Node::new([x, 0.]));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, n)
                .with_material(0, Material::Elastic { e: 1e10 })
                .with_material(1, anchor.clone()),
        ));
        n
    });
    let c2 = d.add_node(Node::new([L, H]));
    let c3 = d.add_node(Node::new([0., H]));
    let layout = StudWallLayout {
        length: L,
        height: H,
        studs: vec![WallStud { x: 0., ea: POST }, WallStud { x: L, ea: POST }],
        top_ea: POST * 1e6,
        top_ei: Some(1.2e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea: POST * 1e6,
        bottom_ei: 6e9,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: GA },
        }),
        foundation: false,
        floor: None,
    };
    let mut loads: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    let line = |wx: f64, wy: f64| WallLoad::TopLine {
        x0: 0.,
        x1: L,
        wx0: wx,
        wx1: wx,
        wy0: wy,
        wy1: wy,
    };
    loads[0].push(line(0., -W_DEAD));
    loads[1].push(line(v_total / L, 0.));
    let wall = StudWall::new([base[0], base[1], c2, c3], vec![], layout, loads).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    // Dead load and lateral shear are two patterns, ramped together.
    for case in 0..2 {
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        d.add_element_load(
            p,
            id,
            ElementLoad {
                wall: Some(WallCaseLoad { case, factor: 1. }),
                ..Default::default()
            },
        );
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        // Absolute displacement-increment norm (mm): 1e-7 is near roundoff with the 1e10 stand-ins.
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-7,
            max_iter: 500,
        })
        .build(d);
    for _ in 0..10 {
        a.step().expect("anchored wall step");
    }
    let dom = a.domain();
    let ux = |n: NodeId| dom.node(n).displacement[0];
    Out {
        drift: (ux(c2) + ux(c3) - ux(base[0]) - ux(base[1])) / 2.,
        moment: dom.element_stud_wall_response(id, 1.).unwrap().piers[0].moment_base,
    }
}

fn near(a: f64, b: f64, tol: f64, what: &str) {
    let rel = (a - b).abs() / a.abs().max(b.abs());
    assert!(rel <= tol, "{what}: {a} vs {b} (rel {rel:.2e})");
}

/// Uplift engaged (`V H/L` above the 3600 N dead load on the post), both load directions: drift equals the
/// O86 11.7.1.2 hand total with `d_a = slack + T/k + crush(C)` (1e-5; observed 2e-8), and the pier moment is `V H` (1e-6).
/// Also no uplift (the overturning couple below the dead load): the tension end stays in bearing and the
/// rotation is the difference of the two crushes.
#[test]
fn anchorage_drift_matches_o86_hand_calc() {
    for v in [12_000., -12_000., 2_000., -2_000.] {
        let out = solve(v, 0.5);
        near(out.drift, hand_drift(v, 0.5), 1e-5, &format!("drift V={v}"));
        near(
            out.moment.abs(),
            v.abs() * H,
            1e-6,
            &format!("moment V={v}"),
        );
    }
}

/// Shrinkage slack is a constant offset once the hold-down is engaged: drift(slack) − drift(0) = `H·slack/L`
/// exactly (1e-6 relative to the offset), which is the part a designer keeps out of the `RdRo/Ie`
/// amplification.
#[test]
fn slack_adds_constant_drift() {
    let s = 6.3;
    let diff = solve(12_000., s).drift - solve(12_000., 0.).drift;
    near(diff, H * s / L, 1e-6, "slack drift");
}
