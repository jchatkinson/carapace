//! Segmented behaviour of walls with openings: a wall with a window, king studs at the opening edges and
//! hold-downs at every pier end must act laterally as separate pier walls joined only by an axial drag strut
//! (segmented shear-wall design, SDPWS 2021 4.3 / CSA O86:19 11.7.1).
//!
//! The plates are hinged at the opening edges, so they tie the piers axially only; the header sits on the
//! jacks and transfers no lateral force. The reference is two separate `StudWall`s (one per pier) whose top
//! corners are joined by an axial `Truss` with the top plate's `EA` over the opening. Both sit on a foundation
//! (`foundation: true`) with the same hold-down springs at every pier end. The match is exact: 1e-10.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    AnchorageBearing, ChordDistribution, Domain, Element, ElementId, LoadSeries, Material, Node,
    NodeId, StudWall, StudWallLayout, Truss, WallChord, WallLoad, WallOpening, WallShearLaw,
    WallSheathing, WallStation, WallStud, ZeroLength,
};

type D2 = Domain<2, 3, NodeId, Element>;
const H: f64 = 2440.;
const PLY: f64 = 9500. * 38. * 140.;
const K_HD: f64 = 20_000.;
const V: f64 = 14_400.;
const X0: f64 = 1600.;
const X1: f64 = 2400.;
const TOL: f64 = 1e-10;

fn layout(
    length: f64,
    xs: &[f64],
    openings: Vec<WallOpening>,
    top: ChordDistribution,
) -> StudWallLayout {
    let studs = xs
        .iter()
        .map(|&x| WallStud {
            x,
            ea: if x == 0. || x == length || x == X0 || x == X1 {
                2. * PLY
            } else {
                PLY
            },
        })
        .collect();
    StudWallLayout {
        length,
        height: H,
        studs,
        top_ea: 2. * PLY,
        // A rim-size EI: with continuous plates this would couple the piers in bending.
        top_ei: (top == ChordDistribution::Beam).then_some(9500. * 38. * 235f64.powi(3) / 12.),
        top_distribution: top,
        bottom_ea: PLY,
        bottom_ei: 9500. * 140. * 38f64.powi(3) / 12.,
        openings,
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: true,
        bearing: None,
        floor: None,
    }
}

/// A bottom node on a hold-down spring (`uy`), fixed sideways.
fn holddown(d: &mut D2, ground: NodeId, x: f64, material: &Material) -> NodeId {
    let n = d.add_node(Node::new([x, 0.]).fix(0));
    d.add_element(Element::ZeroLength(
        ZeroLength::new(ground, n).with_material(1, material.clone()),
    ));
    n
}

fn nonlinear(m: &Material) -> bool {
    !matches!(m, Material::Elastic { .. })
}

struct Out {
    drift: f64,
    piers: Vec<(f64, f64, f64)>, // (q, chords[0], chords[1])
    uplift: Vec<f64>,            // hold-down node uy, by pier end
}

fn solve(
    mut d: D2,
    walls: &[ElementId],
    load_node: NodeId,
    holddowns: &[NodeId],
    nonlinear: bool,
) -> Out {
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(p, load_node, 0, V);
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(if nonlinear {
            Algorithm::Newton {
                tangent: carapace_core::analysis::TangentStrategy::Current,
                line_search: None,
            }
        } else {
            Algorithm::Linear
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 30,
        })
        .build(d);
    a.step().expect("linear solve");
    let dom = a.domain();
    Out {
        drift: dom.node(load_node).displacement[0],
        piers: walls
            .iter()
            .flat_map(|&w| dom.element_stud_wall_response(w, 1.).unwrap().piers)
            .map(|p| (p.q, p.chords[0], p.chords[1]))
            .collect(),
        uplift: holddowns
            .iter()
            .map(|&n| dom.node(n).displacement[1])
            .collect(),
    }
}

fn one_wall(top: ChordDistribution, anchorage: &Material) -> Out {
    let xs = [0., 400., 800., 1200., X0, X1, 2800., 3200., 3600.];
    let window = WallOpening {
        x0: X0,
        x1: X1,
        head_height: 2000.,
        header_ei: 9500. * 76. * 235f64.powi(3) / 12.,
        jack_ea: PLY,
        door: false,
    };
    let l = layout(3600., &xs, vec![window], top);
    assert_eq!(l.pier_ends(), vec![0., X0, X1, 3600.]);
    let mut d = D2::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let c0 = holddown(&mut d, ground, 0., anchorage);
    let c1 = holddown(&mut d, ground, 3600., anchorage);
    let c2 = d.add_node(Node::new([3600., H]));
    let c3 = d.add_node(Node::new([0., H]));
    let s0 = holddown(&mut d, ground, X0, anchorage);
    let s1 = holddown(&mut d, ground, X1, anchorage);
    let stations = [(X0, s0), (X1, s1)]
        .map(|(x, node)| WallStation {
            chord: WallChord::Bottom,
            x,
            node,
        })
        .to_vec();
    let empty: Vec<Vec<WallLoad>> = vec![vec![]; 8];
    let wall = StudWall::new([c0, c1, c2, c3], stations, l, empty).unwrap();
    let id = d.add_element(Element::StudWall(wall));
    solve(d, &[id], c3, &[c0, s0, s1, c1], nonlinear(anchorage))
}

fn separate_walls(top: ChordDistribution, anchorage: &Material) -> Out {
    let mut d = D2::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut ids = Vec::new();
    let mut tops = Vec::new();
    let mut hd = Vec::new();
    for (xa, xb, xs) in [
        (0., X0, vec![0., 400., 800., 1200., X0]),
        (X1, 3600., vec![X1, 2800., 3200., 3600.]),
    ] {
        let c0 = holddown(&mut d, ground, xa, anchorage);
        let c1 = holddown(&mut d, ground, xb, anchorage);
        let c2 = d.add_node(Node::new([xb, H]));
        let c3 = d.add_node(Node::new([xa, H]));
        let local: Vec<f64> = xs.iter().map(|x| x - xa).collect();
        let mut l = layout(xb - xa, &local, vec![], top);
        // End posts keep their 2-ply EA in the shifted frame.
        l.studs.first_mut().unwrap().ea = 2. * PLY;
        l.studs.last_mut().unwrap().ea = 2. * PLY;
        let empty: Vec<Vec<WallLoad>> = vec![vec![]; 8];
        let wall = StudWall::new([c0, c1, c2, c3], vec![], l, empty).unwrap();
        ids.push(d.add_element(Element::StudWall(wall)));
        tops.push((c3, c2));
        hd.extend([c0, c1]);
    }
    // Drag strut over the opening: the top plate's axial stiffness, 2 plies over the 800 mm span.
    d.add_element(Element::Truss(Truss::new(
        tops[0].1,
        tops[1].0,
        1.,
        Material::Elastic { e: 2. * PLY },
    )));
    solve(d, &ids, tops[0].0, &hd, nonlinear(anchorage))
}

fn near(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()),
        "{what}: {a} vs {b} (rel {:.2e})",
        (a - b).abs() / a.abs().max(b.abs())
    );
}

/// A window wall with king studs at the opening edges and hold-downs at all four pier ends equals two separate
/// pier walls joined by an axial drag strut, for both a `Direct` and a rim-stiff `Beam` top plate: drift, pier
/// shear flows, M/jd chord forces and hold-down uplift. Reference: segmented shear-wall idealization; 1e-10.
#[test]
fn window_wall_is_segmented() {
    for top in [ChordDistribution::Direct, ChordDistribution::Beam] {
        let elastic = Material::Elastic { e: K_HD };
        let one = one_wall(top, &elastic);
        let two = separate_walls(top, &elastic);
        near(one.drift, two.drift, TOL, "drift");
        assert_eq!(one.piers.len(), 2);
        for (a, b) in one.piers.iter().zip(&two.piers) {
            near(a.0, b.0, TOL, "q");
            near(a.1, b.1, TOL, "chord left");
            near(a.2, b.2, TOL, "chord right");
        }
        for (a, b) in one.uplift.iter().zip(&two.uplift) {
            near(*a, *b, TOL, "hold-down uy");
        }
    }
}

/// The nonlinear window wall with an O86 crushing/slack anchorage at all four pier ends equals two separately
/// solved pier walls joined by the same drag strut. Reference: segmented shear-wall idealization; 1e-10.
#[test]
fn nonlinear_anchorage_window_wall_is_segmented() {
    let anchor =
        Material::anchorage(K_HD, 0.15, AnchorageBearing::O86Crush { qr: 100_000.0 }).unwrap();
    let one = one_wall(ChordDistribution::Direct, &anchor);
    let two = separate_walls(ChordDistribution::Direct, &anchor);
    near(one.drift, two.drift, TOL, "drift");
    for (a, b) in one.piers.iter().zip(&two.piers) {
        near(a.0, b.0, TOL, "q");
        near(a.1, b.1, TOL, "chord left");
        near(a.2, b.2, TOL, "chord right");
    }
    for (a, b) in one.uplift.iter().zip(&two.uplift) {
        near(*a, *b, TOL, "anchorage displacement");
    }
}
