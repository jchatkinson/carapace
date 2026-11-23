//! Floor zone of a stud wall (`StudWallLayout::floor`): the floor or roof depth between the framing top and
//! the exterior top row, condensed inside the wall (it replaces the former `FloorBand` element, whose bay and
//! rim math it reuses).
//!
//! Equivalence with the former wall + `FloorBand` + rigid vertical links model: the references below were
//! recorded from that model, run from the pre-change code (`git show 41e02b5`: tests `floor_band.rs`), with
//! ZeroLength links between each band bottom node and its top node. Stacks (the links are the same finite
//! springs in both models, `RIGID = 1e10`): the old model's own statics residual is 3e-8, so they agree to
//! about 1e-7 and are asserted at 1e-6. One storey and the window wall: the old model's rigid links were
//! penalties, re-recorded at `1e12` N/mm (the old results drift by about 1e-9 relative between `1e10` and
//! `1e13` from conditioning: the reference's own noise); the zone's posts are exact; 1e-8 relative (observed
//! 5e-9).
//!
//! Hand checks: one storey on a zone, base moment `V(H + d)` (O86:19 Fig. 11.7 statics), 1e-9; zone shear flow
//! integrates to the load, 1e-9. Window wall: pier shears sum to the load (1e-8) and each pier's moment is within
//! 1% of `V_p(H + d)` (the rim at `H + d` and the top plate at `H` both carry drag over the opening and couple the
//! piers slightly, about 0.5% here). 2- and 5-storey walls with rigid 254 mm zones against an exact hand walk
//! with true lever arms and the O86 A.11.7.1 storey terms: level displacements 1e-5, chords and base moments
//! 1e-6. Rigid-body modes give zero force (1e-12), bays over an opening carry no shear, layout errors are named.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    ChordDistribution, Domain, Element, FloorZone, LoadSeries, Material, Node, NodeId, StudWall,
    StudWallLayout, WallCaseLoad, WallChord, WallLoad, WallOpening, WallShearLaw, WallSheathing,
    WallStation, WallStud, ZeroLength,
};

type D2 = Domain<2, 3, NodeId, Element>;

fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b} (relative {:e}, tol {tol:e})",
        (a - b).abs() / a.abs().max(b.abs()).max(1e-20)
    );
}
fn near_all(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        near(*x, *y, tol);
    }
}

const PLY_EA: f64 = 9500. * 38. * 140.;
const BV: f64 = 11_000.;
const D_F: f64 = 3.25;
const RIGID: f64 = 1e10;
const FLOOR_DEPTH: f64 = 254.;

struct StackCase {
    name: &'static str,
    length: f64,
    height: f64,
    forces: Vec<f64>,
    plies: Vec<f64>,
    nail_spacing: Vec<f64>,
    k_hd: f64,
}

fn two_storey_stack() -> StackCase {
    StackCase {
        name: "2-storey",
        length: 2400.,
        height: 3000.,
        forces: vec![6000., 8000.],
        plies: vec![3., 2.],
        nail_spacing: vec![100., 150.],
        k_hd: 20_000.,
    }
}

fn five_storey_stack() -> StackCase {
    StackCase {
        name: "5-storey",
        length: 4800.,
        height: 3000.,
        forces: vec![3000., 6000., 9000., 12000., 12000.],
        plies: vec![6., 5., 4., 3., 2.],
        nail_spacing: vec![75., 100., 100., 150., 150.],
        k_hd: 50_000.,
    }
}

impl StackCase {
    fn n(&self) -> usize {
        self.forces.len()
    }

    fn shear(&self, i: usize) -> f64 {
        self.forces[i..].iter().sum()
    }

    fn force_y(&self, i: usize, depth: f64) -> f64 {
        if i + 1 < self.n() {
            i as f64 * (self.height + depth) + self.height + depth
        } else {
            i as f64 * (self.height + depth) + self.height
        }
    }

    fn wall_base_y(&self, i: usize, depth: f64) -> f64 {
        i as f64 * (self.height + depth)
    }

    fn moment_base(&self, i: usize, depth: f64) -> f64 {
        (i..self.n())
            .map(|j| self.forces[j] * (self.force_y(j, depth) - self.wall_base_y(i, depth)))
            .sum()
    }

    fn moment_top(&self, i: usize, depth: f64) -> f64 {
        let top_y = self.wall_base_y(i, depth) + self.height;
        (i..self.n())
            .map(|j| self.forces[j] * (self.force_y(j, depth) - top_y).max(0.))
            .sum()
    }

    fn nail_slip(&self, i: usize) -> f64 {
        let v = self.shear(i) / self.length;
        (0.013 * v * self.nail_spacing[i] / D_F.powi(2)).powi(2)
    }

    fn secant_ga(&self, i: usize) -> f64 {
        let v = self.shear(i) / self.length;
        1. / (1. / BV + 0.0025 * self.nail_slip(i) / v)
    }
}

struct StackResult {
    level_u: Vec<f64>,
    chord: Vec<f64>,
    base_moment: Vec<f64>,
}

/// Hand walk for the rigid-band model. At each wall base, add the hold-down rotation from the true static
/// base moment (all force heights include the intervening `d`). Across a wall, add `phi H`, the local O86
/// bending, shear and nail slip, then the wall rotation `theta`; across a rigid band, add `phi d`. The
/// one-flexible-storey O86 comparison instead uses `H_i = H+d` for every wall and is printed only: it assigns
/// floor depth to flexible wall deformation, while this model treats the band as rigid.
fn stack_hand_calc(c: &StackCase, depth: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut u, mut phi) = (0., 0.);
    let mut levels = Vec::new();
    let mut base_moments = Vec::new();
    let mut chord = Vec::new();
    for i in 0..c.n() {
        let h = c.height;
        let ei = c.plies[i] * PLY_EA * c.length.powi(2) / 2.;
        let v = c.shear(i);
        let m_top = c.moment_top(i, depth);
        let m_base = c.moment_base(i, depth);
        let theta = m_top * h / ei + v * h.powi(2) / (2. * ei);
        let bend = m_top * h.powi(2) / (2. * ei) + v * h.powi(3) / (3. * ei);
        let shear = v / c.length * h / BV;
        let slip = 0.0025 * h * c.nail_slip(i);
        phi += (m_base / c.length) / c.k_hd / c.length;
        u += phi * h + bend + shear + slip;
        phi += theta;
        base_moments.push(m_base);
        chord.push(m_base / c.length);
        if i + 1 < c.n() {
            u += phi * depth;
        }
        levels.push(u);
    }
    (levels, chord, base_moments)
}

fn stack_o86_comparison(c: &StackCase, storey_height: f64) -> Vec<f64> {
    let mut rotation = 0.;
    let mut anchorage_rotation = 0.;
    let mut total = 0.;
    let mut out = Vec::new();
    for i in 0..c.n() {
        let h = storey_height;
        let ei = c.plies[i] * PLY_EA * c.length.powi(2) / 2.;
        let v = c.shear(i);
        let m = c.moment_top(i, 0.);
        let theta = m * h / ei + v * h.powi(2) / (2. * ei);
        let bending = m * h.powi(2) / (2. * ei) + v * h.powi(3) / (3. * ei) + h * rotation;
        let shear = v / c.length * h / BV;
        let slip = 0.0025 * h * c.nail_slip(i);
        anchorage_rotation += c.moment_base(i, 0.) / c.length / c.k_hd / c.length;
        total += bending + shear + slip + h * anchorage_rotation;
        rotation += theta;
        out.push(total);
    }
    out
}

fn solve(d: D2, tol: f64) -> carapace_core::analysis::Analysis {
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance { tol, max_iter: 4 })
        .build(d);
    a.step().expect("linear floor-zone solve");
    a
}
fn empty_loads() -> Vec<Vec<WallLoad>> {
    vec![Vec::new(); 8]
}

// --- One storey on a zone -------------------------------------------------------------------------------

const ONE_LENGTH: f64 = 3600.;
const ONE_HEIGHT: f64 = 2440.;
const ONE_DEPTH: f64 = 254.;
const ONE_FORCE: f64 = 14_400.;

fn one_storey_layout(sheathing: WallShearLaw) -> StudWallLayout {
    StudWallLayout {
        length: ONE_LENGTH,
        height: ONE_HEIGHT,
        studs: vec![
            WallStud { x: 0., ea: 1e10 },
            WallStud {
                x: ONE_LENGTH,
                ea: 1e10,
            },
        ],
        top_ea: 1e10,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: 1e10,
        bottom_ei: 1e10,
        openings: vec![],
        sheathing: Some(WallSheathing { law: sheathing }),
        foundation: false,
        floor: Some(FloorZone {
            depth: ONE_DEPTH,
            ga: 1e10,
            rim_ea: 1e10,
        }),
    }
}

/// Wall on two ground springs (rigid in `ux`, 20 kN/mm in `uy`), zone on top, `ONE_FORCE` at the exterior top
/// left. Returns the solved analysis, the wall and the exterior top nodes.
fn one_storey(
    law: WallShearLaw,
    newton: bool,
) -> (
    carapace_core::analysis::Analysis,
    carapace_core::model::ElementId,
    [NodeId; 2],
) {
    let mut d: D2 = D2::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut lower = Vec::new();
    for x in [0., ONE_LENGTH] {
        let node = d.add_node(Node::new([x, 0.]));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, node)
                .with_material(0, Material::Elastic { e: 1e10 })
                .with_material(1, Material::Elastic { e: 20_000. }),
        ));
        lower.push(node);
    }
    let top = ONE_HEIGHT + ONE_DEPTH;
    let c2 = d.add_node(Node::new([ONE_LENGTH, top]));
    let c3 = d.add_node(Node::new([0., top]));
    let wall = d.add_element(Element::StudWall(
        StudWall::new(
            [lower[0], lower[1], c2, c3],
            vec![],
            one_storey_layout(law),
            empty_loads(),
        )
        .unwrap(),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(pattern, c3, 0, ONE_FORCE);
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(if newton {
            Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            }
        } else {
            Algorithm::Linear
        })
        .test(if newton {
            ConvergenceTest::NormDispIncr {
                tol: 1e-10,
                max_iter: 100,
            }
        } else {
            ConvergenceTest::NormUnbalance {
                tol: 1e-6,
                max_iter: 4,
            }
        })
        .build(d);
    a.step().expect("one-storey floor-zone solve");
    (a, wall, [c3, c2])
}

/// A lateral load on the exterior top of a wall with a floor zone gives the wall below the floor-depth lever
/// `V(H + d)`, and the zone's shear flow integrates to the load. Reference: O86 Fig. 11.7 statics; moment and
/// reaction equilibrium tolerance 1e-9. The recorded numbers of the former wall + `FloorBand` + rigid links
/// model (ground springs and `ga = rim = 1e10` as here, links `1e12`) agree to 1e-8. They were produced by
/// printing `u(c3)`, `u(c2)`, pier `q`, `moment_base`, `moment_top`, `chords[0]` and `base` from the old test.
#[test]
fn floor_depth_adds_to_wall_overturning_lever() {
    let (a, wall, [c3, c2]) = one_storey(WallShearLaw::Linear { ga: 1e5 }, false);
    let dom = a.domain();
    let r = dom.element_stud_wall_response(wall, 1.).unwrap();
    near(
        r.piers[0].moment_base.abs(),
        ONE_FORCE * (ONE_HEIGHT + ONE_DEPTH),
        1e-9,
    );
    near(r.base[0].abs(), ONE_FORCE, 1e-9);
    // Moment at the framing top: the force still acts a depth `d` above it.
    near(
        r.piers[0].moment_top.abs(),
        ONE_FORCE * ONE_DEPTH,
        1e-9 * (ONE_HEIGHT / ONE_DEPTH),
    );
    assert_eq!(r.floor_bays.len(), 1);
    near(r.floor_bays[0].q * ONE_LENGTH, ONE_FORCE, 1e-9);
    assert_eq!(r.floor_rim.len(), 1);
    // Recorded references of the former band model.
    let (tol, node) = (1e-8, |n: NodeId, i: usize| dom.node(n).displacement[i]);
    near(node(c3, 0), 9.076119085883412e-1, tol);
    near(node(c3, 1), 5.402386235829583e-1, tol);
    near(node(c2, 0), 9.050199085893841e-1, tol);
    near(node(c2, 1), -5.402386234872552e-1, tol);
    let p = &r.piers[0];
    near(p.q, 3.999999998220368e0, tol);
    near(p.moment_base, -3.879359996665017e7, tol);
    near(p.moment_top, -3.657599982282464e6, tol);
    near(p.chords[0], 1.077599999073616e4, tol);
    near(r.base[0], -1.4399999993593334e4, tol);
    near(r.base[2], 3.879359996665016e7, tol);
}

/// The nonlinear `O86` and `Material` sheathing laws keep the zone: the zone is linear and condensed with the
/// other internals. Newton solution; statics hold exactly (pier `q = V/L`, base moment `V(H + d)`, zone flow
/// `V/L`): 1e-8.
#[test]
fn nonlinear_sheathing_keeps_the_floor_zone() {
    for law in [
        WallShearLaw::O86 {
            bv: 11_000.,
            nail_spacing: 100.,
            nail_diameter: 3.25,
        },
        WallShearLaw::Material(Material::Elastic { e: 1e5 }),
    ] {
        let (a, wall, _) = one_storey(law, true);
        let r = a.domain().element_stud_wall_response(wall, 1.).unwrap();
        near(r.piers[0].q * ONE_LENGTH, ONE_FORCE, 1e-8);
        near(
            r.piers[0].moment_base.abs(),
            ONE_FORCE * (ONE_HEIGHT + ONE_DEPTH),
            1e-8,
        );
        near(r.floor_bays[0].q * ONE_LENGTH, ONE_FORCE, 1e-8);
    }
}

/// Top slot loads act on the exterior top row when the wall has a zone: a uniform line load `w` over the wall
/// gives the base moment `w L (H + d)` and equals nodal loads `w a/2` at the top positions (lever rule),
/// while a vertical line load still goes down the posts. Reference: statics, 1e-9; the two load paths agree to
/// 1e-10.
#[test]
fn top_slot_loads_act_on_the_exterior_row() {
    let build = |slot: bool| {
        let mut d: D2 = D2::new();
        let top = ONE_HEIGHT + ONE_DEPTH;
        let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
        let mid = ONE_LENGTH / 2.;
        let mut bottom = Vec::new();
        for x in [0., ONE_LENGTH] {
            let node = d.add_node(Node::new([x, 0.]));
            d.add_element(Element::ZeroLength(
                ZeroLength::new(ground, node)
                    .with_material(0, Material::Elastic { e: 1e10 })
                    .with_material(1, Material::Elastic { e: 2e4 }),
            ));
            bottom.push(node);
        }
        let c2 = d.add_node(Node::new([ONE_LENGTH, top]));
        let c3 = d.add_node(Node::new([0., top]));
        let m = d.add_node(Node::new([mid, top]));
        let mut layout = one_storey_layout(WallShearLaw::Linear { ga: 1e5 });
        layout.studs.insert(1, WallStud { x: mid, ea: 1e8 });
        layout.top_distribution = ChordDistribution::Direct;
        let w = 4.;
        let mut loads = empty_loads();
        if slot {
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: ONE_LENGTH,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
        }
        let wall = d.add_element(Element::StudWall(
            StudWall::new(
                [bottom[0], bottom[1], c2, c3],
                vec![WallStation {
                    chord: WallChord::Top,
                    x: mid,
                    node: m,
                }],
                layout,
                loads,
            )
            .unwrap(),
        ));
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        if slot {
            d.add_element_load(
                p,
                wall,
                carapace_core::model::ElementLoad {
                    wall: Some(WallCaseLoad {
                        case: 0,
                        factor: 1.,
                    }),
                    ..Default::default()
                },
            );
        } else {
            for (n, f) in [(c3, w * mid / 2.), (m, w * mid), (c2, w * mid / 2.)] {
                d.add_nodal_load(p, n, 0, f);
            }
        }
        let a = solve(d, 1e-6);
        let r = a.domain().element_stud_wall_response(wall, 1.).unwrap();
        let u = [c3, m, c2].map(|n| a.domain().node(n).displacement);
        (r, u)
    };
    let (slot, us) = build(true);
    let (nodal, un) = build(false);
    near(
        slot.piers[0].moment_base.abs(),
        4. * ONE_LENGTH * (ONE_HEIGHT + ONE_DEPTH),
        1e-9,
    );
    for (a, b) in us.iter().zip(&un) {
        assert!((a[0] - b[0]).abs() <= 1e-10 * a[0].abs());
        assert!((a[1] - b[1]).abs() <= 1e-10 * a[0].abs());
    }
    near(slot.piers[0].moment_base, nodal.piers[0].moment_base, 1e-9);
}

// --- Window wall -----------------------------------------------------------------------------------------

/// Window wall: king studs at both opening edges, four pier-end hold-downs, top stations at every stud
/// position and a real rim (one 38x235 ply), solved under a diaphragm load spread along the exterior top with
/// an optional rigid-body displacement field `[tx, ty, phi]` imposed on every node instead (no load).
struct Window {
    analysis: carapace_core::analysis::Analysis,
    wall: carapace_core::model::ElementId,
    top: Vec<NodeId>,
}
const WIN_XS: [f64; 10] = [
    0., 400., 800., 1200., 1600., 2000., 2400., 2800., 3200., 3600.,
];
const WIN_FORCE: f64 = 14_400.;
fn window(rigid: Option<[f64; 3]>) -> Window {
    let (length, depth) = (3600., FLOOR_DEPTH);
    let wall_layout = StudWallLayout {
        length,
        height: 2440.,
        studs: WIN_XS
            .iter()
            .map(|&x| WallStud {
                x,
                ea: if x == 0. || x == length || x == 1600. || x == 2400. {
                    2. * PLY_EA
                } else {
                    PLY_EA
                },
            })
            .collect(),
        top_ea: 2. * PLY_EA,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: PLY_EA,
        bottom_ei: PLY_EA,
        openings: vec![WallOpening {
            x0: 1600.,
            x1: 2400.,
            head_height: 2000.,
            header_ei: 1e9,
            jack_ea: PLY_EA,
            door: false,
        }],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: rigid.is_none(),
        floor: Some(FloorZone {
            depth,
            ga: RIGID,
            rim_ea: 9500. * 38. * 235.,
        }),
    };
    let mut d: D2 = Domain::new();
    let node = |x: f64, y: f64| {
        let mut n = Node::new([x, y]);
        if let Some([tx, ty, phi]) = rigid {
            n.displacement = [tx - phi * y, ty + phi * x, 0.];
        }
        n
    };
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut bottom = Vec::new();
    for &x in &[0., 1600., 2400., length] {
        let n = d.add_node(node(x, 0.).fix(0));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, n).with_material(1, Material::Elastic { e: 20_000. }),
        ));
        bottom.push((x, n));
    }
    let ytop = 2440. + depth;
    let top: Vec<NodeId> = WIN_XS.iter().map(|&x| d.add_node(node(x, ytop))).collect();
    let stations: Vec<WallStation> = WIN_XS
        .iter()
        .zip(&top)
        .filter(|(x, _)| **x != 0. && **x != length)
        .map(|(&x, &n)| WallStation {
            chord: WallChord::Top,
            x,
            node: n,
        })
        .chain(bottom.iter().skip(1).take(2).map(|&(x, n)| WallStation {
            chord: WallChord::Bottom,
            x,
            node: n,
        }))
        .collect();
    let wall = d.add_element(Element::StudWall(
        StudWall::new(
            [bottom[0].1, bottom[3].1, top[9], top[0]],
            stations,
            wall_layout,
            empty_loads(),
        )
        .unwrap(),
    ));
    if rigid.is_none() {
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        for (i, &n) in top.iter().enumerate() {
            let lo = if i == 0 {
                0.
            } else {
                (WIN_XS[i] - WIN_XS[i - 1]) / 2.
            };
            let hi = if i == 9 {
                0.
            } else {
                (WIN_XS[i + 1] - WIN_XS[i]) / 2.
            };
            d.add_nodal_load(p, n, 0, WIN_FORCE * (lo + hi) / length);
        }
    }
    let analysis = if rigid.is_none() {
        solve(d, 1e-6)
    } else {
        AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Transformation)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-6,
                max_iter: 4,
            })
            .build(d)
    };
    Window {
        analysis,
        wall,
        top,
    }
}

/// A window wall with king studs at both opening edges, four pier-end hold-downs, top stations at every stud
/// position and inactive zone bays over the opening, under a diaphragm load spread uniformly along the
/// exterior top with a real rim: the pier shears sum to the load (1e-8, roundoff of the `1e10` rigid stand-ins) and each pier's base moment is within 1%
/// of its shear times `H + d` (O86 Fig. 11.7). It is not exact: over the opening the drag runs both in the rim
/// (at `H + d`) and in the wall's top plate (at `H`), and two axial struts a depth `d` apart couple the piers
/// slightly (about 0.5% here). The recorded results of the former band model (links `1e12`) agree to 1e-8:
/// pier `q`, moments and chords, exterior top `ux`/`uy`, base, bay shear flows and rim forces.
#[test]
fn window_wall_zone_preserves_pier_levers() {
    let w = window(None);
    let dom = w.analysis.domain();
    let r = dom.element_stud_wall_response(w.wall, 1.).unwrap();
    assert_eq!(r.piers.len(), 2);
    let shear_total: f64 = r.piers.iter().map(|p| p.q * (p.xb - p.xa)).sum();
    near(shear_total, WIN_FORCE, 1e-8);
    for p in &r.piers {
        near(
            p.moment_base.abs(),
            p.q.abs() * (p.xb - p.xa) * (2440. + FLOOR_DEPTH),
            1e-2,
        );
    }
    assert_eq!(r.floor_bays.len(), 9);
    assert_eq!((r.floor_bays[4].q, r.floor_bays[5].q), (0., 0.));
    let tol = 1e-8;
    let col = |f: fn(&carapace_core::model::WallPierResponse) -> f64| {
        r.piers.iter().map(f).collect::<Vec<_>>()
    };
    near_all(&col(|p| p.q), &[5.614034527203792, 4.514620621105979], tol);
    near_all(
        &col(|p| p.moment_base),
        &[-24309306.346879236, -14484293.638193673],
        tol,
    );
    near_all(
        &col(|p| p.moment_top),
        &[-2392115.552675632, -1265484.459595367],
        // A difference of two moments ten times larger: roundoff of the reference shows at 3e-8.
        1e-7,
    );
    near_all(
        &col(|p| p.chords[0]),
        &[15193.316466799522, 12070.244698494727],
        tol,
    );
    near(r.base[0], -14399.999988853247, 1e-7);
    near(r.base[2], 38793599.96576235, tol);
    let ux = [
        4.572363346962544,
        4.56085000702879,
        4.552061879390479,
        4.572331851490887,
        4.596991074722253,
        4.607447973535971,
        4.61036081624077,
        4.604122145778995,
        4.603656113533042,
        4.60899915907085,
    ];
    let uy = [
        0.6956207216849616,
        0.33487484603424633,
        0.007555200294129211,
        -0.3283462645661408,
        -0.7172475805707766,
        -0.03305644815660516,
        0.651597112693942,
        0.21497526454361843,
        -0.2145014080077428,
        -0.6523247736185178,
    ];
    near_all(
        &w.top
            .iter()
            .map(|&n| dom.node(n).displacement[0])
            .collect::<Vec<_>>(),
        &ux,
        tol,
    );
    for (n, u) in w.top.iter().zip(uy) {
        // The vertical values pass through zero: compare on the scale of the largest.
        assert!((dom.node(*n).displacement[1] - u).abs() <= tol * 0.72);
    }
    let q: Vec<f64> = r.floor_bays.iter().map(|b| b.q).collect();
    let q_ref = [
        -8.209177391679813,
        19.099094789287165,
        19.71520392904033,
        -7.06067743551711,
        0.0,
        0.0,
        -1.704607736430566,
        15.82612983480787,
        -1.665965828831495,
    ];
    for (a, b) in q.iter().zip(q_ref) {
        assert!((a - b).abs() <= tol * 19.8, "bay q {a} != {b}");
    }
    let rim_ref = [
        -2441.835483200034,
        -1863.852020490389,
        4299.0077078452905,
        5229.91300708244,
        2217.7775271542578,
        617.7775271541689,
        -1323.1440215616483,
        -98.83961396371859,
        1133.1931704998683,
    ];
    for (a, b) in r.floor_rim.iter().zip(rim_ref) {
        assert!((a - b).abs() <= tol * 5230., "rim {a} != {b}");
    }
}

/// Rigid translation and infinitesimal rotation of every node of a window wall with a zone give zero element
/// force, bay shear strain and rim force (rigid-body invariance of the wall including the zone); 1e-12 of the
/// force scale `GA a`.
#[test]
fn rigid_body_modes_have_zero_force() {
    let phi = 1e-3;
    let w = window(Some([3., 4., phi]));
    let dom = w.analysis.domain();
    let scale = 6000. * 3600. * phi * 3000.;
    for f in dom.element_end_force(w.wall, 1.).as_slice() {
        assert!(f.abs() <= 1e-9 * scale, "rigid-body force residual {f}");
    }
    let r = dom.element_stud_wall_response(w.wall, 1.).unwrap();
    for b in &r.floor_bays {
        assert!(b.gamma.abs() <= 1e-12, "bay strain {}", b.gamma);
    }
    for f in &r.floor_rim {
        assert!(f.abs() <= 1e-9 * scale, "rim force {f}");
    }
}

/// Invalid zones and top stations are named errors, never silently dropped. Reference: the `StudWallError`
/// contract; exact message fragments.
#[test]
fn validates_the_zone_and_its_stations() {
    let corners = {
        let mut d: D2 = D2::new();
        [(0., 0.), (3600., 0.), (3600., 2694.), (0., 2694.)]
            .map(|(x, y)| d.add_node(Node::new([x, y])))
    };
    let build = |zone: FloorZone, stations: Vec<WallStation>| {
        let mut l = one_storey_layout(WallShearLaw::Linear { ga: 1e5 });
        l.floor = Some(zone);
        l.studs.insert(1, WallStud { x: 1800., ea: 1e8 });
        StudWall::new(corners, stations, l, empty_loads())
            .err()
            .map(|e| e.to_string())
    };
    let ok = FloorZone {
        depth: 254.,
        ga: 1e5,
        rim_ea: 0.,
    };
    assert_eq!(build(ok, vec![]), None);
    for (zone, text) in [
        (FloorZone { depth: 0., ..ok }, "floor.depth"),
        (
            FloorZone {
                depth: f64::NAN,
                ..ok
            },
            "floor.depth",
        ),
        (FloorZone { ga: 0., ..ok }, "floor.ga"),
        (
            FloorZone {
                ga: f64::INFINITY,
                ..ok
            },
            "floor.ga",
        ),
        (FloorZone { rim_ea: -1., ..ok }, "floor.rim_ea"),
    ] {
        let e = build(zone, vec![]).expect("invalid zone must be rejected");
        assert!(e.contains(text), "{e}");
    }
    // A top station off every stud, jack and cripple position cannot take a rigid post.
    let off = {
        let mut d: D2 = D2::new();
        d.add_node(Node::new([1000., 2694.]))
    };
    let e = build(
        ok,
        vec![WallStation {
            chord: WallChord::Top,
            x: 1000.,
            node: off,
        }],
    )
    .expect("off-stud station must be rejected");
    assert!(e.contains("floor-zone top station"), "{e}");
}

// --- Stacks ----------------------------------------------------------------------------------------------

fn zone_stack(c: &StackCase) -> StackResult {
    let mut d: D2 = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let left0 = d.add_node(Node::new([0., 0.]));
    let right0 = d.add_node(Node::new([c.length, 0.]).fix(0).fix(1));
    d.add_element(Element::ZeroLength(
        ZeroLength::new(ground, left0)
            .with_material(0, Material::Elastic { e: RIGID })
            .with_material(1, Material::Elastic { e: c.k_hd }),
    ));
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    let mut bottom = [left0, right0];
    let mut walls = Vec::new();
    let mut force_nodes: Vec<[NodeId; 2]> = Vec::new();
    for i in 0..c.n() {
        let y = c.wall_base_y(i, FLOOR_DEPTH);
        let last = i + 1 == c.n();
        let zone = (!last).then_some(FloorZone {
            depth: FLOOR_DEPTH,
            ga: RIGID,
            rim_ea: RIGID,
        });
        let ytop = y + c.height + zone.map_or(0., |z| z.depth);
        let top_left = d.add_node(Node::new([0., ytop]));
        let top_right = d.add_node(Node::new([c.length, ytop]));
        let post_ea = c.plies[i] * PLY_EA;
        let layout = StudWallLayout {
            length: c.length,
            height: c.height,
            studs: vec![
                WallStud { x: 0., ea: post_ea },
                WallStud {
                    x: c.length,
                    ea: post_ea,
                },
            ],
            top_ea: PLY_EA * 1e6,
            top_ei: Some(1.2e10),
            top_distribution: ChordDistribution::Beam,
            bottom_ea: PLY_EA * 1e6,
            bottom_ei: 6.0e9,
            openings: vec![],
            sheathing: Some(WallSheathing {
                law: WallShearLaw::Linear { ga: c.secant_ga(i) },
            }),
            foundation: false,
            floor: zone,
        };
        let mut loads = empty_loads();
        if last {
            let w = c.forces[i] / c.length;
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: c.length,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
        }
        let wall = StudWall::new(
            [bottom[0], bottom[1], top_right, top_left],
            vec![],
            layout,
            loads,
        )
        .unwrap();
        let wall_id = d.add_element(Element::StudWall(wall));
        if last {
            d.add_element_load(
                pattern,
                wall_id,
                carapace_core::model::ElementLoad {
                    wall: Some(WallCaseLoad {
                        case: 0,
                        factor: 1.,
                    }),
                    ..Default::default()
                },
            );
            force_nodes.push([top_left, top_right]);
        } else {
            // The wall above starts on new nodes at the exterior top, tied horizontally and joined by
            // vertical links (the hold-down at the left, a rigid bearing at the right).
            let bl = d.add_node(Node::new([0., ytop]));
            let br = d.add_node(Node::new([c.length, ytop]));
            for (a, b, stiffness) in [(top_left, bl, c.k_hd), (top_right, br, RIGID)] {
                d.add_element(Element::ZeroLength(
                    ZeroLength::new(a, b).with_material(1, Material::Elastic { e: stiffness }),
                ));
                d.equal_dof(a, b, &[0]);
            }
            d.add_nodal_load(pattern, top_left, 0, c.forces[i]);
            force_nodes.push([top_left, top_left]);
            bottom = [bl, br];
        }
        walls.push(wall_id);
    }
    let a = solve(d, 1e-13);
    let dom = a.domain();
    StackResult {
        level_u: force_nodes
            .iter()
            .map(|n| (dom.node(n[0]).displacement[0] + dom.node(n[1]).displacement[0]) / 2.)
            .collect(),
        chord: walls
            .iter()
            .map(|&id| dom.element_stud_wall_response(id, 1.).unwrap().piers[0].chords[0])
            .collect(),
        base_moment: walls
            .iter()
            .map(|&id| dom.element_stud_wall_response(id, 1.).unwrap().piers[0].moment_base)
            .collect(),
    }
}

fn compare_zone_stack(c: &StackCase) -> Vec<String> {
    let actual = zone_stack(c);
    let (hand, hand_chord, hand_moment) = stack_hand_calc(c, FLOOR_DEPTH);
    let no_zone = stack_hand_calc(c, 0.).0;
    let o86 = stack_o86_comparison(c, c.height + FLOOR_DEPTH);
    eprintln!("\n{} with rigid floor zones (d = {FLOOR_DEPTH} mm)", c.name);
    eprintln!(" level | FE zone (mm) | hand rigid zone | O86 H+d | no zone");
    let mut failures = Vec::new();
    for i in 0..c.n() {
        eprintln!(
            " {:>5} | {:>12.6} | {:>15.6} | {:>8.6} | {:>8.6}",
            i + 1,
            actual.level_u[i],
            hand[i],
            o86[i],
            no_zone[i]
        );
    }
    for i in 0..c.n() {
        for (label, got, expected, tolerance) in [
            ("level displacement", actual.level_u[i], hand[i], 1e-5),
            (
                "M/jd chord force",
                actual.chord[i].abs(),
                hand_chord[i],
                1e-6,
            ),
            (
                "base moment",
                actual.base_moment[i].abs(),
                hand_moment[i],
                1e-6,
            ),
        ] {
            if (got - expected).abs() > tolerance * got.abs().max(expected.abs()) {
                failures.push(format!(
                    "{} level {} {label}: FE {got:.9}, hand {expected:.9} (tol {tolerance:e})",
                    c.name,
                    i + 1
                ));
            }
        }
    }
    failures
}

/// Two- and five-storey code-idealized walls with rigid `254 mm` floor zones match the exact hand walk for
/// level deflections and static M/jd chord forces. The printed O86 `H+d` and no-zone values show the
/// floor-depth effect; O86 treats that depth as flexible wall height, while the zone keeps it rigid.
/// Tolerances 1e-5 for displacement and 1e-6 for chord forces and moments.
#[test]
fn rigid_floor_zones_match_multistorey_hand_calc() {
    let mut failures = compare_zone_stack(&two_storey_stack());
    failures.extend(compare_zone_stack(&five_storey_stack()));
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The same stacks equal the former wall + `FloorBand` + links model, whose results were recorded from the
/// pre-change code (`git show 41e02b5`, test `floor_band.rs`, `{:?}` of the stack results; the links there are
/// the same finite springs as here, so the networks are identical): level displacements, chord forces and base
/// moments, 1e-6 relative. They agree to about 1e-7: the recorded results carry the roundoff of that model's
/// conditioning (its own statics residual is 2.5e-8 on the 2-storey chord and 3e-8 on the 5-storey base
/// moment; the zone model's is 3e-9), so 1e-9 is not attainable against them.
#[test]
fn stacks_equal_the_former_band_model() {
    let two = zone_stack(&two_storey_stack());
    near_all(&two.level_u, &[8.0787600444487, 15.51254671690047], 1e-6);
    near_all(&two.chord, &[28981.666389495196, 9999.99975409579], 1e-6);
    near_all(
        &two.base_moment,
        &[-69555999.33478847, -23999999.409829896],
        1e-6,
    );
    let five = zone_stack(&five_storey_stack());
    near_all(
        &five.level_u,
        &[
            9.312008671719154,
            23.11562924041682,
            35.91218391339955,
            50.106979642611485,
            58.66780912436458,
        ],
        1e-6,
    );
    near_all(
        &five.chord,
        &[
            101052.49420739675,
            72579.99570664734,
            46141.24716821181,
            23769.9986803108,
            7499.999782750062,
        ],
        1e-6,
    );
    near_all(
        &five.base_moment,
        &[
            -485051972.1955044,
            -348383979.3919072,
            -221477986.4074167,
            -114095993.66549185,
            -35999998.957200296,
        ],
        1e-6,
    );
}

/// Zone assembly gives bit-identical serial and threaded displacements and reactions. Reference: the native
/// parallel scatter contract; exact `f64::to_bits` equality.
#[cfg(feature = "parallel")]
#[test]
fn serial_and_parallel_zone_results_are_bit_identical() {
    use carapace_core::analysis::{set_parallelism, Parallelism};
    set_parallelism(Parallelism::Serial);
    let serial = window(None);
    set_parallelism(Parallelism::Threads(2));
    let threaded = window(None);
    set_parallelism(Parallelism::Serial);
    for (a, b) in serial.top.iter().zip(&threaded.top) {
        for dof in 0..2 {
            assert_eq!(
                serial.analysis.domain().node(*a).displacement[dof].to_bits(),
                threaded.analysis.domain().node(*b).displacement[dof].to_bits()
            );
            assert_eq!(
                serial.analysis.domain().reaction(*a, dof, 1.).to_bits(),
                threaded.analysis.domain().reaction(*b, dof, 1.).to_bits()
            );
        }
    }
}
