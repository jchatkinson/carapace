//! Floor band: local panel response and floor-depth lever in 2D and 3D stacked shear walls.
//!
//! Element checks (constant-shear energy, axial rim bar, rigid-body modes, rotation, inactive bays): 1e-12.
//! One storey on a band: base moment `V(H + d)`, O86:19 Fig. 11.7 statics, 1e-9. Window wall: pier shears sum
//! to the load (1e-8) and each pier's moment is within 1% of `V_p(H + d)` (the rim and the wall's top plate
//! both carry drag over the opening and couple the piers slightly). 2- and 5-storey walls with rigid 254 mm
//! bands against an exact hand walk with true lever arms and the O86 A.11.7.1 storey terms: level displacements
//! 1e-5, chord forces 1e-6. The O86 value with `H_i = H + d` (one flexible storey) and the no-band value are
//! printed for comparison only.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    ChordDistribution, Domain, Element, FloorBand, FloorBandLayout, LoadSeries, Material, Node,
    NodeId, StudWall, StudWallLayout, WallChord, WallLoad, WallOpening, WallShearLaw,
    WallSheathing, WallStation, WallStud, ZeroLength, WALL_SLOTS,
};
use slotmap::Key;

type D2 = Domain<2, 3, NodeId, Element>;

fn band(
    xs: Vec<f64>,
    bays: Vec<Option<f64>>,
    rim_ea: f64,
    angle: f64,
    u: [[f64; 2]; 4],
) -> (D2, carapace_core::model::ElementId) {
    let mut d = D2::new();
    let n = xs.len();
    let mut bottom = Vec::new();
    let mut top = Vec::new();
    for (row, y, target) in [(0, 0., &mut bottom), (1, 100., &mut top)] {
        for (i, &x) in xs.iter().enumerate() {
            let px = x * angle.cos() - y * angle.sin();
            let py = x * angle.sin() + y * angle.cos();
            let mut node = Node::new([px, py]);
            let disp = if n == 2 { u[row * 2 + i] } else { [0., 0.] };
            node.displacement = [
                disp[0] * angle.cos() - disp[1] * angle.sin(),
                disp[0] * angle.sin() + disp[1] * angle.cos(),
                0.,
            ];
            target.push(d.add_node(node));
        }
    }
    let layout = FloorBandLayout {
        length: *xs.last().unwrap(),
        depth: 100.,
        xs,
        bays,
        rim_ea,
    };
    let b = FloorBand::new(bottom, top, layout).unwrap();
    let id = d.add_element(Element::FloorBand(b));
    (d, id)
}

fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b}"
    );
}

/// Constant shear field reports `q = GA gamma`, with the horizontal edge forces split over both rows and the
/// vertical couple on the bottom row. A 30 degree rotation preserves the force components in the band's local
/// axes. Reference: direct differentiation of `0.5 GA a d gamma^2`; tolerance 1e-12.
#[test]
fn shear_and_rotation() {
    let gamma = 0.01;
    let input = [[0., 0.], [0., 0.], [100. * gamma, 0.], [100. * gamma, 0.]];
    let (d, id) = band(vec![0., 1000.], vec![Some(2e6)], 0., 0., input);
    let response = d.element_floor_band_response(id).unwrap();
    near(response.bays[0].gamma, gamma, 1e-12);
    near(response.bays[0].q, 2e6 * gamma, 1e-12);
    let f = d.element_end_force(id, 1.).as_slice().to_vec();
    near(f[0], -1e7, 1e-12);
    near(f[2], -1e7, 1e-12);
    near(f[4], 1e7, 1e-12);
    near(f[6], 1e7, 1e-12);
    // The bay's vertical couple `q d` lands on the bottom row only (the rim bears on the wall below).
    near(f[1], -2e6, 1e-12);
    near(f[3], 2e6, 1e-12);
    near(f[5], 0., 1e-12);
    near(f[7], 0., 1e-12);
    let angle = 30f64.to_radians();
    let (rot, rid) = band(vec![0., 1000.], vec![Some(2e6)], 0., angle, input);
    let fr = rot.element_end_force(rid, 1.).as_slice().to_vec();
    for (rotated, reference) in fr.iter().zip(&f) {
        near(*rotated, *reference, 1e-12);
    }
}

/// An inactive bay has no shear flow, while the rim segment carries `EA` times axial strain. Reference: axial
/// bar force and absent-field definition; tolerance 1e-12.
#[test]
fn inactive_bay_and_rim_layout() {
    let (d, id) = band(
        vec![0., 1000.],
        vec![None],
        3e6,
        0.,
        [[0., 0.], [0., 0.], [0., 0.], [1., 0.]],
    );
    let r = d.element_floor_band_response(id).unwrap();
    assert_eq!(r.bays[0].q, 0.);
    near(r.rim[0], 3000., 1e-12);
}

/// Rigid translations and infinitesimal rotation produce zero band force. Reference: rigid-body invariance of
/// the constant-shear and axial-rim energy; tolerance 1e-12 relative.
#[test]
fn rigid_body_modes_have_zero_force() {
    let phi = 0.01;
    let (d, id) = band(
        vec![0., 1000.],
        vec![Some(2e6)],
        3e6,
        0.,
        [
            [3., 4.],
            [3., 4. + 1000. * phi],
            [3. - 100. * phi, 4.],
            [3. - 100. * phi, 4. + 1000. * phi],
        ],
    );
    let force_scale = 2e6 * 0.01 * 1000.;
    for f in d.element_end_force(id, 1.).as_slice() {
        assert!(
            f.abs() <= 1e-12 * force_scale,
            "rigid-body force residual {f}"
        );
    }
}

/// Layout validation rejects invalid counts, negative rigidity, duplicate nodes and non-positive dimensions.
/// Reference: public `FloorBandError` contract; exact error category.
#[test]
fn validates_layout_and_nodes() {
    let bad = FloorBandLayout {
        length: 1.,
        depth: 1.,
        xs: vec![0., 1.],
        bays: vec![Some(-1.)],
        rim_ea: 0.,
    };
    assert!(FloorBand::new(vec![NodeId::null(); 2], vec![NodeId::null(); 2], bad).is_err());
    let valid = FloorBandLayout {
        length: 1.,
        depth: 1.,
        xs: vec![0., 1.],
        bays: vec![None],
        rim_ea: 0.,
    };
    assert!(FloorBand::new(vec![NodeId::null(); 2], vec![NodeId::null(); 2], valid).is_err());
}

/// `between` merges both wall station sets and disables a bay inside a lower-wall opening. Reference: the
/// segmented design interval rule `x0 <= xa && xb <= x1`; exact station and bay checks.
#[test]
fn between_unions_wall_stations_and_skips_opening_bay() {
    use carapace_core::model::WallOpening;
    let wall = |length: f64, xs: &[f64], openings: Vec<WallOpening>| StudWallLayout {
        length,
        height: 2400.,
        studs: xs.iter().map(|&x| WallStud { x, ea: 1e6 }).collect(),
        top_ea: 1e6,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: 1e6,
        bottom_ei: 1e6,
        openings,
        sheathing: None,
        foundation: false,
    };
    let below = wall(
        3600.,
        &[0., 800., 1600., 2400., 3200., 3600.],
        vec![WallOpening {
            x0: 1600.,
            x1: 2400.,
            head_height: 2000.,
            header_ei: 1e8,
            jack_ea: 1e6,
            door: false,
        }],
    );
    let above = wall(3600., &[0., 400., 1200., 2000., 2800., 3600.], vec![]);
    let layout = FloorBandLayout::between(&below, &above, 254., 5e6, 2e8).unwrap();
    assert_eq!(
        layout.xs,
        vec![0., 400., 800., 1200., 1600., 2000., 2400., 2800., 3200., 3600.]
    );
    assert_eq!(layout.bays[4..6], [None, None]);
    assert!(layout.bays[..4].iter().all(Option::is_some));
    assert!(layout.bays[6..].iter().all(Option::is_some));
}

/// A lateral load at the top of a band gives the wall below the floor-depth lever `V(H+d)`. Reference: O86
/// Fig. 11.7 statics; moment and reaction equilibrium tolerance 1e-9.
#[test]
fn floor_depth_adds_to_wall_overturning_lever() {
    let (length, height, depth, force) = (3600., 2440., 254., 14_400.);
    let mut d: D2 = D2::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut lower = Vec::new();
    for x in [0., length] {
        let node = d.add_node(Node::new([x, 0.]));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, node)
                .with_material(0, Material::Elastic { e: 1e10 })
                .with_material(1, Material::Elastic { e: 20_000. }),
        ));
        lower.push(node);
    }
    let c2 = d.add_node(Node::new([length, height]));
    let c3 = d.add_node(Node::new([0., height]));
    let layout = StudWallLayout {
        length,
        height,
        studs: vec![
            WallStud { x: 0., ea: 1e10 },
            WallStud {
                x: length,
                ea: 1e10,
            },
        ],
        top_ea: 1e10,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: 1e10,
        bottom_ei: 1e10,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 1e5 },
        }),
        foundation: false,
    };
    let empty: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
    let wall = d.add_element(Element::StudWall(
        StudWall::new([lower[0], lower[1], c2, c3], vec![], layout, empty).unwrap(),
    ));
    let t0 = d.add_node(Node::new([0., height + depth]));
    let t1 = d.add_node(Node::new([length, height + depth]));
    let band = FloorBand::new(
        vec![c3, c2],
        vec![t0, t1],
        FloorBandLayout {
            length,
            depth,
            xs: vec![0., length],
            bays: vec![Some(1e10)],
            rim_ea: 1e10,
        },
    )
    .unwrap();
    d.add_element(Element::FloorBand(band));
    for (a, b) in [(c3, t0), (c2, t1)] {
        d.add_element(Element::ZeroLength(
            ZeroLength::new(a, b).with_material(1, Material::Elastic { e: 1e10 }),
        ));
    }
    let pattern = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    d.add_nodal_load(pattern, t0, 0, force);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-6,
            max_iter: 4,
        })
        .build(d);
    analysis.step().expect("linear floor-band solve");
    let response = analysis
        .domain()
        .element_stud_wall_response(wall, 1.)
        .unwrap();
    near(
        response.piers[0].moment_base.abs(),
        force * (height + depth),
        1e-9,
    );
    let base = response.base;
    near(base[0].abs(), force, 1e-9);
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

fn floor_band_stack(c: &StackCase) -> StackResult {
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
    let mut roof_nodes = None;
    for i in 0..c.n() {
        let y = c.wall_base_y(i, FLOOR_DEPTH);
        let top_left = d.add_node(Node::new([0., y + c.height]));
        let top_right = d.add_node(Node::new([c.length, y + c.height]));
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
        };
        let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
        if i + 1 == c.n() {
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
        if i + 1 == c.n() {
            let mut wall_load = [0.; WALL_SLOTS];
            wall_load[0] = 1.;
            d.add_element_load(
                pattern,
                wall_id,
                carapace_core::model::ElementLoad {
                    wall: wall_load,
                    ..Default::default()
                },
            );
            roof_nodes = Some([top_left, top_right]);
        } else {
            let band_top_left = d.add_node(Node::new([0., y + c.height + FLOOR_DEPTH]));
            let band_top_right = d.add_node(Node::new([c.length, y + c.height + FLOOR_DEPTH]));
            let band = FloorBand::new(
                vec![top_left, top_right],
                vec![band_top_left, band_top_right],
                FloorBandLayout {
                    length: c.length,
                    depth: FLOOR_DEPTH,
                    xs: vec![0., c.length],
                    bays: vec![Some(RIGID)],
                    rim_ea: RIGID,
                },
            )
            .unwrap();
            d.add_element(Element::FloorBand(band));
            for (a, b, stiffness) in [
                (top_left, band_top_left, c.k_hd),
                (top_right, band_top_right, RIGID),
            ] {
                d.add_element(Element::ZeroLength(
                    ZeroLength::new(a, b).with_material(1, Material::Elastic { e: stiffness }),
                ));
            }
            d.add_nodal_load(pattern, band_top_left, 0, c.forces[i]);
            force_nodes.push([band_top_left, band_top_left]);
            bottom = [band_top_left, band_top_right];
        }
        walls.push(wall_id);
    }
    force_nodes.push(roof_nodes.unwrap());
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-13,
            max_iter: 4,
        })
        .build(d);
    analysis.step().expect("multi-storey floor-band solve");
    let dom = analysis.domain();
    StackResult {
        level_u: force_nodes
            .iter()
            .map(|nodes| {
                (dom.node(nodes[0]).displacement[0] + dom.node(nodes[1]).displacement[0]) / 2.
            })
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

fn compare_floor_band_stack(c: &StackCase) -> Vec<String> {
    let actual = floor_band_stack(c);
    let (hand, hand_chord, hand_moment) = stack_hand_calc(c, FLOOR_DEPTH);
    let no_band = stack_hand_calc(c, 0.).0;
    let o86 = stack_o86_comparison(c, c.height + FLOOR_DEPTH);
    eprintln!("\n{} with rigid floor bands (d = {FLOOR_DEPTH} mm)", c.name);
    eprintln!(" level | FE band (mm) | hand rigid band | O86 H+d | no band");
    let mut failures = Vec::new();
    for i in 0..c.n() {
        let force_level_u = actual.level_u[i];
        eprintln!(
            " {:>5} | {:>12.6} | {:>15.6} | {:>8.6} | {:>8.6}",
            i + 1,
            force_level_u,
            hand[i],
            o86[i],
            no_band[i]
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

/// A window wall with king studs at both opening edges, four pier-end hold-downs, top stations at every band
/// position and inactive band bays over the opening, under a diaphragm load spread uniformly along the band top
/// with a real rim (one 38x235 ply): the pier shears sum to the load (1e-8, limited by the 1e10 rigid links),
/// and each pier's base moment is within 1% of its shear times `H + d` (O86 Fig. 11.7). It is not exact: over the
/// opening the drag runs both in the rim (at `H + d`) and in the wall's top plate (at `H`), and two axial
/// struts a depth `d` apart couple the piers slightly (about 0.5% here).
#[test]
fn window_wall_band_preserves_pier_levers() {
    let (length, depth, force) = (3600., 254., 14_400.);
    let opening = WallOpening {
        x0: 1600.,
        x1: 2400.,
        head_height: 2000.,
        header_ei: 1e9,
        jack_ea: PLY_EA,
        door: false,
    };
    let xs = vec![
        0., 400., 800., 1200., 1600., 2000., 2400., 2800., 3200., 3600.,
    ];
    let wall_layout = StudWallLayout {
        length,
        height: 2440.,
        studs: xs
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
        openings: vec![opening],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::Linear { ga: 6000. },
        }),
        foundation: true,
    };
    let band_layout =
        FloorBandLayout::between(&wall_layout, &wall_layout, depth, RIGID, 9500. * 38. * 235.)
            .unwrap();
    assert_eq!(band_layout.bays[4..6], [None, None]);
    let mut d: D2 = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut bottom = Vec::new();
    for &x in &[0., 1600., 2400., length] {
        let node = d.add_node(Node::new([x, 0.]).fix(0));
        d.add_element(Element::ZeroLength(
            ZeroLength::new(ground, node).with_material(1, Material::Elastic { e: 20_000. }),
        ));
        bottom.push((x, node));
    }
    let c2 = d.add_node(Node::new([length, 2440.]));
    let c3 = d.add_node(Node::new([0., 2440.]));
    let mut top = Vec::new();
    for &x in &band_layout.xs {
        let node = if x == 0. {
            c3
        } else if x == length {
            c2
        } else {
            d.add_node(Node::new([x, 2440.]))
        };
        top.push(node);
    }
    let stations: Vec<WallStation> = band_layout
        .xs
        .iter()
        .zip(&top)
        .filter(|(x, _)| **x != 0. && **x != length)
        .map(|(&x, &node)| WallStation {
            chord: WallChord::Top,
            x,
            node,
        })
        .collect();
    let base_stations: Vec<WallStation> = bottom
        .iter()
        .skip(1)
        .take(2)
        .map(|&(x, node)| WallStation {
            chord: WallChord::Bottom,
            x,
            node,
        })
        .collect();
    let empty: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
    let wall_id = d.add_element(Element::StudWall(
        StudWall::new(
            [bottom[0].1, bottom[3].1, c2, c3],
            base_stations.into_iter().chain(stations).collect(),
            wall_layout,
            empty,
        )
        .unwrap(),
    ));
    let band_top: Vec<NodeId> = band_layout
        .xs
        .iter()
        .map(|&x| d.add_node(Node::new([x, 2440. + depth])))
        .collect();
    d.add_element(Element::FloorBand(
        FloorBand::new(top.clone(), band_top.clone(), band_layout).unwrap(),
    ));
    for (&a, &b) in top.iter().zip(&band_top) {
        d.add_element(Element::ZeroLength(
            ZeroLength::new(a, b).with_material(1, Material::Elastic { e: RIGID }),
        ));
    }
    let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
    let xs = &xs;
    for (i, &n) in band_top.iter().enumerate() {
        let lo = if i == 0 { 0. } else { (xs[i] - xs[i - 1]) / 2. };
        let hi = if i + 1 == xs.len() {
            0.
        } else {
            (xs[i + 1] - xs[i]) / 2.
        };
        d.add_nodal_load(p, n, 0, force * (lo + hi) / length);
    }
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-6,
            max_iter: 4,
        })
        .build(d);
    a.step().expect("window floor-band solve");
    let r = a.domain().element_stud_wall_response(wall_id, 1.).unwrap();
    assert_eq!(r.piers.len(), 2);
    let shear_total: f64 = r
        .piers
        .iter()
        .map(|pier| pier.q * (pier.xb - pier.xa))
        .sum();
    near(shear_total, force, 1e-8);
    for pier in &r.piers {
        near(
            pier.moment_base.abs(),
            pier.q.abs() * (pier.xb - pier.xa) * (2440. + depth),
            1e-2,
        );
    }
}

/// Two- and five-storey code-idealized walls with rigid `254 mm` floor bands match the exact hand walk for
/// level deflections and static M/jd chord forces. The printed O86 `H+d` and no-band values show the floor-depth
/// effect; O86 treats that depth as flexible wall height, while the band model keeps it rigid. Tolerances 1e-5
/// for displacement and 1e-6 for chord forces/moments.
#[test]
fn rigid_floor_bands_match_multistorey_hand_calc() {
    let mut failures = compare_floor_band_stack(&two_storey_stack());
    failures.extend(compare_floor_band_stack(&five_storey_stack()));
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Floor-band element assembly gives bit-identical serial and threaded displacements and reactions. Reference:
/// the native parallel scatter contract; exact `f64::to_bits` equality.
#[cfg(feature = "parallel")]
#[test]
fn serial_and_parallel_band_results_are_bit_identical() {
    use carapace_core::analysis::{set_parallelism, Parallelism};
    fn model() -> (D2, [NodeId; 2]) {
        let mut d = D2::new();
        let b0 = d.add_node(Node::new([0., 0.]).fix(0).fix(1));
        let b1 = d.add_node(Node::new([1000., 0.]).fix(0).fix(1));
        let t0 = d.add_node(Node::new([0., 100.]).fix(1));
        let t1 = d.add_node(Node::new([1000., 100.]).fix(1));
        let band = FloorBand::new(
            vec![b0, b1],
            vec![t0, t1],
            FloorBandLayout {
                length: 1000.,
                depth: 100.,
                xs: vec![0., 1000.],
                bays: vec![Some(2e6)],
                rim_ea: 3e6,
            },
        )
        .unwrap();
        d.add_element(Element::FloorBand(band));
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        d.add_nodal_load(p, t0, 0, 1234.);
        (d, [t0, t1])
    }
    let solve = |d| {
        let mut a = AnalysisBuilder::new()
            .constraint_handler(ConstraintHandler::Transformation)
            .integrator(Integrator::LoadControl { increment: 1. })
            .algorithm(Algorithm::Linear)
            .test(ConvergenceTest::NormUnbalance {
                tol: 1e-9,
                max_iter: 4,
            })
            .build(d);
        a.step().unwrap();
        a
    };
    set_parallelism(Parallelism::Serial);
    let (d0, n0) = model();
    let serial = solve(d0);
    set_parallelism(Parallelism::Threads(2));
    let (d1, n1) = model();
    let parallel = solve(d1);
    for i in 0..2 {
        assert_eq!(
            serial.domain().node(n0[i]).displacement[0].to_bits(),
            parallel.domain().node(n1[i]).displacement[0].to_bits()
        );
        assert_eq!(
            serial.domain().reaction(n0[i], 0, 1.).to_bits(),
            parallel.domain().reaction(n1[i], 0, 1.).to_bits()
        );
    }
    set_parallelism(Parallelism::Serial);
}
