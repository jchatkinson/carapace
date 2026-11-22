//! Multi-storey stacked stud-wall shear walls with hold-downs, compared with the CSA O86:19 Annex A.11.7.1
//! four-term multi-storey deflection method (bending with cumulative rotation, panel shear, nail slip, and
//! cumulative anchorage elongation), for a 2-storey and a 5-storey wall line.
//!
//! Each storey is one `StudWall` (a single full-height pier, no openings). The leeward (compression) corner is
//! shared between storeys; the windward corner is joined by a `ZeroLength` hold-down (rigid in `ux`, `k_hd` in
//! `uy`). The base has the same arrangement against ground. The element is linear, so nail slip enters through
//! a per-storey secant shear stiffness at the storey's design shear, `1/ga = 1/Bv + 0.0025·e_n/v` with
//! `e_n = (0.013·v·s/d_F²)²` (O86 A.11.7): exact at that load level, the only one the code equation evaluates.
//! Every variant is also solved with the nonlinear `O86` shear law (Newton, displacement-increment test, since
//! the rigid stand-ins make an unbalance norm ill-conditioned) and must equal its secant run: 1e-8 relative
//! with real plates, 1e-6 where the stand-ins limit conditioning.
//!
//! The code-idealized wall (interior studs negligible, plates axially rigid, no bearing between the corners,
//! storey force delivered uniformly along the top plate) is the code's own idealization and must reproduce
//! the hand calculation: tolerance 1e-5 relative. The residual, about 2e-6, comes from the finite stand-in
//! plate stiffness `EA × 1e6`; stiffer stand-ins only add roundoff. The other variants switch on one real
//! effect at a time and are characterizations, printed with `--nocapture`. In every variant the pier's M/jd
//! chord force equals the O86 hold-down force by statics (1e-6 relative, limited by the conditioning of the
//! rigid stand-ins; `RIGID = 1e10` rather than larger for the same reason).

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    ChordDistribution, Domain, Element, ElementLoad, Material, Node, NodeId, StudWall,
    StudWallLayout, WallChord, WallLoad, WallSheathing, WallStation, WallStud, ZeroLength,
    WALL_SLOTS,
};

type D2 = Domain<2, 3, NodeId, Element>;

/// One 38×140 SPF stud, E = 9500 MPa: EA in N.
const EA_PLY: f64 = 9500. * 38. * 140.;
/// 12.5 mm OSB shear-through-thickness rigidity (N/mm), assumed for the comparison.
const BV: f64 = 11_000.;
/// Sheathing nail diameter (mm).
const D_F: f64 = 3.25;
/// Stand-in for a rigid link (N/mm).
const RIGID: f64 = 1e10;

struct Case {
    name: &'static str,
    length: f64,
    height: f64,
    /// Specified lateral force at each level, bottom to top (N).
    forces: Vec<f64>,
    /// End-post plies per storey, bottom to top.
    plies: Vec<f64>,
    /// Edge nail spacing per storey (mm).
    nail_spacing: Vec<f64>,
    /// Hold-down stiffness, the same at every level (N/mm).
    k_hd: f64,
}

/// One wall idealization: which real effects are switched on.
#[derive(Clone, Copy, Debug)]
struct Variant {
    name: &'static str,
    /// Interior studs have one ply's EA (otherwise negligible).
    interior_studs: bool,
    /// Plates have their real EA (otherwise axially rigid).
    real_plates: bool,
    /// Interior stud bottoms bear (bilaterally) on the storey below, or on the foundation with sill anchor
    /// bolts (fixed `ux`, `uy`). Otherwise the bottom plate connects only at its corners.
    bearing: bool,
    /// Storey force applied at the two top corners (otherwise uniformly along the top plate).
    corner_load: bool,
    sideways_only: bool,
}

const VARIANTS: [Variant; 6] = [
    Variant {
        name: "code idealized",
        interior_studs: false,
        real_plates: false,
        bearing: false,
        corner_load: false,
        sideways_only: false,
    },
    Variant {
        name: "+ interior studs",
        interior_studs: true,
        real_plates: false,
        bearing: false,
        corner_load: false,
        sideways_only: false,
    },
    Variant {
        name: "+ real plates, corner anchors",
        interior_studs: true,
        real_plates: true,
        bearing: false,
        corner_load: false,
        sideways_only: false,
    },
    Variant {
        name: "+ real plates, corner load",
        interior_studs: true,
        real_plates: true,
        bearing: false,
        corner_load: true,
        sideways_only: false,
    },
    Variant {
        name: "+ real plates, stud bearing",
        interior_studs: true,
        real_plates: true,
        bearing: true,
        corner_load: false,
        sideways_only: false,
    },
    Variant {
        name: "+ real plates, sideways only",
        interior_studs: true,
        real_plates: true,
        bearing: false,
        corner_load: false,
        sideways_only: true,
    },
];

fn two_storey() -> Case {
    Case {
        name: "2-storey",
        length: 2400.,
        height: 3000.,
        forces: vec![6000., 8000.],
        plies: vec![3., 2.],
        nail_spacing: vec![100., 150.],
        k_hd: 20_000.,
    }
}

fn five_storey() -> Case {
    Case {
        name: "5-storey",
        length: 4800.,
        height: 3000.,
        forces: vec![3000., 6000., 9000., 12000., 12000.],
        plies: vec![6., 5., 4., 3., 2.],
        nail_spacing: vec![75., 100., 100., 150., 150.],
        k_hd: 50_000.,
    }
}

impl Case {
    fn n(&self) -> usize {
        self.forces.len()
    }
    /// Storey shear `V_i = Σ_{j≥i} F_j`.
    fn shear(&self, i: usize) -> f64 {
        self.forces[i..].iter().sum()
    }
    /// Overturning moment at the top of storey `i`, `M_i = Σ_{j>i} V_j H_j` (O86 Figure A.6).
    fn moment_top(&self, i: usize) -> f64 {
        (i + 1..self.n()).map(|j| self.shear(j) * self.height).sum()
    }
    /// Hold-down force at the bottom of storey `i`: `M_{i-1} / L_s`.
    fn holddown(&self, i: usize) -> f64 {
        (self.moment_top(i) + self.shear(i) * self.height) / self.length
    }
    fn nail_slip(&self, i: usize) -> f64 {
        let v = self.shear(i) / self.length;
        (0.013 * v * self.nail_spacing[i] / (D_F * D_F)).powi(2)
    }
    /// Secant shear stiffness reproducing `vH/Bv + 0.0025·H·e_n` at the storey's design shear.
    fn secant_ga(&self, i: usize) -> f64 {
        let v = self.shear(i) / self.length;
        1. / (1. / BV + 0.0025 * self.nail_slip(i) / v)
    }
    fn interior_xs(&self) -> Vec<f64> {
        (1..)
            .map(|k| k as f64 * 400.)
            .take_while(|&x| x < self.length - 1e-9)
            .collect()
    }
}

/// O86:19 A.11.7.1 level deflections (bottom to top): `[total, bending, shear, nail slip, anchorage]`, each
/// the sum of the storey drifts up to that level.
fn o86_hand_calc(c: &Case) -> Vec<[f64; 5]> {
    let h = c.height;
    let (mut theta_sum, mut alpha_sum) = (0., 0.);
    let mut acc = [0.; 5];
    let mut out = Vec::new();
    for i in 0..c.n() {
        let ei = c.plies[i] * EA_PLY * c.length * c.length / 2.;
        let (m, v) = (c.moment_top(i), c.shear(i));
        let theta = m * h / ei + v * h * h / (2. * ei);
        let bending = m * h * h / (2. * ei) + v * h.powi(3) / (3. * ei) + h * theta_sum;
        let shear = v / c.length * h / BV;
        let slip = 0.0025 * h * c.nail_slip(i);
        alpha_sum += c.holddown(i) / c.k_hd / c.length;
        let anchorage = h * alpha_sum;
        theta_sum += theta;
        for (a, d) in acc[1..].iter_mut().zip([bending, shear, slip, anchorage]) {
            *a += d;
        }
        acc[0] = acc[1..].iter().sum();
        out.push(acc);
    }
    out
}

fn layout(c: &Case, i: usize, v: Variant, foundation: bool, nonlinear: bool) -> StudWallLayout {
    let post = c.plies[i] * EA_PLY;
    let interior = if v.interior_studs { EA_PLY } else { 1. };
    let mut studs = vec![WallStud { x: 0., ea: post }];
    studs.extend(
        c.interior_xs()
            .into_iter()
            .map(|x| WallStud { x, ea: interior }),
    );
    studs.push(WallStud {
        x: c.length,
        ea: post,
    });
    let (top_ea, bottom_ea) = if v.real_plates {
        (2. * EA_PLY, EA_PLY)
    } else {
        (EA_PLY * 1e6, EA_PLY * 1e6)
    };
    StudWallLayout {
        length: c.length,
        height: c.height,
        studs,
        top_ea,
        top_ei: Some(1.2e10),
        top_distribution: ChordDistribution::Beam,
        bottom_ea,
        bottom_ei: 6.0e9,
        openings: vec![],
        sheathing: Some(WallSheathing {
            law: if nonlinear {
                carapace_core::model::WallShearLaw::O86 {
                    bv: BV,
                    nail_spacing: c.nail_spacing[i],
                    nail_diameter: D_F,
                }
            } else {
                carapace_core::model::WallShearLaw::Linear { ga: c.secant_ga(i) }
            },
        }),
        foundation: foundation && v.bearing && i == 0,
    }
}

fn holddown(d: &mut D2, below: NodeId, above: NodeId, k_hd: f64) {
    d.add_element(Element::ZeroLength(
        ZeroLength::new(below, above)
            .with_material(0, Material::Elastic { e: RIGID })
            .with_material(1, Material::Elastic { e: k_hd }),
    ));
}

struct Result {
    /// Level deflections, bottom to top.
    drift: Vec<f64>,
    /// Hold-down forces below each storey, bottom to top.
    holddown: Vec<f64>,
    /// Left end-post axial force at the base and M/jd design force, by storey.
    left_stud_and_chord: Vec<[f64; 2]>,
    /// (base, top) stud-force moments, by storey.
    pier_moments: Vec<[f64; 2]>,
}

fn element_model(c: &Case, v: Variant) -> Result {
    element_model_mode(c, v, true, false)
}

fn element_model_mode(c: &Case, v: Variant, use_foundation: bool, nonlinear: bool) -> Result {
    let (l, h) = (c.length, c.height);
    let mut d: D2 = Domain::new();
    let ground = d.add_node(Node::new([0., 0.]).fix(0).fix(1).fix(2));
    let mut bottom_left = d.add_node(Node::new([0., 0.]));
    let mut bottom_right = d.add_node(Node::new([l, 0.]).fix(0).fix(1));
    holddown(&mut d, ground, bottom_left, c.k_hd);
    let mut links = vec![(ground, bottom_left)];
    let mut bottom_stations: Vec<(f64, NodeId)> = if v.sideways_only {
        c.interior_xs()
            .into_iter()
            .map(|x| (x, d.add_node(Node::new([x, 0.]).fix(0))))
            .collect()
    } else if v.bearing && !use_foundation {
        c.interior_xs()
            .into_iter()
            .map(|x| (x, d.add_node(Node::new([x, 0.]).fix(0).fix(1))))
            .collect()
    } else {
        vec![]
    };
    let pattern = d.default_pattern();
    let mut levels = Vec::new();
    let mut wall_ids = Vec::new();
    for i in 0..c.n() {
        let y = (i + 1) as f64 * h;
        let top_right = d.add_node(Node::new([l, y]));
        let top_left = d.add_node(Node::new([0., y]));
        let top_stations: Vec<(f64, NodeId)> = if (v.bearing || v.sideways_only) && i + 1 < c.n() {
            let xs = c.interior_xs();
            xs.into_iter()
                .map(|x| (x, d.add_node(Node::new([x, y]))))
                .collect()
        } else {
            vec![]
        };
        let stations = bottom_stations
            .iter()
            .map(|&(x, node)| (WallChord::Bottom, x, node))
            .chain(
                top_stations
                    .iter()
                    .map(|&(x, node)| (WallChord::Top, x, node)),
            )
            .map(|(chord, x, node)| WallStation { chord, x, node })
            .collect();
        let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| vec![]);
        if v.corner_load {
            d.add_nodal_load(pattern, top_left, 0, c.forces[i] / 2.);
            d.add_nodal_load(pattern, top_right, 0, c.forces[i] / 2.);
        } else {
            let w = c.forces[i] / l;
            loads[0].push(WallLoad::TopLine {
                x0: 0.,
                x1: l,
                wx0: w,
                wx1: w,
                wy0: 0.,
                wy1: 0.,
            });
        }
        let corners = [bottom_left, bottom_right, top_right, top_left];
        let wall = StudWall::new(
            corners,
            stations,
            layout(c, i, v, use_foundation, nonlinear),
            loads,
        )
        .unwrap();
        let id = d.add_element(Element::StudWall(wall));
        wall_ids.push(id);
        if !v.corner_load {
            let mut wall = [0.; WALL_SLOTS];
            wall[0] = 1.;
            d.add_element_load(
                pattern,
                id,
                ElementLoad {
                    wall,
                    ..Default::default()
                },
            );
        }
        levels.push((top_left, top_right));
        if i + 1 < c.n() {
            let next_left = d.add_node(Node::new([0., y]));
            holddown(&mut d, top_left, next_left, c.k_hd);
            links.push((top_left, next_left));
            bottom_left = next_left;
            bottom_right = top_right;
            if v.sideways_only {
                bottom_stations = top_stations
                    .iter()
                    .map(|&(x, upper)| {
                        let lower = d.add_node(Node::new([x, y]));
                        d.add_element(Element::ZeroLength(
                            ZeroLength::new(upper, lower)
                                .with_material(0, Material::Elastic { e: RIGID }),
                        ));
                        (x, lower)
                    })
                    .collect();
            } else {
                bottom_stations = top_stations;
            }
        }
    }
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
        .test(if nonlinear {
            ConvergenceTest::NormDispIncr {
                tol: 1e-9,
                max_iter: 100,
            }
        } else {
            ConvergenceTest::NormUnbalance {
                tol: 1e-13,
                max_iter: 4,
            }
        })
        .build(d);
    a.step().expect("stacked shear wall solve");
    let dom = a.domain();
    let u = |n: NodeId, k: usize| dom.node(n).displacement[k];
    Result {
        drift: levels
            .iter()
            .map(|&(tl, tr)| (u(tl, 0) + u(tr, 0)) / 2.)
            .collect(),
        holddown: links
            .iter()
            .map(|&(below, above)| c.k_hd * (u(above, 1) - u(below, 1)))
            .collect(),
        left_stud_and_chord: wall_ids
            .iter()
            .map(|&id| {
                let r = dom.element_stud_wall_response(id, 1.).unwrap();
                [r.studs.first().unwrap().bottom, r.piers[0].chords[0]]
            })
            .collect(),
        pier_moments: wall_ids
            .iter()
            .map(|&id| {
                let r = dom.element_stud_wall_response(id, 1.).unwrap();
                [r.piers[0].moment_base, r.piers[0].moment_top]
            })
            .collect(),
    }
}

fn compare(c: &Case) {
    let hand = o86_hand_calc(c);
    eprintln!(
        "\n{} — O86 A.11.7.1 level deflection (mm) = bending + shear + nail slip + anchorage",
        c.name
    );
    for (i, r) in hand.iter().enumerate() {
        eprintln!(
            "  level {}: {:7.3} = {:6.3} + {:6.3} + {:6.3} + {:6.3}   hold-down below: {:6.0} N",
            i + 1,
            r[0],
            r[1],
            r[2],
            r[3],
            r[4],
            c.holddown(i)
        );
    }
    eprintln!("  element / O86 (deflection per level; hold-down per level):");
    for v in VARIANTS {
        let r = element_model(c, v);
        let nonlinear = element_model_mode(c, v, true, true);
        if v.bearing {
            let station = element_model_mode(c, v, false, false);
            for (a, b) in r.drift.iter().zip(&station.drift) {
                near(*a, *b, 1e-9);
            }
            for (a, b) in r.holddown.iter().zip(&station.holddown) {
                near(*a, *b, 1e-9);
            }
            for (a, b) in r
                .left_stud_and_chord
                .iter()
                .flatten()
                .zip(station.left_stud_and_chord.iter().flatten())
            {
                near(*a, *b, 1e-9);
            }
        }
        let fmt = |xs: Vec<f64>| {
            xs.iter()
                .map(|x| format!("{x:6.3}"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let fmt_precise = |xs: Vec<f64>| {
            xs.iter()
                .map(|x| format!("{x:.9}"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        eprintln!(
            "  {:30} Δ {}   T {}   left stud / M/jd {}",
            v.name,
            fmt((0..c.n()).map(|i| r.drift[i] / hand[i][0]).collect()),
            fmt((0..c.n()).map(|i| r.holddown[i] / c.holddown(i)).collect()),
            fmt(r.left_stud_and_chord.iter().map(|f| f[0] / f[1]).collect()),
        );
        // The O86 law equals the secant at the design shear, so the two paths differ only by conditioning:
        // ~1e-9 with real plates, ~1e-7 with the `EA × 1e6` plate stand-ins, ~3e-8 with `RIGID` station links.
        let tol = if v.real_plates && !v.sideways_only {
            1e-8
        } else {
            1e-6
        };
        for (a, b) in r.drift.iter().zip(&nonlinear.drift) {
            near(*a, *b, tol);
        }
        for (a, b) in r.holddown.iter().zip(&nonlinear.holddown) {
            near(*a, *b, tol);
        }
        eprintln!(
            "    O86 law / secant Δ {}   T {}",
            fmt_precise(
                (0..c.n())
                    .map(|i| nonlinear.drift[i] / r.drift[i])
                    .collect()
            ),
            fmt_precise(
                (0..c.n())
                    .map(|i| nonlinear.holddown[i] / r.holddown[i])
                    .collect()
            ),
        );
        if v.name == "code idealized" {
            for (i, drift) in nonlinear.drift.iter().enumerate() {
                near(*drift, hand[i][0], 1e-5);
                near(nonlinear.holddown[i], c.holddown(i), 1e-5);
            }
        }
        for i in 0..c.n() {
            eprintln!(
                "    storey {}: M/jd {:.9} × O86; Mbase {:.6e}, Mtop {:.6e} (O86 Mtop {:.6e})",
                i + 1,
                r.left_stud_and_chord[i][1] / c.holddown(i),
                r.pier_moments[i][0],
                r.pier_moments[i][1],
                -c.moment_top(i),
            );
        }
        // Statics: the pier M/jd chord force is the O86 hold-down force `M_{i-1}/L_s`, and the moment
        // entering the top of each storey is the overturning from the storeys above, in every variant.
        // Exact up to the conditioning of the rigid stand-ins (`RIGID`, plates `EA × 1e6`): 1e-6 relative.
        for i in 0..c.n() {
            near(r.left_stud_and_chord[i][1], c.holddown(i), 1e-6);
            let scale = r.pier_moments[i][0].abs();
            assert!(
                (r.pier_moments[i][1] + c.moment_top(i)).abs() <= 1e-6 * scale,
                "{} {} storey {}: moment_top {} != {}",
                c.name,
                v.name,
                i + 1,
                r.pier_moments[i][1],
                -c.moment_top(i)
            );
        }
        if v.name == "code idealized" {
            for (i, (drift, t)) in r.drift.iter().zip(&r.holddown).enumerate() {
                near(*drift, hand[i][0], 1e-5);
                near(*t, c.holddown(i), 1e-5);
            }
        }
    }
}

fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()),
        "{a} != {b} (relative tolerance {tol})"
    );
}

/// 2-storey wall line (L = 2400, H = 3000, 3-ply/2-ply end posts, hold-downs at both levels) against the
/// O86:19 A.11.7.1 four-term multi-storey method. The code-idealized wall matches level deflections and
/// hold-down forces to 1e-5 relative; the other variants are printed characterizations.
#[test]
fn two_storey_shear_wall_matches_o86_multistorey_deflection() {
    compare(&two_storey());
}

/// 5-storey wall line (L = 4800, H = 3000, 6- to 2-ply end posts, hold-downs at every level) against the
/// O86:19 A.11.7.1 four-term multi-storey method, including cumulative bending and anchorage rotation. The
/// code-idealized wall matches to 1e-5 relative; the other variants are printed characterizations.
#[test]
fn five_storey_shear_wall_matches_o86_multistorey_deflection() {
    compare(&five_storey());
}
