//! Stud-wall loads belong to load patterns: a wall owns a list of load cases and each pattern selects one with a
//! factor (`ElementLoad::wall`), so any number of patterns may load one wall. Verified against the same wall
//! loaded by one pattern carrying the combined factored loads (linear algebra, 1e-12 relative to the largest
//! response of each kind; nonlinear Newton, 1e-9), against hand sums of the applied load, and for the validation
//! errors (case out of range, two cases from one pattern). 2D and Z-up 3D walls, with a window opening so that
//! headers, jacks and cripples are loaded, and with the O86 pier law so the state depends on the load.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    ChordDistribution, Domain, Domain3, Element, Element3, ElementLoad, ElementLoad3, LoadSeries,
    ModelError, Node, Node3, StudWall, StudWallLayout, WallCaseLoad, WallLoad, WallOpening,
    WallShearLaw, WallSheathing, WallStud,
};

const L: f64 = 3600.;
const H: f64 = 2440.;

fn layout(law: WallShearLaw) -> StudWallLayout {
    StudWallLayout {
        length: L,
        height: H,
        studs: (0..=9)
            .map(|i| WallStud {
                x: 400. * i as f64,
                ea: 1e8,
            })
            .collect(),
        top_ea: 1e12,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: 1e12,
        bottom_ei: 1e12,
        openings: vec![WallOpening {
            x0: 1200.,
            x1: 2400.,
            head_height: 2000.,
            header_ei: 1e12,
            jack_ea: 1e8,
            door: false,
        }],
        sheathing: Some(WallSheathing { law }),
        foundation: true,
        bearing: None,
        floor: None,
    }
}

fn linear() -> WallShearLaw {
    WallShearLaw::Linear { ga: 1e7 }
}

fn o86() -> WallShearLaw {
    WallShearLaw::O86 {
        bv: 100.,
        nail_spacing: 100.,
        nail_diameter: 3.25,
    }
}

/// Load case `k` of the family used below: a vertical line load, a point load with a horizontal part, and a
/// horizontal line load, each different per case.
fn case_loads(k: usize) -> Vec<WallLoad> {
    let f = 1. + k as f64;
    vec![
        WallLoad::TopLine {
            x0: 100. * f,
            x1: 1500. + 150. * f,
            wx0: 0.,
            wx1: 0.,
            wy0: -2. * f,
            wy1: -3. * f,
        },
        WallLoad::TopPoint {
            x: 2900. - 130. * f,
            fx: 50. * f,
            fy: -400. * f,
        },
        WallLoad::TopLine {
            x0: 0.,
            x1: L,
            wx0: 0.1 * f,
            wx1: 0.2 * f,
            wy0: 0.,
            wy1: 0.,
        },
    ]
}

fn scaled(loads: &[WallLoad], s: f64) -> Vec<WallLoad> {
    loads
        .iter()
        .map(|l| match *l {
            WallLoad::TopPoint { x, fx, fy } => WallLoad::TopPoint {
                x,
                fx: fx * s,
                fy: fy * s,
            },
            WallLoad::TopLine {
                x0,
                x1,
                wx0,
                wx1,
                wy0,
                wy1,
            } => WallLoad::TopLine {
                x0,
                x1,
                wx0: wx0 * s,
                wx1: wx1 * s,
                wy0: wy0 * s,
                wy1: wy1 * s,
            },
        })
        .collect()
}

/// Net applied `(horizontal, vertical)` force of a case list, from the load definitions.
fn applied(loads: &[WallLoad]) -> (f64, f64) {
    loads.iter().fold((0., 0.), |(h, v), l| match *l {
        WallLoad::TopPoint { fx, fy, .. } => (h + fx, v + fy),
        WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => (
            h + 0.5 * (wx0 + wx1) * (x1 - x0),
            v + 0.5 * (wy0 + wy1) * (x1 - x0),
        ),
    })
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Profile {
    TwoD,
    ZUp3D,
}

/// Quantities compared between a patterned and a combined model.
#[derive(Debug)]
struct Out {
    /// Top-corner displacements, in-plane `[x, up]` of both top corners.
    top: Vec<f64>,
    /// Base reaction `(fx, fy)` per base node.
    base: Vec<f64>,
    /// Pier shear flows, header shears and moments, top-chord axial forces.
    internal: Vec<f64>,
    /// `Domain::element_end_force` (the wall's load-inclusive local force).
    end_force: Vec<f64>,
}

/// Solves the wall of `cases` loaded by `patterns` = `(pattern scale factor, case)`, in `steps` equal steps.
fn solve(
    profile: Profile,
    law: WallShearLaw,
    cases: Vec<Vec<WallLoad>>,
    patterns: &[(f64, usize)],
    steps: usize,
) -> Out {
    let algorithm = if matches!(law, WallShearLaw::Linear { .. }) {
        Algorithm::Linear
    } else {
        Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        }
    };
    // A macro, not a closure: the two profiles' domains have different types.
    macro_rules! analysis {
        ($d:expr) => {
            AnalysisBuilder::new()
                .constraint_handler(ConstraintHandler::Transformation)
                .integrator(Integrator::LoadControl {
                    increment: 1. / steps as f64,
                })
                .algorithm(algorithm.clone())
                .test(ConvergenceTest::NormUnbalance {
                    tol: 1e-6,
                    max_iter: 60,
                })
                .build($d)
        };
    }
    let case = |case: usize| {
        Some(WallCaseLoad {
            case: case as u32,
            factor: 1.,
        })
    };
    let (mut top, response, end_force);
    match profile {
        Profile::TwoD => {
            let mut d: Domain = Domain::new();
            let c = [
                d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
                d.add_node(Node::new([L, 0.]).fix(0).fix(1)),
                d.add_node(Node::new([L, H])),
                d.add_node(Node::new([0., H])),
            ];
            let wall = StudWall::new(c, vec![], layout(law), cases).unwrap();
            let id = d.add_element(Element::StudWall(wall));
            for &(scale, k) in patterns {
                let p = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, scale);
                d.add_element_load(
                    p,
                    id,
                    ElementLoad {
                        wall: case(k),
                        ..Default::default()
                    },
                );
            }
            let mut a = analysis!(d);
            for _ in 0..steps {
                a.step().unwrap();
            }
            top = [c[2], c[3]]
                .iter()
                .flat_map(|&n| [0, 1].map(|dof| a.domain().node(n).displacement[dof]))
                .collect::<Vec<_>>();
            response = a.domain().element_stud_wall_response(id, 1.).unwrap();
            end_force = a.domain().element_end_force(id, 1.).as_slice().to_vec();
        }
        Profile::ZUp3D => {
            let mut d = Domain3::new();
            // The wall lies in the Y-Z plane with Z up: local x is +Y, local y is +Z.
            let c = [
                d.add_node(Node3::new([0., 0., 0.]).fix(0).fix(1).fix(2)),
                d.add_node(Node3::new([0., L, 0.]).fix(0).fix(1).fix(2)),
                d.add_node(Node3::new([0., L, H])),
                d.add_node(Node3::new([0., 0., H])),
            ];
            let wall = StudWall::new3(c, vec![], layout(law), cases).unwrap();
            let id = d.add_element(Element3::StudWall(wall));
            for &(scale, k) in patterns {
                let p = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1. }, scale);
                d.add_element_load(
                    p,
                    id,
                    ElementLoad3 {
                        wall: case(k),
                        ..Default::default()
                    },
                );
            }
            let mut a = analysis!(d);
            for _ in 0..steps {
                a.step().unwrap();
            }
            top = [c[2], c[3]]
                .iter()
                .flat_map(|&n| [1, 2].map(|dof| a.domain().node(n).displacement[dof]))
                .collect::<Vec<_>>();
            response = a.domain().element_stud_wall_response(id, 1.).unwrap();
            end_force = a.domain().element_end_force(id, 1.).as_slice().to_vec();
        }
    }
    let base = response
        .base_nodes
        .iter()
        .flat_map(|n| [n.fx, n.fy])
        .collect();
    let mut internal: Vec<f64> = response.piers.iter().map(|p| p.q).collect();
    internal.extend(response.header.iter().flat_map(|h| [h.moment, h.shear]));
    internal.extend(response.top_chord.iter().map(|c| c.axial));
    top.shrink_to_fit();
    Out {
        top,
        base,
        internal,
        end_force,
    }
}

/// Asserts `a` and `b` agree to `tol` relative to the largest magnitude in the vector.
fn agree(a: &[f64], b: &[f64], tol: f64, what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: length");
    let scale = a.iter().chain(b).fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(scale > 0., "{what}: all zero");
    let worst = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    println!(
        "{what}: worst difference {:.3e} of scale {scale:.3e}",
        worst / scale
    );
    assert!(
        worst <= tol * scale,
        "{what}: {worst:e} > {tol:e} * {scale:e}"
    );
}

fn compare(x: &Out, y: &Out, tol: f64, what: &str) {
    agree(&x.top, &y.top, tol, &format!("{what} top displacement"));
    agree(&x.base, &y.base, tol, &format!("{what} base reaction"));
    agree(&x.internal, &y.internal, tol, &format!("{what} internal"));
    agree(
        &x.end_force,
        &y.end_force,
        tol,
        &format!("{what} end force"),
    );
}

/// Dead x 1.25 and live x 1.5 as two patterns equal one pattern with the combined factored loads, in 2D and Z-up
/// 3D (linear sheathing, 1e-12 relative): the effective load of a wall sums over patterns that select different
/// cases. The base vertical reaction also equals the hand sum of the applied vertical load.
#[test]
fn two_patterns_with_different_factors_equal_the_combined_pattern() {
    for profile in [Profile::TwoD, Profile::ZUp3D] {
        let (dead, live) = (case_loads(0), case_loads(1));
        let mut combined = scaled(&dead, 1.25);
        combined.extend(scaled(&live, 1.5));
        let split = solve(
            profile,
            linear(),
            vec![dead.clone(), live.clone()],
            &[(1.25, 0), (1.5, 1)],
            1,
        );
        let one = solve(profile, linear(), vec![combined], &[(1., 0)], 1);
        compare(&split, &one, 1e-12, &format!("{profile:?}"));
        let (_, vd) = applied(&dead);
        let (_, vl) = applied(&live);
        let vertical: f64 = split.base.iter().skip(1).step_by(2).sum();
        let expected = -(1.25 * vd + 1.5 * vl);
        assert!(
            (vertical - expected).abs() <= 1e-9 * expected.abs(),
            "{profile:?}: base vertical {vertical} vs applied {expected}"
        );
    }
}

/// More than 8 patterns on one wall: twelve cases, each in its own pattern with its own factor, equal the single
/// combined pattern (1e-12), and the base vertical and horizontal reactions equal the hand sums of the applied
/// loads (1e-9; the wall is statically supported in the plane).
#[test]
fn twelve_patterns_load_one_wall() {
    let factors: Vec<f64> = (0..12).map(|k| 0.5 + 0.15 * k as f64).collect();
    let cases: Vec<_> = (0..12).map(case_loads).collect();
    let mut combined = Vec::new();
    for (k, f) in factors.iter().enumerate() {
        combined.extend(scaled(&cases[k], *f));
    }
    let patterns: Vec<_> = factors.iter().copied().zip(0..).collect();
    for profile in [Profile::TwoD, Profile::ZUp3D] {
        let many = solve(profile, linear(), cases.clone(), &patterns, 1);
        let one = solve(profile, linear(), vec![combined.clone()], &[(1., 0)], 1);
        compare(&many, &one, 1e-12, &format!("{profile:?} 12 patterns"));
        let (h, v) = cases.iter().zip(&factors).fold((0., 0.), |(h, v), (c, f)| {
            let (ch, cv) = applied(c);
            (h + f * ch, v + f * cv)
        });
        let fx: f64 = many.base.iter().step_by(2).sum();
        let fy: f64 = many.base.iter().skip(1).step_by(2).sum();
        assert!(
            (fx + h).abs() <= 1e-9 * h.abs() && (fy + v).abs() <= 1e-9 * v.abs(),
            "{profile:?}: base ({fx}, {fy}) vs applied ({h}, {v})"
        );
    }
}

/// The O86 pier law makes the wall state depend on the total load, so the effective cases enter the nonlinear
/// pier solve: dead x 1.25 and live x 1.5 in two patterns equal one combined pattern through four Newton steps
/// (1e-9 relative).
#[test]
fn nonlinear_walls_sum_the_cases_of_all_patterns() {
    for profile in [Profile::TwoD, Profile::ZUp3D] {
        let (dead, live) = (case_loads(0), case_loads(1));
        let mut combined = scaled(&dead, 1.25);
        combined.extend(scaled(&live, 1.5));
        let split = solve(profile, o86(), vec![dead, live], &[(1.25, 0), (1.5, 1)], 4);
        let one = solve(profile, o86(), vec![combined], &[(1., 0)], 4);
        compare(&split, &one, 1e-9, &format!("{profile:?} O86"));
    }
}

/// Several patterns may select the same case: the factors add.
#[test]
fn patterns_selecting_one_case_add_their_factors() {
    let loads = case_loads(2);
    let two = solve(
        Profile::TwoD,
        linear(),
        vec![loads.clone()],
        &[(0.75, 0), (1.5, 0)],
        1,
    );
    let one = solve(Profile::TwoD, linear(), vec![loads], &[(2.25, 0)], 1);
    compare(&two, &one, 1e-12, "same case");
}

/// A case index past the wall's cases, and two different cases added to one element in one pattern (which the
/// pattern cannot represent), are model errors rather than silently dropped or merged.
#[test]
fn invalid_wall_cases_are_validation_errors() {
    let build = |first: u32, second: Option<u32>| {
        let mut d: Domain = Domain::new();
        let c = [
            d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
            d.add_node(Node::new([L, 0.]).fix(0).fix(1)),
            d.add_node(Node::new([L, H])),
            d.add_node(Node::new([0., H])),
        ];
        let wall = StudWall::new(
            c,
            vec![],
            layout(linear()),
            vec![case_loads(0), case_loads(1)],
        )
        .unwrap();
        let id = d.add_element(Element::StudWall(wall));
        let p = d.add_load_pattern(LoadSeries::Linear { slope: 1. });
        for case in std::iter::once(first).chain(second) {
            d.add_element_load(
                p,
                id,
                ElementLoad {
                    wall: Some(WallCaseLoad { case, factor: 1. }),
                    ..Default::default()
                },
            );
        }
        d.validate()
    };
    assert!(build(1, None).is_ok());
    assert!(build(1, Some(1)).is_ok());
    assert!(matches!(
        build(2, None),
        Err(ModelError::IncompatibleElementLoad { .. })
    ));
    assert!(matches!(
        build(0, Some(1)),
        Err(ModelError::IncompatibleElementLoad { .. })
    ));
}
