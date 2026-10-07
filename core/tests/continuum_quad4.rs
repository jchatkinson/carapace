//! Plain `Quad4` against beam theory and the thick-cylinder (Lame) solution.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{Domain, Element, Node, NodeId, PlaneMaterial, Quad4};

fn solve(domain: Domain) -> Analysis {
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis.step().expect("linear solve");
    analysis
}

const L: f64 = 10.0;
const H: f64 = 1.0;
const T: f64 = 1.0;
const E: f64 = 1000.0;
const NU: f64 = 0.3;
const P: f64 = 1.0;

/// A cantilever of `nx` by `ny` elements, clamped at x = 0 (both translations),
/// unit tip shear shared over the tip nodes. Returns the mean tip deflection.
fn cantilever_tip(nx: usize, ny: usize) -> f64 {
    let mut domain = Domain::new();
    let mut ids = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let mut node = Node::new([L * i as f64 / nx as f64, H * j as f64 / ny as f64]);
            if i == 0 {
                node = node.fix(0).fix(1);
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |i: usize, j: usize| ids[j * (nx + 1) + i];
    let material = PlaneMaterial::plane_stress(E, NU).unwrap();
    for j in 0..ny {
        for i in 0..nx {
            let nodes = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            domain.add_element(Element::Quad4(Quad4::new(nodes, T, material.clone())));
        }
    }
    for j in 0..=ny {
        let share = if j == 0 || j == ny { 0.5 } else { 1.0 } / ny as f64;
        domain.load_node(at(nx, j), 1, -P * share);
    }
    let analysis = solve(domain);
    -(0..=ny)
        .map(|j| analysis.domain().node(at(nx, j)).displacement[1])
        .sum::<f64>()
        / (ny + 1) as f64
}

fn timoshenko_tip() -> f64 {
    let (i, a) = (T * H.powi(3) / 12.0, T * H);
    let g = E / (2.0 * (1.0 + NU));
    P * L.powi(3) / (3.0 * E * i) + P * L / (5.0 / 6.0 * g * a)
}

/// The plain element locks in bending: with one element through the depth the tip
/// deflection is far too small (documented, not hidden), and the ratio converges to
/// the Timoshenko value at roughly second order under refinement.
#[test]
fn cantilever_tip_deflection_converges_to_timoshenko_beam_theory() {
    let exact = timoshenko_tip();
    let locked = cantilever_tip(5, 1) / exact;
    assert!(locked < 0.45, "5x1 should lock badly, ratio {locked}");
    assert!(cantilever_tip(10, 1) / exact < 0.75, "10x1 still locks");

    // Refinement keeping the element aspect ratio: errors shrink by about four per halving.
    let ratios: Vec<f64> = [(10, 2), (20, 4), (40, 8), (80, 16)]
        .into_iter()
        .map(|(nx, ny)| cantilever_tip(nx, ny) / exact)
        .collect();
    let errors: Vec<f64> = ratios.iter().map(|r| 1.0 - r).collect();
    for pair in errors.windows(2) {
        assert!(
            pair[1] > 0.0 && pair[0] / pair[1] > 3.0,
            "errors {errors:?}"
        );
    }
    assert!(ratios[3] > 0.99 && ratios[3] < 1.0, "{ratios:?}");
}

/// Quarter of a thick cylinder, inner radius 1, outer radius 2, internal pressure
/// `p` (the outer surface free), symmetry on both axes. The pressure enters as
/// the consistent radial nodal forces of the inner arc (equal angular spacing:
/// `p t a dtheta` per node, half at the two ends). Returns the relative error of
/// the radial displacement at the inner and outer radius against Lame's
/// plane-stress solution `u_r = a^2 p / (E (b^2 - a^2)) ((1 - nu) r + (1 + nu) b^2 / r)`.
fn lame_errors(nr: usize, nt: usize, material: PlaneMaterial, plane_strain: bool) -> (f64, f64) {
    let (a, b, p) = (1.0, 2.0, 1.0);
    let mut domain = Domain::new();
    let mut ids: Vec<NodeId> = Vec::new();
    for k in 0..=nr {
        for l in 0..=nt {
            let r = a + (b - a) * k as f64 / nr as f64;
            let theta = std::f64::consts::FRAC_PI_2 * l as f64 / nt as f64;
            let mut node = Node::new([r * theta.cos(), r * theta.sin()]);
            if l == 0 {
                node = node.fix(1); // y = 0: no vertical motion
            }
            if l == nt {
                node = node.fix(0); // x = 0: no horizontal motion
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |k: usize, l: usize| ids[k * (nt + 1) + l];
    for k in 0..nr {
        for l in 0..nt {
            let nodes = [at(k, l), at(k + 1, l), at(k + 1, l + 1), at(k, l + 1)];
            domain.add_element(Element::Quad4(Quad4::new(nodes, T, material.clone())));
        }
    }
    let dtheta = std::f64::consts::FRAC_PI_2 / nt as f64;
    for l in 0..=nt {
        let share = if l == 0 || l == nt { 0.5 } else { 1.0 };
        let f = p * T * a * dtheta * share;
        let theta = dtheta * l as f64;
        domain.load_node(at(0, l), 0, f * theta.cos());
        domain.load_node(at(0, l), 1, f * theta.sin());
    }
    let analysis = solve(domain);
    let (e, nu) = match material {
        PlaneMaterial::PlaneStress { e, nu } | PlaneMaterial::PlaneStrain { e, nu } => (e, nu),
        _ => unreachable!(),
    };
    // Plane strain is plane stress with E' = E / (1 - nu^2), nu' = nu / (1 - nu).
    let (e, nu) = if plane_strain {
        (e / (1.0 - nu * nu), nu / (1.0 - nu))
    } else {
        (e, nu)
    };
    let exact =
        |r: f64| a * a * p / (e * (b * b - a * a)) * ((1.0 - nu) * r + (1.0 + nu) * b * b / r);
    let error_at = |k: usize, r: f64| {
        let worst = (0..=nt)
            .map(|l| {
                let theta = dtheta * l as f64;
                let u = analysis.domain().node(at(k, l)).displacement;
                let radial = u[0] * theta.cos() + u[1] * theta.sin();
                (radial / exact(r) - 1.0).abs()
            })
            .fold(0.0, f64::max);
        worst
    };
    (error_at(0, a), error_at(nr, b))
}

#[test]
fn thick_cylinder_matches_lame_on_a_coarse_mesh_and_converges() {
    let stress = PlaneMaterial::plane_stress(E, NU).unwrap();
    let (inner, outer) = lame_errors(2, 4, stress.clone(), false);
    assert!(
        inner < 0.02 && outer < 0.01,
        "coarse mesh: {inner}, {outer}"
    );
    let (inner_fine, outer_fine) = lame_errors(8, 16, stress, false);
    assert!(
        inner_fine < inner / 8.0 && outer_fine < outer / 6.0,
        "{inner_fine}, {outer_fine}"
    );

    let strain = PlaneMaterial::plane_strain(E, NU).unwrap();
    let (inner, outer) = lame_errors(2, 4, strain.clone(), true);
    assert!(
        inner < 0.03 && outer < 0.02,
        "plane strain: {inner}, {outer}"
    );
    let (inner_fine, _) = lame_errors(8, 16, strain, true);
    assert!(inner_fine < inner / 8.0);
}
