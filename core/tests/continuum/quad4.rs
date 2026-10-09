//! Plain `Quad4` against beam theory and the thick-cylinder (Lame) solution.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{Domain, Element, Node, NodeId, PlaneMaterial, Quad4, Quad4Formulation};

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
    plane_stress_cantilever(nx, ny, Quad4Formulation::Full)
}

fn plane_stress_cantilever(nx: usize, ny: usize, formulation: Quad4Formulation) -> f64 {
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
            domain.add_element(Element::Quad4(
                Quad4::new(nodes, T, material.clone()).with_formulation(formulation),
            ));
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
fn lame_errors(
    nr: usize,
    nt: usize,
    material: PlaneMaterial,
    plane_strain: bool,
    formulation: Quad4Formulation,
) -> (f64, f64) {
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
            domain.add_element(Element::Quad4(
                Quad4::new(nodes, T, material.clone()).with_formulation(formulation),
            ));
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

/// Thick-walled cylinder against the Lame closed-form radial and hoop stresses, for plane stress
/// and plane strain. A coarse mesh is already within a few percent, and refining the mesh reduces
/// the error by the expected factors.
#[test]
fn thick_cylinder_matches_lame_on_a_coarse_mesh_and_converges() {
    let stress = PlaneMaterial::plane_stress(E, NU).unwrap();
    let (inner, outer) = lame_errors(2, 4, stress.clone(), false, Quad4Formulation::Full);
    assert!(
        inner < 0.02 && outer < 0.01,
        "coarse mesh: {inner}, {outer}"
    );
    let (inner_fine, outer_fine) = lame_errors(8, 16, stress, false, Quad4Formulation::Full);
    assert!(
        inner_fine < inner / 8.0 && outer_fine < outer / 6.0,
        "{inner_fine}, {outer_fine}"
    );

    let strain = PlaneMaterial::plane_strain(E, NU).unwrap();
    let (inner, outer) = lame_errors(2, 4, strain.clone(), true, Quad4Formulation::Full);
    assert!(
        inner < 0.03 && outer < 0.02,
        "plane strain: {inner}, {outer}"
    );
    let (inner_fine, _) = lame_errors(8, 16, strain, true, Quad4Formulation::Full);
    assert!(inner_fine < inner / 8.0);
}

/// MacNeal-Harder-style cantilever (L = 6, h = 0.2, t = 0.1, E = 1e7, nu = 0.3, unit tip
/// shear) of `n` elements in a row, one through the depth, with the interior element
/// edges shifted by `shape` to make parallelograms or trapezoids. Returns the tip
/// deflection (positive in the load direction).
fn macneal_harder(n: usize, shape: fn(usize) -> (f64, f64), formulation: Quad4Formulation) -> f64 {
    let (len, depth, t, e, nu) = (6.0, 0.2, 0.1, 1e7, 0.3);
    let mut domain = Domain::new();
    let mut bottom = Vec::new();
    let mut top = Vec::new();
    for i in 0..=n {
        let (db, dt) = shape(i);
        let x = len * i as f64 / n as f64;
        let fixed = |node: Node| if i == 0 { node.fix(0).fix(1) } else { node };
        bottom.push(domain.add_node(fixed(Node::new([x + db, 0.0]))));
        top.push(domain.add_node(fixed(Node::new([x + dt, depth]))));
    }
    let material = PlaneMaterial::plane_stress(e, nu).unwrap();
    for i in 0..n {
        let nodes = [bottom[i], bottom[i + 1], top[i + 1], top[i]];
        domain.add_element(Element::Quad4(
            Quad4::new(nodes, t, material.clone()).with_formulation(formulation),
        ));
    }
    domain.load_node(bottom[n], 1, -0.5);
    domain.load_node(top[n], 1, -0.5);
    let analysis = solve(domain);
    -0.5 * (analysis.domain().node(bottom[n]).displacement[1]
        + analysis.domain().node(top[n]).displacement[1])
}

fn rectangles(_: usize) -> (f64, f64) {
    (0.0, 0.0)
}

/// Every interior top edge is shifted right by 0.3 of an element: parallelograms
/// (the clamped and tip edges stay vertical).
fn parallelograms(i: usize) -> (f64, f64) {
    if i == 0 || i == 6 {
        (0.0, 0.0)
    } else {
        (0.0, 0.3)
    }
}

/// Interior edges shifted in opposite directions top and bottom: trapezoids.
fn trapezoids(i: usize) -> (f64, f64) {
    if i == 0 || i == 6 {
        (0.0, 0.0)
    } else {
        (-0.2, 0.2)
    }
}

/// Cook's membrane: corners (0,0), (48,44), (48,60), (0,44); clamped at x = 0, unit
/// shear spread uniformly along the free edge x = 48. Returns the vertical
/// displacement of the top right corner.
fn cook(n: usize, material: PlaneMaterial, formulation: Quad4Formulation) -> f64 {
    let mut domain = Domain::new();
    let mut ids = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (s, t) = (i as f64 / n as f64, j as f64 / n as f64);
            let x = 48.0 * s;
            let y = 44.0 * s + t * (44.0 + 16.0 * s - 44.0 * s);
            let mut node = Node::new([x, y]);
            if i == 0 {
                node = node.fix(0).fix(1);
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |i: usize, j: usize| ids[j * (n + 1) + i];
    for j in 0..n {
        for i in 0..n {
            let nodes = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            domain.add_element(Element::Quad4(
                Quad4::new(nodes, 1.0, material.clone()).with_formulation(formulation),
            ));
        }
    }
    for j in 0..=n {
        let share = if j == 0 || j == n { 0.5 } else { 1.0 } / n as f64;
        domain.load_node(at(n, j), 1, share);
    }
    solve(domain).domain().node(at(n, n)).displacement[1]
}

/// Pure bending by an end couple (consistent nodal forces of a linear end traction), nu = 0 so
/// the clamped cantilever's exact field `u = (k x y, -k x^2 / 2)` also satisfies the
/// clamp. The enhanced element reproduces it exactly with one element through the
/// depth on rectangles; the plain element locks.
fn pure_bending_tip(nx: usize, formulation: Quad4Formulation) -> (f64, f64) {
    let moment = 1.0;
    let mut domain = Domain::new();
    let mut ids = Vec::new();
    for j in 0..=1 {
        for i in 0..=nx {
            let mut node = Node::new([L * i as f64 / nx as f64, H * j as f64 - H / 2.0]);
            if i == 0 {
                node = node.fix(0).fix(1);
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |i: usize, j: usize| ids[j * (nx + 1) + i];
    let material = PlaneMaterial::plane_stress(E, 0.0).unwrap();
    for i in 0..nx {
        let nodes = [at(i, 0), at(i + 1, 0), at(i + 1, 1), at(i, 1)];
        domain.add_element(Element::Quad4(
            Quad4::new(nodes, T, material.clone()).with_formulation(formulation),
        ));
    }
    // Couple: +F at the top node, -F at the bottom, F h = M.
    domain.load_node(at(nx, 1), 0, moment / H);
    domain.load_node(at(nx, 0), 0, -moment / H);
    let analysis = solve(domain);
    let i = T * H.powi(3) / 12.0;
    // Tension on top (+F there) curves the cantilever so its tip moves down: M L^2 / (2 E I).
    let exact = -moment * L * L / (2.0 * E * i);
    (analysis.domain().node(at(nx, 0)).displacement[1], exact)
}

/// Pure bending on a rectangular mesh: the enhanced Quad4 gives the exact beam-theory tip
/// deflection for any mesh (1, 4 or 10 elements along the span), while the plain element is too
/// stiff, and a single plain element locks (less than 70% of the exact answer).
#[test]
fn enhanced_pure_bending_is_exact_on_rectangles_and_beats_the_plain_element() {
    for nx in [1, 4, 10] {
        let (got, exact) = pure_bending_tip(nx, Quad4Formulation::Enhanced);
        assert!(
            ((got - exact) / exact).abs() < 1e-10,
            "nx {nx}: {got} vs {exact}"
        );
        let (plain, _) = pure_bending_tip(nx, Quad4Formulation::Full);
        assert!(
            plain / exact < 1.0 - 1e-6,
            "plain is too stiff: {}",
            plain / exact
        );
        if nx == 1 {
            assert!(
                plain / exact < 0.7,
                "one plain element locks: {}",
                plain / exact
            );
        }
    }
}

/// Cantilever tip load against the Timoshenko tip deflection on a 10x2 mesh: the enhanced element
/// is within 3% of the exact answer and at least 0.2 closer than the plain element.
#[test]
fn enhanced_cantilever_converges_faster_than_the_plain_element() {
    let exact = timoshenko_tip();
    let enhanced = enhanced_cantilever_tip(10, 2) / exact;
    let plain = cantilever_tip(10, 2) / exact;
    assert!(
        enhanced > 0.97 && enhanced > plain + 0.2,
        "{enhanced} vs {plain}"
    );
}

fn enhanced_cantilever_tip(nx: usize, ny: usize) -> f64 {
    plane_stress_cantilever(nx, ny, Quad4Formulation::Enhanced)
}

/// MacNeal-Harder cantilever results, recorded (not asserted to be better than the plain
/// element): rectangles are near beam theory, distorted meshes are not, for either
/// formulation. The numbers are regression values and are compared with OpenSees'
/// `enhancedQuad` in `comparison/`.
#[test]
fn macneal_harder_results_are_recorded() {
    let (len, depth, t, e, nu): (f64, f64, f64, f64, f64) = (6.0, 0.2, 0.1, 1e7, 0.3);
    let i = t * depth.powi(3) / 12.0;
    let beam = len.powi(3) / (3.0 * e * i) + len / (5.0 / 6.0 * e / (2.0 * (1.0 + nu)) * t * depth);
    let rect = macneal_harder(6, rectangles, Quad4Formulation::Enhanced);
    assert!(((rect - beam) / beam).abs() < 0.01, "{rect} vs {beam}");
    for (name, shape, full, enhanced) in [
        (
            "rectangles",
            rectangles as fn(usize) -> (f64, f64),
            0.010088,
            0.107328,
        ),
        ("parallelograms", parallelograms, 0.002330, 0.061063),
        ("trapezoids", trapezoids, 0.002135, 0.066750),
    ] {
        let f = macneal_harder(6, shape, Quad4Formulation::Full);
        let q = macneal_harder(6, shape, Quad4Formulation::Enhanced);
        assert!(((f - full) / full).abs() < 1e-3, "{name} full: {f}");
        assert!(
            ((q - enhanced) / enhanced).abs() < 1e-3,
            "{name} enhanced: {q}"
        );
    }
}

/// Nearly incompressible plane strain (nu = 0.4999) on Cook's membrane. The plain element
/// locks completely; the enhanced element, which is *not* a mixed formulation and carries
/// no volumetric-locking guarantee, is observed to behave much better here and in the
/// thick-cylinder test. Recorded as observations (the 64x64 enhanced result is the
/// reference), not as a property of the formulation.
#[test]
fn cook_membrane_at_nu_4999_records_the_plain_elements_locking() {
    let material = PlaneMaterial::plane_strain(1.0, 0.4999).unwrap();
    let reference = cook(64, material.clone(), Quad4Formulation::Enhanced);
    let plain = cook(4, material.clone(), Quad4Formulation::Full) / reference;
    let enhanced = cook(4, material, Quad4Formulation::Enhanced) / reference;
    assert!(plain < 0.4, "plain {plain}");
    assert!(enhanced > 0.8 && enhanced < 1.0, "enhanced {enhanced}");
}

/// Thick cylinder under internal pressure in plane strain with Poisson's ratio 0.4999 (nearly
/// incompressible): the enhanced element stays within 1% of the Lame solution, while the plain
/// element locks (more than 50% error).
#[test]
fn enhanced_thick_cylinder_is_accurate_even_at_nu_4999() {
    let (inner, outer) = lame_errors(
        4,
        8,
        PlaneMaterial::plane_strain(E, 0.4999).unwrap(),
        true,
        Quad4Formulation::Enhanced,
    );
    assert!(inner < 0.01 && outer < 0.01, "{inner}, {outer}");
    let (inner, _) = lame_errors(
        4,
        8,
        PlaneMaterial::plane_strain(E, 0.4999).unwrap(),
        true,
        Quad4Formulation::Full,
    );
    assert!(inner > 0.5, "the plain element locks: {inner}");
}
