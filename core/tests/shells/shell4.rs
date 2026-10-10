//! `Shell4` (the MITC4 flat shell) against OpenSees `ShellMITC4` and closed-form plate theory.
//!
//! The OpenSees references were produced with OpenSeesPy (`ShellMITC4` + `ElasticMembranePlateSection`)
//! on the exact models built here; the tolerance is 1e-8 relative to the largest value.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest,
    Integrator,
};
use carapace_core::model::{
    Domain3, Element3, Element3Id, ElementLoad3, Node3, Node3Id, Shell4, ShellSection,
};

fn solve(domain: Domain3) -> Analysis<3, 6, Node3Id, Element3> {
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

fn assert_close(name: &str, got: &[f64], want: &[f64], scale: f64) {
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(
            (g - w).abs() <= 1e-8 * scale,
            "{name}[{i}]: got {g}, OpenSees {w}"
        );
    }
}

/// Two distorted quads on the tilted plane `z = 0.3 x + 0.1 y`, nodes 1 and 4 clamped,
/// point loads on nodes 5 and 6.
fn tilted_patch(self_weight: bool) -> (Analysis<3, 6, Node3Id, Element3>, Vec<Node3Id>) {
    let (analysis, ids, _) = tilted_patch_elements(self_weight);
    (analysis, ids)
}

fn tilted_patch_elements(
    self_weight: bool,
) -> (
    Analysis<3, 6, Node3Id, Element3>,
    Vec<Node3Id>,
    Vec<Element3Id>,
) {
    let mut domain = Domain3::new();
    let xy = [
        (0.0, 0.0),
        (2.0, 0.2),
        (2.3, 1.7),
        (-0.1, 1.4),
        (4.2, 0.1),
        (4.5, 1.6),
    ];
    let ids: Vec<Node3Id> = xy
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| {
            let mut node = Node3::new([x, y, 0.3 * x + 0.1 * y]);
            if i == 0 || i == 3 {
                node = (0..6).fold(node, |n, d| n.fix(d));
            }
            domain.add_node(node)
        })
        .collect();
    let section = ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 2e-3).unwrap();
    let a = domain.add_element(Element3::Shell4(Shell4::new(
        [ids[0], ids[1], ids[2], ids[3]],
        section.clone(),
    )));
    let b = domain.add_element(Element3::Shell4(Shell4::new(
        [ids[1], ids[4], ids[5], ids[2]],
        section,
    )));
    for (dof, v) in [1.0, -2.0, 3.0, 0.5, -0.7, 0.2].into_iter().enumerate() {
        domain.load_node(ids[4], dof, v);
    }
    for (dof, v) in [-0.5, 1.5, -4.0, 0.1, 0.3, -0.9].into_iter().enumerate() {
        domain.load_node(ids[5], dof, v);
    }
    if self_weight {
        let pattern = domain.default_pattern();
        for e in [a, b] {
            domain.add_element_load(pattern, e, ElementLoad3::body(0.0, 0.0, 9.81));
        }
    }
    (solve(domain), ids, vec![a, b])
}

fn displacements(
    analysis: &Analysis<3, 6, Node3Id, Element3>,
    ids: &[Node3Id],
    nodes: &[usize],
) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|&n| analysis.domain().node(ids[n - 1]).displacement)
        .collect()
}

fn reactions(analysis: &Analysis<3, 6, Node3Id, Element3>, ids: &[Node3Id]) -> Vec<f64> {
    [1, 4]
        .iter()
        .flat_map(|&n| (0..6).map(move |d| (n, d)))
        .map(|(n, d)| analysis.domain().reaction(ids[n - 1], d, 1.0))
        .collect()
}

/// Displacements and clamp reactions of two distorted, tilted quads under point loads.
#[test]
fn tilted_distorted_patch_matches_opensees() {
    let (analysis, ids) = tilted_patch(false);
    #[rustfmt::skip]
    let u = [
        0.007990161495, 0.002571574432, -0.026470016599, -0.025269653013, 0.027840741278, -0.004873468159,
        0.022971754648, 0.007582000056, -0.07665193686, -0.025250773291, 0.056701209378, -0.001639744718,
        0.032863282016, 0.011036245391, -0.108117982242, -0.054959707678, 0.05198534551, -0.010401219149,
        0.063758791612, 0.021889234115, -0.213901519408, -0.060055311292, 0.069112918456, -0.011504155255,
    ];
    assert_close(
        "u",
        &displacements(&analysis, &ids, &[2, 3, 5, 6]),
        &u,
        0.22,
    );
    #[rustfmt::skip]
    let r = [
        0.61050538, 0.377650149, -1.9379159604, 0.8819314137, -1.191351685, 0.2349107973,
        -1.11050538, 0.122349851, 2.9379159604, 0.2434447253, -4.4952843192, -0.1273833442,
    ];
    assert_close("reaction", &reactions(&analysis, &ids), &r, 4.5);
}

/// Self-weight loads each node with `rho h` times its tributary area times the body vector.
/// OpenSees' `ShellMITC4` applies the opposite sign (`-rho h b`), so the reference run's
/// `-selfWeight 0 0 -9.81` is `body(0, 0, 9.81)` here.
#[test]
fn self_weight_matches_opensees() {
    let (analysis, ids) = tilted_patch(true);
    #[rustfmt::skip]
    let u = [
        0.007792643317, 0.002507963778, -0.025803241067, -0.025178725878, 0.027194792958, -0.004908751132,
        0.022685408794, 0.007489908113, -0.075697553691, -0.025244815065, 0.055921694505, -0.001713965308,
        0.032155942744, 0.010806864944, -0.105747978579, -0.054893757487, 0.05108379486, -0.010469481701,
        0.062939225536, 0.021624192929, -0.21116828581, -0.059986904599, 0.068273037741, -0.011565387538,
    ];
    assert_close(
        "u",
        &displacements(&analysis, &ids, &[2, 3, 5, 6]),
        &u,
        0.22,
    );
    #[rustfmt::skip]
    let r = [
        0.6047371292, 0.3754064199, -1.9599166596, 0.8710044742, -1.1392640098, 0.2354102215,
        -1.1047371292, 0.1245935801, 2.9064559779, 0.2519141916, -4.4282059216, -0.1195828443,
    ];
    assert_close("reaction", &reactions(&analysis, &ids), &r, 4.5);
}

/// The first element's committed strain and resultants at its four Gauss points (OpenSees'
/// `eleResponse(1, 'stresses')` order and sign convention).
#[test]
fn gauss_point_resultants_match_opensees() {
    let (analysis, _, elements) = tilted_patch_elements(false);
    let responses = analysis
        .domain()
        .element_shell_responses(elements[0])
        .unwrap();
    #[rustfmt::skip]
    let want = [
        [0.15300997729714821, 0.0974009509245366, -0.23308159393440123, 2.7744668677482176, 0.3973853405481096, 0.7672074855918711, 3.8671655973183383, -0.5842252743485605],
        [0.17242928940233437, 0.1851868064060666, -0.37431818354232427, 2.645743301963054, -0.15656032296915945, 1.3122547152760908, 3.9950292388576747, 0.008179010645220794],
        [-0.05086409883474502, 0.11513846224336377, -0.33075103842361564, 3.517578161013844, 0.15089568368938155, 1.0970782600071913, -5.1923860413437835, 1.9555771235774841],
        [-0.07303903767993183, 0.03464760731300698, -0.20270484297288047, 3.651777120998986, 0.6559790378035427, 0.6019449432183839, -5.014740841576701, 1.2955148590313055],
    ];
    assert_eq!(responses.len(), 4);
    for (g, (response, want)) in responses.iter().zip(&want).enumerate() {
        assert_close(
            &format!("resultant at point {g}"),
            &response.resultant,
            want,
            5.0,
        );
    }
}

/// A pressure on a distorted quad is a consistent load: the clamp reactions sum to `-p A e3`,
/// with `e3` the normal given by the node order (here up the tilted plane).
#[test]
fn pressure_resultant_is_p_area_along_the_normal() {
    let mut domain = Domain3::new();
    let xy = [(0.0, 0.0), (2.0, 0.2), (2.3, 1.7), (-0.1, 1.4)];
    let z = |x: f64, y: f64| 0.3 * x + 0.1 * y;
    let ids: Vec<Node3Id> = xy
        .iter()
        .map(|&(x, y)| {
            let node = (0..6).fold(Node3::new([x, y, z(x, y)]), |n, d| n.fix(d));
            domain.add_node(node)
        })
        .collect();
    let section = ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 0.0).unwrap();
    let shell = domain.add_element(Element3::Shell4(Shell4::new(
        [ids[0], ids[1], ids[2], ids[3]],
        section,
    )));
    let pattern = domain.default_pattern();
    domain.add_element_load(pattern, shell, ElementLoad3::pressure(2.5));
    let analysis = solve(domain);

    let p = |i: usize| {
        let (x, y) = xy[i];
        nalgebra::Vector3::new(x, y, z(x, y))
    };
    let cross = (p(2) - p(0)).cross(&(p(3) - p(1)));
    let area = 0.5 * cross.norm();
    let normal = cross.normalize();
    let mut total = nalgebra::Vector3::zeros();
    for &id in &ids {
        for d in 0..3 {
            total[d] += analysis.domain().reaction(id, d, 1.0);
        }
    }
    let want = -2.5 * area * normal;
    assert!(normal.z > 0.0);
    assert_close("total reaction", total.as_slice(), want.as_slice(), area);
}

/// A simply supported square plate under uniform pressure: the centre deflection is
/// `0.004062 p a^4 / D` (thin-plate Navier series).
#[test]
fn simply_supported_square_plate_centre_deflection() {
    let (a, h, e, nu, p) = (10.0, 0.2, 2.0e5, 0.3, 0.01);
    let n = 8;
    let mut domain = Domain3::new();
    let mut at = vec![];
    for j in 0..=n {
        for i in 0..=n {
            let (x, y) = (a * i as f64 / n as f64, a * j as f64 / n as f64);
            let mut node = Node3::new([x, y, 0.0]);
            if i == 0 || j == 0 || i == n || j == n {
                node = node.fix(2);
            }
            if i == 0 && j == 0 {
                node = node.fix(0).fix(1);
            }
            if i == n && j == 0 {
                node = node.fix(1);
            }
            at.push(domain.add_node(node));
        }
    }
    let section = ShellSection::elastic_membrane_plate(e, nu, h, 0.0).unwrap();
    let pattern = domain.default_pattern();
    for j in 0..n {
        for i in 0..n {
            let ids = [
                at[j * (n + 1) + i],
                at[j * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i],
            ];
            let shell = domain.add_element(Element3::Shell4(Shell4::new(ids, section.clone())));
            domain.add_element_load(pattern, shell, ElementLoad3::pressure(p));
        }
    }
    let analysis = solve(domain);
    let centre = analysis
        .domain()
        .node(at[(n / 2) * (n + 1) + n / 2])
        .displacement[2];
    let d = e * h.powi(3) / (12.0 * (1.0 - nu * nu));
    let exact = 0.004062 * p * a.powi(4) / d;
    assert!(
        (centre / exact - 1.0).abs() < 0.01,
        "centre deflection {centre}, plate theory {exact}"
    );
}

/// The fundamental frequency of a simply supported square plate is
/// `2 pi^2 / a^2 sqrt(D / (rho h))`. Carapace lumps the shell mass (OpenSees uses the
/// consistent mass), so this is a closed-form check, not an OpenSees match.
#[test]
fn simply_supported_square_plate_fundamental_frequency() {
    let (a, h, e, nu, rho) = (10.0, 0.2, 2.0e5, 0.3, 1.0e-3);
    let n = 8;
    let mut domain = Domain3::new();
    let mut at = vec![];
    for j in 0..=n {
        for i in 0..=n {
            let (x, y) = (a * i as f64 / n as f64, a * j as f64 / n as f64);
            let mut node = Node3::new([x, y, 0.0]);
            if i == 0 || j == 0 || i == n || j == n {
                node = node.fix(2);
            }
            if i == 0 && j == 0 {
                node = node.fix(0).fix(1);
            }
            if i == n && j == 0 {
                node = node.fix(1);
            }
            at.push(domain.add_node(node));
        }
    }
    let section = ShellSection::elastic_membrane_plate(e, nu, h, rho).unwrap();
    for j in 0..n {
        for i in 0..n {
            let ids = [
                at[j * (n + 1) + i],
                at[j * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i],
            ];
            domain.add_element(Element3::Shell4(Shell4::new(ids, section.clone())));
        }
    }
    let modes = modal_analysis(&mut domain, 1).expect("modal analysis");
    let d = e * h.powi(3) / (12.0 * (1.0 - nu * nu));
    let exact = 2.0 * std::f64::consts::PI.powi(2) / (a * a) * (d / (rho * h)).sqrt();
    assert!(
        (modes[0].frequency / exact - 1.0).abs() < 0.02,
        "omega {}, plate theory {exact}",
        modes[0].frequency
    );
}

/// A quad with one node lifted out of the plane: the shell is treated as flat in its mean frame, exactly
/// as OpenSees does, so the two agree.
#[test]
fn warped_quad_matches_opensees() {
    let mut domain = Domain3::new();
    let coords = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.1, 1.6, 0.35],
        [0.0, 1.5, 0.1],
    ];
    let ids: Vec<Node3Id> = coords
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let node = Node3::new(c);
            domain.add_node(if i == 0 || i == 3 {
                (0..6).fold(node, |n, d| n.fix(d))
            } else {
                node
            })
        })
        .collect();
    let section = ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 0.0).unwrap();
    domain.add_element(Element3::Shell4(Shell4::new(
        [ids[0], ids[1], ids[2], ids[3]],
        section,
    )));
    for (dof, v) in [0.7, -0.3, 1.1, 0.2, -0.1, 0.05].into_iter().enumerate() {
        domain.load_node(ids[1], dof, v);
    }
    for (dof, v) in [-0.4, 0.9, -2.0, 0.1, 0.3, -0.2].into_iter().enumerate() {
        domain.load_node(ids[2], dof, v);
    }
    let analysis = solve(domain);
    #[rustfmt::skip]
    let want = [
        0.0003951192081206001, 0.0006565489175840232, -0.001354850918113836, -0.012718400126295543, 0.002622879541392172, 0.0001565639841421267,
        0.0010162814698049458, 0.003843794827473608, -0.023071703428927916, -0.012086555364308228, 0.01950868983916031, 0.002276775183219031,
    ];
    let got = displacements(&analysis, &ids, &[2, 3]);
    assert_close("warped", &got, &want, 0.025);
}

/// A column framing into a shell corner: the shell stiffens all six DOFs, so the beam's rotations
/// and the shell's share one node with no extra constraint.
#[test]
fn shell_and_beam_sharing_a_node_match_opensees() {
    use carapace_core::model::{ElasticBeamColumn3, GeomTransf3};
    let mut domain = Domain3::new();
    let coords = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.5, 0.0],
        [0.0, 1.5, 0.0],
        [2.0, 1.5, 1.2],
    ];
    let ids: Vec<Node3Id> = coords
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let node = Node3::new(c);
            domain.add_node(if i == 0 || i == 3 {
                (0..6).fold(node, |n, d| n.fix(d))
            } else {
                node
            })
        })
        .collect();
    let section = ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 0.0).unwrap();
    domain.add_element(Element3::Shell4(Shell4::new(
        [ids[0], ids[1], ids[2], ids[3]],
        section,
    )));
    domain.add_element(Element3::ElasticBeamColumn3(ElasticBeamColumn3::new(
        ids[2],
        ids[4],
        3e4,
        1.2e4,
        0.05,
        2e-3,
        1e-4,
        2e-4,
        GeomTransf3::linear([1.0, 0.0, 0.0]),
    )));
    for (dof, v) in [0.5, -0.25, 0.0, 0.0, 0.1, 0.0].into_iter().enumerate() {
        domain.load_node(ids[4], dof, v);
    }
    let analysis = solve(domain);
    #[rustfmt::skip]
    let want = [
        0.000200645212817101, -0.00026234467386033843, -0.006506458818532393, 0.0012373139414665803, 0.007021946044397192, -0.0001933884102067844,
        0.1286269804660937, -0.02574712140362023, -0.006506458818532393, 0.031237313941466577, 0.16702194604439718, -0.0001933884102067844,
    ];
    let got = displacements(&analysis, &ids, &[3, 5]);
    assert_close("shell+beam", &got, &want, 0.17);
}

/// The Scordelis-Lo barrel vault roof (quarter model, 8x8 mesh) against OpenSees, and against the
/// published vertical deflection of 0.3024 at the midpoint of the free edge: a mesh-converged value
/// of about 0.30 (0.2965 at 12x12 in OpenSees).
#[test]
fn scordelis_lo_roof_matches_opensees() {
    let (r, half, t, e) = (25.0_f64, 25.0, 0.25, 4.32e8);
    let theta = 40.0_f64.to_radians();
    let n = 8;
    let mut domain = Domain3::new();
    let mut ids = vec![];
    for j in 0..=n {
        for i in 0..=n {
            let a = theta * i as f64 / n as f64;
            let y = half * (1.0 - j as f64 / n as f64);
            let mut node = Node3::new([r * a.sin(), y, r * a.cos()]);
            // Symmetry at the crown (x = 0) and at mid-length (y = 0); a rigid diaphragm at the end.
            if i == 0 {
                node = node.fix(0).fix(4).fix(5);
            }
            if j == n {
                node = node.fix(1).fix(3).fix(5);
            }
            if j == 0 {
                node = node.fix(0).fix(2);
            }
            ids.push(domain.add_node(node));
        }
    }
    let section = ShellSection::elastic_membrane_plate(e, 0.0, t, 360.0).unwrap();
    let pattern = domain.default_pattern();
    for j in 0..n {
        for i in 0..n {
            let at = |i: usize, j: usize| ids[j * (n + 1) + i];
            let shell = domain.add_element(Element3::Shell4(Shell4::new(
                [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)],
                section.clone(),
            )));
            domain.add_element_load(pattern, shell, ElementLoad3::body(0.0, 0.0, -1.0));
        }
    }
    let analysis = solve(domain);
    let corner = analysis.domain().node(ids[n * (n + 1) + n]).displacement;
    assert_close(
        "free-edge midpoint",
        &[corner[0], corner[2]],
        &[-0.15373219605989855, -0.29150057494553583],
        0.3,
    );
    assert!((corner[2] / -0.3024 - 1.0).abs() < 0.05, "uz {}", corner[2]);
}
