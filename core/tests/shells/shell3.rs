//! `Shell3` (the DKT/Allman flat triangle) against OpenSees `ShellDKGT`, on three distorted triangles on the tilted
//! plane `z = 0.3 x + 0.1 y`. Reference values are from OpenSeesPy (`ShellDKGT` + `ElasticMembranePlateSection`).

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain3, Element3, Element3Id, ElementLoad3, Node3, Node3Id, Shell3, ShellSection,
};

type A = Analysis<3, 6, Node3Id, Element3>;

fn patch(body: Option<f64>) -> (A, Vec<Node3Id>, Vec<Element3Id>) {
    let mut domain = Domain3::new();
    let xy = [(0.0, 0.0), (2.0, 0.2), (2.3, 1.7), (-0.1, 1.4), (4.2, 0.1)];
    let ids: Vec<Node3Id> = xy
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| {
            let node = Node3::new([x, y, 0.3 * x + 0.1 * y]);
            domain.add_node(if i == 0 || i == 3 {
                (0..6).fold(node, |n, d| n.fix(d))
            } else {
                node
            })
        })
        .collect();
    let section = ShellSection::elastic_membrane_plate(3e4, 0.25, 0.4, 2e-3).unwrap();
    let elements: Vec<Element3Id> = [[0, 1, 2], [0, 2, 3], [1, 4, 2]]
        .iter()
        .map(|t| {
            domain.add_element(Element3::Shell3(Shell3::new(
                [ids[t[0]], ids[t[1]], ids[t[2]]],
                section.clone(),
            )))
        })
        .collect();
    for (dof, v) in [1.0, -2.0, 3.0, 0.5, -0.7, 0.2].into_iter().enumerate() {
        domain.load_node(ids[4], dof, v);
    }
    for (dof, v) in [-0.5, 1.5, -4.0, 0.1, 0.3, -0.9].into_iter().enumerate() {
        domain.load_node(ids[2], dof, v);
    }
    if let Some(bz) = body {
        let pattern = domain.default_pattern();
        for &e in &elements {
            domain.add_element_load(pattern, e, ElementLoad3::body(0.0, 0.0, bz));
        }
    }
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
    (analysis, ids, elements)
}

fn check(name: &str, got: &[f64], want: &[f64], scale: f64) {
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(
            (g - w).abs() <= 1e-8 * scale,
            "{name}[{i}]: got {g}, OpenSees {w}"
        );
    }
}

fn displacements(analysis: &A, ids: &[Node3Id], nodes: &[usize]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|&n| analysis.domain().node(ids[n - 1]).displacement)
        .collect()
}

fn reactions(analysis: &A, ids: &[Node3Id]) -> Vec<f64> {
    [1, 4]
        .iter()
        .flat_map(|&n| (0..6).map(move |d| (n, d)))
        .map(|(n, d)| analysis.domain().reaction(ids[n - 1], d, 1.0))
        .collect()
}

/// Displacements, clamp reactions and the first triangle's four Gauss-point resultants.
#[test]
fn triangle_patch_matches_opensees() {
    let (analysis, ids, elements) = patch(None);
    #[rustfmt::skip]
    let u = [
        -0.014252395946010287, -0.005802946837076639, 0.04514517192159253, -0.015761165796629638, -0.04761939116729468, -0.011503508018656981,
        -0.005639334645083282, -0.003939420483684128, 0.022462710212293367, -0.03165095499114575, -0.03451615688302404, -0.015083441363777745,
        -0.06623758349489814, -0.028278257772653304, 0.21662847499794943, -0.03716998390490014, -0.11139964702985417, -0.02700235006241788,
    ];
    check("u", &displacements(&analysis, &ids, &[2, 3, 5]), &u, 0.22);
    #[rustfmt::skip]
    let r = [
        2.4122577017030284, 0.3587065279484115, -3.0249510569195497, -1.2009890580275089, 2.4138861416797686, -0.05546332918383601,
        -2.912257701703038, 0.14129347205162568, 4.024951056919511, 0.23159986026585244, 0.46396709981562356, 0.8924318940048726,
    ];
    check("reaction", &reactions(&analysis, &ids), &r, 4.1);
    #[rustfmt::skip]
    let stresses = [
        [-1.80855385059053, 0.23958202870492318, -0.5865540146456433, -3.1185429322185767, -0.1705207998253595, 1.52975528028147, 0.0, 0.0],
        [-2.8713173084320607, 0.26941733144207836, -0.9472147801478176, -3.364477023143204, -0.21127108011858053, 1.6314261177932294, 0.0, 0.0],
        [-3.181194301949967, -0.9700906426295406, 0.8955332094339158, -3.303498285713952, -0.4885915305077063, 1.2825703501507058, 0.0, 0.0],
        [0.6268500586104436, 1.419419397302228, -1.7079804732230155, -2.6876534877985727, 0.188300211150203, 1.6752693729004784, 0.0, 0.0],
    ];
    let responses = analysis
        .domain()
        .element_shell_responses(elements[0])
        .unwrap();
    assert_eq!(responses.len(), 4);
    for (g, (response, want)) in responses.iter().zip(&stresses).enumerate() {
        check(&format!("resultant {g}"), &response.resultant, want, 4.0);
    }
}

/// OpenSees' `ShellDKGT` integrates its self-weight with `w * xsj` instead of `0.5 * w * xsj` (the mass and
/// stiffness use the latter), applying twice the weight. Carapace applies the correct weight, so the OpenSees
/// run with `-selfWeight 0 0 -9.81` equals Carapace with twice the body acceleration `(0, 0, 2 * 9.81)`.
#[test]
fn opensees_doubles_the_triangle_self_weight() {
    let (analysis, ids, _) = patch(Some(2.0 * 9.81));
    #[rustfmt::skip]
    let u = [
        -0.014510281694702342, -0.005885996522735786, 0.04601619051137522, -0.015591272894839537, -0.048313806842992614, -0.011518125256484817,
        -0.006016647395908907, -0.004060488635656864, 0.023720329524298898, -0.03158459504137322, -0.03533928157637824, -0.01514212171137059,
        -0.06704918086467386, -0.028536861541990973, 0.21935423740512028, -0.03710214260875464, -0.11233972419343093, -0.027066158264137996,
    ];
    check("u", &displacements(&analysis, &ids, &[2, 3, 5]), &u, 0.22);
}
