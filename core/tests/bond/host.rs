//! Bonded bars in a continuum host (bond plan, WP2): the native builder `bond_polyline` against today's perfect-bond
//! `EmbeddedBar`. A strip of Quad4 elements is pulled along x with an elastic bar embedded off the nodes.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
    TangentStrategy,
};
use carapace_core::model::continuum::{bond_polyline, embed_polyline, BondSpec};
use carapace_core::model::{
    Domain, Element, ElementId, EmbeddedBar, Material, Node, NodeId, PlaneMaterial, Quad4,
};

const LENGTH: f64 = 200.0;
const HEIGHT: f64 = 40.0;
const NX: usize = 8;
const NY: usize = 2;
const BAR_AREA: f64 = 100.0;

struct Strip {
    domain: Domain,
    nodes: Vec<NodeId>,
    segments: Vec<EmbeddedBar>,
}

fn strip(polyline: &[[f64; 2]]) -> Strip {
    let mut domain: Domain = Domain::new();
    let mut nodes = Vec::new();
    for j in 0..=NY {
        for i in 0..=NX {
            let mut node =
                Node::new([LENGTH * i as f64 / NX as f64, HEIGHT * j as f64 / NY as f64]);
            if i == 0 {
                node = node.fix(0);
            }
            if i == 0 && j == 0 {
                node = node.fix(1);
            }
            nodes.push(domain.add_node(node));
        }
    }
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    for j in 0..NY {
        for i in 0..NX {
            let n = |a: usize, b: usize| nodes[b * (NX + 1) + a];
            domain.add_element(Element::Quad4(Quad4::new(
                [n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1)],
                100.0,
                material.clone(),
            )));
        }
    }
    for j in 0..=NY {
        domain.load_node(
            nodes[j * (NX + 1) + NX],
            0,
            if j == 1 { 2000.0 } else { 1000.0 },
        );
    }
    let bar = Material::Elastic { e: 200_000.0 };
    let segments = embed_polyline(&domain, polyline, BAR_AREA, &bar).unwrap();
    Strip {
        domain,
        nodes,
        segments,
    }
}

fn analysis(domain: Domain) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: 1e-6,
            moment_tol: 1e-6,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 30,
        })
        .build(domain)
}

fn solve(mut analysis: Analysis, load: f64) -> Analysis {
    // Reference load = 4 N total at factor 1; scale with the factor.
    let mut remaining = load;
    while remaining > 1e-12 {
        analysis.step().expect("converged");
        remaining -= 1.0;
    }
    analysis
}

fn bonded_spec(stiff: f64) -> BondSpec {
    BondSpec {
        perimeter: 31.4,
        fbd: 1e9,
        gb: stiff,
        rf: 1e-5,
        slip_limit: 10.0,
        lock_angle: 15f64.to_radians(),
        ends: [None, None],
    }
}

/// B2: with a bond stiffness far above the concrete's, the bonded bar reproduces the perfect-bond `EmbeddedBar`:
/// the same mesh displacements and the same bar force, for a bar along x, a bar inclined across several elements, and
/// a bar bent at a (faceted) vertex.
#[test]
fn b2_a_very_stiff_bond_reproduces_the_perfect_bond_bar() {
    let cases: Vec<(&str, Vec<[f64; 2]>)> = vec![
        ("along x", vec![[0.0, 15.0], [200.0, 15.0]]),
        ("inclined", vec![[5.0, 5.0], [195.0, 35.0]]),
        (
            "faceted bend",
            vec![[0.0, 12.0], [100.0, 12.0], [140.0, 20.0], [200.0, 20.0]],
        ),
    ];
    for (name, polyline) in cases {
        let plain = strip(&polyline);
        let mut domain = plain.domain.clone();
        let plain_ids: Vec<ElementId> = plain
            .segments
            .iter()
            .cloned()
            .map(|s| domain.add_element(Element::EmbeddedBar(s)))
            .collect();
        let reference = solve(analysis(domain), 1.0);

        let bonded = strip(&polyline);
        let mut domain = bonded.domain.clone();
        let spec = bonded_spec(1e9);
        let built = bond_polyline(
            &mut domain,
            &bonded.segments,
            BAR_AREA,
            &Material::Elastic { e: 200_000.0 },
            &spec,
        )
        .unwrap();
        let rod_ids: Vec<ElementId> = built
            .rods
            .into_iter()
            .map(|s| domain.add_element(Element::EmbeddedBar(s)))
            .collect();
        for spring in built.springs {
            domain.add_element(Element::ZeroLength(spring.element));
        }
        let result = solve(analysis(domain), 1.0);

        let scale = plain
            .nodes
            .iter()
            .map(|&n| reference.domain().node(n).displacement[0].abs())
            .fold(0.0, f64::max);
        for (k, &node) in plain.nodes.iter().enumerate() {
            let (a, b) = (
                reference.domain().node(node).displacement,
                result.domain().node(bonded.nodes[k]).displacement,
            );
            for dof in 0..2 {
                assert!(
                    (a[dof] - b[dof]).abs() <= 1e-5 * scale,
                    "{name}: node {k} dof {dof}: {} vs {}",
                    a[dof],
                    b[dof]
                );
            }
        }
        assert_eq!(
            plain_ids.len(),
            rod_ids.len(),
            "{name}: one rod segment per host segment"
        );
        for (p, r) in plain_ids.iter().zip(&rod_ids) {
            let (fp, fr) = (
                reference
                    .domain()
                    .element_embedded_bar_response(*p)
                    .unwrap()
                    .0,
                result.domain().element_embedded_bar_response(*r).unwrap().0,
            );
            assert!(
                (fp - fr).abs() <= 1e-4 * fp.abs().max(1e-6),
                "{name}: bar force {fp} vs {fr}"
            );
        }
    }
}

/// With a soft bond the bar carries less of the load than the perfect-bond bar (it slips), and the concrete more.
#[test]
fn a_soft_bond_lets_the_bar_slip() {
    let polyline = [[0.0, 15.0], [200.0, 15.0]];
    let force_at_mid = |gb: f64| {
        let model = strip(&polyline);
        let mut domain = model.domain.clone();
        let spec = bonded_spec(gb);
        let built = bond_polyline(
            &mut domain,
            &model.segments,
            BAR_AREA,
            &Material::Elastic { e: 200_000.0 },
            &spec,
        )
        .unwrap();
        let mid = built.rods.len() / 2;
        let rod_ids: Vec<ElementId> = built
            .rods
            .into_iter()
            .map(|s| domain.add_element(Element::EmbeddedBar(s)))
            .collect();
        for spring in built.springs {
            domain.add_element(Element::ZeroLength(spring.element));
        }
        let result = solve(analysis(domain), 1.0);
        result
            .domain()
            .element_embedded_bar_response(rod_ids[mid])
            .unwrap()
            .0
    };
    let (stiff, soft) = (force_at_mid(1e9), force_at_mid(1.0));
    assert!(soft < stiff, "{soft} vs {stiff}");
}
