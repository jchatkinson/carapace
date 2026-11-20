//! `embed_polyline` (plan `docs/csfm-plan.md` §2.3, §4 D): segment counts
//! against the CalcsApp golden fixtures (crossing, along/shared edges, a
//! T3 diagonal and through-node path, a distorted Q4, the deep beam and
//! span beam), zero-length vertices dropped, an outside-mesh error, and
//! an ambiguous-host error on deliberately overlapping elements.

use carapace_core::model::continuum::{embed_polyline, embed_polyline_in};
use carapace_core::model::{
    Domain, Element, Material, ModelError, Node, NodeId, PlaneMaterial, Quad4, Tri3,
};
use serde_json::Value;

fn fixture() -> Value {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/csfm/fixtures/embedding.json"
    );
    let text = std::fs::read_to_string(path).expect("embedding.json fixture");
    serde_json::from_str(&text).expect("valid JSON")
}

fn domain_from_mesh(mesh: &Value) -> (Domain, Vec<NodeId>) {
    let mut domain: Domain = Domain::new();
    let node_x = mesh["nodeX"].as_array().unwrap();
    let node_y = mesh["nodeY"].as_array().unwrap();
    let node_ids: Vec<NodeId> = (0..node_x.len())
        .map(|i| {
            domain.add_node(Node::new([
                node_x[i].as_f64().unwrap(),
                node_y[i].as_f64().unwrap(),
            ]))
        })
        .collect();

    let element_types = mesh["elementTypes"].as_array().unwrap();
    let element_offsets: Vec<usize> = mesh["elementOffsets"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    let element_nodes: Vec<usize> = mesh["elementNodes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();

    for (e, kind) in element_types.iter().enumerate() {
        let start = element_offsets[e];
        match kind.as_str().unwrap() {
            "T3" => {
                let n = [
                    node_ids[element_nodes[start]],
                    node_ids[element_nodes[start + 1]],
                    node_ids[element_nodes[start + 2]],
                ];
                domain.add_element(Element::Tri3(Tri3::new(
                    n[0],
                    n[1],
                    n[2],
                    1.0,
                    material.clone(),
                )));
            }
            "Q4" => {
                let n: [NodeId; 4] = std::array::from_fn(|i| node_ids[element_nodes[start + i]]);
                domain.add_element(Element::Quad4(Quad4::new(n, 1.0, material.clone())));
            }
            other => panic!("unknown element type {other}"),
        }
    }
    (domain, node_ids)
}

fn polyline_from_path(path: &Value) -> Vec<[f64; 2]> {
    path["polyline"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| [p["x"].as_f64().unwrap(), p["y"].as_f64().unwrap()])
        .collect()
}

fn dummy_material() -> Material {
    Material::Elastic { e: 200_000.0 }
}

/// Segment counts match the fixtures (`rods.len()` per path), including
/// the deep beam's 36 total segments and the span beam's 30 (plan §4 D,
/// ported from `rod.test.ts`).
#[test]
fn segment_counts_match_the_fixtures() {
    let fixture = fixture();
    let mut deep_beam_total = 0;
    let mut span_beam_total = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        if case.get("expectError").and_then(Value::as_bool) == Some(true) {
            continue;
        }
        let (domain, _) = domain_from_mesh(&case["mesh"]);
        for resolved in case["resolved"].as_array().unwrap() {
            let path_id = resolved["id"].as_str().unwrap();
            let path = case["paths"]
                .as_array()
                .unwrap()
                .iter()
                .find(|p| p["id"] == path_id)
                .unwrap();
            let polyline = polyline_from_path(path);
            let expected_rods = resolved["rods"].as_array().unwrap().len();
            let bars = embed_polyline(&domain, &polyline, 1.0, &dummy_material())
                .unwrap_or_else(|e| panic!("case {name} path {path_id}: {e}"));
            assert_eq!(
                bars.len(),
                expected_rods,
                "case {name} path {path_id}: got {} segments, fixture has {expected_rods}",
                bars.len()
            );
            if name == "deep-beam" {
                deep_beam_total += bars.len();
            }
            if name == "span-beam" {
                span_beam_total += bars.len();
            }
        }
    }
    assert_eq!(deep_beam_total, 36, "deep beam total segments");
    assert_eq!(span_beam_total, 30, "span beam total segments");
}

/// A point outside the mesh is an error.
#[test]
fn outside_mesh_is_an_error() {
    let fixture = fixture();
    let case = fixture["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["name"] == "outside-mesh")
        .unwrap();
    let (domain, _) = domain_from_mesh(&case["mesh"]);
    let path = &case["paths"].as_array().unwrap()[0];
    let polyline = polyline_from_path(path);
    let result = embed_polyline(&domain, &polyline, 1.0, &dummy_material());
    assert!(matches!(
        result,
        Err(ModelError::EmbeddedBarOutsideMesh { .. })
    ));
}

/// Zero-length (duplicate-vertex) spans are dropped rather than producing
/// a degenerate segment.
#[test]
fn zero_length_vertices_are_dropped() {
    let mut domain: Domain = Domain::new();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    let n0 = domain.add_node(Node::new([0.0, 0.0]));
    let n1 = domain.add_node(Node::new([2.0, 0.0]));
    let n2 = domain.add_node(Node::new([2.0, 2.0]));
    let n3 = domain.add_node(Node::new([0.0, 2.0]));
    domain.add_element(Element::Quad4(Quad4::new([n0, n1, n2, n3], 1.0, material)));

    let clean = embed_polyline(&domain, &[[0.2, 0.2], [1.8, 1.8]], 1.0, &dummy_material()).unwrap();
    let with_duplicate = embed_polyline(
        &domain,
        &[[0.2, 0.2], [0.2, 0.2], [1.8, 1.8]],
        1.0,
        &dummy_material(),
    )
    .unwrap();
    assert_eq!(clean.len(), with_duplicate.len());
}

/// An empty or single-vertex polyline is a degenerate-path error.
#[test]
fn empty_or_single_vertex_polyline_is_an_error() {
    let mut domain: Domain = Domain::new();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    let n0 = domain.add_node(Node::new([0.0, 0.0]));
    let n1 = domain.add_node(Node::new([2.0, 0.0]));
    let n2 = domain.add_node(Node::new([2.0, 2.0]));
    let n3 = domain.add_node(Node::new([0.0, 2.0]));
    domain.add_element(Element::Quad4(Quad4::new([n0, n1, n2, n3], 1.0, material)));

    assert!(matches!(
        embed_polyline(&domain, &[], 1.0, &dummy_material()),
        Err(ModelError::EmbeddedBarDegeneratePath)
    ));
    assert!(matches!(
        embed_polyline(&domain, &[[1.0, 1.0]], 1.0, &dummy_material()),
        Err(ModelError::EmbeddedBarDegeneratePath)
    ));
    assert!(matches!(
        embed_polyline(&domain, &[[1.0, 1.0], [1.0, 1.0]], 1.0, &dummy_material()),
        Err(ModelError::EmbeddedBarDegeneratePath)
    ));
}

/// Two exactly overlapping quads make a segment's midpoint ambiguous.
#[test]
fn overlapping_hosts_are_an_ambiguous_host_error() {
    let mut domain: Domain = Domain::new();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    let n0 = domain.add_node(Node::new([0.0, 0.0]));
    let n1 = domain.add_node(Node::new([2.0, 0.0]));
    let n2 = domain.add_node(Node::new([2.0, 2.0]));
    let n3 = domain.add_node(Node::new([0.0, 2.0]));
    domain.add_element(Element::Quad4(Quad4::new(
        [n0, n1, n2, n3],
        1.0,
        material.clone(),
    )));
    let m0 = domain.add_node(Node::new([0.0, 0.0]));
    let m1 = domain.add_node(Node::new([2.0, 0.0]));
    let m2 = domain.add_node(Node::new([2.0, 2.0]));
    let m3 = domain.add_node(Node::new([0.0, 2.0]));
    domain.add_element(Element::Quad4(Quad4::new([m0, m1, m2, m3], 1.0, material)));

    let result = embed_polyline(&domain, &[[0.2, 0.2], [1.8, 1.8]], 1.0, &dummy_material());
    assert!(matches!(
        result,
        Err(ModelError::EmbeddedBarAmbiguousHost { .. })
    ));
}

/// Two adjacent unit-height quads; returns the domain and both element ids.
fn two_quads() -> (Domain, [carapace_core::model::ElementId; 2]) {
    let mut domain: Domain = Domain::new();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    let n: Vec<NodeId> = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [2.0, 1.0],
    ]
    .iter()
    .map(|&c| domain.add_node(Node::new(c)))
    .collect();
    let q0 = domain.add_element(Element::Quad4(Quad4::new(
        [n[0], n[1], n[4], n[3]],
        1.0,
        material.clone(),
    )));
    let q1 = domain.add_element(Element::Quad4(Quad4::new(
        [n[1], n[2], n[5], n[4]],
        1.0,
        material,
    )));
    (domain, [q0, q1])
}

/// `hosts` restricts both the edge-crossing split and the host search: a path that runs into the
/// excluded element is outside the mesh, and one confined to the allowed element embeds in it.
#[test]
fn restricting_hosts_excludes_elements_from_the_search() {
    let (domain, [q0, q1]) = two_quads();
    let path = [[0.2, 0.5], [1.8, 0.5]];

    let all = embed_polyline_in(&domain, &path, 1.0, &dummy_material(), None).unwrap();
    assert_eq!(all.len(), 2, "split at the shared edge");

    let only_q0 = embed_polyline_in(&domain, &path, 1.0, &dummy_material(), Some(&[q0]));
    assert!(matches!(
        only_q0,
        Err(ModelError::EmbeddedBarOutsideMesh { .. })
    ));

    let inside_q1 = embed_polyline_in(
        &domain,
        &[[1.2, 0.5], [1.8, 0.5]],
        1.0,
        &dummy_material(),
        Some(&[q1]),
    )
    .unwrap();
    assert_eq!(inside_q1.len(), 1);
}

/// Overlapping elements are ambiguous unless `hosts` names only one of them.
#[test]
fn restricting_hosts_resolves_an_ambiguous_overlap() {
    let mut domain: Domain = Domain::new();
    let material = PlaneMaterial::plane_stress(30_000.0, 0.2).unwrap();
    let mut ids = Vec::new();
    for _ in 0..2 {
        let n: Vec<NodeId> = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]
            .iter()
            .map(|&c| domain.add_node(Node::new(c)))
            .collect();
        ids.push(domain.add_element(Element::Quad4(Quad4::new(
            [n[0], n[1], n[2], n[3]],
            1.0,
            material.clone(),
        ))));
    }
    let path = [[0.2, 0.2], [1.8, 1.8]];
    assert!(embed_polyline_in(&domain, &path, 1.0, &dummy_material(), None).is_err());
    let bars = embed_polyline_in(&domain, &path, 1.0, &dummy_material(), Some(&[ids[1]])).unwrap();
    assert_eq!(bars.len(), 1);
}

/// Segments remember their end points, in drawing order, for the wire's segment table.
#[test]
fn segments_report_their_end_points() {
    let (domain, _) = two_quads();
    let bars = embed_polyline(&domain, &[[0.2, 0.5], [1.8, 0.5]], 1.0, &dummy_material()).unwrap();
    assert_eq!(bars[0].end_points(), Some([[0.2, 0.5], [1.0, 0.5]]));
    assert_eq!(bars[1].end_points(), Some([[1.0, 0.5], [1.8, 0.5]]));
}
