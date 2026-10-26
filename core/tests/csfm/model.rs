//! Builds a native `Domain` from one of the CSFM golden fixture's `model` objects
//! (`core/tests/csfm/fixtures/models/*.json`), for the native end-to-end runs in `nonlinear.rs`
//! (plan `docs/csfm-plan.md` §4 F/E/H). Mirrors `core/tests/continuum/embed_polyline.rs`'s
//! `domain_from_mesh` (mesh/element construction) and extends it with the CSFM materials, exact
//! supports, loads and embedded rebar the model fixtures carry.
//!
//! Deliberately reads the fixture as a generic `serde_json::Value` rather than through typed
//! structs — the fixture's own schema is CalcsApp's export, not a Carapace wire format, and there
//! is no reuse to be had from strongly typing a read-only, test-only consumer of it.

use carapace_core::model::continuum::embed_polyline;
use carapace_core::model::{
    CsfmConcreteLaw, Domain, Element, ElementId, Material, Node, NodeId, PlaneMaterial, Quad4, Tri3,
};
use serde_json::Value;

/// One `rebarPaths[i]`'s embedded-bar elements, in the fixture's segment (`rodId`) order, for
/// matching `locations` entries tagged `{ kind: "rod", pathId, rodId }`.
pub struct RodPath {
    pub path_id: String,
    pub bars: Vec<ElementId>,
}

pub struct BuiltModel {
    pub domain: Domain,
    /// Node-table order, matching the fixture's `[ux0, uy0, ux1, uy1, ...]` displacement layout.
    pub nodes: Vec<NodeId>,
    pub rod_paths: Vec<RodPath>,
}

fn f64_of(v: &Value) -> f64 {
    v.as_f64()
        .unwrap_or_else(|| panic!("expected a number, got {v}"))
}

fn usize_of(v: &Value) -> usize {
    v.as_u64()
        .unwrap_or_else(|| panic!("expected a non-negative integer, got {v}")) as usize
}

/// Builds the `Domain` for one `model` object: mesh and `CsfmConcrete` elements, exact fixed-DOF
/// supports (plan §2.4 — Carapace has no penalty mode to opt into), nodal loads, and `CsfmRebar`
/// embedded bars via `embed_polyline`. Panics (via the fixture's own invariants) rather than
/// returning `Result` — a fixture that doesn't build is a test bug, not a case under test.
pub fn build(model: &Value) -> BuiltModel {
    assert!(
        model["nodeSprings"].as_array().unwrap().is_empty(),
        "no current CSFM fixture uses node springs; add support to model.rs before relying on this"
    );

    let mesh = &model["mesh"];
    let node_x = mesh["nodeX"].as_array().unwrap();
    let node_y = mesh["nodeY"].as_array().unwrap();
    let thickness = f64_of(&mesh["thickness"]);

    let mut fixed: Vec<Vec<usize>> = vec![Vec::new(); node_x.len()];
    for bc in model["boundaryConditions"].as_array().unwrap() {
        let value = f64_of(&bc["value"]);
        assert_eq!(
            value, 0.0,
            "only zero-displacement supports are modeled here (prescribed-displacement fixtures need step_prescribed)"
        );
        fixed[usize_of(&bc["node"])].push(usize_of(&bc["dof"]));
    }

    let mut domain = Domain::new();
    let node_ids: Vec<NodeId> = (0..node_x.len())
        .map(|i| {
            let mut node = Node::new([f64_of(&node_x[i]), f64_of(&node_y[i])]);
            for &dof in &fixed[i] {
                node = node.fix(dof);
            }
            domain.add_node(node)
        })
        .collect();

    let concrete = &model["concrete"];
    let law = match concrete["law"].as_str().unwrap() {
        "parabola-rectangle" => CsfmConcreteLaw::ParabolaRectangle,
        "elastic-ideally-plastic" => CsfmConcreteLaw::ElasticIdeallyPlastic,
        other => panic!("unknown concrete law {other}"),
    };
    let softening = concrete["softeningOn"].as_bool().unwrap();
    let gamma = concrete
        .get("tensileRegularization")
        .and_then(|r| r.get("gamma"))
        .map(f64_of);
    let fck = f64_of(&model["fck"]);
    let plane_material = PlaneMaterial::csfm_concrete(fck, None, law, softening, gamma)
        .expect("fixture's concrete parameters must be valid");
    let (fctm, ecm) = match &plane_material {
        PlaneMaterial::CsfmConcrete(fields) => (fields.fctm, fields.ecm),
        _ => unreachable!(),
    };

    let element_types = mesh["elementTypes"].as_array().unwrap();
    let element_offsets: Vec<usize> = mesh["elementOffsets"]
        .as_array()
        .unwrap()
        .iter()
        .map(usize_of)
        .collect();
    let element_nodes: Vec<usize> = mesh["elementNodes"]
        .as_array()
        .unwrap()
        .iter()
        .map(usize_of)
        .collect();

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
                    thickness,
                    plane_material.clone(),
                )));
            }
            "Q4" => {
                let n: [NodeId; 4] = std::array::from_fn(|i| node_ids[element_nodes[start + i]]);
                domain.add_element(Element::Quad4(Quad4::new(
                    n,
                    thickness,
                    plane_material.clone(),
                )));
            }
            other => panic!("unknown element type {other}"),
        }
    }

    for load in model["loads"].as_array().unwrap() {
        domain.load_node(
            node_ids[usize_of(&load["node"])],
            usize_of(&load["dof"]),
            f64_of(&load["value"]),
        );
    }

    let lambda = f64_of(&model["lambda"]);
    let tension_stiffening = model["tensionStiffening"].as_bool().unwrap();
    let rod_paths: Vec<RodPath> = model["rebarPaths"]
        .as_array()
        .unwrap()
        .iter()
        .map(|path| {
            let polyline: Vec<[f64; 2]> = path["polyline"]
                .as_array()
                .unwrap()
                .iter()
                .map(|p| [f64_of(&p["x"]), f64_of(&p["y"])])
                .collect();
            let bar_area = f64_of(&path["barAreaMm2"]) * f64_of(&path["barCount"]);
            let rho_eff = bar_area / f64_of(&path["effectiveAreaMm2"]);
            let material = Material::csfm_rebar(
                f64_of(&path["elasticModulusMpa"]),
                f64_of(&path["fyMpa"]),
                f64_of(&path["ftMpa"]),
                f64_of(&path["ultimateStrain"]),
                fctm,
                ecm,
                f64_of(&path["barDiameterMm"]),
                rho_eff,
                lambda,
                tension_stiffening,
            )
            .expect("fixture's rebar path parameters must be valid");
            let bars = embed_polyline(&domain, &polyline, bar_area, &material)
                .expect("fixture's rebar path must embed cleanly in its mesh");
            let ids: Vec<ElementId> = bars
                .into_iter()
                .map(|bar| domain.add_element(Element::EmbeddedBar(bar)))
                .collect();
            RodPath {
                path_id: path["id"].as_str().unwrap().to_string(),
                bars: ids,
            }
        })
        .collect();

    BuiltModel {
        domain,
        nodes: node_ids,
        rod_paths,
    }
}
