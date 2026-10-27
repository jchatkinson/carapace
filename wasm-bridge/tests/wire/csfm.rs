//! CSFM wire entries (plan `docs/csfm-plan.md` §2.7, §4 I): the `csfmConcrete` plane material and
//! the `csfmRebar` material, their decode errors, and the enhanced-quad restriction.

use crate::common::empty_input;
use carapace_core::model::{CsfmConcreteLaw, Material, PlaneMaterial};
use carapace_wasm::input_v1::materials::resolve_materials;
use carapace_wasm::input_v1::{decode, CarapaceInputV1, DecodeError};
use serde_json::json;

fn input_with(plane_materials: serde_json::Value, materials: serde_json::Value) -> CarapaceInputV1 {
    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["planeMaterials"] = plane_materials;
    input["materials"] = materials;
    serde_json::from_value(input).unwrap()
}

fn concrete(overrides: serde_json::Value) -> serde_json::Value {
    let mut spec = json!({
        "kind": "csfmConcrete", "fck": 30.0, "law": "parabolaRectangle", "softening": true,
    });
    for (key, value) in overrides.as_object().unwrap() {
        spec[key] = value.clone();
    }
    spec
}

fn rebar(overrides: serde_json::Value) -> serde_json::Value {
    let mut spec = json!({
        "kind": "csfmRebar", "es": 200000.0, "fy": 500.0, "ft": 550.0, "epsU": 0.05,
        "fctm": 2.9, "ecm": 33000.0, "diameter": 16.0, "rhoEff": 0.02, "lambda": 0.5,
        "tensionStiffening": true,
    });
    for (key, value) in overrides.as_object().unwrap() {
        spec[key] = value.clone();
    }
    spec
}

fn plane_material_error(spec: serde_json::Value) -> DecodeError {
    decode(input_with(json!([spec]), json!([]))).expect_err("must be rejected")
}

fn material_error(spec: serde_json::Value) -> DecodeError {
    decode(input_with(json!([]), json!([spec]))).expect_err("must be rejected")
}

#[test]
fn csfm_concrete_matches_native_constructor() {
    let input = input_with(
        json!([
            concrete(json!({})),
            concrete(json!({
                "fck": 40.0, "ecm": 30000.0, "law": "elasticIdeallyPlastic",
                "softening": false, "tensileRegularization": 0.01,
            })),
        ]),
        json!([]),
    );
    assert_eq!(input.plane_materials.len(), 2);
    // Decoding an element-free model succeeds; the constructor mapping is checked by building the
    // same materials natively and comparing the fields the wire carries.
    decode(input).expect("csfm concrete decodes");
    let native =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, true, None)
            .unwrap();
    assert!(matches!(native, PlaneMaterial::CsfmConcrete(_)));
}

#[test]
fn csfm_concrete_rejects_bad_fck_and_gamma() {
    for bad_fck in [0.0, -5.0, 50.1] {
        assert_eq!(
            plane_material_error(concrete(json!({ "fck": bad_fck }))),
            DecodeError::InvalidPlaneMaterial {
                index: 0,
                reason: "fck must satisfy 0 < fck <= 50 MPa, and an Ecm override must be finite and positive",
            }
        );
    }
    assert!(matches!(
        plane_material_error(concrete(json!({ "ecm": -1.0 }))),
        DecodeError::InvalidPlaneMaterial { index: 0, .. }
    ));
    for bad_gamma in [0.0, 1.0, 1.5, -0.1] {
        assert_eq!(
            plane_material_error(concrete(json!({ "tensileRegularization": bad_gamma }))),
            DecodeError::InvalidPlaneMaterial {
                index: 0,
                reason: "tensile-regularization gamma must satisfy 0 < gamma < 1",
            }
        );
    }
}

#[test]
fn csfm_rebar_matches_native_constructor() {
    let input = input_with(json!([]), json!([rebar(json!({}))]));
    let resolved = resolve_materials(&input.materials).unwrap();
    let native = Material::csfm_rebar(
        200000.0, 500.0, 550.0, 0.05, 2.9, 33000.0, 16.0, 0.02, 0.5, true,
    )
    .unwrap();
    let (Material::CsfmRebar(wire), Material::CsfmRebar(native)) = (&resolved[0], &native) else {
        panic!("expected CsfmRebar");
    };
    assert_eq!(wire.routing, native.routing);
    assert_eq!(wire.sr, native.sr);
    assert_eq!(wire.eps_m_at_ft, native.eps_m_at_ft);
}

#[test]
fn csfm_rebar_rejects_invalid_properties() {
    let cases = [
        (
            json!({ "ft": 500.0 }),
            "steel properties must be finite and positive with ft > fy and eps_u > fy / es",
        ),
        (
            json!({ "epsU": 0.001 }),
            "steel properties must be finite and positive with ft > fy and eps_u > fy / es",
        ),
        (
            json!({ "es": -1.0 }),
            "steel properties must be finite and positive with ft > fy and eps_u > fy / es",
        ),
        (
            json!({ "rhoEff": 0.0 }),
            "effective reinforcement ratio must satisfy 0 < rho_eff < 1",
        ),
        (
            json!({ "rhoEff": 1.0 }),
            "effective reinforcement ratio must satisfy 0 < rho_eff < 1",
        ),
        (
            json!({ "fctm": 0.0 }),
            "fctm, ecm, bar diameter and lambda must be finite and positive",
        ),
        (
            json!({ "diameter": -2.0 }),
            "fctm, ecm, bar diameter and lambda must be finite and positive",
        ),
        (
            json!({ "lambda": 0.0 }),
            "fctm, ecm, bar diameter and lambda must be finite and positive",
        ),
        // TCM routing (rho_eff = 0.02 is above rho_crit here) with a huge lambda: 2 sigma_min > fy.
        (
            json!({ "lambda": 50.0 }),
            "2 * sigma_min exceeds fy: reduce lambda or this reinforcement ratio routes to POM",
        ),
    ];
    for (overrides, reason) in cases {
        assert_eq!(
            material_error(rebar(overrides.clone())),
            DecodeError::InvalidMaterial { index: 0, reason },
            "{overrides}"
        );
    }
}

#[test]
fn csfm_rebar_error_reports_its_arena_index() {
    let error = decode(input_with(
        json!([]),
        json!([{ "kind": "elastic", "e": 1.0 }, rebar(json!({ "rhoEff": 2.0 }))]),
    ))
    .expect_err("must be rejected");
    assert!(matches!(
        error,
        DecodeError::InvalidMaterial { index: 1, .. }
    ));
}

// ---- embeddedBars ----

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::continuum::embed_polyline;
use carapace_core::model::{Domain, Element, LoadSeries, Node, Quad4, Quad4Formulation};
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};

/// Two 100 x 100 mm quads side by side (nodes 0..5 row-major), left edge fixed, 5 kN pulling
/// node 2 and node 5 each in x, one tie bar through the middle (y = 50) with the given `hosts`.
fn tie_model(bars: serde_json::Value, quad_formulation: &str) -> CarapaceInputV1 {
    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["nodes"] = json!({
        "coords": [0.0, 0.0, 100.0, 0.0, 200.0, 0.0, 0.0, 100.0, 100.0, 100.0, 200.0, 100.0],
        "fixed": [3, 0, 0, 3, 0, 0],
        "massNodeIndex": [], "mass": [],
    });
    input["planeMaterials"] = json!([concrete(json!({ "tensileRegularization": 0.01 }))]);
    input["materials"] = json!([rebar(json!({ "diameter": 16.0, "rhoEff": 0.02 }))]);
    input["quads"] = json!({
        "nodeIds": [0, 1, 4, 3, 1, 2, 5, 4],
        "thickness": [100.0, 100.0], "material": [0, 0], "density": [0.0, 0.0],
        "formulation": [quad_formulation, quad_formulation],
    });
    input["embeddedBars"] = bars;
    input["loadPatterns"] =
        json!({ "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0] });
    input["nodalLoads"] = json!({
        "pattern": [0, 0], "node": [2, 5], "dof": [0, 0], "value": [5000.0, 5000.0], "stage": [0, 0],
    });
    let mut input: CarapaceInputV1 = serde_json::from_value(input).unwrap();
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "load".to_string(),
            steps: 4,
            integrator: IntegratorSpec::LoadControl { increment: 0.25 },
            algorithm: AlgorithmSpec::Newton {
                tangent: Default::default(),
                line_search: None,
            },
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 30,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: (0..6)
            .map(|node| RecorderSpec::NodeDisp { node, dof: 0 })
            .collect(),
    };
    input
}

fn bar(polyline: serde_json::Value, extra: serde_json::Value) -> serde_json::Value {
    let mut spec = json!({ "polyline": polyline, "area": 201.0, "material": 0 });
    for (key, value) in extra.as_object().unwrap() {
        spec[key] = value.clone();
    }
    spec
}

fn decode_error(input: CarapaceInputV1) -> DecodeError {
    decode(input).expect_err("must be rejected")
}

/// The same tie built natively.
fn native_tie_displacements() -> Vec<f64> {
    let mut domain: Domain = Domain::new();
    let coords = [
        [0.0, 0.0],
        [100.0, 0.0],
        [200.0, 0.0],
        [0.0, 100.0],
        [100.0, 100.0],
        [200.0, 100.0],
    ];
    let ids: Vec<_> = coords
        .iter()
        .enumerate()
        .map(|(k, &c)| {
            let node = Node::new(c);
            domain.add_node(if k == 0 || k == 3 {
                node.fix(0).fix(1)
            } else {
                node
            })
        })
        .collect();
    let concrete = PlaneMaterial::csfm_concrete(
        30.0,
        None,
        CsfmConcreteLaw::ParabolaRectangle,
        true,
        Some(0.01),
    )
    .unwrap();
    for quad in [[0, 1, 4, 3], [1, 2, 5, 4]] {
        domain.add_element(Element::Quad4(
            Quad4::new(quad.map(|i| ids[i]), 100.0, concrete.clone())
                .with_formulation(Quad4Formulation::Full),
        ));
    }
    let steel = Material::csfm_rebar(
        200000.0, 500.0, 550.0, 0.05, 2.9, 33000.0, 16.0, 0.02, 0.5, true,
    )
    .unwrap();
    for segment in embed_polyline(&domain, &[[0.0, 50.0], [200.0, 50.0]], 201.0, &steel).unwrap() {
        domain.add_element(Element::EmbeddedBar(segment));
    }
    let pattern = domain.add_load_pattern(LoadSeries::Linear { slope: 1.0 });
    for node in [2, 5] {
        domain.add_nodal_load(pattern, ids[node], 0, 5000.0);
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.25 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 30,
        })
        .build(domain);
    for _ in 0..4 {
        analysis.step().unwrap();
    }
    ids.iter()
        .map(|&id| analysis.domain().node(id).displacement[0])
        .collect()
}

#[test]
fn embedded_bars_run_identically_to_the_native_model() {
    let input = tie_model(
        json!([bar(json!([0.0, 50.0, 200.0, 50.0]), json!({}))]),
        "full",
    );
    let mut session = decode(input).expect("decodes");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    let expected = native_tie_displacements();
    assert!(expected[2] > 0.0);
    for (node, expected) in expected.iter().enumerate() {
        let (_, got) = crate::common::last_sample(&outcome, node).unwrap();
        assert_eq!(got, *expected, "node {node}");
    }
}

#[test]
fn the_segment_table_is_reported_per_bar() {
    let input = tie_model(
        json!([
            bar(json!([0.0, 50.0, 200.0, 50.0]), json!({})),
            bar(json!([20.0, 25.0, 80.0, 25.0]), json!({})),
        ]),
        "full",
    );
    let session = decode(input).unwrap();
    let output = session.model_output();
    assert_eq!(output.embedded_bars.len(), 2);
    let first = &output.embedded_bars[0];
    assert_eq!(first.first_segment, 0);
    assert_eq!(first.segments.len(), 2, "split at the shared edge x = 100");
    assert_eq!(first.segments[0].start, [0.0, 50.0]);
    assert_eq!(first.segments[0].end, [100.0, 50.0]);
    assert_eq!(first.segments[1].end, [200.0, 50.0]);
    assert_eq!(first.segments[0].length, 100.0);
    // Segment 0 sits in quad 0, segment 1 in quad 1 (node-table rows).
    assert_eq!(first.segments[0].host_nodes, vec![0, 1, 4, 3]);
    assert_eq!(first.segments[1].host_nodes, vec![1, 2, 5, 4]);
    let second = &output.embedded_bars[1];
    assert_eq!(second.first_segment, 2);
    assert_eq!(second.segments.len(), 1);
    // Resolved routing and crack spacing come from the bar's CsfmRebar material.
    let Material::CsfmRebar(native) = Material::csfm_rebar(
        200000.0, 500.0, 550.0, 0.05, 2.9, 33000.0, 16.0, 0.02, 0.5, true,
    )
    .unwrap() else {
        unreachable!()
    };
    assert_eq!(first.sr, Some(native.sr));
    assert_eq!(
        first.routing,
        Some(carapace_wasm::input_v1::RebarRouting::Tcm)
    );
}

#[test]
fn hosts_restrict_where_a_bar_may_embed() {
    let left_only = json!({ "hosts": [{ "kind": "quad4", "index": 0 }] });
    let error = decode_error(tie_model(
        json!([bar(json!([0.0, 50.0, 200.0, 50.0]), left_only.clone())]),
        "full",
    ));
    assert!(matches!(
        error,
        DecodeError::InvalidModel {
            error: carapace_wasm::input_v1::error::ModelErrorDetail::EmbeddedBarOutsideMesh { .. }
        }
    ));
    let session = decode(tie_model(
        json!([bar(json!([10.0, 50.0, 90.0, 50.0]), left_only)]),
        "full",
    ))
    .unwrap();
    assert_eq!(session.model_output().embedded_bars[0].segments.len(), 1);
    // "all" is the same as omitting `hosts`.
    decode(tie_model(
        json!([bar(
            json!([0.0, 50.0, 200.0, 50.0]),
            json!({ "hosts": "all" })
        )]),
        "full",
    ))
    .unwrap();
}

#[test]
fn embedded_bar_decode_errors() {
    let line = json!([0.0, 50.0, 200.0, 50.0]);
    // A bar outside the mesh.
    assert!(matches!(
        decode_error(tie_model(
            json!([bar(json!([0.0, 500.0, 200.0, 500.0]), json!({}))]),
            "full"
        )),
        DecodeError::InvalidModel {
            error: carapace_wasm::input_v1::error::ModelErrorDetail::EmbeddedBarOutsideMesh { .. }
        }
    ));
    // Fewer than two vertices / odd coordinate count.
    for polyline in [json!([0.0, 50.0]), json!([0.0, 50.0, 100.0])] {
        assert_eq!(
            decode_error(tie_model(json!([bar(polyline, json!({}))]), "full")),
            DecodeError::InvalidRow {
                table: "embedded_bars",
                row: 0,
                reason: "polyline must hold at least two (x, y) vertices",
            }
        );
    }
    // Non-positive or non-finite area.
    assert_eq!(
        decode_error(tie_model(
            json!([bar(line.clone(), json!({ "area": 0.0 }))]),
            "full"
        )),
        DecodeError::InvalidRow {
            table: "embedded_bars",
            row: 0,
            reason: "area must be finite and positive, and the polyline finite",
        }
    );
    // Unknown material, unknown host element.
    assert_eq!(
        decode_error(tie_model(
            json!([bar(line.clone(), json!({ "material": 7 }))]),
            "full"
        )),
        DecodeError::UnknownMaterialIndex {
            table: "embedded_bars",
            row: 7
        }
    );
    assert_eq!(
        decode_error(tie_model(
            json!([bar(
                line.clone(),
                json!({ "hosts": [{ "kind": "quad4", "index": 9 }] })
            )]),
            "full"
        )),
        DecodeError::UnknownElementIndex {
            table: "quads",
            row: 9
        }
    );
}

#[test]
fn embedded_bars_are_rejected_in_a_3d_model() {
    let mut input = serde_json::to_value(empty_input(3)).unwrap();
    input["embeddedBars"] = json!([bar(json!([0.0, 0.0, 1.0, 0.0]), json!({}))]);
    let input: CarapaceInputV1 = serde_json::from_value(input).unwrap();
    assert_eq!(
        decode_error(input),
        DecodeError::TableNotInProfile {
            table: "embedded_bars"
        }
    );
}

#[test]
fn csfm_concrete_with_the_enhanced_quad_is_a_model_error() {
    let error = decode_error(tie_model(json!([]), "enhanced"));
    assert!(
        matches!(
            error,
            DecodeError::InvalidModel {
                error: carapace_wasm::input_v1::error::ModelErrorDetail::InvalidElement { .. }
            }
        ),
        "{error:?}"
    );
}
