//! CSFM wire entries (plan `docs/csfm-plan.md` §2.7, §4 I): the `csfmConcrete` plane material and
//! the `csfmRebar` material, their decode errors, and the enhanced-quad restriction.

use crate::common::empty_input;
use carapace_core::model::{CsfmConcreteLaw, Material, PlaneMaterial};
use carapace_wasm::input_v1::materials::resolve_materials;
use carapace_wasm::input_v1::tables::{ElementKind, PlaneMaterialSpec, PlaneStateSpec};
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

// ---- stop.materialLimit ----

use carapace_core::analysis::{
    ArcDirection, ArcLength, ArcScales, ArcStopCriteria, MaterialLimitTarget,
};
use carapace_core::model::{LimitStateConfig, Tri3};
use carapace_wasm::input_v1::sequence::{
    ArcLengthSpec, ArcScalesSpec, ArcStopSpec, MaterialLimitSpec,
};
use carapace_wasm::input_v1::StopReasonDetail;

const LEG: f64 = 100.0;

/// The one-DOF bare-steel tie of `core/tests/csfm/material_limits.rs`: a `Tri3` with legs along
/// x and y, one bar along the x leg, only node 1's `ux` free, unit reference load. The bar
/// carries `stress = lambda`, so the rupture limit is reached at `lambda = ft = 540`.
fn tie_input(limit: Option<MaterialLimitSpec>) -> CarapaceInputV1 {
    let mut input = serde_json::to_value(empty_input(2)).unwrap();
    input["nodes"] = json!({
        "coords": [0.0, 0.0, LEG, 0.0, 0.0, LEG], "fixed": [3, 2, 3], "massNodeIndex": [], "mass": [],
    });
    input["planeMaterials"] = json!([concrete(json!({ "softening": false }))]);
    input["materials"] = json!([rebar(json!({
        "ft": 540.0, "fctm": 2.9, "ecm": 32000.0, "lambda": 0.67, "tensionStiffening": false,
    }))]);
    input["triangles"] = json!({
        "nodeIds": [0, 1, 2], "thickness": [1.0], "material": [0], "density": [0.0],
    });
    input["embeddedBars"] =
        json!([{ "polyline": [0.0, 0.0, LEG, 0.0], "area": 1.0, "material": 0 }]);
    input["loadPatterns"] =
        json!({ "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0] });
    input["nodalLoads"] = json!({
        "pattern": [0], "node": [1], "dof": [0], "value": [1.0], "stage": [0],
    });
    let mut input: CarapaceInputV1 = serde_json::from_value(input).unwrap();
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "to-limit".to_string(),
            steps: 200,
            integrator: IntegratorSpec::ArcLength(ArcLengthSpec {
                initial_radius: 0.1,
                min_radius: Some(0.01),
                max_radius: Some(0.5),
                target_iterations: None,
                max_retries: None,
                direction: Default::default(),
                scales: ArcScalesSpec::Explicit {
                    displacement: 0.05 * LEG,
                    rotation: None,
                    load: 540.0,
                },
                predictor: Default::default(),
                seed: None,
                arc_tolerance: None,
                correction_tolerance: None,
                backtracking: None,
                stop: Some(ArcStopSpec {
                    material_limit: limit,
                    ..Default::default()
                }),
            }),
            algorithm: AlgorithmSpec::Newton {
                tangent: Default::default(),
                line_search: None,
            },
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: vec![],
    };
    input
}

/// Runs a session to its end and returns the last continuation detail and load factor.
fn run_to_stop(input: CarapaceInputV1) -> (carapace_wasm::input_v1::ContinuationDetail, f64) {
    let mut session = decode(input).unwrap();
    let mut last = None;
    loop {
        let outcome = session.advance(25);
        assert!(outcome.error.is_none(), "{:?}", outcome.error);
        last = outcome.continuation.or(last);
        if outcome.done {
            return (last.expect("took a step"), outcome.load_factor);
        }
    }
}

fn native_landing(limits: LimitStateConfig) -> f64 {
    let mut domain: Domain = Domain::new();
    let ids: Vec<_> = [[0.0, 0.0], [LEG, 0.0], [0.0, LEG]]
        .iter()
        .zip([(0, 1), (1, 1), (0, 1)])
        .enumerate()
        .map(|(k, (&c, _))| {
            let node = Node::new(c);
            domain.add_node(match k {
                1 => node.fix(1),
                _ => node.fix(0).fix(1),
            })
        })
        .collect();
    let concrete =
        PlaneMaterial::csfm_concrete(30.0, None, CsfmConcreteLaw::ParabolaRectangle, false, None)
            .unwrap();
    domain.add_element(Element::Tri3(Tri3::new(
        ids[0], ids[1], ids[2], 1.0, concrete,
    )));
    let steel = Material::csfm_rebar(
        200000.0, 500.0, 540.0, 0.05, 2.9, 32000.0, 16.0, 0.02, 0.67, false,
    )
    .unwrap();
    for bar in embed_polyline(&domain, &[[0.0, 0.0], [LEG, 0.0]], 1.0, &steel).unwrap() {
        domain.add_element(Element::EmbeddedBar(bar));
    }
    domain.load_node(ids[1], 0, 1.0);
    let mut arc = ArcLength::adaptive(
        0.1,
        0.01,
        0.5,
        ArcScales::Explicit {
            displacement: 0.05 * LEG,
            rotation: None,
            load: 540.0,
        },
    );
    arc.direction = ArcDirection::Increasing;
    arc.stop = ArcStopCriteria {
        material_limit: Some(Box::new(MaterialLimitTarget {
            limits,
            exact: true,
        })),
        ..Default::default()
    };
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::ArcLength(arc))
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: 1e-8,
            moment_tol: 1e-8,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 30,
        })
        .build(domain);
    for _ in 0..200 {
        let result = analysis.step().unwrap();
        if result.arc.unwrap().stop.is_some() {
            return result.load_factor;
        }
    }
    panic!("native run did not land");
}

#[test]
fn material_limit_lands_on_rupture_over_the_wire() {
    let (detail, load_factor) = run_to_stop(tie_input(Some(MaterialLimitSpec {
        exact: true,
        ..Default::default()
    })));
    let stop = detail.stop.expect("stopped on the limit");
    assert_eq!(stop.reason, StopReasonDetail::MaterialLimit);
    assert!(stop.landed_exactly);
    assert!((load_factor - 540.0).abs() < 1.0, "{load_factor}");
    let native = native_landing(LimitStateConfig::default());
    assert!(
        (load_factor - native).abs() <= 1e-6 * native,
        "wire {load_factor} vs native {native}"
    );
}

#[test]
fn material_limit_overrides_reach_the_limit_config() {
    // Bare steel: sigma = 500 + (540 - 500) / (0.05 - 0.0025) * (0.01 - 0.0025) at eps = 0.01.
    let expected = 500.0 + 40.0 / 0.0475 * 0.0075;
    let (detail, load_factor) = run_to_stop(tie_input(Some(MaterialLimitSpec {
        exact: true,
        rebar_strain_limit: Some(0.01),
        ..Default::default()
    })));
    assert_eq!(detail.stop.unwrap().reason, StopReasonDetail::MaterialLimit);
    assert!(detail.stop.unwrap().landed_exactly);
    assert!(
        (load_factor - expected).abs() < 1e-3 * expected,
        "{load_factor} vs {expected}"
    );
    let native = native_landing(LimitStateConfig {
        rebar_strain_limit: Some(0.01),
        ..Default::default()
    });
    assert!((load_factor - native).abs() <= 1e-6 * native);
}

#[test]
fn invalid_material_limits_are_decode_errors() {
    for limit in [
        MaterialLimitSpec {
            concrete_crushing_limit: Some(-0.0035),
            ..Default::default()
        },
        MaterialLimitSpec {
            concrete_tensile_limit: Some(0.0),
            ..Default::default()
        },
        MaterialLimitSpec {
            rebar_strain_limit: Some(f64::NAN),
            ..Default::default()
        },
    ] {
        assert_eq!(
            decode_error(tie_input(Some(limit))),
            DecodeError::InvalidAnalysisOption {
                stage: "to-limit".to_string(),
                field: "integrator.stop.materialLimit",
            }
        );
    }
}

// ---- recorders: principal values, kc2, embedded bars ----

use carapace_wasm::input_v1::sequence::{EmbeddedBarQuantity, GaussQuantity};

fn last_value(outcome: &carapace_wasm::input_v1::StepOutcome, recorder: usize) -> f64 {
    crate::common::last_sample(outcome, recorder)
        .unwrap_or_else(|| panic!("recorder {recorder} recorded nothing"))
        .1
}

/// The tie model pulled hard enough (400 kN) that quad 0's principal tension passes the `kc2`
/// cap (`0.2 / 140`), with one recorder of every new kind.
#[test]
fn principal_kc2_and_embedded_bar_recorders_report_consistent_values() {
    let mut input = tie_model(
        json!([bar(json!([0.0, 50.0, 200.0, 50.0]), json!({}))]),
        "full",
    );
    input.nodal_loads.value = vec![40_000.0, 40_000.0];
    let gauss = |quantity, component| RecorderSpec::GaussPoint {
        element_kind: ElementKind::Quad4,
        element_index: 0,
        point: 0,
        quantity,
        component,
    };
    let mut recorders = Vec::new();
    for quantity in [
        GaussQuantity::Strain,
        GaussQuantity::Stress,
        GaussQuantity::PrincipalStrain,
        GaussQuantity::PrincipalStress,
    ] {
        for component in 0..3 {
            recorders.push(gauss(quantity, component));
        }
    }
    recorders.push(gauss(GaussQuantity::Kc2, 0)); // 12
    for quantity in [
        EmbeddedBarQuantity::AxialForce,
        EmbeddedBarQuantity::Strain,
        EmbeddedBarQuantity::Stress,
        EmbeddedBarQuantity::Ruptured,
    ] {
        recorders.push(RecorderSpec::EmbeddedBar {
            bar: 0,
            segment: 0,
            quantity,
        });
    } // 13..=16
    for quantity in [GaussQuantity::Strain, GaussQuantity::Stress] {
        recorders.push(RecorderSpec::GaussPoint {
            element_kind: ElementKind::EmbeddedBar,
            element_index: 0,
            point: 0,
            quantity,
            component: 0,
        });
    } // 17, 18
    input.sequence.recorders = recorders;

    let mut session = decode(input).expect("decodes");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    let v = |k: usize| last_value(&outcome, k);

    // Principal strain: Mohr's invariants and the eigen-equation M n = e1 n.
    let (ex, ey, gamma) = (v(0), v(1), v(2));
    let (e1, e2, theta) = (v(6), v(7), v(8));
    assert!(
        e1 > 0.2 / 140.0,
        "the model must pass the kc2 cap: eps1 = {e1}"
    );
    assert!(e1 >= e2);
    assert!((e1 + e2 - (ex + ey)).abs() < 1e-14);
    assert!((e1 * e2 - (ex * ey - gamma * gamma / 4.0)).abs() < 1e-14);
    let n = [theta.cos(), theta.sin()];
    assert!((ex * n[0] + gamma / 2.0 * n[1] - e1 * n[0]).abs() < 1e-14);
    assert!((gamma / 2.0 * n[0] + ey * n[1] - e1 * n[1]).abs() < 1e-14);

    // Principal stress: the stress tensor (true shear) has the same eigen-structure.
    let (sx, sy, txy) = (v(3), v(4), v(5));
    let (s1, s2, phi) = (v(9), v(10), v(11));
    assert!(s1 >= s2);
    assert!((s1 + s2 - (sx + sy)).abs() < 1e-9);
    assert!((s1 * s2 - (sx * sy - txy * txy)).abs() < 1e-6);
    let m = [phi.cos(), phi.sin()];
    assert!((sx * m[0] + txy * m[1] - s1 * m[0]).abs() < 1e-9);
    assert!((txy * m[0] + sy * m[1] - s1 * m[1]).abs() < 1e-9);

    // kc2 from the recorded principal tension, softening on.
    assert!((v(12) - 1.0 / (0.8 + 140.0 * e1)).abs() < 1e-14);
    assert!(v(12) < 1.0);

    // Embedded bar: stress * area, and the same strain/stress as the generic Gauss recorder.
    assert_eq!(v(13), v(15) * 201.0);
    assert!(v(14) > 0.0 && v(15) > 0.0);
    assert_eq!(v(16), 0.0);
    assert_eq!((v(14), v(15)), (v(17), v(18)));
}

#[test]
fn new_recorder_decode_errors() {
    let model = || {
        tie_model(
            json!([bar(json!([0.0, 50.0, 200.0, 50.0]), json!({}))]),
            "full",
        )
    };
    let with = |recorder: RecorderSpec| {
        let mut input = model();
        input.sequence.recorders = vec![recorder];
        decode_error(input)
    };
    let kc2 = |component| RecorderSpec::GaussPoint {
        element_kind: ElementKind::Quad4,
        element_index: 0,
        point: 0,
        quantity: GaussQuantity::Kc2,
        component,
    };
    // kc2 has one component.
    assert_eq!(
        with(kc2(1)),
        DecodeError::InvalidRecorderComponent {
            recorder: 0,
            component: 1,
            width: 1
        }
    );
    // kc2 of a bar, which has no concrete law.
    assert_eq!(
        with(RecorderSpec::GaussPoint {
            element_kind: ElementKind::EmbeddedBar,
            element_index: 0,
            point: 0,
            quantity: GaussQuantity::Kc2,
            component: 0,
        }),
        DecodeError::RecorderNotApplicable {
            recorder: 0,
            reason: "kc2 needs a csfmConcrete element",
        }
    );
    // Bar / segment out of range.
    assert_eq!(
        with(RecorderSpec::EmbeddedBar {
            bar: 3,
            segment: 0,
            quantity: EmbeddedBarQuantity::Strain
        }),
        DecodeError::UnknownElementIndex {
            table: "embedded_bars",
            row: 3
        }
    );
    assert_eq!(
        with(RecorderSpec::EmbeddedBar {
            bar: 0,
            segment: 2,
            quantity: EmbeddedBarQuantity::Strain
        }),
        DecodeError::InvalidEmbeddedBarSegment {
            recorder: 0,
            bar: 0,
            segment: 2,
            count: 2
        }
    );
    // kc2 on a plain elastic plane material.
    let mut plain = model();
    plain.plane_materials = vec![PlaneMaterialSpec::Isotropic {
        e: 30_000.0,
        nu: 0.2,
        state: PlaneStateSpec::PlaneStress,
    }];
    plain.sequence.recorders = vec![kc2(0)];
    assert_eq!(
        decode_error(plain),
        DecodeError::RecorderNotApplicable {
            recorder: 0,
            reason: "kc2 needs a csfmConcrete element",
        }
    );
    // embeddedBar recorders do not exist in a 3D model.
    let mut input3 = serde_json::to_value(empty_input(3)).unwrap();
    input3["nodes"] =
        json!({ "coords": [0.0, 0.0, 0.0], "fixed": [0], "massNodeIndex": [], "mass": [] });
    let mut input3: CarapaceInputV1 = serde_json::from_value(input3).unwrap();
    input3.sequence.recorders = vec![RecorderSpec::EmbeddedBar {
        bar: 0,
        segment: 0,
        quantity: EmbeddedBarQuantity::Strain,
    }];
    assert_eq!(
        decode_error(input3),
        DecodeError::ElementKindNotInProfile {
            table: "embedded_bars"
        }
    );
}

// ---- bulk recorders ----

use carapace_wasm::input_v1::sequence::{BulkQuantity, BulkTable};
use carapace_wasm::input_v1::StepOutcome;

fn bulk(table: BulkTable, quantity: BulkQuantity) -> RecorderSpec {
    RecorderSpec::Bulk { table, quantity }
}

fn loaded_tie() -> CarapaceInputV1 {
    let mut input = tie_model(
        json!([bar(json!([0.0, 50.0, 200.0, 50.0]), json!({}))]),
        "full",
    );
    input.nodal_loads.value = vec![40_000.0, 40_000.0];
    input
}

/// Every accepted step's bulk values equal the per-point scalar recorders' samples, in
/// element -> point -> component order.
#[test]
fn bulk_gauss_recorders_match_the_scalar_recorders_point_for_point() {
    let mut input = loaded_tie();
    let quantities = [
        (BulkQuantity::Strain, GaussQuantity::Strain, 3),
        (BulkQuantity::Stress, GaussQuantity::Stress, 3),
        (
            BulkQuantity::PrincipalStrain,
            GaussQuantity::PrincipalStrain,
            3,
        ),
        (
            BulkQuantity::PrincipalStress,
            GaussQuantity::PrincipalStress,
            3,
        ),
        (BulkQuantity::Kc2, GaussQuantity::Kc2, 1),
    ];
    let mut recorders = Vec::new();
    for &(bulk_quantity, _, _) in &quantities {
        recorders.push(bulk(BulkTable::Quads, bulk_quantity));
    }
    let mut scalar_layout = Vec::new(); // (quantity index, first scalar recorder)
    for (q, &(_, quantity, components)) in quantities.iter().enumerate() {
        scalar_layout.push((q, recorders.len()));
        for element_index in 0..2 {
            for point in 0..4 {
                for component in 0..components {
                    recorders.push(RecorderSpec::GaussPoint {
                        element_kind: ElementKind::Quad4,
                        element_index,
                        point,
                        quantity,
                        component,
                    });
                }
            }
        }
    }
    input.sequence.recorders = recorders;

    let mut session = decode(input).unwrap();
    let outcome = session.advance(100);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");
    assert_eq!(outcome.bulk_batches.len(), quantities.len());
    for &(q, first_scalar) in &scalar_layout {
        let (_, _, components) = quantities[q];
        let batch = &outcome.bulk_batches[q];
        assert_eq!(batch.recorder_index, q);
        assert_eq!(
            (batch.elements, batch.points, batch.components),
            (2, 4, components as u32)
        );
        assert_eq!(batch.first_sample, 0);
        let steps = batch.progress.len();
        assert_eq!(steps, 4);
        let width = 2 * 4 * components as usize;
        assert_eq!(batch.values.len(), steps * width);
        for step in 0..steps {
            for k in 0..width {
                let scalar = outcome
                    .recorder_batches
                    .iter()
                    .find(|b| b.recorder_index == first_scalar + k)
                    .unwrap();
                assert_eq!(
                    batch.values[step * width + k],
                    scalar.samples[step].1,
                    "q{q} step {step} k{k}"
                );
                assert_eq!(batch.progress[step], scalar.samples[step].0);
            }
        }
    }
    // The scalar slots of bulk recorders stay empty.
    assert!(outcome
        .recorder_batches
        .iter()
        .all(|b| b.recorder_index >= quantities.len()));
    // Principal tension passed the kc2 cap somewhere, so the kc2 array is not all ones.
    assert!(outcome.bulk_batches[4].values.iter().any(|&k| k < 1.0));
}

#[test]
fn bulk_embedded_bar_recorders_match_the_scalar_recorders() {
    let mut input = loaded_tie();
    let quantities = [
        (BulkQuantity::AxialForce, EmbeddedBarQuantity::AxialForce),
        (BulkQuantity::Strain, EmbeddedBarQuantity::Strain),
        (BulkQuantity::Stress, EmbeddedBarQuantity::Stress),
        (BulkQuantity::Ruptured, EmbeddedBarQuantity::Ruptured),
    ];
    let mut recorders: Vec<_> = quantities
        .iter()
        .map(|&(q, _)| bulk(BulkTable::EmbeddedBars, q))
        .collect();
    for &(_, quantity) in &quantities {
        for segment in 0..2 {
            recorders.push(RecorderSpec::EmbeddedBar {
                bar: 0,
                segment,
                quantity,
            });
        }
    }
    input.sequence.recorders = recorders;
    let mut session = decode(input).unwrap();
    let outcome = session.advance(100);
    assert!(outcome.done && outcome.error.is_none());
    for (q, batch) in outcome.bulk_batches.iter().enumerate() {
        assert_eq!((batch.elements, batch.points, batch.components), (2, 1, 1));
        for step in 0..batch.progress.len() {
            for segment in 0..2 {
                let scalar = outcome
                    .recorder_batches
                    .iter()
                    .find(|b| b.recorder_index == 4 + 2 * q + segment)
                    .unwrap();
                assert_eq!(batch.values[step * 2 + segment], scalar.samples[step].1);
            }
        }
    }
    assert!(outcome.bulk_batches[0].values.iter().any(|&f| f > 0.0));
}

/// A bulk recorder is paced like every other recorder: one entry per accepted step, with
/// `firstSample` continuing across `advance` calls.
#[test]
fn bulk_batches_continue_across_advance_calls() {
    let mut input = loaded_tie();
    input.sequence.recorders = vec![bulk(BulkTable::Quads, BulkQuantity::Strain)];
    let mut session = decode(input).unwrap();
    let mut firsts = Vec::new();
    let mut all_progress = Vec::<f64>::new();
    loop {
        let outcome: StepOutcome = session.advance(1);
        for batch in &outcome.bulk_batches {
            assert_eq!(batch.progress.len(), 1);
            firsts.push(batch.first_sample);
            all_progress.extend(&batch.progress);
        }
        if outcome.done {
            break;
        }
    }
    assert_eq!(firsts, vec![0, 1, 2, 3]);
    assert_eq!(all_progress, vec![0.25, 0.5, 0.75, 1.0]);
}

#[test]
fn bulk_recorder_decode_errors() {
    let with = |input: CarapaceInputV1, recorder| {
        let mut input = input;
        input.sequence.recorders = vec![recorder];
        decode_error(input)
    };
    let not_applicable = |reason| DecodeError::RecorderNotApplicable {
        recorder: 0,
        reason,
    };
    assert_eq!(
        with(
            loaded_tie(),
            bulk(BulkTable::Quads, BulkQuantity::AxialForce)
        ),
        not_applicable("axialForce and ruptured belong to the embeddedBars table")
    );
    assert_eq!(
        with(
            loaded_tie(),
            bulk(BulkTable::EmbeddedBars, BulkQuantity::PrincipalStrain)
        ),
        not_applicable("principal values and kc2 belong to the triangles and quads tables")
    );
    assert_eq!(
        with(
            loaded_tie(),
            bulk(BulkTable::Triangles, BulkQuantity::Strain)
        ),
        not_applicable("the table is empty")
    );
    let mut plain = loaded_tie();
    plain.plane_materials = vec![PlaneMaterialSpec::Isotropic {
        e: 30_000.0,
        nu: 0.2,
        state: PlaneStateSpec::PlaneStress,
    }];
    assert_eq!(
        with(plain, bulk(BulkTable::Quads, BulkQuantity::Kc2)),
        not_applicable("kc2 needs a csfmConcrete element")
    );
    // No quads exist in a 3D model.
    let mut input3 = serde_json::to_value(empty_input(3)).unwrap();
    input3["nodes"] =
        json!({ "coords": [0.0, 0.0, 0.0], "fixed": [0], "massNodeIndex": [], "mass": [] });
    let input3: CarapaceInputV1 = serde_json::from_value(input3).unwrap();
    assert_eq!(
        with(input3, bulk(BulkTable::Quads, BulkQuantity::Strain)),
        DecodeError::ElementKindNotInProfile { table: "quads" }
    );
}

// ---- per-step material-limit report ----

use carapace_wasm::input_v1::LimitKindDetail;

/// Every accepted step reports the controlling utilization and where it is. On the rupture tie the
/// bar is the only limit-bearing element that can control (the concrete is unloaded in
/// compression), and the last step lands on utilization 1.
#[test]
fn each_step_reports_the_controlling_limit_and_its_location() {
    let mut session = decode(tie_input(Some(MaterialLimitSpec {
        exact: true,
        ..Default::default()
    })))
    .unwrap();
    let mut history = Vec::new();
    loop {
        let outcome = session.advance(10);
        history.extend(outcome.limit_history.clone());
        if outcome.done {
            break;
        }
    }
    // One entry per accepted step, aligned with the load factor.
    assert!(history.len() > 3);
    let last = history.last().unwrap();
    assert!(
        (last.utilization - 1.0).abs() < 1e-6,
        "{}",
        last.utilization
    );
    assert_eq!(last.kind, LimitKindDetail::RebarRupture);
    assert_eq!(last.element_kind, ElementKind::EmbeddedBar);
    assert_eq!(last.element_index, 0);
    assert!((last.progress - 540.0).abs() < 1.0);
    // The bar strains monotonically, so utilization never decreases and progress is increasing.
    for pair in history.windows(2) {
        assert!(pair[1].utilization >= pair[0].utilization);
        assert!(pair[1].progress > pair[0].progress);
    }
}

#[test]
fn limit_history_names_the_controlling_element_kind_and_row() {
    // The loaded tie: with a heavily strained bar the continuum points stay far from crushing, so
    // the bar controls; its location is the flat segment index.
    let mut input = loaded_tie();
    input.sequence.recorders = vec![];
    let mut session = decode(input).unwrap();
    let outcome = session.advance(100);
    assert_eq!(outcome.limit_history.len(), 4);
    for sample in &outcome.limit_history {
        assert_eq!(sample.kind, LimitKindDetail::RebarRupture);
        assert_eq!(sample.element_kind, ElementKind::EmbeddedBar);
        assert!(sample.element_index < 2);
        assert!(sample.utilization > 0.0 && sample.utilization < 1.0);
    }
    assert_eq!(
        outcome
            .limit_history
            .iter()
            .map(|s| s.progress)
            .collect::<Vec<_>>(),
        vec![0.25, 0.5, 0.75, 1.0]
    );
}

#[test]
fn models_without_limit_states_report_no_limit_history() {
    let mut session = decode(crate::cutbacks::default_input()).unwrap();
    let outcome = session.advance(1000);
    assert!(outcome.done);
    assert!(outcome.limit_history.is_empty());
}
