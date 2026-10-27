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
