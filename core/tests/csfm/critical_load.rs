//! `findCriticalLoad` parity (plan `docs/csfm-plan.md` §4 G, third bullet):
//! Carapace's own bracketing, via the generic `material_limit` stop
//! criterion generalized to the concrete tensile-strain `LimitKind`,
//! reproduces TS's `findCriticalLoad` bisection on the deep-beam fixture
//! (`core/tests/csfm/fixtures/critical-load.json`, phase 0) within its own
//! `loadFractionTolerance`. This is the mechanism plan §2.5 point 3 means
//! by "this replaces `findCriticalLoad`'s full re-solve bisection": TS
//! re-solves from zero at each trial fraction (a genuine bisection over
//! load fraction), while Carapace retakes one step with an interpolated
//! radius from the last *accepted* state (`docs/arclength.md` §7's
//! existing bracket-and-retake). These are different procedures that
//! should still agree, because `CsfmConcrete`/`CsfmRebar` are stateless
//! and nonlinear-elastic (`commit` is a no-op for both): the equilibrium
//! state at a given total load fraction is a function of the load alone
//! on this model's single smooth pre-peak branch, not of how it got
//! there, so TS's from-zero resolve and Carapace's incremental retake
//! land on the same point.
//!
//! **Fixture note**: `exportFindCriticalLoad`'s `model` object omits
//! `loads` (`findCriticalLoad` reads load through `resolved.model.loadMap`
//! internally in TS, not a static list) — not a defect, just a narrower
//! export than `models/*.json`'s. Its mesh, boundary conditions,
//! `rebarPaths`, `fck`, `concrete` and `lambda`/`tensionStiffening` are
//! confirmed byte-identical to `models/deep-beam.json`'s (both come from
//! calcs.app's `deepBeam()`), so this test borrows that fixture's `loads`
//! array rather than hand-deriving one or widening the export script for
//! a single field.

use std::fs;

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ArcDirection, ArcLength, ArcScales, ArcStopCriteria,
    ConstraintHandler, ConvergenceTest, Integrator, MaterialLimitTarget, StopReason,
    TangentStrategy,
};
use carapace_core::model::{LimitKind, LimitStateConfig};
use serde_json::Value;

use crate::model::{build, BuiltModel};

fn fixture(name: &str) -> Value {
    let path = format!("{}/tests/csfm/fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: {e}"))
}

#[test]
fn arc_length_bracketing_matches_ts_find_critical_load() {
    let critical_load = fixture("critical-load.json");
    let deep_beam = fixture("models/deep-beam.json");

    let mut model = critical_load["model"].clone();
    assert!(
        model.get("loads").is_none(),
        "critical-load.json's model unexpectedly carries its own loads; drop the borrowed-loads fixture note and use them directly"
    );
    for key in [
        "mesh",
        "boundaryConditions",
        "rebarPaths",
        "fck",
        "lambda",
        "tensionStiffening",
    ] {
        assert_eq!(
            model[key], deep_beam["model"][key],
            "critical-load.json and models/deep-beam.json must share the same underlying deepBeam() model for {key} to be borrowed safely"
        );
    }
    model["loads"] = deep_beam["model"]["loads"].clone();

    let BuiltModel { domain, .. } = build(&model);

    let stop_criteria = &critical_load["stopCriteria"];
    let max_compressive = stop_criteria["maxConcreteCompressiveStrain"]
        .as_f64()
        .unwrap();
    let max_tensile = stop_criteria["maxConcreteTensileStrain"].as_f64().unwrap();
    let load_fraction_tolerance = critical_load["loadFractionTolerance"].as_f64().unwrap();
    let outcome = &critical_load["outcome"];
    assert_eq!(outcome["kind"], "criterion_reached");
    let expected_load_fraction = outcome["loadFraction"].as_f64().unwrap();
    assert_eq!(
        outcome["violatedAt"]["kind"], "concrete",
        "this test only generalizes the concrete tensile-strain criterion"
    );

    let limits = LimitStateConfig {
        concrete_crushing_limit: Some(max_compressive),
        concrete_tensile_limit: Some(max_tensile),
        rebar_strain_limit: None,
    };
    let mut arc = ArcLength::adaptive(0.5, 0.05, 8.0, ArcScales::Auto { load: 1.0 });
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
            force_tol: 1e-6,
            moment_tol: 1e-6,
            relative_tol: 1e-9,
            displacement_tol: None,
            max_iter: 30,
        })
        .try_build(domain)
        .unwrap_or_else(|e| panic!("invalid model: {e:?}"));

    let mut landed = None;
    for _ in 0..200 {
        let result = analysis
            .step()
            .unwrap_or_else(|e| panic!("arc-length step failed before landing: {e:?}"));
        let arc_info = result.arc.expect("ArcLength reports continuation info");
        if let Some(stop) = arc_info.stop {
            landed = Some((stop, result.load_factor));
            break;
        }
    }
    let (stop, load_fraction) = landed
        .unwrap_or_else(|| panic!("did not land on the material limit within the step budget"));

    assert_eq!(stop.reason, StopReason::MaterialLimit);
    assert!(stop.landed_exactly, "exact landing was requested");

    let (element, utilization, kind) = analysis
        .domain()
        .controlling_limit_utilization(&LimitStateConfig {
            concrete_crushing_limit: Some(max_compressive),
            concrete_tensile_limit: Some(max_tensile),
            rebar_strain_limit: None,
        })
        .expect("a limit state must be present at the landed step");
    assert_eq!(kind, LimitKind::ConcreteTensileStrain);
    assert!(
        (utilization - 1.0).abs() < 1e-6,
        "landed utilization {utilization}"
    );

    eprintln!(
        "carapace load_fraction={load_fraction} ts expected={expected_load_fraction} diff={} tol={load_fraction_tolerance}",
        (load_fraction - expected_load_fraction).abs()
    );
    assert!(
        (load_fraction - expected_load_fraction).abs() <= load_fraction_tolerance,
        "Carapace's bracketing ({load_fraction}) must agree with TS's findCriticalLoad outcome \
         ({expected_load_fraction}) within its own loadFractionTolerance ({load_fraction_tolerance}); \
         controlling element {element:?}"
    );
}
