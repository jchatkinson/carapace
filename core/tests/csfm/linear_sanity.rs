//! Linear sanity checks on the CSFM fixtures (plan `docs/csfm-plan.md` §4 E): these don't exercise
//! `CsfmConcrete`/`CsfmRebar`'s nonlinear branches at all, and are meant to catch a topology or
//! assembly bug (free-DOF bookkeeping, embedded-bar coupling, support elimination) independently of
//! whether the nonlinear materials are right.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use serde_json::Value;

use crate::model::build_linear;

fn load(name: &str) -> Value {
    let path = format!("{}/tests/csfm/fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: invalid JSON: {e}"))
}

/// Ported from CalcsApp's `expectLinearDeepBeamToBalance` (`test/fem/csfm-nonlinear.test.ts`):
/// "proves free residual accounting on the deep-beam topology before testing material
/// nonlinearity." A genuinely linear system (isotropic concrete, plain elastic rods through the
/// same embedded-bar topology as the real deep-beam fixture) must balance its free-DOF residual in
/// at most 2 Newton iterations — one to apply the first (exact, for a linear system) correction,
/// one to confirm the residual is now zero. More than that means the embedded-bar coupling or the
/// free/fixed DOF bookkeeping is wrong, independently of any CSFM material question.
#[test]
fn one_newton_correction_balances_the_linear_deep_beam() {
    let fixture = load("models/deep-beam.json");
    let built = build_linear(&fixture["model"]);

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.05 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-10,
            max_iter: 5,
        })
        .try_build(built.domain)
        .expect("the linear deep-beam topology is a valid model");

    let result = analysis
        .step()
        .expect("a linear system must balance well within the iteration budget");
    assert!(
        result.iterations <= 2,
        "expected at most 2 iterations (one correction, one confirmation), got {}",
        result.iterations
    );
}

/// `|a - b| <= abs + rel * max(|a|, |b|)`.
fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

/// `linear-elastic.json` (TS `solveLinearElastic` on the deep-beam topology, plan §4 E) against
/// Carapace's exact-support `Algorithm::Linear` solve of the equivalent native model. TS imposes
/// its boundary conditions by penalty; `summary.json`'s measured deep-beam penalty sensitivity
/// (`penalty1e11.displacement`, about 4.2e-12 relative) bounds how much of any difference here is
/// attributable to that penalty rather than to a real discrepancy, so the tolerance here is set
/// comfortably — several orders of magnitude — above it rather than reusing it exactly, since this
/// is a different comparison (two different solves, not two penalty scales of the same one).
#[test]
fn linear_elastic_displacements_match_ts_within_the_measured_penalty_sensitivity() {
    let fixture = load("linear-elastic.json");
    let built = build_linear(&fixture["model"]);
    let nodes = built.nodes;

    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 1,
        })
        .try_build(built.domain)
        .expect("the linear-elastic fixture's model is valid");
    analysis
        .step()
        .expect("a one-shot linear solve never fails");

    let expected = fixture["displacements"].as_array().unwrap();
    let max_u = expected
        .iter()
        .map(|v| v.as_f64().unwrap().abs())
        .fold(0.0_f64, f64::max);

    // Several orders of magnitude above the measured ~4.2e-12 penalty sensitivity (see doc
    // comment) - this is checking real agreement, not just "within penalty noise."
    let (abs, rel) = (1e-12 * max_u.max(1.0), 1e-6);
    for (n, &node) in nodes.iter().enumerate() {
        for dof in 0..2 {
            let got = analysis.domain().node(node).displacement[dof];
            let want = expected[2 * n + dof].as_f64().unwrap();
            assert!(
                close(got, want, abs, rel),
                "node {n} dof {dof}: got {got}, TS {want}"
            );
        }
    }
}
