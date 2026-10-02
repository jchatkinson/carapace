use nalgebra::DVector;

use super::AnalysisError;

/// Convergence criterion for `Algorithm::Newton`/`KrylovNewton`. Closed enum
/// (§2.1). Unused by `Algorithm::Linear` (which never iterates).
///
/// Every variant is evaluated against the **accepted** state of an
/// iteration (`IterationContext`): the residual reassembled *after* the
/// iteration's correction (and any line-search rescaling) was applied, not
/// the residual that correction was solved against. A force test can
/// therefore never pass on a stale pre-correction residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConvergenceTest {
    /// Converged once the Euclidean norm of the post-update unbalanced
    /// force residual drops below `tol`.
    NormUnbalance { tol: f64, max_iter: usize },
    /// Converged once the displacement correction's norm drops below `tol`.
    /// Does not certify force equilibrium on its own (a stiff or stalled
    /// iteration can take a tiny step at a large residual).
    NormDispIncr { tol: f64, max_iter: usize },
    /// Converged once the (absolute) incremental work `|du . residual|`
    /// (accepted correction against the accepted residual) drops below
    /// `tol`. Pairs force with translation and moment with rotation
    /// consistently, but its value changes with the work unit and can be
    /// small through cancellation, so it does not certify force
    /// equilibrium on its own either.
    EnergyIncr { tol: f64, max_iter: usize },
    /// Every configured criterion must pass:
    ///
    /// - **Force equilibrium** (always): on each reduced equation `i`,
    ///   `|r_i| <= abs_i + relative_tol * F_ref_i`, where `abs_i` is
    ///   `force_tol` for translational equations and `moment_tol` for
    ///   rotational ones, and `F_ref_i` is the largest of the committed and
    ///   predicted external force and the committed internal force on that
    ///   equation, frozen at the start of the step. A per-equation maximum
    ///   is used rather than a norm so many quiet DOFs can't dilute one
    ///   unbalanced one, and the reference is never the initial or current
    ///   residual (either can be zero) or the current load alone (which
    ///   can pass through zero during continuation).
    /// - **Displacement correction** (when `displacement_tol` is set): the
    ///   largest component of the last correction is at most
    ///   `displacement_tol`. Under arc length this is the last *undamped*
    ///   correction, so a tiny accepted backtracking step can't fake it.
    Combined {
        force_tol: f64,
        moment_tol: f64,
        relative_tol: f64,
        displacement_tol: Option<f64>,
        max_iter: usize,
    },
}

/// What a convergence test sees at the end of an iteration: the residual
/// at the accepted (post-update) state, the correction actually applied to
/// reach it, and — for `ConvergenceTest::Combined` — the per-equation
/// force tolerance frozen at the start of the step.
pub(crate) struct IterationContext<'a> {
    pub(crate) residual: &'a DVector<f64>,
    pub(crate) du: &'a DVector<f64>,
    pub(crate) force_tolerance: Option<&'a ForceTolerance>,
}

impl ConvergenceTest {
    pub(crate) fn max_iter(&self) -> usize {
        match self {
            ConvergenceTest::NormUnbalance { max_iter, .. }
            | ConvergenceTest::NormDispIncr { max_iter, .. }
            | ConvergenceTest::EnergyIncr { max_iter, .. }
            | ConvergenceTest::Combined { max_iter, .. } => *max_iter,
        }
    }

    /// Rejects nonfinite or nonpositive tolerances and a zero iteration
    /// budget. `relative_tol` may be zero; the absolute tolerances may not.
    /// Every Newton-family step runs this; it is public so wire decoders
    /// can reject bad input up front.
    pub fn validate(&self) -> Result<(), AnalysisError> {
        let positive = |value: f64| value.is_finite() && value > 0.0;
        let invalid = |field| Err(AnalysisError::InvalidOption { field });
        match *self {
            ConvergenceTest::NormUnbalance { tol, max_iter }
            | ConvergenceTest::NormDispIncr { tol, max_iter }
            | ConvergenceTest::EnergyIncr { tol, max_iter } => {
                if !positive(tol) {
                    return invalid("convergence.tol");
                }
                if max_iter == 0 {
                    return invalid("convergence.maxIter");
                }
            }
            ConvergenceTest::Combined {
                force_tol,
                moment_tol,
                relative_tol,
                displacement_tol,
                max_iter,
            } => {
                if !positive(force_tol) {
                    return invalid("convergence.forceTol");
                }
                if !positive(moment_tol) {
                    return invalid("convergence.momentTol");
                }
                if !relative_tol.is_finite() || relative_tol < 0.0 {
                    return invalid("convergence.relativeTol");
                }
                if displacement_tol.is_some_and(|tol| !positive(tol)) {
                    return invalid("convergence.displacementTol");
                }
                if max_iter == 0 {
                    return invalid("convergence.maxIter");
                }
            }
        }
        Ok(())
    }

    /// Whether this test can certify a state reached without any
    /// correction (the step's predictor) — true for the force-based tests,
    /// which then treat the correction as zero. Displacement/work tests
    /// need a correction to measure, so they always take at least one.
    pub(crate) fn accepts_predictor(&self) -> bool {
        matches!(
            self,
            ConvergenceTest::NormUnbalance { .. } | ConvergenceTest::Combined { .. }
        )
    }

    /// Whether the iteration described by `ctx` has converged.
    pub(crate) fn check(&self, ctx: &IterationContext<'_>) -> bool {
        if !ctx.residual.iter().all(|value| value.is_finite())
            || !ctx.du.iter().all(|value| value.is_finite())
        {
            return false;
        }
        match self {
            ConvergenceTest::NormUnbalance { tol, .. } => ctx.residual.norm() < *tol,
            ConvergenceTest::NormDispIncr { tol, .. } => ctx.du.norm() < *tol,
            ConvergenceTest::EnergyIncr { tol, .. } => ctx.du.dot(ctx.residual).abs() < *tol,
            ConvergenceTest::Combined {
                displacement_tol, ..
            } => {
                let force_ok = ctx
                    .force_tolerance
                    .expect("Combined convergence always receives a force tolerance")
                    .scaled_max(ctx.residual)
                    <= 1.0;
                force_ok && displacement_tol.is_none_or(|tol| ctx.du.amax() <= tol)
            }
        }
    }
}

/// Per-equation force tolerance for `ConvergenceTest::Combined` (and the
/// force gate arc-length continuation applies under either force test),
/// frozen for one step. `thresholds[i]` is the largest admissible
/// `|r_i|`; dividing a residual by it gives a dimensionless measure where
/// `<= 1` means equilibrium on that equation.
#[derive(Debug, Clone)]
pub(crate) struct ForceTolerance {
    thresholds: DVector<f64>,
}

impl ForceTolerance {
    /// `rotational[i]` marks moment equations. `magnitudes` are the force
    /// vectors whose componentwise maximum forms `F_ref` (committed and
    /// predicted external load, committed internal force).
    pub(crate) fn combined(
        force_tol: f64,
        moment_tol: f64,
        relative_tol: f64,
        rotational: &[bool],
        magnitudes: &[&DVector<f64>],
    ) -> Self {
        let thresholds = DVector::from_fn(rotational.len(), |i, _| {
            let reference = magnitudes
                .iter()
                .map(|vector| vector[i].abs())
                .fold(0.0, f64::max);
            let absolute = if rotational[i] { moment_tol } else { force_tol };
            absolute + relative_tol * reference
        });
        Self { thresholds }
    }

    /// The tolerance `ConvergenceTest::NormUnbalance` implies for arc
    /// length's merit function: a uniform per-equation scale, so the
    /// scaled squared norm is `||r||^2 / tol^2`.
    pub(crate) fn uniform(n: usize, tol: f64) -> Self {
        Self {
            thresholds: DVector::from_element(n, tol),
        }
    }

    /// Builds the tolerance a step needs for `test`, or `None` when `test`
    /// doesn't use one.
    pub(crate) fn for_test(
        test: &ConvergenceTest,
        rotational: &[bool],
        magnitudes: &[&DVector<f64>],
    ) -> Option<Self> {
        match *test {
            ConvergenceTest::Combined {
                force_tol,
                moment_tol,
                relative_tol,
                ..
            } => Some(Self::combined(
                force_tol,
                moment_tol,
                relative_tol,
                rotational,
                magnitudes,
            )),
            _ => None,
        }
    }

    /// `max_i |r_i| / threshold_i`.
    pub(crate) fn scaled_max(&self, residual: &DVector<f64>) -> f64 {
        residual
            .iter()
            .zip(self.thresholds.iter())
            .map(|(r, t)| r.abs() / t)
            .fold(0.0, f64::max)
    }

    /// `sum_i (r_i / threshold_i)^2`.
    pub(crate) fn scaled_norm_squared(&self, residual: &DVector<f64>) -> f64 {
        residual
            .iter()
            .zip(self.thresholds.iter())
            .map(|(r, t)| (r / t).powi(2))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn combined(displacement_tol: Option<f64>) -> ConvergenceTest {
        ConvergenceTest::Combined {
            force_tol: 1e-3,
            moment_tol: 1e-1,
            relative_tol: 1e-6,
            displacement_tol,
            max_iter: 10,
        }
    }

    #[test]
    fn combined_uses_moment_tolerance_on_rotational_equations_only() {
        let load = DVector::from_column_slice(&[1000.0, 0.0]);
        let tol = ForceTolerance::combined(1e-3, 1e-1, 1e-6, &[false, true], &[&load]);
        let test = combined(None);
        let du = DVector::zeros(2);
        // Translational threshold 1e-3 + 1e-6 * 1000 = 2e-3; rotational 1e-1.
        let pass = DVector::from_column_slice(&[1.9e-3, 0.09]);
        let fail_force = DVector::from_column_slice(&[2.1e-3, 0.0]);
        let fail_moment = DVector::from_column_slice(&[0.0, 0.11]);
        let ctx = |r| IterationContext {
            residual: r,
            du: &du,
            force_tolerance: Some(&tol),
        };
        assert!(test.check(&ctx(&pass)));
        assert!(!test.check(&ctx(&fail_force)));
        assert!(!test.check(&ctx(&fail_moment)));
    }

    #[test]
    fn small_correction_never_hides_large_residual_under_combined() {
        let tol = ForceTolerance::combined(1e-6, 1e-6, 0.0, &[false], &[]);
        let test = combined(Some(1.0));
        let residual = DVector::from_column_slice(&[10.0]);
        let du = DVector::from_column_slice(&[1e-12]);
        assert!(!test.check(&IterationContext {
            residual: &residual,
            du: &du,
            force_tolerance: Some(&tol),
        }));
    }

    #[test]
    fn nonfinite_state_never_converges() {
        let test = ConvergenceTest::NormDispIncr {
            tol: 1.0,
            max_iter: 1,
        };
        let residual = DVector::from_column_slice(&[f64::NAN]);
        let du = DVector::from_column_slice(&[0.0]);
        assert!(!test.check(&IterationContext {
            residual: &residual,
            du: &du,
            force_tolerance: None,
        }));
    }

    #[test]
    fn validation_rejects_bad_combined_options() {
        assert!(combined(None).validate().is_ok());
        let bad = ConvergenceTest::Combined {
            force_tol: 0.0,
            moment_tol: 1.0,
            relative_tol: 0.0,
            displacement_tol: None,
            max_iter: 1,
        };
        assert_eq!(
            bad.validate(),
            Err(AnalysisError::InvalidOption {
                field: "convergence.forceTol"
            })
        );
    }
}
