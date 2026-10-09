/// Post-processes a raw Newton correction to actually reduce the residual,
/// rather than accepting it outright. Closed enum. Ported from
/// Xara's `SRC/analysis/algorithm/equiSolnAlgo/search/` (`docs/
/// algorithms.md` §4): both variants bracket a sign change in the 1-D
/// "energy" `s(eta) = du . residual(eta)` by geometric expansion of `eta`,
/// then refine within the bracket — bisection, or false-position
/// (`RegulaFalsi`, usually fewer iterations than bisection once bracketed).
/// Only these two for now — `Secant`/`InitialInterpolatedLineSearch` can
/// diverge without a bracket and aren't needed until a specific model
/// demands the extra speed (`docs/algorithms.md` §4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineSearch {
    Bisection {
        tol: f64,
        max_iter: usize,
        max_eta: f64,
    },
    RegulaFalsi {
        tol: f64,
        max_iter: usize,
        max_eta: f64,
    },
}

impl LineSearch {
    fn params(&self) -> (f64, usize, f64) {
        match *self {
            LineSearch::Bisection {
                tol,
                max_iter,
                max_eta,
            }
            | LineSearch::RegulaFalsi {
                tol,
                max_iter,
                max_eta,
            } => (tol, max_iter, max_eta),
        }
    }

    /// Given `s0 = du . residual` at `eta = 0` and `s1 = du . residual` at
    /// `eta = 1` (the raw, un-rescaled correction), and a closure to
    /// evaluate `s(eta)` at any other trial `eta`, returns the `eta` this
    /// strategy converges to. `s0` and `s1` must have opposite signs
    /// (checked by the caller — `Analysis::try_step` only calls this once
    /// the ordinary convergence test at `eta = 1` has already failed, per
    /// `docs/algorithms.md` §4 step 3; a same-sign `s0`/`s1` means the raw
    /// step already reduced the energy monotonically, so bracketing
    /// expansion below is what finds the far side of the sign change).
    pub(crate) fn resolve(&self, s0: f64, s1: f64, mut eval: impl FnMut(f64) -> f64) -> f64 {
        let (tol, max_iter, max_eta) = self.params();

        // Bracket a sign change by geometric expansion, starting from the
        // raw step's own [0, 1] interval (Xara's `BisectionLineSearch`/
        // `RegulaFalsiLineSearch` both start this way).
        let (mut eta_lo, mut s_lo) = (0.0, s0);
        let (mut eta_hi, mut s_hi) = (1.0, s1);
        let mut expand_iter = 0;
        while s_lo.signum() == s_hi.signum() && eta_hi < max_eta && expand_iter < max_iter {
            eta_lo = eta_hi;
            s_lo = s_hi;
            eta_hi = (eta_hi * 4.0).min(max_eta);
            s_hi = eval(eta_hi);
            expand_iter += 1;
        }
        if s_lo.signum() == s_hi.signum() {
            // Never bracketed (e.g. `s` monotonically one-signed out to
            // `max_eta`) — the best available point is whichever endpoint
            // has the smaller `|s|`, not a fabricated root.
            return if s_lo.abs() < s_hi.abs() {
                eta_lo
            } else {
                eta_hi
            };
        }

        for _ in 0..max_iter {
            let eta_trial = match self {
                LineSearch::Bisection { .. } => 0.5 * (eta_lo + eta_hi),
                // False position: linear-interpolate the root between the
                // bracket's two endpoints instead of always bisecting.
                LineSearch::RegulaFalsi { .. } => {
                    eta_lo + (eta_hi - eta_lo) * (-s_lo) / (s_hi - s_lo)
                }
            };
            let s_trial = eval(eta_trial);
            if s_trial.abs() < tol {
                return eta_trial;
            }
            if s_trial.signum() == s_lo.signum() {
                eta_lo = eta_trial;
                s_lo = s_trial;
            } else {
                eta_hi = eta_trial;
                s_hi = s_trial;
            }
        }
        0.5 * (eta_lo + eta_hi)
    }
}
