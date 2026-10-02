use nalgebra::{DMatrix, DVector};

/// Krylov-subspace acceleration of a stale-tangent Newton correction —
/// ported from Xara's `KrylovAccelerator.cpp` (`docs/algorithms.md` §5). A
/// concrete struct, not a trait object (§2.1/§2.6): exactly one accelerator
/// implementation exists, same reasoning as `SparseSolver` committing to
/// one solver.
///
/// Keeps the last `max_dimension` pairs of (raw correction `v_i`, the
/// change in residual it produced, `av_i = r_{i-1} - r_i`). Given a new raw
/// correction `v_k` (solved against the stale tangent) and its residual
/// `r`, solves the small `(n x k)` least-squares problem `minimize || r -
/// Av * c ||` over `c`, then corrects `v_k += sum(c_i * v_i)` — using the
/// history of how past corrections actually moved the residual to
/// extrapolate a better correction than the stale tangent alone would
/// give. `n` (the full system size) can be large, but `k <= max_dimension`
/// is small (Xara defaults `maxDimension = 3`), so this is a tall-skinny
/// dense least-squares problem — `nalgebra`'s QR decomposition covers it
/// exactly, no need for a LAPACK binding.
pub(crate) struct KrylovAccelerator {
    max_dimension: usize,
    v: Vec<DVector<f64>>,
    av: Vec<DVector<f64>>,
    last_residual: Option<DVector<f64>>,
}

impl KrylovAccelerator {
    pub(crate) fn new(max_dimension: usize) -> Self {
        KrylovAccelerator {
            max_dimension,
            v: Vec::new(),
            av: Vec::new(),
            last_residual: None,
        }
    }

    /// Whether the subspace has grown past `max_dimension` — the caller's
    /// signal to re-form (and re-factor) the tangent and start a fresh
    /// `KrylovAccelerator`, rather than keep accelerating against an
    /// increasingly stale basis. A fresh `Analysis::step`/`TransientAnalysis::
    /// step` also starts from a new `KrylovAccelerator::new` (the driver
    /// constructs one locally per call — a new step's stale tangent makes
    /// the previous step's correction history invalid as a subspace for
    /// it, so there's nothing to reset, only to not carry over).
    pub(crate) fn should_refactor(&self) -> bool {
        self.v.len() >= self.max_dimension
    }

    /// Accelerate one raw correction. `residual` is the unbalance the raw
    /// solve (`Op^-1 * residual`, already computed by the caller as
    /// `v_raw`) was solved against.
    pub(crate) fn accelerate(
        &mut self,
        v_raw: DVector<f64>,
        residual: &DVector<f64>,
    ) -> DVector<f64> {
        // This correction's contribution to the residual history is only
        // known once the *next* residual arrives — record the previous
        // (v, av) pair (`v` was already pushed by the previous call) now
        // that we have both ends of it.
        if let Some(last_residual) = self.last_residual.take() {
            self.av.push(last_residual - residual);
            debug_assert_eq!(self.v.len(), self.av.len());
        }
        self.last_residual = Some(residual.clone());

        let accelerated = if self.av.is_empty() {
            v_raw.clone()
        } else {
            // Small (k <= max_dimension) dense least-squares via the
            // normal equations — `k` is tiny (Xara defaults
            // `maxDimension = 3`), so this k x k solve is cheap and well
            // enough conditioned; no need for a rank-revealing
            // decomposition the way a large or ill-conditioned system
            // would (nalgebra's `QR::solve` only handles square systems
            // anyway, so a tall-skinny `av_mat` can't use it directly).
            let n = residual.len();
            let k = self.av.len();
            let av_mat = DMatrix::from_fn(n, k, |i, j| self.av[j][i]);
            let ata = av_mat.transpose() * &av_mat;
            let atb = av_mat.transpose() * residual;
            let c = ata.lu().solve(&atb).unwrap_or_else(|| DVector::zeros(k));
            let mut corrected = v_raw.clone();
            for (i, c_i) in c.iter().enumerate() {
                corrected += &self.v[i] * *c_i;
            }
            corrected
        };

        self.v.push(v_raw);
        accelerated
    }
}
