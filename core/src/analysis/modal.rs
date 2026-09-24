use nalgebra::{DMatrix, DVector, SymmetricEigen};
use slotmap::Key;

use crate::model::{Domain, ElementOps};

use super::{AnalysisError, SparseSolver};

/// One computed mode: natural circular frequency (rad/time) and the
/// corresponding free-DOF mode shape, M-normalized (`shape^T * M * shape = 1`)
/// — the normalization modal damping (M6+) needs to apply a per-mode
/// damping ratio consistently.
#[derive(Debug, Clone)]
pub struct Mode {
    pub frequency: f64,
    pub shape: DVector<f64>,
}

/// Resolves implementation-plan §5 open decision #2 for real: a
/// shift-invert Lanczos eigensolver (shift = 0, so the operator is simply
/// `K^-1 * M`), not a dense full-spectrum `nalgebra::SymmetricEigen` and
/// not an ARPACK FFI binding either. Two things made the earlier dense
/// choice wrong beyond just performance, once `SparseSolver` itself went
/// sparse (see §5 decision #1):
///
/// - **Memory.** A dense N×N matrix costs O(n²) regardless of how sparse
///   the model actually is — on a model `SparseSolver` handles comfortably,
///   the dense eigensolve path could exhaust a wasm worker's memory before
///   it's even "slow". `Domain` no longer builds one for this at all.
/// - **Shape of the problem.** A full-spectrum solve computes *all* n
///   eigenpairs even when only the lowest `num_modes` are wanted — the
///   usual case, and the one modal-superposition damping (M6+) actually
///   needs. Lanczos naturally converges the extreme eigenvalues of the
///   shift-inverted operator first, so asking for just `num_modes` does
///   the proportionate amount of work.
///
/// This isn't "roll our own ARPACK": the hard numerical pieces are both
/// already library code — `SparseSolver`'s sparse LU for each shift-invert
/// solve `K*w = M*v`, and `nalgebra::SymmetricEigen` for the small
/// `m*m` (`m` = Lanczos subspace size, not `n`) projected tridiagonal
/// eigenproblem. What's ours is the well-understood outer Lanczos
/// recurrence and (full, since `m` is always small) M-orthogonal
/// re-orthogonalization loop.
///
/// Only `shift = 0` is supported (find the `num_modes` lowest natural
/// frequencies) — by far the standard structural-dynamics case. A nonzero
/// shift (targeting modes near a chosen frequency, or handling a
/// near-singular `K` for buckling analysis) is a natural extension, not
/// implemented until something actually needs it.
///
/// `M` (lumped, `Node::mass`) must be diagonal and strictly positive on
/// every free DOF — a real physical requirement of modal analysis, not an
/// implementation gap.
///
/// Generic over the same kinematic profile as `Domain`/`Analysis` (see
/// `Domain`'s doc comment) — the Lanczos recurrence and mass-normalization
/// above only ever touch `domain` through its generic free-DOF assembly
/// interface (`assemble_mass_diagonal`, `num_free_dofs`,
/// `assemble_tangent_and_resistance`), never anything DOF-count-specific,
/// so one implementation covers both `Domain`/`Domain3` — inferred from
/// `domain`'s concrete type at each call site, same as
/// `AnalysisBuilder::build`.
pub fn modal_analysis<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
    num_modes: usize,
) -> Result<Vec<Mode>, AnalysisError>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
{
    let mass = domain.assemble_mass_diagonal();
    let n = domain.num_free_dofs();

    if num_modes == 0 || num_modes > n {
        return Err(AnalysisError::InvalidModeCount {
            requested: num_modes,
            free_dofs: n,
        });
    }
    for i in 0..n {
        if mass[i] <= 0.0 {
            return Err(AnalysisError::SingularSystem);
        }
    }

    let (k, _resistance) = domain.assemble_tangent_and_resistance();
    let solver = SparseSolver::new();

    let m_dot = |x: &DVector<f64>, y: &DVector<f64>| -> f64 {
        (0..n).map(|i| x[i] * mass[i] * y[i]).sum()
    };
    let m_norm = |x: &DVector<f64>| -> f64 { m_dot(x, x).sqrt() };

    // Lanczos subspace size: bigger than `num_modes` for accuracy (the
    // extra dimensions are "scratch" that improve convergence of the
    // modes we keep), capped at `n` (can't build a larger Krylov subspace
    // than the space itself).
    let subspace_size = (2 * num_modes + 8).min(n);

    // `vectors[0]` is an unused zero placeholder so `vectors[j]`/`beta[j]`
    // line up with the standard 1-indexed Lanczos recurrence; the real
    // basis is `vectors[1..]`.
    let mut vectors: Vec<DVector<f64>> = Vec::with_capacity(subspace_size + 1);
    vectors.push(DVector::zeros(n));
    let seed = DVector::from_element(n, 1.0);
    vectors.push(&seed / m_norm(&seed));

    let mut alpha = Vec::with_capacity(subspace_size);
    let mut beta = vec![0.0];

    let mut built = 0;
    for j in 1..=subspace_size {
        let mv = mass.component_mul(&vectors[j]);
        let mut w = solver.solve(&k, &mv)?;
        w -= beta[j - 1] * &vectors[j - 1];

        let alpha_j = m_dot(&w, &vectors[j]);
        w -= alpha_j * &vectors[j];

        // Full M-orthogonal re-orthogonalization against every prior
        // basis vector — cheap since `subspace_size` is always small,
        // and necessary since three-term recurrence alone loses
        // orthogonality to rounding error within a handful of steps.
        for v in vectors.iter().skip(1) {
            let proj = m_dot(&w, v);
            w -= proj * v;
        }

        alpha.push(alpha_j);
        built = j;

        let beta_j = m_norm(&w);
        if beta_j < 1e-10 || j == subspace_size {
            break;
        }
        beta.push(beta_j);
        vectors.push(&w / beta_j);
    }

    let mut t = DMatrix::<f64>::zeros(built, built);
    for i in 0..built {
        t[(i, i)] = alpha[i];
    }
    for i in 0..built.saturating_sub(1) {
        t[(i, i + 1)] = beta[i + 1];
        t[(i + 1, i)] = beta[i + 1];
    }

    let eigen = SymmetricEigen::new(t);

    // Ritz values `theta` of `K^-1 * M` relate to the generalized problem's
    // `lambda = omega^2` via `lambda = 1/theta` (shift = 0). The largest
    // `theta` are the best-converged Lanczos results and correspond to the
    // smallest `lambda` — the lowest frequencies, which is what we want.
    let mut modes: Vec<Mode> = (0..built)
        .filter(|&i| eigen.eigenvalues[i] > 1e-12)
        .map(|i| {
            let theta = eigen.eigenvalues[i];
            let y = eigen.eigenvectors.column(i);
            let mut shape = DVector::<f64>::zeros(n);
            for l in 0..built {
                shape += y[l] * &vectors[l + 1];
            }
            Mode {
                frequency: (1.0 / theta).sqrt(),
                shape,
            }
        })
        .collect();
    modes.sort_by(|a, b| a.frequency.partial_cmp(&b.frequency).unwrap());
    modes.truncate(num_modes);

    Ok(modes)
}
