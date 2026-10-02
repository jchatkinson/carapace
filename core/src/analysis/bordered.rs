//! The `(n+1)`-unknown bordered system arc-length continuation solves for
//! its tangent and its coupled corrections (`docs/arclength.md` §4–§5):
//!
//! ```text
//! [ K      column ] [x_q     ]   [rhs_q     ]
//! [ row^T  corner ] [x_lambda] = [rhs_lambda]
//! ```
//!
//! **Block elimination by default.** Only `K` is factored; the border is
//! eliminated through its Schur complement `corner - row^T K^-1 column`.
//! Factoring the bordered matrix directly lets partial pivoting pull the
//! dense constraint row in early, and on a 40-story, 8-bay frame that fills
//! `L+U` to ~10x what `K` alone needs and costs ~17x the factor time
//! (measured by the ignored `bordered_fill_and_factor_time_on_a_frame`
//! test), so the direct factorization is only the **fallback**: used when
//! `K`'s factorization breaks down (an exactly singular `K` at a smooth
//! load limit, where the border keeps the system nonsingular) or when the
//! Schur complement suffers cancellation. Every solve is then checked
//! against the full bordered system's backward error either way.
//!
//! Uses `faer`'s low-level simplicial LU rather than `SparseSolver`'s
//! high-level one because the determinant sign — the bifurcation
//! diagnostic in `docs/arclength.md` §7 — needs the `U` diagonal and both
//! permutations, which the high-level type keeps private. With block
//! elimination, `det A = det K * schur`. COLAMD orderings are cached per
//! pattern (`K`'s pattern only changes with the model; the border is stored
//! with every structural entry, including zeros).

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::linalg::colamd;
use faer::sparse::linalg::lu::simplicial::{
    factorize_simplicial_numeric_lu, factorize_simplicial_numeric_lu_scratch,
    solve_in_place_scratch, SimplicialLu,
};
use faer::sparse::{SymbolicSparseColMat, Triplet};
use faer::{Conj, Mat, Par};
use nalgebra::DVector;

use crate::model::SparseMatrix;

use super::ArcFailure;

/// Scaled backward-error target for an accepted bordered solve
/// (`docs/arclength.md` §6). One step of iterative refinement is tried
/// before a solve is reported inaccurate.
const BACKWARD_ERROR_TARGET: f64 = 1e-10;

/// Block elimination falls back to the direct bordered factorization when
/// the Schur complement keeps less than this fraction of its terms'
/// magnitude (`|schur| <= SCHUR_CANCELLATION * (|corner| + |row|.|b|)`).
const SCHUR_CANCELLATION: f64 = 1e-8;

/// Owns the COLAMD ordering caches for one continuation phase. Kept apart
/// from `SparseSolver`'s `n`-DOF symbolic cache so neither evicts the other.
#[derive(Default)]
pub(crate) struct BorderedSolver {
    k_ordering: Option<Box<Ordering>>,
    full_ordering: Option<Box<Ordering>>,
}

struct Ordering {
    pattern: SymbolicSparseColMat<usize>,
    col_fwd: Vec<usize>,
    col_inv: Vec<usize>,
}

impl Ordering {
    fn for_pattern<'a>(
        cache: &'a mut Option<Box<Ordering>>,
        matrix: &SparseMatrix,
    ) -> Result<&'a Ordering, ArcFailure> {
        let pattern = matrix.symbolic();
        let matches = cache.as_ref().is_some_and(|cached| {
            cached.pattern.col_ptr() == pattern.col_ptr()
                && cached.pattern.row_idx() == pattern.row_idx()
        });
        if !matches {
            let size = matrix.nrows();
            let mut col_fwd = vec![0usize; size];
            let mut col_inv = vec![0usize; size];
            let mut buffer = MemBuffer::new(colamd::order_scratch::<usize>(
                size,
                size,
                pattern.compute_nnz(),
            ));
            colamd::order(
                &mut col_fwd,
                &mut col_inv,
                pattern,
                colamd::Control::default(),
                MemStack::new(&mut buffer),
            )
            .map_err(|_| ArcFailure::SingularSystem)?;
            *cache = Some(Box::new(Ordering {
                pattern: pattern.to_owned().map_err(|_| ArcFailure::SingularSystem)?,
                col_fwd,
                col_inv,
            }));
        }
        Ok(cache.as_ref().expect("populated above"))
    }
}

/// The border of one bordered system. `column` multiplies the scalar
/// unknown in the first `n` equations; `row`/`corner` form the last one.
pub(crate) struct Border<'a> {
    pub(crate) column: &'a DVector<f64>,
    pub(crate) row: &'a DVector<f64>,
    pub(crate) corner: f64,
}

/// A simplicial LU (`P A Q = L U`) with its permutations.
struct Lu {
    lu: SimplicialLu<usize, f64>,
    row_fwd: Vec<usize>,
    row_inv: Vec<usize>,
    col_fwd: Vec<usize>,
    col_inv: Vec<usize>,
}

impl Lu {
    fn factor(matrix: &SparseMatrix, ordering: &Ordering) -> Result<Self, ArcFailure> {
        let size = matrix.nrows();
        let mut lu = SimplicialLu::new();
        let mut row_fwd = vec![0usize; size];
        let mut row_inv = vec![0usize; size];
        let mut buffer = MemBuffer::new(factorize_simplicial_numeric_lu_scratch::<usize, f64>(
            size, size,
        ));
        factorize_simplicial_numeric_lu(
            &mut row_fwd,
            &mut row_inv,
            &mut lu,
            matrix.as_ref(),
            PermRef::new_checked(&ordering.col_fwd, &ordering.col_inv, size),
            MemStack::new(&mut buffer),
        )
        .map_err(|_| ArcFailure::SingularSystem)?;
        Ok(Lu {
            lu,
            row_fwd,
            row_inv,
            col_fwd: ordering.col_fwd.clone(),
            col_inv: ordering.col_inv.clone(),
        })
    }

    /// Sign of the factored matrix's determinant; a zero or nonfinite
    /// pivot means the matrix is (numerically) singular.
    fn det_sign(&self) -> Result<i8, ArcFailure> {
        let mut sign = permutation_sign(&self.row_fwd) * permutation_sign(&self.col_fwd);
        let u = self.lu.u_factor_unsorted();
        for j in 0..self.row_fwd.len() {
            let diagonal = u
                .row_idx_of_col(j)
                .zip(u.val_of_col(j))
                .find_map(|(i, &value)| (i == j).then_some(value))
                .unwrap_or(0.0);
            if diagonal == 0.0 || !diagonal.is_finite() {
                return Err(ArcFailure::SingularSystem);
            }
            if diagonal < 0.0 {
                sign = -sign;
            }
        }
        Ok(sign)
    }

    fn solve(&self, b: &DVector<f64>) -> DVector<f64> {
        let size = b.len();
        let mut x = Mat::<f64>::from_fn(size, 1, |i, _| b[i]);
        let mut buffer = MemBuffer::new(solve_in_place_scratch::<usize, f64>(size, 1, Par::Seq));
        self.lu.solve_in_place_with_conj(
            PermRef::new_checked(&self.row_fwd, &self.row_inv, size),
            PermRef::new_checked(&self.col_fwd, &self.col_inv, size),
            Conj::No,
            x.as_mut(),
            Par::Seq,
            MemStack::new(&mut buffer),
        );
        DVector::from_fn(size, |i, _| x[(i, 0)])
    }

    #[cfg(test)]
    fn fill(&self) -> usize {
        self.lu.l_factor_unsorted().compute_nnz() + self.lu.u_factor_unsorted().compute_nnz()
    }
}

enum Method {
    /// `K_hat` factored; `b = K_hat^-1 column_hat`, `schur = corner_hat -
    /// row_hat . b`.
    Block {
        lu: Lu,
        row: DVector<f64>,
        b: DVector<f64>,
        schur: f64,
    },
    /// The whole bordered matrix factored directly.
    Full(Lu),
}

/// A factorization of the equilibrated bordered matrix `A_hat = R * A * C`,
/// with `R` the row equilibration and `C` the caller's unknown scales.
pub(crate) struct BorderedFactorization {
    method: Method,
    matrix: SparseMatrix,
    row_scale: DVector<f64>,
    col_scale: DVector<f64>,
    det_sign: i8,
}

impl BorderedSolver {
    /// Factors the bordered matrix built from `k` and `border`.
    /// `unknown_scale[i]` is the characteristic size of unknown `i` (the
    /// displacement/rotation scale for equations, the load scale for the
    /// scalar), so the factored matrix is dimensionless; rows are then
    /// equilibrated by their largest entry.
    pub(crate) fn factor(
        &mut self,
        k: &SparseMatrix,
        border: Border<'_>,
        unknown_scale: &DVector<f64>,
    ) -> Result<BorderedFactorization, ArcFailure> {
        let n = k.nrows();
        let size = n + 1;

        let mut triplets = Vec::with_capacity(k.compute_nnz() + 2 * n + 1);
        let k_ref = k.as_ref();
        for j in 0..n {
            for (i, &value) in k_ref.row_idx_of_col(j).zip(k_ref.val_of_col(j)) {
                triplets.push(Triplet::new(i, j, value * unknown_scale[j]));
            }
        }
        for i in 0..n {
            triplets.push(Triplet::new(i, n, border.column[i] * unknown_scale[n]));
            triplets.push(Triplet::new(n, i, border.row[i] * unknown_scale[i]));
        }
        triplets.push(Triplet::new(n, n, border.corner * unknown_scale[n]));
        if !triplets.iter().all(|t| t.val.is_finite()) {
            return Err(ArcFailure::NonFinite);
        }

        let mut row_max = DVector::<f64>::zeros(size);
        for t in &triplets {
            row_max[t.row] = row_max[t.row].max(t.val.abs());
        }
        if row_max.iter().any(|&m| m == 0.0) {
            return Err(ArcFailure::SingularSystem);
        }
        let row_scale = row_max.map(|m| 1.0 / m);
        for t in &mut triplets {
            t.val *= row_scale[t.row];
        }
        let matrix = SparseMatrix::try_new_from_triplets(size, size, &triplets)
            .expect("bordered indices are always in [0, n]");

        let mut column = DVector::<f64>::zeros(n);
        let mut row = DVector::<f64>::zeros(n);
        let mut corner = 0.0;
        let k_triplets: Vec<_> = triplets
            .iter()
            .filter(|t| {
                if t.row == n && t.col == n {
                    corner += t.val;
                } else if t.col == n {
                    column[t.row] += t.val;
                } else if t.row == n {
                    row[t.col] += t.val;
                }
                t.row < n && t.col < n
            })
            .copied()
            .collect();
        let k_hat = SparseMatrix::try_new_from_triplets(n, n, &k_triplets)
            .expect("K indices are always in [0, n)");

        let (method, det_sign) = match self.block(&k_hat, &column, row, corner)? {
            Some(block) => block,
            None => {
                let ordering = Ordering::for_pattern(&mut self.full_ordering, &matrix)?;
                let lu = Lu::factor(&matrix, ordering)?;
                let sign = lu.det_sign()?;
                (Method::Full(lu), sign)
            }
        };

        Ok(BorderedFactorization {
            method,
            matrix,
            row_scale,
            col_scale: unknown_scale.clone(),
            det_sign,
        })
    }

    /// Block elimination, or `None` when `K_hat` is singular or the Schur
    /// complement cancels (the caller then factors the bordered matrix).
    fn block(
        &mut self,
        k_hat: &SparseMatrix,
        column: &DVector<f64>,
        row: DVector<f64>,
        corner: f64,
    ) -> Result<Option<(Method, i8)>, ArcFailure> {
        let ordering = Ordering::for_pattern(&mut self.k_ordering, k_hat)?;
        let Ok(lu) = Lu::factor(k_hat, ordering) else {
            return Ok(None);
        };
        let Ok(k_sign) = lu.det_sign() else {
            return Ok(None);
        };
        let b = lu.solve(column);
        if !b.iter().all(|v| v.is_finite()) {
            return Ok(None);
        }
        let coupling = row.dot(&b);
        let schur = corner - coupling;
        let magnitude = corner.abs()
            + row
                .iter()
                .zip(b.iter())
                .map(|(r, b)| (r * b).abs())
                .sum::<f64>();
        if !schur.is_finite() || schur.abs() <= SCHUR_CANCELLATION * magnitude {
            return Ok(None);
        }
        let sign = if schur > 0.0 { k_sign } else { -k_sign };
        Ok(Some((Method::Block { lu, row, b, schur }, sign)))
    }
}

impl BorderedFactorization {
    /// Sign of `det(A)` for the unscaled bordered matrix (the positive
    /// diagonal scalings don't change it).
    pub(crate) fn det_sign(&self) -> i8 {
        self.det_sign
    }

    /// Solves `A x = rhs` (unscaled quantities in, unscaled solution out),
    /// with one step of iterative refinement when the scaled backward error
    /// misses its target, and rejects nonfinite or still-inaccurate results.
    pub(crate) fn solve(&self, rhs: &DVector<f64>) -> Result<DVector<f64>, ArcFailure> {
        let b = rhs.component_mul(&self.row_scale);
        if !b.iter().all(|v| v.is_finite()) {
            return Err(ArcFailure::NonFinite);
        }
        let mut y = self.solve_scaled(&b);
        let mut error = self.backward_error(&y, &b);
        if error > BACKWARD_ERROR_TARGET && error.is_finite() {
            let residual = &b - self.multiply(&y);
            y += self.solve_scaled(&residual);
            error = self.backward_error(&y, &b);
        }
        if !y.iter().all(|v| v.is_finite()) || !error.is_finite() {
            return Err(ArcFailure::NonFinite);
        }
        if error > BACKWARD_ERROR_TARGET {
            return Err(ArcFailure::InaccurateSolve);
        }
        Ok(y.component_mul(&self.col_scale))
    }

    fn solve_scaled(&self, rhs: &DVector<f64>) -> DVector<f64> {
        match &self.method {
            Method::Full(lu) => lu.solve(rhs),
            Method::Block { lu, row, b, schur } => {
                let n = b.len();
                let a = lu.solve(&rhs.rows(0, n).into_owned());
                let y = (rhs[n] - row.dot(&a)) / schur;
                let x = a - b * y;
                DVector::from_fn(n + 1, |i, _| if i < n { x[i] } else { y })
            }
        }
    }

    fn multiply(&self, y: &DVector<f64>) -> DVector<f64> {
        let a = self.matrix.as_ref();
        let mut result = DVector::<f64>::zeros(y.len());
        for j in 0..y.len() {
            for (i, &value) in a.row_idx_of_col(j).zip(a.val_of_col(j)) {
                result[i] += value * y[j];
            }
        }
        result
    }

    /// `||A_hat y - b||_inf / (||A_hat||_inf ||y||_inf + ||b||_inf)`.
    fn backward_error(&self, y: &DVector<f64>, b: &DVector<f64>) -> f64 {
        let a = self.matrix.as_ref();
        let mut row_sums = DVector::<f64>::zeros(y.len());
        for j in 0..y.len() {
            for (i, &value) in a.row_idx_of_col(j).zip(a.val_of_col(j)) {
                row_sums[i] += value.abs();
            }
        }
        let norm_a = row_sums.amax();
        let residual = (b - self.multiply(y)).amax();
        let denominator = norm_a * y.amax() + b.amax();
        if denominator == 0.0 {
            0.0
        } else {
            residual / denominator
        }
    }

    #[cfg(test)]
    fn uses_block_elimination(&self) -> bool {
        matches!(self.method, Method::Block { .. })
    }
}

/// `+1` for an even permutation, `-1` for an odd one (by cycle count).
fn permutation_sign(forward: &[usize]) -> i8 {
    let mut visited = vec![false; forward.len()];
    let mut sign = 1i8;
    for start in 0..forward.len() {
        if visited[start] {
            continue;
        }
        let mut length = 0;
        let mut i = start;
        while !visited[i] {
            visited[i] = true;
            i = forward[i];
            length += 1;
        }
        if length % 2 == 0 {
            sign = -sign;
        }
    }
    sign
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(n: usize, entries: &[(usize, usize, f64)]) -> SparseMatrix {
        let triplets: Vec<_> = entries
            .iter()
            .map(|&(i, j, v)| Triplet::new(i, j, v))
            .collect();
        SparseMatrix::try_new_from_triplets(n, n, &triplets).unwrap()
    }

    fn ones(n: usize) -> DVector<f64> {
        DVector::from_element(n, 1.0)
    }

    #[test]
    fn singular_k_with_nonsingular_border_solves_exactly() {
        // K = [[1, 1], [1, 1]] is singular; the border makes A nonsingular.
        let k = matrix(2, &[(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]);
        let column = DVector::from_column_slice(&[-1.0, 0.0]);
        let row = DVector::from_column_slice(&[1.0, -1.0]);
        let mut solver = BorderedSolver::default();
        let factor = solver
            .factor(
                &k,
                Border {
                    column: &column,
                    row: &row,
                    corner: 0.0,
                },
                &ones(3),
            )
            .unwrap();
        assert!(
            !factor.uses_block_elimination(),
            "singular K needs the direct fallback"
        );
        // A = [[1,1,-1],[1,1,0],[1,-1,0]], x = [1,2,3] -> b = [0,3,-1].
        let x = factor
            .solve(&DVector::from_column_slice(&[0.0, 3.0, -1.0]))
            .unwrap();
        assert!((x - DVector::from_column_slice(&[1.0, 2.0, 3.0])).amax() < 1e-12);
        // det(A) = 1*(0-0) - 1*(0-0) + (-1)*(-1-1) = 2.
        assert_eq!(factor.det_sign(), 1);
    }

    #[test]
    fn determinant_sign_matches_dense_value_and_unknown_scales_are_inert() {
        // K = diag(-2, 3), border column 0, row 0, corner 5: det = -30.
        let k = matrix(2, &[(0, 0, -2.0), (1, 1, 3.0)]);
        let zeros = DVector::zeros(2);
        let mut solver = BorderedSolver::default();
        for scale in [ones(3), DVector::from_column_slice(&[1e-3, 10.0, 1e4])] {
            let factor = solver
                .factor(
                    &k,
                    Border {
                        column: &zeros,
                        row: &zeros,
                        corner: 5.0,
                    },
                    &scale,
                )
                .unwrap();
            assert!(factor.uses_block_elimination());
            assert_eq!(factor.det_sign(), -1);
            let x = factor
                .solve(&DVector::from_column_slice(&[-2.0, 6.0, 10.0]))
                .unwrap();
            assert!((x - DVector::from_column_slice(&[1.0, 2.0, 2.0])).amax() < 1e-12);
        }
    }

    #[test]
    fn truly_singular_bordered_matrix_is_reported() {
        let k = matrix(2, &[(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]);
        let column = DVector::from_column_slice(&[1.0, 1.0]);
        let row = DVector::from_column_slice(&[1.0, 1.0]);
        let mut solver = BorderedSolver::default();
        let result = solver.factor(
            &k,
            Border {
                column: &column,
                row: &row,
                corner: 1.0,
            },
            &ones(3),
        );
        assert!(matches!(
            result.map(|f| f.solve(&ones(3))),
            Err(ArcFailure::SingularSystem) | Ok(Err(_))
        ));
    }

    /// `docs/arclength.md` §5/§10 step 2: bordered fill and factor time on
    /// a representative frame, against the `n`-DOF supernodal LU the other
    /// integrators use. Not a pass/fail test — run with
    /// `cargo test --release -p carapace-core bordered_fill -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn bordered_fill_and_factor_time_on_a_frame() {
        use crate::analysis::SparseSolver;
        use crate::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Node};
        use std::time::Instant;

        for (stories, bays) in [(10, 3), (20, 5), (40, 8)] {
            let mut domain = Domain::new();
            let mut grid = vec![vec![]; stories + 1];
            for (level, row) in grid.iter_mut().enumerate() {
                for bay in 0..=bays {
                    let node = Node::new([bay as f64 * 6.0, level as f64 * 3.5]);
                    let node = if level == 0 {
                        node.fix(0).fix(1).fix(2)
                    } else {
                        node
                    };
                    row.push(domain.add_node(node));
                }
            }
            for level in 1..=stories {
                for (&below, &above) in grid[level - 1].iter().zip(&grid[level]) {
                    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                        below,
                        above,
                        2e8,
                        0.02,
                        4e-4,
                        GeomTransf::Corotational,
                    )));
                }
                for bay in 0..bays {
                    domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                        grid[level][bay],
                        grid[level][bay + 1],
                        2e8,
                        0.015,
                        3e-4,
                        GeomTransf::Corotational,
                    )));
                }
                domain.load_node(grid[level][0], 0, level as f64);
            }
            domain.number_dofs();
            let (k, _) = domain.assemble_tangent_and_resistance();
            let p = domain.assemble_reference_load_sensitivity(0.0);
            let n = p.len();
            let row = DVector::from_element(n, 1.0);
            let minus_p = -&p;
            let scale = DVector::from_element(n + 1, 1.0);
            let repeats = 20;

            let solver = SparseSolver::new();
            solver.factor(&k).unwrap();
            let start = Instant::now();
            for _ in 0..repeats {
                solver.factor(&k).unwrap();
            }
            let plain = start.elapsed().as_secs_f64() / repeats as f64;

            let mut bordered = BorderedSolver::default();
            let border = || Border {
                column: &minus_p,
                row: &row,
                corner: 1.0,
            };
            let factor = bordered.factor(&k, border(), &scale).unwrap();
            let start = Instant::now();
            for _ in 0..repeats {
                bordered.factor(&k, border(), &scale).unwrap();
            }
            let with_border = start.elapsed().as_secs_f64() / repeats as f64;
            let block_fill = match &factor.method {
                Method::Block { lu, .. } => lu.fill(),
                Method::Full(lu) => lu.fill(),
            };
            let full = {
                let ordering =
                    Ordering::for_pattern(&mut bordered.full_ordering, &factor.matrix).unwrap();
                let start = Instant::now();
                let lu = Lu::factor(&factor.matrix, ordering).unwrap();
                (lu.fill(), start.elapsed().as_secs_f64())
            };
            println!(
                "{stories}x{bays}: n={n} nnz(K)={} | n-DOF LU {:.3} ms | block elimination \
                 nnz(L+U)={block_fill} {:.3} ms ({:.2}x) | direct bordered nnz(L+U)={} {:.3} ms",
                k.compute_nnz(),
                plain * 1e3,
                with_border * 1e3,
                with_border / plain,
                full.0,
                full.1 * 1e3,
            );
            assert!(factor.uses_block_elimination());
        }
    }

    #[test]
    fn permutation_parity() {
        assert_eq!(permutation_sign(&[0, 1, 2]), 1);
        assert_eq!(permutation_sign(&[1, 0, 2]), -1);
        assert_eq!(permutation_sign(&[1, 2, 0]), 1);
    }
}
