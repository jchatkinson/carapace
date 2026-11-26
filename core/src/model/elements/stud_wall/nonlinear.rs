//! Nonlinear sheathing laws: the condensed per-pier system solved at each trial state.

use super::*;

type Comps = ArrayVec<f64, MAX_WALL_PIERS>;

impl<K: NodeKey> StudWall<K> {
    /// The components' responses `q`, tangents `d` and residuals `gamma - y + C q` at the trial strains `gamma`.
    fn local_residual(
        &self,
        cache: &NonlinearCache,
        y: &Comps,
        gamma: &Comps,
    ) -> Option<(Comps, Comps, Comps)> {
        let count = gamma.len();
        let (mut q, mut d) = (Comps::new(), Comps::new());
        for (i, &g) in gamma.iter().enumerate() {
            let (stress, tangent) = if i < self.pier_comps {
                shear_law_response(
                    &self.layout.sheathing.as_ref()?.law,
                    self.material_states[i].as_ref(),
                    g,
                )
            } else {
                self.material_states[i].as_ref()?.trial_stress_tangent(g)
            };
            q.push(stress);
            d.push(tangent);
        }
        let r = (0..count)
            .map(|i| {
                gamma[i] - y[i]
                    + (0..count)
                        .map(|j| cache.compliance[(i, j)] * q[j])
                        .sum::<f64>()
            })
            .collect();
        Some((q, d, r))
    }

    /// Solves the components' local system `gamma + C q(gamma) = y` (`y` the strains the exterior displacements and loads impose, `C`
    /// the condensed frame's compliance) by Newton with a backtracking line search on `|r|`. The laws are monotone and `C` is positive
    /// definite, so the Jacobian `I + C diag(d)` is never singular and its Newton step always reduces `|r|`; full steps alone can cycle
    /// across a kinked law (a hold-down between tension and bearing), which the line search breaks. Returns the strains, the responses
    /// and `dq/dy`.
    pub(super) fn nonlinear_trial(
        &self,
        ue: &[f64],
        cases: &[WallCaseLoad],
    ) -> Option<(
        Comps,
        Comps,
        ArrayVec<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>,
    )> {
        let cache = self.nonlinear.as_ref()?;
        let count = self.comps.len();
        let mut y = Comps::new();
        for i in 0..count {
            let mut yi = (0..ue.len())
                .map(|j| cache.strain_map[(i, j)] * ue[j])
                .sum::<f64>();
            for c in cases {
                yi += cache.load_strain[c.case as usize][i] * c.factor;
            }
            y.push(yi);
        }
        let jacobian = |d: &Comps| {
            let mut jac = ArrayVec::<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>::new();
            jac.extend((0..count * count).map(|k| {
                f64::from(k / count == k % count)
                    + cache.compliance[(k / count, k % count)] * d[k % count]
            }));
            jac
        };
        let norm2 = |r: &Comps| r.iter().map(|v| v * v).sum::<f64>().sqrt();
        let scale = 1. + y.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut gamma = y.clone();
        let (mut q, mut d, mut r) = self.local_residual(cache, &y, &gamma)?;
        for _ in 0..60 {
            if r.iter().fold(0.0f64, |m, v| m.max(v.abs())) <= 1e-14 * scale {
                let jac = jacobian(&d);
                let mut dqdy = ArrayVec::<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>::new();
                dqdy.extend((0..count * count).map(|_| 0.));
                for col in 0..count {
                    let mut a = jac.clone();
                    let mut rhs = Comps::new();
                    rhs.extend((0..count).map(|i| f64::from(i == col)));
                    if !small_lu_solve(&mut a, &mut rhs, count) {
                        return None;
                    }
                    for row in 0..count {
                        dqdy[row * count + col] = d[row] * rhs[row];
                    }
                }
                return Some((gamma, q, dqdy));
            }
            let mut a = jacobian(&d);
            let mut step = r.clone();
            if !small_lu_solve(&mut a, &mut step, count) {
                return None;
            }
            // Backtrack until |r| decreases (Armijo); a vanishing step is taken anyway rather than stalling.
            let n0 = norm2(&r);
            let mut t = 1.;
            loop {
                let trial: Comps = (0..count).map(|i| gamma[i] - t * step[i]).collect();
                if !trial.iter().all(|v| v.is_finite()) {
                    return None;
                }
                let next = self.local_residual(cache, &y, &trial)?;
                if norm2(&next.2) <= (1. - 1e-4 * t) * n0 || t < 1e-8 {
                    gamma = trial;
                    (q, d, r) = next;
                    break;
                }
                t *= 0.5;
            }
        }
        None
    }
}

pub(super) fn make_nonlinear_cache(
    model: &FrameBuild,
    ne: usize,
    piers: &[PierRec],
) -> Result<NonlinearCache, StudWallError> {
    let free: Vec<_> = (ne..model.ndof)
        .filter(|d| !model.fixed_internal.contains(d))
        .collect();
    let p = piers.len();
    let kii = DMatrix::from_fn(free.len(), free.len(), |i, j| model.k[(free[i], free[j])]);
    let lu = kii.lu();
    let beta_i = DMatrix::from_fn(free.len(), p, |i, j| piers[j].beta[free[i]]);
    let mut force_i = beta_i.clone();
    for (j, pier) in piers.iter().enumerate() {
        force_i.column_mut(j).scale_mut(pier.area);
    }
    let solved = if free.is_empty() {
        DMatrix::zeros(0, p)
    } else {
        lu.solve(&force_i)
            .ok_or_else(|| StudWallError("unsheathed internal stiffness K_ii is singular".into()))?
    };
    let mut internal_recovery = DMatrix::zeros(model.ndof - ne, p);
    for (i, &d) in free.iter().enumerate() {
        for j in 0..p {
            internal_recovery[(d - ne, j)] = solved[(i, j)];
        }
    }
    let kei = DMatrix::from_fn(ne, free.len(), |i, j| model.k[(i, free[j])]);
    let mut coupling = DMatrix::zeros(ne, p);
    for j in 0..p {
        for i in 0..ne {
            coupling[(i, j)] = piers[j].area * piers[j].beta[i]
                - (0..free.len())
                    .map(|k| kei[(i, k)] * solved[(k, j)])
                    .sum::<f64>();
        }
    }
    let strain_map = DMatrix::from_fn(p, ne, |i, j| coupling[(j, i)] / piers[i].area);
    let load_solved = if free.is_empty() {
        DMatrix::zeros(0, p)
    } else {
        lu.solve(&beta_i)
            .ok_or_else(|| StudWallError("unsheathed internal stiffness K_ii is singular".into()))?
    };
    let mut load_projection = DMatrix::zeros(p, model.ndof);
    for i in 0..p {
        for (k, &d) in free.iter().enumerate() {
            load_projection[(i, d)] = load_solved[(k, i)];
        }
    }
    let mut compliance: DMatrix<f64> = DMatrix::zeros(p, p);
    for i in 0..p {
        for j in 0..p {
            if !free.is_empty() {
                compliance[(i, j)] = (0..free.len())
                    .map(|k| beta_i[(k, i)] * solved[(k, j)])
                    .sum();
            }
            if i == j {
                compliance[(i, j)] += piers[i].c / piers[i].area;
            }
        }
    }
    if !compliance.iter().all(|x| x.is_finite()) {
        return Err(StudWallError(
            "nonlinear pier condensation is non-finite".into(),
        ));
    }
    Ok(NonlinearCache {
        coupling,
        strain_map,
        compliance,
        load_projection,
        load_strain: Vec::new(),
        internal_recovery,
    })
}

pub(super) fn pier_load_strain(
    projection: &DMatrix<f64>,
    ne: usize,
    load: &DVector<f64>,
) -> DVector<f64> {
    DVector::from_fn(projection.nrows(), |i, _| {
        (ne..load.len()).map(|d| projection[(i, d)] * load[d]).sum()
    })
}

pub(super) fn small_lu_solve(
    a: &mut ArrayVec<f64, { MAX_WALL_PIERS * MAX_WALL_PIERS }>,
    rhs: &mut ArrayVec<f64, MAX_WALL_PIERS>,
    n: usize,
) -> bool {
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        if !a[pivot * n + k].is_finite() || a[pivot * n + k].abs() <= f64::MIN_POSITIVE {
            return false;
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            rhs.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for i in (0..n).rev() {
        let tail: f64 = (i + 1..n).map(|j| a[i * n + j] * rhs[j]).sum();
        rhs[i] = (rhs[i] - tail) / a[i * n + i];
    }
    rhs.iter().all(|v| v.is_finite())
}
