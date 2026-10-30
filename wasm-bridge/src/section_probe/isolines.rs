//! Direct fixed-axial-load / fixed-moment-angle solve on the crushing-limited PMM boundary.
//!
//! The boundary is parametrized by `(theta, s)`: `theta` in `[0, pi]` is the profile angle and `s` a continuous
//! position along the closed loop of strain profiles `(et, eb)` (linear interpolation between consecutive
//! profiles). The loop must be symmetric under `(et, eb) -> (eb, et)` with reversed order, so that
//! `(theta + pi, s)` and `(theta, last - s)` are the same strain field; this is the seam the solver wraps over.
//!
//! For each requested axial force `P` and moment angle `alpha` the solver finds `(theta, s)` with `N = P` and the
//! moment vector `(m2, m3)` along `alpha`, by finite-difference Newton with halving line search and a
//! pattern-search polish (a port of CalcsApp's `solvePMMStateFromSeed`). The root with the largest moment along
//! `alpha` is returned (the outer branch of the interaction curve).
//!
//! Seeding, after the investigation in `docs/CARAPACE_MIGRATION_PLAN.md` (probe seeding):
//! - walk `P` from high (tension) to low (compression) and continue from the previous level's root, which tracks
//!   the outer branch toward the axial apex; a continuation root is accepted only if its moment is close to the
//!   previous one;
//! - otherwise seed from the sampled grid: axial-force crossings along every `seed_angle_stride`-th sampled angle
//!   (linear interpolation of the sampled `N`, no extra evaluations), keeping the seeds with the smallest
//!   moment-alignment residual and the largest moment along `alpha`;
//! - if no sampled line crosses `P` (it equals the sampled extreme), the nearest sampled states.
//!
//! Known limitation: within roughly 10% of an angle's axial apex the interaction curve has an inner and an outer
//! root that merge at a fold, and the outer root's basin is narrow; those few points may return the inner root.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use super::{Factoring, SectionProbe, SectionProbeError};

#[derive(Debug, Clone)]
pub(super) struct GridState {
    pub theta: f64,
    pub profile: f64,
    pub n: f64,
    /// Factored axial force (equal to `n` without factoring).
    pub nf: f64,
}

#[derive(Debug, Clone)]
pub(super) struct Grid {
    profiles: Vec<[f64; 2]>,
    states: Vec<GridState>,
    /// Per sampled angle: indices into `states`, in profile order.
    by_angle: Vec<Vec<usize>>,
    n_min: f64,
    n_max: f64,
    m_max: f64,
    nf_min: f64,
    nf_max: f64,
    mf_max: f64,
    factoring: Option<Factoring>,
}

#[derive(Debug, Clone, Copy, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct IsolineOptions {
    #[serde(default = "d_seed_stride")]
    pub seed_angle_stride: usize,
    #[serde(default = "d_three")]
    pub residual_seeds: usize,
    #[serde(default = "d_three")]
    pub moment_seeds: usize,
    #[serde(default = "d_eight")]
    pub nearest_seeds: usize,
    #[serde(default = "d_true")]
    pub continuation: bool,
    /// Largest accepted change in moment along `alpha` between consecutive levels for a continuation root, as a
    /// fraction of the largest sampled moment.
    #[serde(default = "d_jump")]
    pub continuation_jump: f64,
    #[serde(default = "d_accept")]
    pub accept_tolerance: f64,
    /// Fraction of the sampled axial range, at the compression end, where a root is not continued from the
    /// neighbouring moment angle (inner and outer roots merge at a fold there; see the module docs).
    #[serde(default = "d_apex")]
    pub apex_fraction: f64,
}

fn d_seed_stride() -> usize {
    5
}
fn d_three() -> usize {
    3
}
fn d_eight() -> usize {
    8
}
fn d_true() -> bool {
    true
}
fn d_jump() -> f64 {
    0.3
}
fn d_apex() -> f64 {
    0.1
}
fn d_accept() -> f64 {
    2e-6
}

impl Default for IsolineOptions {
    fn default() -> Self {
        IsolineOptions {
            seed_angle_stride: d_seed_stride(),
            residual_seeds: 3,
            moment_seeds: 3,
            nearest_seeds: 8,
            continuation: true,
            continuation_jump: d_jump(),
            accept_tolerance: d_accept(),
            apex_fraction: d_apex(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct IsolineRequest {
    /// Axial forces to solve at (tension positive), in any order.
    pub levels: Vec<f64>,
    /// Moment-vector directions `alpha` in radians: `(m2, m3) = M (cos alpha, sin alpha)`.
    pub moment_angles: Vec<f64>,
    #[serde(default)]
    pub options: IsolineOptions,
    /// Solve in factored space (needs a factoring in `build_grid`): `levels` are factored axial forces and the
    /// moments are factored.
    #[serde(default)]
    pub factored: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct IsolineResult {
    pub points: Vec<IsolinePoint>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct IsolinePoint {
    pub moment_angle_index: u32,
    pub level_index: u32,
    pub found: bool,
    pub n: f64,
    pub mz: f64,
    pub my: f64,
    pub m2: f64,
    pub m3: f64,
    /// Moment magnitude along the requested angle.
    pub moment: f64,
    /// Profile angle and position along the profile loop.
    pub theta: f64,
    pub profile: f64,
    pub top_strain: f64,
    pub bottom_strain: f64,
    /// Neutral-axis depth from the extreme compression edge; `None` for uniform strain.
    pub c: Option<f64>,
    /// ACI `phi` and the largest bar tensile strain behind it (see `Factoring::Aci`).
    pub phi: Option<f64>,
    pub eps_t: Option<f64>,
    /// Section evaluations spent on this point.
    pub evals: u32,
}

/// Normalized residual below which Newton and the pattern search stop (the acceptance tolerance is 2e-6).
const NEWTON_STOP: f64 = 1e-6;

struct Target {
    p: f64,
    cos: f64,
    sin: f64,
    sp: f64,
    sm: f64,
}

impl Target {
    fn residual(&self, n: f64, m2: f64, m3: f64) -> (f64, f64) {
        (
            (n - self.p) / self.sp,
            (m2 * self.sin - m3 * self.cos) / self.sm,
        )
    }
    fn signed(&self, m2: f64, m3: f64) -> f64 {
        m2 * self.cos + m3 * self.sin
    }
}

#[derive(Clone, Copy)]
struct Root {
    value: f64,
    theta: f64,
    profile: f64,
    n: f64,
    m2: f64,
    m3: f64,
}

struct Solver<'a> {
    probe: &'a SectionProbe,
    grid: &'a Grid,
    factored: bool,
    evals: u32,
}

impl Grid {
    fn axial(&self, state: &GridState, factored: bool) -> f64 {
        if factored {
            state.nf
        } else {
            state.n
        }
    }

    fn last_profile(&self) -> f64 {
        (self.profiles.len() - 1) as f64
    }

    fn interpolate(&self, s: f64) -> (f64, f64) {
        let b = s.clamp(0.0, self.last_profile());
        let lo = b.floor() as usize;
        let hi = (lo + 1).min(self.profiles.len() - 1);
        let t = b - lo as f64;
        let (a, c) = (self.profiles[lo], self.profiles[hi]);
        (a[0] + (c[0] - a[0]) * t, a[1] + (c[1] - a[1]) * t)
    }

    /// Wraps `(theta, s)` into `theta in [0, pi]`, mirroring `s` across the seam.
    fn wrap(&self, mut theta: f64, mut s: f64) -> (f64, f64) {
        let last = self.last_profile();
        while theta < 0.0 {
            theta += PI;
            s = last - s;
        }
        while theta > PI {
            theta -= PI;
            s = last - s;
        }
        (theta, s.clamp(0.0, last))
    }
}

impl SectionProbe {
    /// Samples the PMM boundary at every `(angle, profile)` pair and keeps it for [`Self::solve_isolines`].
    /// `profiles` holds `(et, eb)` pairs forming the closed, seam-symmetric loop described in the module docs.
    pub fn build_grid(
        &mut self,
        angles: &[f64],
        profiles: &[f64],
        factoring: Option<Factoring>,
    ) -> Result<usize, SectionProbeError> {
        if let Some(f) = &factoring {
            f.validate(self.group_count)
                .map_err(|reason| SectionProbeError::InvalidRequest { reason })?;
        }
        let pairs: Vec<[f64; 2]> = profiles.chunks_exact(2).map(|p| [p[0], p[1]]).collect();
        if pairs.len() < 3 || angles.is_empty() {
            return Err(SectionProbeError::InvalidRequest {
                reason: "a grid needs at least one angle and three profiles".into(),
            });
        }
        let mut states = Vec::with_capacity(angles.len() * pairs.len());
        let mut by_angle = Vec::with_capacity(angles.len());
        let (mut n_min, mut n_max, mut m_max) = (f64::INFINITY, f64::NEG_INFINITY, 1.0f64);
        let (mut nf_min, mut nf_max, mut mf_max) = (f64::INFINITY, f64::NEG_INFINITY, 1.0f64);
        for &theta in angles {
            let mut row = Vec::with_capacity(pairs.len());
            for (pi, p) in pairs.iter().enumerate() {
                let (e0, kz, ky) = self.profile_state(theta, p[0], p[1])?;
                let r = self.evaluate(e0, kz, ky);
                let rf = match &factoring {
                    Some(f) => self.evaluate_factored(f, e0, kz, ky).response,
                    None => r,
                };
                n_min = n_min.min(r.n);
                n_max = n_max.max(r.n);
                m_max = m_max.max(r.m2.hypot(r.m3));
                nf_min = nf_min.min(rf.n);
                nf_max = nf_max.max(rf.n);
                mf_max = mf_max.max(rf.m2.hypot(rf.m3));
                row.push(states.len());
                states.push(GridState {
                    theta,
                    profile: pi as f64,
                    n: r.n,
                    nf: rf.n,
                });
            }
            by_angle.push(row);
        }
        let count = states.len();
        self.grid = Some(Grid {
            profiles: pairs,
            states,
            by_angle,
            n_min,
            n_max,
            m_max,
            nf_min,
            nf_max,
            mf_max,
            factoring,
        });
        Ok(count)
    }

    /// `(n_min, n_max)` over the sampled boundary, if a grid is built (factored axial forces with `factored`).
    pub fn grid_axial_range(&self, factored: bool) -> Option<(f64, f64)> {
        self.grid.as_ref().map(|g| {
            if factored {
                (g.nf_min, g.nf_max)
            } else {
                (g.n_min, g.n_max)
            }
        })
    }

    /// Every `(angle, profile)` state, angles outermost, as rows of [`STATE_STRIDE`] numbers: nominal `n, m2, m3`,
    /// factored `n, m2, m3` (equal to nominal without factoring), `phi`, `eps_t` (NaN when not defined), `c`
    /// (NaN for uniform strain) and the profile's `et, eb`.
    pub fn sample_states(
        &self,
        angles: &[f64],
        profiles: &[f64],
        factoring: Option<&Factoring>,
    ) -> Result<Vec<f64>, SectionProbeError> {
        if let Some(f) = factoring {
            f.validate(self.group_count)
                .map_err(|reason| SectionProbeError::InvalidRequest { reason })?;
        }
        let mut out = Vec::with_capacity(angles.len() * profiles.len() / 2 * STATE_STRIDE);
        for &theta in angles {
            let (lo, hi) = self.extents(theta);
            for p in profiles.chunks_exact(2) {
                let (e0, kz, ky) = self.profile_state(theta, p[0], p[1])?;
                let r = self.evaluate(e0, kz, ky);
                let (rf, phi, eps_t) = match factoring {
                    Some(f) => {
                        let x = self.evaluate_factored(f, e0, kz, ky);
                        (x.response, x.phi, x.eps_t)
                    }
                    None => (r, None, None),
                };
                out.extend_from_slice(&[
                    r.n,
                    r.m2,
                    r.m3,
                    rf.n,
                    rf.m2,
                    rf.m3,
                    phi.unwrap_or(f64::NAN),
                    eps_t.unwrap_or(f64::NAN),
                    neutral_axis_depth(hi - lo, p[0], p[1]).unwrap_or(f64::NAN),
                    p[0],
                    p[1],
                ]);
            }
        }
        Ok(out)
    }

    pub fn solve_isolines(
        &self,
        request: &IsolineRequest,
    ) -> Result<Vec<IsolinePoint>, SectionProbeError> {
        let grid = self.grid.as_ref().ok_or(SectionProbeError::NoGrid)?;
        let factored = request.factored;
        if factored && grid.factoring.is_none() {
            return Err(SectionProbeError::InvalidRequest {
                reason: "factored isolines need a factoring in buildGrid".into(),
            });
        }
        let opts = request.options;
        let (sp, sm) = if factored {
            ((grid.nf_max - grid.nf_min).max(1.0), grid.mf_max)
        } else {
            ((grid.n_max - grid.n_min).max(1.0), grid.m_max)
        };
        let (grid_lo, grid_hi) = if factored {
            (grid.nf_min, grid.nf_max)
        } else {
            (grid.n_min, grid.n_max)
        };
        let mut solver = Solver {
            probe: self,
            grid,
            factored,
            evals: 0,
        };
        let mut order: Vec<usize> = (0..request.levels.len()).collect();
        order.sort_by(|a, b| request.levels[*b].partial_cmp(&request.levels[*a]).unwrap());

        let mut points = vec![None; request.levels.len() * request.moment_angles.len()];
        // Roots of the previous moment angle by level: a second continuation candidate for neighbouring angles.
        let mut across: Vec<Option<(Root, f64)>> = vec![None; request.levels.len()];
        for (ai, &alpha) in request.moment_angles.iter().enumerate() {
            let (sin, cos) = alpha.sin_cos();
            let mut previous: Option<(Root, f64)> = None;
            let mut current: Vec<Option<(Root, f64)>> = vec![None; request.levels.len()];
            for &li in &order {
                let target = Target {
                    p: request.levels[li],
                    cos,
                    sin,
                    sp,
                    sm,
                };
                let before = solver.evals;
                let in_apex_zone = (target.p - grid_lo) < opts.apex_fraction * (grid_hi - grid_lo);
                let neighbour = if in_apex_zone { None } else { across[li] };
                // Beyond the sampled axial range there is no boundary point to find.
                let slack = 1e-9 * (grid_hi - grid_lo).abs().max(1.0);
                let root = if target.p < grid_lo - slack || target.p > grid_hi + slack {
                    None
                } else {
                    solver.solve_level(&target, &opts, &[previous, neighbour])
                };
                let evals = solver.evals - before;
                previous = root.map(|r| (r, target.signed(r.m2, r.m3)));
                current[li] = previous;
                points[ai * request.levels.len() + li] = Some(match root {
                    Some(r) => {
                        let (et, eb) = grid.interpolate(r.profile);
                        let moment = target.signed(r.m2, r.m3);
                        let (lo, hi) = self.extents(r.theta);
                        let (phi, eps_t) = match (&grid.factoring, factored) {
                            (Some(f), true) => {
                                let (e0, kz, ky) = self.profile_state(r.theta, et, eb)?;
                                let x = self.evaluate_factored(f, e0, kz, ky);
                                (x.phi, x.eps_t)
                            }
                            _ => (None, None),
                        };
                        IsolinePoint {
                            moment_angle_index: ai as u32,
                            level_index: li as u32,
                            found: true,
                            n: r.n,
                            mz: -r.m2,
                            my: r.m3,
                            m2: r.m2,
                            m3: r.m3,
                            moment,
                            theta: r.theta,
                            profile: r.profile,
                            top_strain: et,
                            bottom_strain: eb,
                            c: neutral_axis_depth(hi - lo, et, eb),
                            phi,
                            eps_t,
                            evals,
                        }
                    }
                    None => IsolinePoint {
                        moment_angle_index: ai as u32,
                        level_index: li as u32,
                        found: false,
                        n: target.p,
                        mz: 0.0,
                        my: 0.0,
                        m2: 0.0,
                        m3: 0.0,
                        moment: 0.0,
                        theta: 0.0,
                        profile: 0.0,
                        top_strain: 0.0,
                        bottom_strain: 0.0,
                        c: None,
                        phi: None,
                        eps_t: None,
                        evals,
                    },
                });
            }
            across = current;
        }
        Ok(points
            .into_iter()
            .map(|p| p.expect("every slot filled"))
            .collect())
    }
}

impl Solver<'_> {
    /// `(n, m2, m3)` at a wrapped `(theta, s)`, or `None` for a degenerate section.
    fn eval(&mut self, theta: f64, s: f64) -> Option<(f64, f64, f64)> {
        self.evals += 1;
        let (theta, s) = self.grid.wrap(theta, s);
        let (et, eb) = self.grid.interpolate(s);
        let (e0, kz, ky) = self.probe.profile_state(theta, et, eb).ok()?;
        let r = match (&self.grid.factoring, self.factored) {
            (Some(f), true) => self.probe.evaluate_factored(f, e0, kz, ky).response,
            _ => self.probe.evaluate(e0, kz, ky),
        };
        Some((r.n, r.m2, r.m3))
    }

    fn residual_at(
        &mut self,
        t: &Target,
        theta: f64,
        s: f64,
    ) -> Option<((f64, f64), (f64, f64, f64))> {
        let state = self.eval(theta, s)?;
        Some((t.residual(state.0, state.1, state.2), state))
    }

    /// Newton from `seed`, then pattern search from the best iterate; the best state seen is returned.
    fn newton(&mut self, t: &Target, seed: (f64, f64)) -> Option<Root> {
        let last = self.grid.last_profile();
        let (mut theta, mut s) = self.grid.wrap(seed.0, seed.1);
        let (mut f, mut state) = self.residual_at(t, theta, s)?;
        let mut best = Root {
            value: f64::INFINITY,
            theta,
            profile: s,
            n: state.0,
            m2: state.1,
            m3: state.2,
        };
        for _ in 0..16 {
            let value = f.0.hypot(f.1);
            if value < best.value {
                best = Root {
                    value,
                    theta,
                    profile: s,
                    n: state.0,
                    m2: state.1,
                    m3: state.2,
                };
            }
            if value < NEWTON_STOP {
                break;
            }
            let (da, ds) = (2e-4, 0.02);
            let Some((fa, _)) = self.residual_at(t, theta + da, s) else {
                break;
            };
            let Some((fp, _)) = self.residual_at(t, theta, (s + ds).clamp(0.0, last)) else {
                break;
            };
            let (j11, j21) = ((fa.0 - f.0) / da, (fa.1 - f.1) / da);
            let (j12, j22) = ((fp.0 - f.0) / ds, (fp.1 - f.1) / ds);
            let det = j11 * j22 - j12 * j21;
            if !det.is_finite() || det.abs() < 1e-11 {
                break;
            }
            let (na, ns) = (
                (-f.0 * j22 + j12 * f.1) / det,
                (j21 * f.0 - j11 * f.1) / det,
            );
            let mut accepted = false;
            let mut damping = 1.0;
            while damping >= 1.0 / 32.0 {
                let (wa, ws) = self.grid.wrap(theta + damping * na, s + damping * ns);
                if let Some((trial_f, trial_state)) = self.residual_at(t, wa, ws) {
                    if trial_f.0.hypot(trial_f.1) < value {
                        theta = wa;
                        s = ws;
                        f = trial_f;
                        state = trial_state;
                        accepted = true;
                        break;
                    }
                }
                damping /= 2.0;
            }
            if !accepted {
                break;
            }
        }
        // Yielding can make the finite-difference Jacobian singular at a profile knot: finish with a bounded
        // pattern search from the best iterate instead of dropping the level.
        // A root this far off is in the wrong basin; polishing it only burns evaluations.
        if best.value > 1e-2 {
            return best.value.is_finite().then_some(best);
        }
        let (mut ra, mut rs) = (0.08, 1.0);
        let (mut ca, mut cs) = (best.theta, best.profile);
        for _ in 0..30 {
            if best.value <= NEWTON_STOP {
                break;
            }
            let mut improved = false;
            for (dx, dy) in [
                (1.0, 0.0),
                (-1.0, 0.0),
                (0.0, 1.0),
                (0.0, -1.0),
                (1.0, 1.0),
                (1.0, -1.0),
                (-1.0, 1.0),
                (-1.0, -1.0),
            ] {
                let (wa, ws) = self.grid.wrap(ca + dx * ra, cs + dy * rs);
                if let Some((trial_f, trial_state)) = self.residual_at(t, wa, ws) {
                    let value = trial_f.0.hypot(trial_f.1);
                    if value < best.value {
                        best = Root {
                            value,
                            theta: wa,
                            profile: ws,
                            n: trial_state.0,
                            m2: trial_state.1,
                            m3: trial_state.2,
                        };
                        ca = wa;
                        cs = ws;
                        improved = true;
                    }
                }
            }
            if !improved {
                ra /= 2.0;
                rs /= 2.0;
                if ra < 1e-6 && rs < 1e-4 {
                    break;
                }
            }
        }
        best.value.is_finite().then_some(best)
    }

    fn accept(&self, t: &Target, root: &Root, tolerance: f64) -> bool {
        root.value < tolerance && t.signed(root.m2, root.m3) >= 0.0
    }

    fn solve_level(
        &mut self,
        t: &Target,
        opts: &IsolineOptions,
        previous: &[Option<(Root, f64)>],
    ) -> Option<Root> {
        if opts.continuation {
            for &(prev, prev_signed) in previous.iter().flatten() {
                if let Some(root) = self.newton(t, (prev.theta, prev.profile)) {
                    let signed = t.signed(root.m2, root.m3);
                    if self.accept(t, &root, opts.accept_tolerance)
                        && (signed - prev_signed).abs() < opts.continuation_jump * t.sm
                    {
                        return Some(root);
                    }
                }
            }
        }
        let grid = self.grid;
        let last = grid.last_profile();
        let mut seeds: Vec<(f64, f64)> = Vec::new();
        for ai in (0..grid.by_angle.len()).step_by(opts.seed_angle_stride.max(1)) {
            let row = &grid.by_angle[ai];
            for w in row.windows(2) {
                let (a, b) = (&grid.states[w[0]], &grid.states[w[1]]);
                let (r0, r1) = (
                    grid.axial(a, self.factored) - t.p,
                    grid.axial(b, self.factored) - t.p,
                );
                if r0 * r1 > 0.0 {
                    continue;
                }
                let f = if (r1 - r0).abs() < 1e-30 {
                    0.0
                } else {
                    -r0 / (r1 - r0)
                };
                let s = a.profile + f * (b.profile - a.profile);
                seeds.push((a.theta, s));
                if ai == 0 {
                    seeds.push((PI, last - s));
                }
            }
        }
        if seeds.is_empty() {
            let mut near: Vec<(f64, usize)> = grid
                .states
                .iter()
                .enumerate()
                .map(|(i, g)| ((grid.axial(g, self.factored) - t.p).abs(), i))
                .collect();
            near.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            for &(_, i) in near.iter().take(opts.nearest_seeds) {
                seeds.push((grid.states[i].theta, grid.states[i].profile));
            }
        }
        // One evaluation per seed: keep the best by alignment residual and by moment along alpha.
        let mut scored: Vec<(f64, f64, (f64, f64))> = Vec::new();
        for seed in seeds {
            if let Some((r, state)) = self.residual_at(t, seed.0, seed.1) {
                let signed = t.signed(state.1, state.2);
                if signed >= 0.0 {
                    scored.push((r.1.abs(), signed, seed));
                }
            }
        }
        let mut picks: Vec<(f64, f64)> = Vec::new();
        scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        picks.extend(scored.iter().take(opts.residual_seeds).map(|q| q.2));
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        for q in scored.iter().take(opts.moment_seeds) {
            if !picks
                .iter()
                .any(|p| (p.0 - q.2 .0).abs() < 1e-9 && (p.1 - q.2 .1).abs() < 1e-9)
            {
                picks.push(q.2);
            }
        }
        let mut found: Option<(f64, Root)> = None;
        for seed in picks {
            if let Some(root) = self.newton(t, seed) {
                if self.accept(t, &root, opts.accept_tolerance) {
                    let signed = t.signed(root.m2, root.m3);
                    if found.as_ref().map_or(true, |(best, _)| signed > *best) {
                        found = Some((signed, root));
                    }
                }
            }
        }
        found.map(|(_, root)| root)
    }
}

/// Numbers per row of [`SectionProbe::sample_states`].
pub const STATE_STRIDE: usize = 11;

/// Depth of the neutral axis from the extreme compression edge for a profile over section height `h`; `None`
/// for (near) uniform strain, where it lies at infinity.
fn neutral_axis_depth(h: f64, et: f64, eb: f64) -> Option<f64> {
    if (et - eb).abs() <= 1e-12 {
        None
    } else if et <= eb {
        Some(h * et / (et - eb))
    } else {
        Some(h * eb / (eb - et))
    }
}
