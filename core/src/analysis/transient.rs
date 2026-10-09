use nalgebra::DVector;
use slotmap::Key;

use crate::model::{
    Domain, Element, Element3, ElementOps, Node3Id, NodeId, NDF, PLANAR_NDIM, SPATIAL_NDF,
    SPATIAL_NDIM,
};

use super::{
    convergence::ForceTolerance, iterate_to_equilibrium, Algorithm, AnalysisError, ConvergenceTest,
    GroundMotion, RayleighDamping, SparseFactorization, SparseSolver,
};

#[derive(Debug, Clone, Copy)]
pub struct TransientStepResult {
    pub step: usize,
    pub time: f64,
    /// `1` for `Algorithm::Linear` (never iterates — today's original,
    /// only behavior before `docs/algorithms.md` §6.1's corrector existed).
    pub iterations: usize,
    /// How many times the effective operator was actually factored this
    /// step — see `StepResult::factorizations`' doc comment for what this
    /// proves about `TangentStrategy`.
    pub factorizations: usize,
}

/// Time-history (transient) analysis via Newmark-beta integration — a
/// distinct type from `Analysis`, not a third `Integrator`/`Algorithm`
/// variant squeezed into it: static analysis is about a scalar load factor
/// and (optionally) Newton-iterating displacement at a fixed time, dynamic
/// analysis is about propagating a very different state (displacement,
/// velocity, acceleration) forward through *actual* time. Trying to unify
/// them into one `Analysis`/`AnalysisBuilder` would blur two genuinely
/// different physical processes for no real benefit — same reasoning as
/// `modal_analysis` being a free function rather than an `Analysis` mode.
///
/// Fixed at the "average acceleration" Newmark parameters
/// (`beta = 1/4, gamma = 1/2`) — unconditionally stable, the standard
/// default (matches OpenSees' own `Newmark` integrator default). Exact for
/// genuinely linear elastic response (`Truss`/
/// `ZeroLength`/`ElasticBeamColumn` with linear materials) since the
/// tangent stiffness formed once per step doesn't change with displacement
/// — one linear solve per step, no Newton iteration, mirroring `Algorithm::
/// Linear`'s same caveat for the static case. A Newton-iterated corrector
/// (needed once nonlinear materials are driven dynamically) is a natural
/// extension, not built until something needs it.
///
/// Generic over the same kinematic profile as `Domain`/`Analysis` (see
/// `Domain`'s doc comment) — Newmark stepping, the initial-acceleration
/// solve, and ground-motion excitation only ever touch `domain` through its
/// generic free-DOF interface (mass diagonal, gather/scatter, tangent
/// assembly, `direction_incidence`), never anything DOF-count-specific, so
/// one implementation covers both `Domain`/`Domain3`. Defaults to the
/// planar profile, so bare `TransientAnalysis` keeps working unchanged;
/// `TransientAnalysis3` (spatial) is a type alias below. A spatial
/// `GroundMotion::direction` of `0..3` selects `ux`/`uy`/`uz`, unchanged
/// from planar beyond the extra valid index.
pub struct TransientAnalysis<
    const NDIM: usize = PLANAR_NDIM,
    const NDOF: usize = NDF,
    NId = NodeId,
    E = Element,
> where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    domain: Domain<NDIM, NDOF, NId, E>,
    mass: DVector<f64>,
    damping: RayleighDamping,
    dt: f64,
    solver: SparseSolver,
    time: f64,
    step_count: usize,
    /// Each active `GroundMotion` paired with its precomputed `mass ⊙
    /// direction_incidence` (`Domain::direction_incidence`) — computed once
    /// (here, and again whenever `with_ground_motion` adds one), not every
    /// step, mirroring how `mass` itself is precomputed once.
    ground_motions: Vec<(GroundMotion, DVector<f64>)>,
    /// `Algorithm::Linear` (today's original, only behavior) unless
    /// `with_algorithm` opts into a real Newton corrector — see
    /// `docs/algorithms.md` §6.1 and `with_algorithm`'s doc comment.
    algorithm: Algorithm,
    test: ConvergenceTest,
    /// Only ever populated/read when `algorithm` is `Newton`/`KrylovNewton`
    /// — see `Analysis::cached_factorization`'s doc comment for what this
    /// caches and why it persists across `step()` calls.
    cached_factorization: Option<SparseFactorization>,
}

/// `TransientAnalysis`'s spatial instantiation — see `Domain3`'s doc
/// comment for why this is a type alias rather than a hand-duplicated
/// struct.
pub type TransientAnalysis3 = TransientAnalysis<SPATIAL_NDIM, SPATIAL_NDF, Node3Id, Element3>;

impl<const NDIM: usize, const NDOF: usize, NId, E> TransientAnalysis<NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    /// Builds a transient analysis from a domain whose nodes may already
    /// carry initial conditions (`Node::with_initial_displacement`/
    /// `with_initial_velocity`). Computes the one thing Newmark needs that
    /// isn't a direct initial condition: a consistent initial acceleration
    /// from equilibrium at t=0, `M*a0 = F0 - K*u0 - C*v0 - M*ι*ag(0)` (the
    /// last term only present once `with_ground_motion` adds one — see
    /// `recompute_initial_acceleration`). `M` is diagonal, so this is a
    /// single elementwise division, no solve needed.
    pub fn new(
        mut domain: Domain<NDIM, NDOF, NId, E>,
        damping: RayleighDamping,
        dt: f64,
    ) -> Result<Self, AnalysisError> {
        assert!(dt > 0.0, "Newmark dt must be positive");

        domain.validate().map_err(AnalysisError::InvalidModel)?;
        let mass = domain.assemble_mass_diagonal();
        let n = domain.num_free_dofs();
        for i in 0..n {
            if mass[i] <= 0.0 {
                return Err(AnalysisError::SingularSystem);
            }
        }

        let mut analysis = TransientAnalysis {
            domain,
            mass,
            damping,
            dt,
            solver: SparseSolver::new(),
            time: 0.0,
            step_count: 0,
            ground_motions: Vec::new(),
            algorithm: Algorithm::Linear,
            // Never read while `algorithm` is `Linear` — a placeholder
            // that would be an obviously-wrong choice if `with_algorithm`
            // were accidentally skipped, rather than a plausible-looking
            // default silently doing the wrong thing.
            test: ConvergenceTest::NormUnbalance {
                tol: f64::NAN,
                max_iter: 0,
            },
            cached_factorization: None,
        };
        analysis.recompute_initial_acceleration();
        Ok(analysis)
    }

    /// Opts into a real Newton corrector (`docs/algorithms.md` §6.1)
    /// instead of the default `Algorithm::Linear` (today's original
    /// behavior: one effective-system solve per step, no iteration, exact
    /// only for linear-material response). Mirrors `with_ground_motion`'s
    /// consuming-builder style; unlike it, doesn't need to re-derive
    /// anything (the initial acceleration doesn't depend on `algorithm`).
    pub fn with_algorithm(mut self, algorithm: Algorithm, test: ConvergenceTest) -> Self {
        self.algorithm = algorithm;
        self.test = test;
        self
    }

    /// Adds a `GroundMotion` and re-derives the initial acceleration to
    /// include its contribution at `t=0` (usually zero — accelerograms
    /// conventionally start at `ag(0) = 0` — but not assumed). Must be
    /// called before any `step()`; re-scatters/re-commits the domain's
    /// state exactly as `new` did, so it's equivalent to having passed this
    /// motion to `new` directly.
    pub fn with_ground_motion(mut self, motion: GroundMotion) -> Self {
        let incidence = self.domain.direction_incidence(motion.direction);
        let mass_incidence = self.mass.component_mul(&incidence);
        self.ground_motions.push((motion, mass_incidence));
        self.recompute_initial_acceleration();
        self
    }

    fn recompute_initial_acceleration(&mut self) {
        let n = self.domain.num_free_dofs();
        let u0 = self.domain.gather_displacement();
        let v0 = self.domain.gather_velocity();
        let (_k, resistance0) = self.domain.assemble_tangent_and_resistance(self.time);
        let load0 = self.domain.assemble_reference_load(self.time);
        let k_v0 = self.domain.multiply_stiffness(&v0, self.time);
        let ground_force0 = self.ground_force(self.time);

        let mut a0 = DVector::<f64>::zeros(n);
        for i in 0..n {
            let c_v0_i =
                self.damping.alpha_m * self.mass[i] * v0[i] + self.damping.beta_k * k_v0[i];
            a0[i] = (load0[i] - resistance0[i] - c_v0_i + ground_force0[i]) / self.mass[i];
        }
        self.domain.scatter_state(&u0, &v0, &a0);
        // Align committed material state with the given initial condition
        // (matters if `with_initial_displacement` was used) before any
        // stepping begins — see `Material`'s doc comment.
        self.domain.commit(self.time);
    }

    /// `-M·ι·ag(t)`, summed over every active `GroundMotion` — the
    /// effective inertial force contribution to `f_eff` (`try_step`) and to
    /// the initial-acceleration equilibrium (`recompute_initial_acceleration`).
    fn ground_force(&self, time: f64) -> DVector<f64> {
        let mut force = DVector::<f64>::zeros(self.domain.num_free_dofs());
        for (motion, mass_incidence) in &self.ground_motions {
            let ag = motion.acceleration(time);
            if ag != 0.0 {
                force -= mass_incidence * ag;
            }
        }
        force
    }

    pub fn domain(&self) -> &Domain<NDIM, NDOF, NId, E> {
        &self.domain
    }

    /// Ends this (dynamic) phase and hands back the `Domain` — symmetric
    /// with `Analysis::into_domain`, so a dynamic phase can itself be
    /// followed by another phase (static or dynamic).
    pub fn into_domain(self) -> Domain<NDIM, NDOF, NId, E> {
        self.domain
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    /// Advance one Newmark step of size `dt`. External load each step is
    /// `Domain::assemble_reference_load(time)` (every `LoadPattern`'s
    /// current/frozen contribution — see its doc comment) plus each active
    /// `GroundMotion`'s effective inertial force at `time` (`ground_force`).
    ///
    /// On failure, the domain and `time`/`step_count` are restored to their
    /// state before this call — same reasoning as `Analysis::step`.
    pub fn step(&mut self) -> Result<TransientStepResult, AnalysisError> {
        let snapshot = self.domain.clone();
        let (snapshot_time, snapshot_step_count) = (self.time, self.step_count);

        let result = self.try_step();

        if result.is_err() {
            self.domain = snapshot;
            self.time = snapshot_time;
            self.step_count = snapshot_step_count;
        } else {
            self.domain.commit(self.time);
        }

        result
    }

    fn try_step(&mut self) -> Result<TransientStepResult, AnalysisError> {
        self.step_count += 1;
        self.time += self.dt;

        let (beta, gamma, dt) = (0.25, 0.5, self.dt);
        let a1 = 1.0 / (beta * dt * dt);
        let a2 = 1.0 / (beta * dt);
        let a3 = 1.0 / (2.0 * beta) - 1.0;
        let a4 = gamma / (beta * dt);
        let a5 = gamma / beta - 1.0;
        let a6 = dt * (gamma / (2.0 * beta) - 1.0);

        let u_n = self.domain.gather_displacement();
        let v_n = self.domain.gather_velocity();
        let a_n = self.domain.gather_acceleration();

        let mass_coeff = a1 + a4 * self.damping.alpha_m;
        let stiffness_coeff = 1.0 + a4 * self.damping.beta_k;
        let time = self.time;
        let external = self.domain.assemble_reference_load(time) + self.ground_force(time);

        let (iterations, factorizations) = match self.algorithm {
            // One effective-system solve at `u_n`, unconditionally
            // accepted — today's original behavior, kept byte-for-byte
            // (see `docs/algorithms.md` §6.1's load-bearing correctness
            // check: the `Newton`/`KrylovNewton` arm below reduces to
            // exactly this for a linear-material model converged in one
            // iteration — verified in `m21_transient_corrector.rs` — so
            // this arm isn't strictly *necessary*, but keeping it removes
            // any risk of this well-tested closed-form path regressing).
            Algorithm::Linear => {
                // Newmark predictor terms — see the module doc comment's
                // citation of the standard (e.g. Chopra, "Dynamics of
                // Structures") formula this whole arm derives from.
                let mass_vec = &u_n * a1 + &v_n * a2 + &a_n * a3;
                let damp_vec = &u_n * a4 + &v_n * a5 + &a_n * a6;

                let (k_eff, _resistance) = self.domain.assemble_newmark_system(
                    &self.mass,
                    mass_coeff,
                    stiffness_coeff,
                    self.time,
                );
                let k_damp_vec = self.domain.multiply_stiffness(&damp_vec, self.time);

                let n = self.domain.num_free_dofs();
                let mut f_eff = external;
                for i in 0..n {
                    f_eff[i] += self.mass[i] * (mass_vec[i] + self.damping.alpha_m * damp_vec[i])
                        + self.damping.beta_k * k_damp_vec[i];
                }

                let u_new = self.solver.solve(&k_eff, &f_eff)?;

                let delta_u = &u_new - &u_n;
                let a_new = &delta_u * a1 - &v_n * a2 - &a_n * a3;
                let v_new = &delta_u * a4 - &v_n * a5 - &a_n * a6;
                self.domain.scatter_state(&u_new, &v_new, &a_new);
                (1, 1)
            }
            // Newton corrector (`docs/algorithms.md` §6.1): iterate the
            // *effective* operator `Op = stiffness_coeff*K(u_trial) +
            // mass_coeff*M` against the residual `external - F_int(u_trial)
            // - C*v_trial - M*a_trial`, `v_trial`/`a_trial` recomputed each
            // iteration from the *cumulative* displacement change since
            // step start (`u_trial - u_n`), via the same Newmark relations
            // the `Linear` arm above uses only once, at the end.
            Algorithm::Newton { .. } | Algorithm::KrylovNewton { .. } => {
                let algorithm = self.algorithm;
                let solver = &self.solver;
                let mass = &self.mass;
                let (alpha_m, beta_k) = (self.damping.alpha_m, self.damping.beta_k);

                // `v_trial`/`a_trial` at the domain's current cumulative
                // displacement change since step start — shared by
                // `form_system` (needs them for the residual) and
                // `after_increment` (needs them to update the domain's
                // velocity/acceleration state), so it's a closure of its
                // own rather than duplicated inline in both places.
                let newmark_state = |u_trial: &DVector<f64>| {
                    let delta_u = u_trial - &u_n;
                    let v_trial = &delta_u * a4 - &v_n * a5 - &a_n * a6;
                    let a_trial = &delta_u * a1 - &v_n * a2 - &a_n * a3;
                    (v_trial, a_trial)
                };

                // `Combined`'s per-equation reference: this step's external
                // force and the committed internal force.
                let force_tolerance = matches!(self.test, ConvergenceTest::Combined { .. })
                    .then(|| {
                        let (_k, resistance) =
                            self.domain.assemble_tangent_and_resistance(self.time);
                        ForceTolerance::for_test(
                            &self.test,
                            &self.domain.rotational_equations(),
                            &[&external, &resistance],
                        )
                    })
                    .flatten();

                let (outcome, _) = iterate_to_equilibrium(
                    self.step_count,
                    &algorithm,
                    &self.test,
                    solver,
                    &mut self.cached_factorization,
                    &mut self.domain,
                    0.0,
                    |domain, _scalar| {
                        let u_trial = domain.gather_displacement();
                        let (v_trial, a_trial) = newmark_state(&u_trial);
                        let (k_eff, resistance) =
                            domain.assemble_newmark_system(mass, mass_coeff, stiffness_coeff, time);
                        let k_v = domain.multiply_stiffness(&v_trial, time);

                        let n = resistance.len();
                        let mut residual = external.clone();
                        for i in 0..n {
                            residual[i] -= resistance[i]
                                + alpha_m * mass[i] * v_trial[i]
                                + beta_k * k_v[i]
                                + mass[i] * a_trial[i];
                        }
                        (k_eff, residual)
                    },
                    // No `Integrator::DisplacementControl`-style corrector
                    // here — `TransientAnalysis` has no solved-for pseudo-
                    // time/load-factor scalar at all (`time` is fixed for
                    // the whole step by `dt`), so this hook is unused.
                    |_domain, _k, _current_factorization, du, _scalar, _iteration| Ok((du, 0.0)),
                    |domain| {
                        let u_trial = domain.gather_displacement();
                        let (v_trial, a_trial) = newmark_state(&u_trial);
                        domain.scatter_state(&u_trial, &v_trial, &a_trial);
                    },
                    force_tolerance.as_ref(),
                )?;
                (outcome.iterations, outcome.factorizations)
            }
        };

        Ok(TransientStepResult {
            step: self.step_count,
            time: self.time,
            iterations,
            factorizations,
        })
    }
}
