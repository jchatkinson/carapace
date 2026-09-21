use nalgebra::DVector;

use crate::model::Domain;

use super::{AnalysisError, RayleighDamping, SparseSolver};

#[derive(Debug, Clone, Copy)]
pub struct TransientStepResult {
    pub step: usize,
    pub time: f64,
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
/// genuinely linear elastic response (this milestone's scope: `Truss`/
/// `ZeroLength`/`ElasticBeamColumn` with linear materials) since the
/// tangent stiffness formed once per step doesn't change with displacement
/// — one linear solve per step, no Newton iteration, mirroring `Algorithm::
/// Linear`'s same caveat for the static case. A Newton-iterated corrector
/// (needed once nonlinear materials are driven dynamically) is a natural
/// extension, not built until something needs it.
pub struct TransientAnalysis {
    domain: Domain,
    mass: DVector<f64>,
    damping: RayleighDamping,
    dt: f64,
    solver: SparseSolver,
    time: f64,
    step_count: usize,
}

impl TransientAnalysis {
    /// Builds a transient analysis from a domain whose nodes may already
    /// carry initial conditions (`Node::with_initial_displacement`/
    /// `with_initial_velocity`). Computes the one thing Newmark needs that
    /// isn't a direct initial condition: a consistent initial acceleration
    /// from equilibrium at t=0, `M*a0 = F0 - K*u0 - C*v0`. `M` is diagonal,
    /// so this is a single elementwise division, no solve needed.
    pub fn new(mut domain: Domain, damping: RayleighDamping, dt: f64) -> Result<Self, AnalysisError> {
        assert!(dt > 0.0, "Newmark dt must be positive");

        let mass = domain.assemble_mass_diagonal();
        let n = domain.num_free_dofs();
        for i in 0..n {
            if mass[i] <= 0.0 {
                return Err(AnalysisError::SingularSystem);
            }
        }

        let u0 = domain.gather_displacement();
        let v0 = domain.gather_velocity();
        let (_k, resistance0) = domain.assemble_tangent_and_resistance();
        let load0 = domain.assemble_reference_load();
        let k_v0 = domain.multiply_stiffness(&v0);

        let mut a0 = DVector::<f64>::zeros(n);
        for i in 0..n {
            let c_v0_i = damping.alpha_m * mass[i] * v0[i] + damping.beta_k * k_v0[i];
            a0[i] = (load0[i] - resistance0[i] - c_v0_i) / mass[i];
        }
        domain.scatter_state(&u0, &v0, &a0);
        // Align committed material state with the given initial condition
        // (matters if `with_initial_displacement` was used) before any
        // stepping begins — see `Material`'s doc comment.
        domain.commit();

        Ok(TransientAnalysis {
            domain,
            mass,
            damping,
            dt,
            solver: SparseSolver::new(),
            time: 0.0,
            step_count: 0,
        })
    }

    pub fn domain(&self) -> &Domain {
        &self.domain
    }

    pub fn time(&self) -> f64 {
        self.time
    }

    /// Advance one Newmark step of size `dt`. Constant external load
    /// (`Domain::assemble_reference_load`) throughout — a time-varying load
    /// (e.g. ground motion) is out of scope until something needs it (§4.3
    /// only asks for Newmark integration + Rayleigh damping at M6).
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
            self.domain.commit();
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

        // Newmark predictor terms — see the module doc comment's citation
        // of the standard (e.g. Chopra, "Dynamics of Structures") formula
        // this whole step derives from.
        let mass_vec = &u_n * a1 + &v_n * a2 + &a_n * a3;
        let damp_vec = &u_n * a4 + &v_n * a5 + &a_n * a6;

        let mass_coeff = a1 + a4 * self.damping.alpha_m;
        let stiffness_coeff = 1.0 + a4 * self.damping.beta_k;
        let (k_eff, k_damp_vec) = self
            .domain
            .assemble_newmark_system(&self.mass, mass_coeff, stiffness_coeff, &damp_vec);

        let n = self.domain.num_free_dofs();
        let mut f_eff = self.domain.assemble_reference_load();
        for i in 0..n {
            f_eff[i] += self.mass[i] * (mass_vec[i] + self.damping.alpha_m * damp_vec[i])
                + self.damping.beta_k * k_damp_vec[i];
        }

        let u_new = self.solver.solve(&k_eff, &f_eff)?;

        let delta_u = &u_new - &u_n;
        let a_new = &delta_u * a1 - &v_n * a2 - &a_n * a3;
        let v_new = &delta_u * a4 - &v_n * a5 - &a_n * a6;
        self.domain.scatter_state(&u_new, &v_new, &a_new);

        Ok(TransientStepResult {
            step: self.step_count,
            time: self.time,
        })
    }
}
