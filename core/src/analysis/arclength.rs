//! Arc-length continuation (`Integrator::ArcLength`) — the design in
//! `docs/arclength.md`: a scaled spherical constraint, an oriented
//! predictor (secant by default, bordered tangent when asked or after a
//! cutback or sharp turn), full Newton correction of displacement *and*
//! load factor through the bordered system, coupled backtracking, and
//! bounded adaptive cutbacks — all transactional, so only an accepted step
//! ever commits material history or advances continuation state.
//!
//! With `q` the independent DOFs, `lambda` the pseudo-time/load parameter,
//! `p` the reference-load sensitivity, `W`/`beta` the metric:
//!
//! ```text
//! r(q, lambda) = P(lambda) - F_int(q)
//! g            = dq^T W dq + beta^2 dlambda^2 - s^2      (dq = q - q_n, ...)
//! [ K              -p              ] [delta q     ]   [ r ]
//! [ 2 dq^T W   2 beta^2 dlambda    ] [delta lambda] = [-g ]
//! ```

use nalgebra::DVector;
use slotmap::Key;

use crate::model::{Domain, ElementOps, NodeId, SparseMatrix};

use super::bordered::{Border, BorderedFactorization, BorderedSolver};
use super::convergence::ForceTolerance;
use super::{
    Algorithm, Analysis, AnalysisError, ArcFailure, ConvergenceTest, StepResult, TangentStrategy,
};

/// Arc-length continuation settings (`Integrator::ArcLength`). Build with
/// [`ArcLength::fixed`] or [`ArcLength::adaptive`] and adjust fields from
/// there; every field is validated when the first step runs.
#[derive(Debug, Clone, PartialEq)]
pub struct ArcLength<NId = NodeId> {
    /// Radius of the first step, in the dimensionless metric.
    pub initial_radius: f64,
    /// Cutbacks never go below this radius.
    pub min_radius: f64,
    /// Growth never goes above this radius. Setting all three radii equal
    /// gives a fixed-radius run with no retries or growth.
    pub max_radius: f64,
    /// Corrector solves the adaptive radius aims for (`docs/arclength.md`
    /// §7: `s_next = s * clamp(sqrt(target / N), 0.5, 1.5)`).
    pub target_iterations: usize,
    /// Cutbacks allowed per step before the step fails.
    pub max_retries: usize,
    /// Initial direction of travel in `lambda`.
    pub direction: ArcDirection,
    /// How the metric is scaled.
    pub scales: ArcScales,
    /// Predictor used after the first step.
    pub predictor: ArcPredictor,
    /// Orients the first tangent when it is singular (or overrides the
    /// default `direction`-based orientation when given).
    pub seed: Option<ArcSeed<NId>>,
    /// Acceptance tolerance on the constraint, `|g| / s^2`.
    pub arc_tolerance: f64,
    /// Optional acceptance tolerance on the last undamped coupled
    /// correction, relative to the radius: `||delta||_M <= tol * s`.
    pub correction_tolerance: Option<f64>,
    pub backtracking: Backtracking,
    pub stop: ArcStopCriteria<NId>,
}

/// Which way `lambda` moves on the first step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ArcDirection {
    #[default]
    Increasing,
    Decreasing,
}

/// The dimensionless continuation metric `W = diag(1/u*^2 | 1/theta*^2)`,
/// `beta = 1/lambda*` (`docs/arclength.md` §3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArcScales {
    /// Characteristic translation `displacement`, rotation `rotation`
    /// (required when the model has independent rotational DOFs) and load
    /// factor `load`, all finite and positive.
    Explicit {
        displacement: f64,
        rotation: Option<f64>,
        load: f64,
    },
    /// Derive the translation/rotation scales from the first elastic
    /// tangent: with `K v = p`, `u* = load * max |v_i|` over independent
    /// translations (`theta*` likewise over rotations), so a load-factor
    /// change of `load` and the displacement it causes weigh equally.
    /// Derived once per phase and reported in `ArcStepInfo`.
    Auto { load: f64 },
}

/// Predictor after the first step (`docs/arclength.md` §4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ArcPredictor {
    /// The normalized last accepted increment: no extra factorization and
    /// oriented by construction. The bordered tangent is still used for
    /// the step after a cutback or a sharp turn.
    #[default]
    Secant,
    /// The bordered tangent at the committed state, every step.
    Tangent,
}

/// A user-supplied first direction, in physical node/DOF terms plus a load
/// component. Normalized in the metric and used to orient the first
/// bordered tangent (not as the predictor itself).
#[derive(Debug, Clone, PartialEq)]
pub struct ArcSeed<NId = NodeId> {
    pub components: Vec<(NId, usize, f64)>,
    pub load: f64,
}

/// Coupled backtracking on the merit function (`docs/arclength.md` §5).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Backtracking {
    /// Armijo sufficient-decrease constant.
    pub armijo: f64,
    /// Smallest step fraction tried before the attempt is rejected.
    pub min_step: f64,
}

impl Default for Backtracking {
    fn default() -> Self {
        Self {
            armijo: 1e-4,
            min_step: 1.0 / 128.0,
        }
    }
}

/// Optional stop criteria (`docs/arclength.md` §7). Any one that is met
/// ends the phase after the step that met it is accepted; later `step()`
/// calls return `AnalysisError::ContinuationComplete`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ArcStopCriteria<NId = NodeId> {
    pub displacement: Option<DisplacementTarget<NId>>,
    pub load_factor: Option<LoadFactorTarget>,
    /// Stop at the first step on which `lambda` crosses zero.
    pub load_factor_zero_crossing: bool,
    pub max_chord_length: Option<f64>,
    pub max_steps: Option<usize>,
}

/// Stop when `node`'s `dof` reaches `value` (in either direction). With
/// `exact`, the bracketing step is retaken with an interpolated radius so
/// the phase ends on the target.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisplacementTarget<NId = NodeId> {
    pub node: NId,
    pub dof: usize,
    pub value: f64,
    pub exact: bool,
}

/// Stop when `lambda` reaches `value` (in either direction); `exact` as
/// for `DisplacementTarget`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LoadFactorTarget {
    pub value: f64,
    pub exact: bool,
}

/// Per-step continuation diagnostics (`StepResult::arc`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArcStepInfo {
    /// Radius of the accepted attempt.
    pub radius: f64,
    /// Radius proposed for the next step.
    pub next_radius: f64,
    /// Cutbacks before the accepted attempt.
    pub retries: usize,
    /// Corrector solves in the accepted attempt.
    pub corrector_solves: usize,
    /// Bordered factorizations across every attempt, predictors included.
    pub factorizations: usize,
    /// Force-test measure at acceptance (`||r||` for `NormUnbalance`,
    /// the largest scaled component for `Combined`).
    pub force_measure: f64,
    /// `|g| / s^2` at acceptance.
    pub arc_error: f64,
    /// Accumulated accepted chord length in the metric.
    pub chord_length: f64,
    /// Sign of `det K` near the accepted state (from the last bordered
    /// factorization of the step), when one was available.
    pub det_sign: Option<i8>,
    /// `det K` changed sign while `lambda` kept its direction: the step
    /// probably crossed a bifurcation and stayed on the primary branch.
    pub bifurcation_suspected: bool,
    /// The metric's translation scale (explicit or derived).
    pub displacement_scale: f64,
    /// The metric's rotation scale, when the model has rotational DOFs.
    pub rotation_scale: Option<f64>,
    /// Set on the step that met a stop criterion.
    pub stop: Option<ArcStop>,
}

/// Which stop criterion ended the phase, and how closely it was met.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArcStop {
    pub reason: StopReason,
    /// Whether an exact-landing retake hit the target within tolerance.
    pub landed_exactly: bool,
    /// Final value minus target (zero for criteria without a target).
    pub overshoot: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    DisplacementTarget,
    LoadFactorTarget,
    LoadFactorZeroCrossing,
    ChordLength,
    StepCount,
}

impl<NId> ArcLength<NId> {
    /// Fixed radius `radius` with every other setting at its default.
    pub fn fixed(radius: f64, scales: ArcScales) -> Self {
        Self::adaptive(radius, radius, radius, scales)
    }

    /// Adaptive radius in `[min_radius, max_radius]`, starting at
    /// `initial_radius`, with every other setting at its default.
    pub fn adaptive(
        initial_radius: f64,
        min_radius: f64,
        max_radius: f64,
        scales: ArcScales,
    ) -> Self {
        Self {
            initial_radius,
            min_radius,
            max_radius,
            target_iterations: 6,
            max_retries: 8,
            direction: ArcDirection::Increasing,
            scales,
            predictor: ArcPredictor::Secant,
            seed: None,
            arc_tolerance: 1e-8,
            correction_tolerance: None,
            backtracking: Backtracking::default(),
            stop: ArcStopCriteria {
                displacement: None,
                load_factor: None,
                load_factor_zero_crossing: false,
                max_chord_length: None,
                max_steps: None,
            },
        }
    }

    /// Checks this configuration together with the algorithm and
    /// convergence test it will run under — the same checks every arc step
    /// makes, exposed so wire decoders can reject bad input up front.
    /// Arc length needs `Newton { tangent: Current, line_search: None }`
    /// (its coupled backtracking replaces line search) and a force test
    /// (`NormUnbalance` or `Combined`); displacement/work tests alone don't
    /// certify equilibrium.
    pub fn validate_with(
        &self,
        algorithm: &Algorithm,
        test: &ConvergenceTest,
    ) -> Result<(), AnalysisError> {
        self.validate()?;
        test.validate()?;
        match algorithm {
            Algorithm::Newton {
                tangent: TangentStrategy::Current,
                line_search: None,
            } => {}
            _ => return Err(AnalysisError::InvalidOption { field: "algorithm" }),
        }
        if !matches!(
            test,
            ConvergenceTest::NormUnbalance { .. } | ConvergenceTest::Combined { .. }
        ) {
            return Err(AnalysisError::InvalidOption {
                field: "convergence",
            });
        }
        Ok(())
    }

    /// Range checks on this configuration alone.
    pub fn validate(&self) -> Result<(), AnalysisError> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        let invalid = |field| Err(AnalysisError::InvalidOption { field });
        if !positive(self.min_radius) {
            return invalid("integrator.minRadius");
        }
        if !positive(self.initial_radius) || self.initial_radius < self.min_radius {
            return invalid("integrator.initialRadius");
        }
        if !positive(self.max_radius) || self.max_radius < self.initial_radius {
            return invalid("integrator.maxRadius");
        }
        if self.target_iterations == 0 {
            return invalid("integrator.targetIterations");
        }
        if self.max_retries > 64 {
            return invalid("integrator.maxRetries");
        }
        match self.scales {
            ArcScales::Explicit {
                displacement,
                rotation,
                load,
            } => {
                if !positive(displacement) {
                    return invalid("integrator.displacementScale");
                }
                if rotation.is_some_and(|r| !positive(r)) {
                    return invalid("integrator.rotationScale");
                }
                if !positive(load) {
                    return invalid("integrator.loadScale");
                }
            }
            ArcScales::Auto { load } => {
                if !positive(load) {
                    return invalid("integrator.loadScale");
                }
            }
        }
        if !positive(self.arc_tolerance) {
            return invalid("integrator.arcTolerance");
        }
        if self.correction_tolerance.is_some_and(|t| !positive(t)) {
            return invalid("integrator.correctionTolerance");
        }
        let Backtracking { armijo, min_step } = self.backtracking;
        if !(armijo.is_finite() && armijo > 0.0 && armijo < 0.5) {
            return invalid("integrator.backtracking.armijo");
        }
        if !(min_step.is_finite() && min_step > 0.0 && min_step <= 1.0) {
            return invalid("integrator.backtracking.minStep");
        }
        if let Some(seed) = &self.seed {
            if !seed.load.is_finite() || seed.components.iter().any(|c| !c.2.is_finite()) {
                return invalid("integrator.seed");
            }
        }
        let stop = &self.stop;
        if stop
            .displacement
            .as_ref()
            .is_some_and(|t| !t.value.is_finite())
        {
            return invalid("integrator.stop.displacement");
        }
        if stop.load_factor.is_some_and(|t| !t.value.is_finite()) {
            return invalid("integrator.stop.loadFactor");
        }
        if stop.max_chord_length.is_some_and(|c| !positive(c)) {
            return invalid("integrator.stop.maxChordLength");
        }
        if stop.max_steps == Some(0) {
            return invalid("integrator.stop.maxSteps");
        }
        Ok(())
    }

    /// Whether switching from `self` to `other` keeps the metric and load
    /// definition, so accepted continuation history stays valid. Radius
    /// bounds, retry/iteration budgets, tolerances and stop criteria may
    /// change freely.
    pub(crate) fn same_path_definition(&self, other: &Self) -> bool
    where
        NId: PartialEq,
    {
        self.scales == other.scales
            && self.direction == other.direction
            && self.seed == other.seed
            && self.predictor == other.predictor
    }
}

/// Accepted continuation history for one phase (`docs/arclength.md` §4/§7),
/// snapshotted with the domain so a failed step restores it exactly.
#[derive(Debug, Clone)]
pub(crate) struct ArcState {
    weights: DVector<f64>,
    beta: f64,
    /// `C` for the bordered system: `u*`/`theta*` per equation, then
    /// `lambda*`.
    unknown_scale: DVector<f64>,
    displacement_scale: f64,
    rotation_scale: Option<f64>,
    rotational: Vec<bool>,
    /// Metric-normalized last accepted increment (`n+1` components).
    direction: Option<DVector<f64>>,
    /// Metric-normalized seed, consumed by the first step.
    seed: Option<DVector<f64>>,
    last_delta_lambda: f64,
    next_radius: f64,
    chord_length: f64,
    accepted_steps: usize,
    force_tangent_next: bool,
    det_sign: Option<i8>,
    stopped: bool,
}

/// One attempt's accepted state, before it is committed.
struct Accepted {
    /// Tangent at the accepted state, for the determinant-sign fallback.
    k: SparseMatrix,
    delta_lambda: f64,
    corrector_solves: usize,
    damped: bool,
    force_measure: f64,
    arc_error: f64,
    det_sign: Option<i8>,
    increment: DVector<f64>,
}

/// Everything an attempt needs that stays fixed for one step.
struct StepContext<'a> {
    q_n: &'a DVector<f64>,
    lambda_n: f64,
    p: &'a DVector<f64>,
    external_n: &'a DVector<f64>,
    internal_n: &'a DVector<f64>,
    predictor: &'a DVector<f64>,
}

impl ArcState {
    /// Called when `set_integrator` installs an `ArcLength` with the same
    /// path definition: keeps the metric, direction, radius and chord
    /// length, but clears a previous stop and restarts the step count, so
    /// new stop criteria can continue the same path.
    pub(crate) fn reopen(&mut self) {
        self.stopped = false;
        self.accepted_steps = 0;
    }

    fn metric_dot(&self, a: &DVector<f64>, b: &DVector<f64>) -> f64 {
        let n = self.weights.len();
        let mut sum = self.beta * self.beta * a[n] * b[n];
        for i in 0..n {
            sum += self.weights[i] * a[i] * b[i];
        }
        sum
    }

    fn metric_norm(&self, a: &DVector<f64>) -> f64 {
        self.metric_dot(a, a).sqrt()
    }
}

/// `true` when continuation history must be discarded because the
/// integrator change alters the path definition (or leaves arc length).
pub(crate) fn integrator_change_invalidates<NId: PartialEq>(
    old: &super::Integrator<NId>,
    new: &super::Integrator<NId>,
) -> bool {
    match (old, new) {
        (super::Integrator::ArcLength(a), super::Integrator::ArcLength(b)) => {
            !a.same_path_definition(b)
        }
        _ => true,
    }
}

/// Dense `(n+1)` concatenation `(q, lambda)`.
fn join(q: &DVector<f64>, lambda: f64) -> DVector<f64> {
    let n = q.len();
    DVector::from_fn(n + 1, |i, _| if i < n { q[i] } else { lambda })
}

fn unit_last(n: usize) -> DVector<f64> {
    DVector::from_fn(n + 1, |i, _| if i == n { 1.0 } else { 0.0 })
}

fn all_finite(v: &DVector<f64>) -> bool {
    v.iter().all(|x| x.is_finite())
}

impl<const NDIM: usize, const NDOF: usize, NId, E> Analysis<NDIM, NDOF, NId, E>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    /// One arc-length step (`Analysis::step` dispatches here). On any error
    /// the domain, continuation history, step count and load factor are
    /// restored exactly; only success commits material history.
    pub(crate) fn arc_step(
        &mut self,
        config: &ArcLength<NId>,
    ) -> Result<StepResult, AnalysisError> {
        if self.arc.as_ref().is_some_and(|state| state.stopped) {
            return Err(AnalysisError::ContinuationComplete);
        }
        config.validate_with(&self.algorithm, &self.test)?;

        let snapshot = self.domain.clone();
        let snapshot_arc = self.arc.clone();
        let (snapshot_step_count, snapshot_load_factor) = (self.step_count, self.load_factor);

        let result = self.try_arc_step(config);
        match result {
            Ok(step) => {
                self.domain.commit(self.load_factor);
                Ok(step)
            }
            Err(error) => {
                self.domain = snapshot;
                self.arc = snapshot_arc;
                self.step_count = snapshot_step_count;
                self.load_factor = snapshot_load_factor;
                Err(error)
            }
        }
    }

    fn try_arc_step(&mut self, config: &ArcLength<NId>) -> Result<StepResult, AnalysisError> {
        let lambda_n = self.load_factor;
        let p = self.domain.assemble_reference_load_sensitivity(lambda_n);
        let (k_n, internal_n) = self.domain.assemble_tangent_and_resistance(lambda_n);
        let external_n = self.domain.assemble_reference_load(lambda_n);
        let q_n = self.domain.gather_displacement();
        let mut factorizations = 0;

        if self.arc.is_none() {
            let state = self.initialize_arc(
                config,
                &p,
                &k_n,
                &internal_n,
                &external_n,
                &mut factorizations,
            )?;
            self.arc = Some(Box::new(state));
        }
        self.step_count += 1;
        let step = self.step_count;

        let state = self.arc.as_ref().expect("initialized above").clone();
        let n = q_n.len();

        // Predictor direction for this step (independent of the radius, so
        // computed once and shared by every attempt).
        let tangent_wanted = state.direction.is_none()
            || config.predictor == ArcPredictor::Tangent
            || state.force_tangent_next;
        let predictor = if tangent_wanted {
            let t_prev = match (&state.direction, &state.seed) {
                (Some(direction), _) => direction.clone(),
                (None, Some(seed)) => seed.clone(),
                (None, None) => {
                    let sign = match config.direction {
                        ArcDirection::Increasing => 1.0,
                        ArcDirection::Decreasing => -1.0,
                    };
                    unit_last(n) * (sign / state.beta)
                }
            };
            factorizations += 1;
            match bordered_tangent(&mut self.bordered, &state, &k_n, &p, &t_prev) {
                Ok((t, _det_at_committed_state)) => t,
                Err(failure) => match &state.direction {
                    Some(direction) => direction.clone(),
                    None if state.seed.is_none() && failure == ArcFailure::SingularSystem => {
                        return Err(AnalysisError::MissingSeedDirection)
                    }
                    None => {
                        return Err(AnalysisError::CutbacksExhausted {
                            step,
                            attempts: 0,
                            radius: state.next_radius,
                            last_failure: failure,
                        })
                    }
                },
            }
        } else {
            state
                .direction
                .clone()
                .expect("secant needs a previous step")
        };

        let context = StepContext {
            q_n: &q_n,
            lambda_n,
            p: &p,
            external_n: &external_n,
            internal_n: &internal_n,
            predictor: &predictor,
        };

        let start = self.domain.clone();
        let mut radius = state
            .next_radius
            .clamp(config.min_radius, config.max_radius);
        let mut retries = 0;
        let mut attempts = 0;
        let accepted = loop {
            attempts += 1;
            let outcome = attempt(
                &mut self.domain,
                &mut self.bordered,
                &state,
                config,
                &self.test,
                &context,
                radius,
                &mut factorizations,
            );
            match outcome {
                Ok(accepted) => break accepted,
                Err(failure) => {
                    self.domain = start.clone();
                    let exhausted = retries >= config.max_retries || radius <= config.min_radius;
                    if exhausted {
                        return Err(AnalysisError::CutbacksExhausted {
                            step,
                            attempts,
                            radius,
                            last_failure: failure,
                        });
                    }
                    radius = (radius * 0.5).max(config.min_radius);
                    retries += 1;
                }
            }
        };

        let mut accepted = accepted;
        let mut radius_used = radius;
        self.load_factor = lambda_n + accepted.delta_lambda;

        // Stop criteria, with optional exact landing on a bracketed target.
        let mut stop = None;
        let stop_criteria = &config.stop;
        if let Some(target) = stop_criteria.displacement {
            let eq_value = |domain: &Domain<NDIM, NDOF, NId, E>| {
                domain.node(target.node).displacement[target.dof]
            };
            let before = start.node(target.node).displacement[target.dof];
            let after = eq_value(&self.domain);
            if brackets(before, after, target.value) {
                let landed = target.exact
                    && self.land_exactly(
                        &start,
                        &state,
                        config,
                        &context,
                        &mut accepted,
                        &mut radius_used,
                        &mut factorizations,
                        before,
                        after,
                        target.value,
                        &eq_value,
                    );
                stop = Some(ArcStop {
                    reason: StopReason::DisplacementTarget,
                    landed_exactly: landed,
                    overshoot: eq_value(&self.domain) - target.value,
                });
            }
        }
        if stop.is_none() {
            if let Some(target) = stop_criteria.load_factor {
                let after = self.load_factor;
                if brackets(lambda_n, after, target.value) {
                    let landed = target.exact
                        && self.land_exactly_on_load(
                            &start,
                            &state,
                            config,
                            &context,
                            &mut accepted,
                            &mut radius_used,
                            &mut factorizations,
                            after,
                            target.value,
                        );
                    stop = Some(ArcStop {
                        reason: StopReason::LoadFactorTarget,
                        landed_exactly: landed,
                        overshoot: self.load_factor - target.value,
                    });
                }
            }
        }
        if stop.is_none()
            && stop_criteria.load_factor_zero_crossing
            && lambda_n != 0.0
            && (self.load_factor == 0.0 || self.load_factor.signum() != lambda_n.signum())
        {
            stop = Some(ArcStop {
                reason: StopReason::LoadFactorZeroCrossing,
                landed_exactly: false,
                overshoot: self.load_factor,
            });
        }

        // Commit continuation history for the accepted increment.
        let norm = state.metric_norm(&accepted.increment);
        let new_direction = &accepted.increment / norm;
        let sharp_turn = state
            .direction
            .as_ref()
            .is_some_and(|old| state.metric_dot(old, &new_direction) < SHARP_TURN_COS);
        // From the step's last corrector factorization; when the predictor
        // was accepted outright, factor once at the accepted state (a sign
        // carried over from earlier steps would hide a crossing).
        let det_sign = accepted.det_sign.or_else(|| {
            factorizations += 1;
            accepted_det_sign(&mut self.bordered, &accepted.k, &p)
        });
        let state = self.arc.as_mut().expect("initialized above");
        let bifurcation_suspected = matches!(
            (state.det_sign, det_sign),
            (Some(old), Some(new)) if old != 0 && new != 0 && old != new
        ) && state.last_delta_lambda != 0.0
            && accepted.delta_lambda.signum() == state.last_delta_lambda.signum();

        let grow = retries == 0 && !accepted.damped && !sharp_turn;
        let factor = (config.target_iterations as f64 / accepted.corrector_solves.max(1) as f64)
            .sqrt()
            .clamp(0.5, 1.5);
        let factor = if grow { factor } else { factor.min(1.0) };
        state.next_radius = (radius_used * factor).clamp(config.min_radius, config.max_radius);
        state.direction = Some(new_direction);
        state.seed = None;
        state.last_delta_lambda = accepted.delta_lambda;
        state.chord_length += norm;
        state.accepted_steps += 1;
        state.force_tangent_next = retries > 0 || sharp_turn;
        state.det_sign = det_sign;

        if stop.is_none()
            && stop_criteria
                .max_chord_length
                .is_some_and(|c| state.chord_length >= c)
        {
            stop = Some(ArcStop {
                reason: StopReason::ChordLength,
                landed_exactly: false,
                overshoot: 0.0,
            });
        }
        if stop.is_none()
            && stop_criteria
                .max_steps
                .is_some_and(|m| state.accepted_steps >= m)
        {
            stop = Some(ArcStop {
                reason: StopReason::StepCount,
                landed_exactly: false,
                overshoot: 0.0,
            });
        }
        state.stopped = stop.is_some();

        let info = super::ArcStepInfo {
            radius: radius_used,
            next_radius: state.next_radius,
            retries,
            corrector_solves: accepted.corrector_solves,
            factorizations,
            force_measure: accepted.force_measure,
            arc_error: accepted.arc_error,
            chord_length: state.chord_length,
            det_sign,
            bifurcation_suspected,
            displacement_scale: state.displacement_scale,
            rotation_scale: state.rotation_scale,
            stop,
        };
        Ok(StepResult {
            step,
            load_factor: self.load_factor,
            iterations: accepted.corrector_solves,
            factorizations,
            arc: Some(info),
        })
    }

    /// Validates the load definition and entry equilibrium, builds the
    /// metric, and resolves the seed (`docs/arclength.md` §3–§4).
    fn initialize_arc(
        &mut self,
        config: &ArcLength<NId>,
        p: &DVector<f64>,
        k_n: &SparseMatrix,
        internal_n: &DVector<f64>,
        external_n: &DVector<f64>,
        factorizations: &mut usize,
    ) -> Result<ArcState, AnalysisError> {
        if self.domain.has_unfrozen_path_series() {
            return Err(AnalysisError::UnsupportedLoadSeries);
        }
        if p.amax() == 0.0 {
            return Err(AnalysisError::ZeroLoadSensitivity);
        }
        let n = p.len();
        let rotational = self.domain.rotational_equations();
        let has_rotation = rotational.iter().any(|&r| r);
        let has_translation = rotational.iter().any(|&r| !r);

        let residual = external_n - internal_n;
        let (balanced, measure) = force_check(
            &self.test,
            &rotational,
            &residual,
            &[external_n, internal_n],
        );
        if !balanced {
            return Err(AnalysisError::InitialStateNotInEquilibrium { measure });
        }

        let (displacement_scale, rotation_scale, load_scale) = match config.scales {
            ArcScales::Explicit {
                displacement,
                rotation,
                load,
            } => {
                if has_rotation && rotation.is_none() {
                    return Err(AnalysisError::InvalidOption {
                        field: "integrator.rotationScale",
                    });
                }
                (displacement, rotation.filter(|_| has_rotation), load)
            }
            ArcScales::Auto { load } => {
                // `K v = p` through the bordered system with a pure-load
                // row: x_lambda = 1, K x_q = p.
                let zeros = DVector::zeros(n);
                let minus_p = -p;
                *factorizations += 1;
                let factor = self
                    .bordered
                    .factor(
                        k_n,
                        Border {
                            column: &minus_p,
                            row: &zeros,
                            corner: 1.0,
                        },
                        &DVector::from_element(n + 1, 1.0),
                    )
                    .map_err(|_| AnalysisError::SingularSystem)?;
                let v = factor
                    .solve(&unit_last(n))
                    .map_err(|_| AnalysisError::SingularSystem)?;
                let largest = |want_rotation: bool| {
                    (0..n)
                        .filter(|&i| rotational[i] == want_rotation)
                        .map(|i| v[i].abs())
                        .fold(0.0, f64::max)
                };
                let displacement = load * largest(false);
                let rotation = has_rotation.then(|| load * largest(true));
                let usable = |s: f64| s.is_finite() && s > 0.0;
                if has_translation && !usable(displacement) {
                    return Err(AnalysisError::InvalidOption {
                        field: "integrator.displacementScale",
                    });
                }
                if rotation.is_some_and(|r| !usable(r)) {
                    return Err(AnalysisError::InvalidOption {
                        field: "integrator.rotationScale",
                    });
                }
                (
                    if has_translation { displacement } else { 1.0 },
                    rotation,
                    load,
                )
            }
        };

        let unknown_scale = DVector::from_fn(n + 1, |i, _| {
            if i == n {
                load_scale
            } else if rotational[i] {
                rotation_scale.expect("rotation scale required above")
            } else {
                displacement_scale
            }
        });
        let weights = DVector::from_fn(n, |i, _| unknown_scale[i].powi(-2));
        let beta = 1.0 / load_scale;

        let mut state = ArcState {
            weights,
            beta,
            unknown_scale,
            displacement_scale,
            rotation_scale,
            rotational,
            direction: None,
            seed: None,
            last_delta_lambda: 0.0,
            next_radius: config.initial_radius,
            chord_length: 0.0,
            accepted_steps: 0,
            force_tangent_next: false,
            det_sign: None,
            stopped: false,
        };

        if let Some(seed) = &config.seed {
            let mut q = DVector::zeros(n);
            for &(node, dof, value) in &seed.components {
                if dof >= NDOF {
                    return Err(AnalysisError::InvalidOption {
                        field: "integrator.seed",
                    });
                }
                let eq =
                    self.domain
                        .equation_of(node, dof)
                        .ok_or(AnalysisError::InvalidOption {
                            field: "integrator.seed",
                        })?;
                q[eq] += value;
            }
            let direction = join(&q, seed.load);
            let norm = state.metric_norm(&direction);
            if !(norm.is_finite() && norm > 0.0) {
                return Err(AnalysisError::InvalidOption {
                    field: "integrator.seed",
                });
            }
            state.seed = Some(direction / norm);
        }
        Ok(state)
    }

    /// Retakes the bracketing step with interpolated radii (regula falsi on
    /// the radius) until `value_of` hits `target`, keeping the original
    /// bracketing step if that doesn't succeed. Returns whether it landed.
    #[allow(clippy::too_many_arguments)]
    fn land_exactly(
        &mut self,
        start: &Domain<NDIM, NDOF, NId, E>,
        state: &ArcState,
        config: &ArcLength<NId>,
        context: &StepContext<'_>,
        accepted: &mut Accepted,
        radius_used: &mut f64,
        factorizations: &mut usize,
        before: f64,
        after: f64,
        target: f64,
        value_of: &impl Fn(&Domain<NDIM, NDOF, NId, E>) -> f64,
    ) -> bool {
        self.retake_until(
            start,
            state,
            config,
            context,
            accepted,
            radius_used,
            factorizations,
            before,
            after,
            target,
            &|domain, _lambda| value_of(domain),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn land_exactly_on_load(
        &mut self,
        start: &Domain<NDIM, NDOF, NId, E>,
        state: &ArcState,
        config: &ArcLength<NId>,
        context: &StepContext<'_>,
        accepted: &mut Accepted,
        radius_used: &mut f64,
        factorizations: &mut usize,
        after: f64,
        target: f64,
    ) -> bool {
        let before = context.lambda_n;
        self.retake_until(
            start,
            state,
            config,
            context,
            accepted,
            radius_used,
            factorizations,
            before,
            after,
            target,
            &|_domain, lambda| lambda,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn retake_until(
        &mut self,
        start: &Domain<NDIM, NDOF, NId, E>,
        state: &ArcState,
        config: &ArcLength<NId>,
        context: &StepContext<'_>,
        accepted: &mut Accepted,
        radius_used: &mut f64,
        factorizations: &mut usize,
        before: f64,
        after: f64,
        target: f64,
        value_of: &impl Fn(&Domain<NDIM, NDOF, NId, E>, f64) -> f64,
    ) -> bool {
        let tolerance = LANDING_TOLERANCE * (after - before).abs().max(f64::MIN_POSITIVE);
        if (after - target).abs() <= tolerance {
            return true;
        }
        let bracket_domain = self.domain.clone();
        let bracket_lambda = self.load_factor;

        // f(s) = value(s) - target, bracketed on [0, s_used].
        let (mut s_lo, mut f_lo) = (0.0, before - target);
        let (mut s_hi, mut f_hi) = (*radius_used, after - target);
        for _ in 0..MAX_LANDING_RETAKES {
            let s = s_lo + (s_hi - s_lo) * (-f_lo) / (f_hi - f_lo);
            if !(s.is_finite() && s > 0.0) {
                break;
            }
            self.domain = start.clone();
            let mut retake_factorizations = 0;
            let Ok(retaken) = attempt(
                &mut self.domain,
                &mut self.bordered,
                state,
                config,
                &self.test,
                context,
                s,
                &mut retake_factorizations,
            ) else {
                *factorizations += retake_factorizations;
                break;
            };
            *factorizations += retake_factorizations;
            let lambda = context.lambda_n + retaken.delta_lambda;
            let f = value_of(&self.domain, lambda) - target;
            if f.abs() <= tolerance {
                self.load_factor = lambda;
                *accepted = retaken;
                *radius_used = s;
                return true;
            }
            if f.signum() == f_lo.signum() {
                s_lo = s;
                f_lo = f;
            } else {
                s_hi = s;
                f_hi = f;
            }
        }
        self.domain = bracket_domain;
        self.load_factor = bracket_lambda;
        false
    }
}

/// Turn (in the metric) beyond which the next step suppresses radius growth
/// and uses the bordered tangent: about 30 degrees.
const SHARP_TURN_COS: f64 = 0.866;
/// Orientation guard: the accepted increment must satisfy
/// `d^T M t >= ORIENTATION_EPS * s^2` (`docs/arclength.md` §4).
const ORIENTATION_EPS: f64 = 1e-6;
/// Exact landing hits the target to this fraction of the bracketing step's
/// change in the controlled quantity.
const LANDING_TOLERANCE: f64 = 1e-8;
const MAX_LANDING_RETAKES: usize = 6;
/// Iterations over which the merit function must at least halve.
const STAGNATION_WINDOW: usize = 5;

fn brackets(before: f64, after: f64, target: f64) -> bool {
    before != target && (before - target) * (after - target) <= 0.0
}

/// The force half of acceptance: `(passes, measure)`.
fn force_check(
    test: &ConvergenceTest,
    rotational: &[bool],
    residual: &DVector<f64>,
    magnitudes: &[&DVector<f64>],
) -> (bool, f64) {
    match *test {
        ConvergenceTest::NormUnbalance { tol, .. } => {
            let norm = residual.norm();
            (norm < tol, norm)
        }
        _ => {
            let tolerance = ForceTolerance::for_test(test, rotational, magnitudes)
                .expect("arc length only accepts force tests");
            let measure = tolerance.scaled_max(residual);
            (measure <= 1.0, measure)
        }
    }
}

/// The bordered tangent at the committed state, oriented by `t_prev`:
/// `[K -p; (M t_prev)^T] v = e_{n+1}`, `t = v / ||v||_M`. Also returns the
/// sign of `det K` (Cramer: `v_lambda = det K / det A`).
fn bordered_tangent(
    bordered: &mut BorderedSolver,
    state: &ArcState,
    k: &SparseMatrix,
    p: &DVector<f64>,
    t_prev: &DVector<f64>,
) -> Result<(DVector<f64>, Option<i8>), ArcFailure> {
    let n = p.len();
    let row = DVector::from_fn(n, |i, _| state.weights[i] * t_prev[i]);
    let minus_p = -p;
    let factor = bordered.factor(
        k,
        Border {
            column: &minus_p,
            row: &row,
            corner: state.beta * state.beta * t_prev[n],
        },
        &state.unknown_scale,
    )?;
    let v = factor.solve(&unit_last(n))?;
    let norm = state.metric_norm(&v);
    if !(norm.is_finite() && norm > 0.0) {
        return Err(ArcFailure::NonFinite);
    }
    let t = v / norm;
    // Positive by construction (= 1/||v||_M); anything else means the solve
    // can't be trusted.
    if state.metric_dot(&t, t_prev) <= 0.0 {
        return Err(ArcFailure::InaccurateSolve);
    }
    let det_k = if t[n] == 0.0 {
        0
    } else {
        factor.det_sign() * t[n].signum() as i8
    };
    Ok((t, Some(det_k)))
}

/// Sign of `det K` at an accepted state reached without any corrector
/// factorization: `[K -p; 0 1]` has `det = det K`.
fn accepted_det_sign(
    bordered: &mut BorderedSolver,
    k: &SparseMatrix,
    p: &DVector<f64>,
) -> Option<i8> {
    let n = p.len();
    let zeros = DVector::zeros(n);
    let minus_p = -p;
    match bordered.factor(
        k,
        Border {
            column: &minus_p,
            row: &zeros,
            corner: 1.0,
        },
        &DVector::from_element(n + 1, 1.0),
    ) {
        Ok(factor) => Some(factor.det_sign()),
        Err(ArcFailure::SingularSystem) => Some(0),
        Err(_) => None,
    }
}

/// Sign of `det K` from a corrector factorization: `x_lambda = det K / det A`
/// for `A x = e_{n+1}`.
fn det_sign_from(factor: &BorderedFactorization, n: usize) -> Option<i8> {
    let x = factor.solve(&unit_last(n)).ok()?;
    let last = x[n];
    if last == 0.0 {
        Some(0)
    } else {
        Some(factor.det_sign() * last.signum() as i8)
    }
}

/// One attempt at radius `s` from the committed state (the caller restores
/// the domain on failure): apply the predictor, then coupled Newton with
/// backtracking until force equilibrium, the constraint, and any
/// correction criteria all hold at the accepted state.
#[allow(clippy::too_many_arguments)]
fn attempt<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, NId, E>,
    bordered: &mut BorderedSolver,
    state: &ArcState,
    config: &ArcLength<NId>,
    test: &ConvergenceTest,
    context: &StepContext<'_>,
    s: f64,
    factorizations: &mut usize,
) -> Result<Accepted, ArcFailure>
where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    let n = context.q_n.len();
    let t = context.predictor;
    let s2 = s * s;

    let predictor_q = DVector::from_fn(n, |i, _| s * t[i]);
    domain.apply_displacement_increment(&predictor_q);
    let mut lambda = context.lambda_n + s * t[n];

    // Frozen for the attempt: committed and predicted external force and
    // committed internal force (`docs/arclength.md` §6).
    let external_predicted = domain.assemble_reference_load(lambda);
    let force_tolerance = match *test {
        ConvergenceTest::NormUnbalance { tol, .. } => ForceTolerance::uniform(n, tol),
        _ => ForceTolerance::for_test(
            test,
            &state.rotational,
            &[context.external_n, &external_predicted, context.internal_n],
        )
        .expect("arc length only accepts force tests"),
    };
    let force_ok = |residual: &DVector<f64>| -> (bool, f64) {
        match *test {
            ConvergenceTest::NormUnbalance { tol, .. } => {
                let norm = residual.norm();
                (norm < tol, norm)
            }
            _ => {
                let measure = force_tolerance.scaled_max(residual);
                (measure <= 1.0, measure)
            }
        }
    };
    let displacement_tol = match *test {
        ConvergenceTest::Combined {
            displacement_tol, ..
        } => displacement_tol,
        _ => None,
    };
    let merit = |residual: &DVector<f64>, g: f64| {
        0.5 * (force_tolerance.scaled_norm_squared(residual)
            + (g / (s2 * config.arc_tolerance)).powi(2))
    };
    let constraint = |dq: &DVector<f64>, dl: f64| {
        let mut sum = state.beta * state.beta * dl * dl - s2;
        for i in 0..n {
            sum += state.weights[i] * dq[i] * dq[i];
        }
        sum
    };
    let evaluate = |domain: &Domain<NDIM, NDOF, NId, E>, lambda: f64| {
        let (k, internal) = domain.assemble_tangent_and_resistance(lambda);
        let residual = domain.assemble_reference_load(lambda) - internal;
        (k, residual)
    };

    let (mut k, mut residual) = evaluate(domain, lambda);
    let mut corrector_solves = 0;
    let mut damped = false;
    let mut last_undamped: Option<DVector<f64>> = None;
    let mut last_factor: Option<BorderedFactorization> = None;
    let mut merit_history: Vec<f64> = Vec::new();
    let max_iter = test.max_iter();

    loop {
        let dq = domain.gather_displacement() - context.q_n;
        let dl = lambda - context.lambda_n;
        if !all_finite(&residual) || !all_finite(&dq) || !dl.is_finite() {
            return Err(ArcFailure::NonFinite);
        }
        let g = constraint(&dq, dl);
        let (force_pass, force_measure) = force_ok(&residual);
        let arc_error = g.abs() / s2;
        let correction_pass = match (&last_undamped, config.correction_tolerance) {
            (Some(delta), Some(tol)) => state.metric_norm(delta) <= tol * s,
            _ => true,
        };
        let displacement_pass = match (&last_undamped, displacement_tol) {
            (Some(delta), Some(tol)) => delta.rows(0, n).amax() <= tol,
            _ => true,
        };

        if force_pass && arc_error <= config.arc_tolerance && correction_pass && displacement_pass {
            let increment = join(&dq, dl);
            if state.metric_dot(&increment, t) < ORIENTATION_EPS * s2 {
                return Err(ArcFailure::BranchOrientation);
            }
            let det_sign = last_factor.as_ref().and_then(|f| det_sign_from(f, n));
            return Ok(Accepted {
                k,
                delta_lambda: dl,
                corrector_solves,
                damped,
                force_measure,
                arc_error,
                det_sign,
                increment,
            });
        }
        if corrector_solves >= max_iter {
            return Err(ArcFailure::NotConverged);
        }
        if dq.amax() == 0.0 && dl == 0.0 {
            return Err(ArcFailure::DegenerateConstraint);
        }

        let minus_p = -context.p;
        let row = DVector::from_fn(n, |i, _| 2.0 * state.weights[i] * dq[i]);
        *factorizations += 1;
        let factor = bordered.factor(
            &k,
            Border {
                column: &minus_p,
                row: &row,
                corner: 2.0 * state.beta * state.beta * dl,
            },
            &state.unknown_scale,
        )?;
        let delta = factor.solve(&join(&residual, -g))?;
        corrector_solves += 1;

        // Coupled backtracking from this iteration's base state. Each
        // candidate `base + eta * delta` is projected radially back onto
        // the sphere, so every trial satisfies the constraint exactly and
        // the merit function measures force imbalance. The correction is
        // tangent to the sphere (its constraint row reads `2 dq^T W dq_c +
        // 2 beta^2 dl dl_c = -g ~ 0`), so the projection only changes the
        // step at second order and `delta` stays a descent direction.
        let phi0 = merit(&residual, g);
        let base = join(&dq, dl);
        let mut eta = 1.0;
        let accepted_candidate = loop {
            let candidate = &base + eta * &delta;
            let norm = state.metric_norm(&candidate);
            if norm.is_finite() && norm > 0.0 {
                let target = candidate * (s / norm);
                let current = domain.gather_displacement() - context.q_n;
                let shift = DVector::from_fn(n, |i, _| target[i] - current[i]);
                domain.apply_displacement_increment(&shift);
                let candidate_lambda = context.lambda_n + target[n];
                let (k_c, r_c) = evaluate(domain, candidate_lambda);
                let dq_c = domain.gather_displacement() - context.q_n;
                let g_c = constraint(&dq_c, candidate_lambda - context.lambda_n);
                let phi = merit(&r_c, g_c);
                let converged_here = force_ok(&r_c).0 && g_c.abs() / s2 <= config.arc_tolerance;
                if phi.is_finite()
                    && (phi <= (1.0 - 2.0 * config.backtracking.armijo * eta) * phi0
                        || converged_here)
                {
                    break Some((k_c, r_c, candidate_lambda, phi));
                }
            }
            eta *= 0.5;
            if eta < config.backtracking.min_step {
                break None;
            }
        };
        let Some((k_c, r_c, candidate_lambda, phi)) = accepted_candidate else {
            return Err(ArcFailure::NoDescent);
        };
        damped |= eta < 1.0;
        k = k_c;
        residual = r_c;
        lambda = candidate_lambda;
        last_undamped = Some(delta);
        last_factor = Some(factor);

        merit_history.push(phi);
        let len = merit_history.len();
        if len > STAGNATION_WINDOW && phi > 0.5 * merit_history[len - 1 - STAGNATION_WINDOW] {
            return Err(ArcFailure::Stagnated);
        }
    }
}
