# Arc-length continuation implementation plan

Status: implemented — `core::analysis::Integrator::ArcLength` (native) and
`integrator: { kind: "arcLength", ... }` (wire, both profiles). This
document remains the design record. §0 lists where the implementation
departs from or narrows the original plan; the sections below have been
updated to match the code.

## 0. Implementation status and deviations

Code: `core/src/analysis/arclength.rs` (configuration, continuation
driver, attempt/corrector), `bordered.rs` (bordered solves),
`convergence.rs` (accepted-state context, `Combined`), `newton_loop.rs`
(post-update checks), `state.rs` (dispatch, history invalidation). Wire:
`wasm-bridge/src/input_v1/sequence.rs` (`ArcLengthSpec`,
`ConvergenceSpec::Combined`), `decode.rs`, `session.rs`
(`StepOutcome::continuation`, structured errors).

Departures from the plan, each for a reason found while implementing:

- **Radial projection in the corrector** (§5). Without it the kink step
  stalled at every radius.
- **Block elimination by default, direct bordered LU as fallback** (§5).
  The direct factorization filled in about 10x on a 40x8 frame.
- **The determinant sign** comes from the step's last corrector
  factorization. When the predictor is accepted with no corrector solve,
  one extra factorization at the accepted state supplies it, because a sign
  carried over from earlier steps would hide a crossing (§7).
- **`Combined` convergence** (§6) folds the "characteristic force/moment
  scales" into the absolute tolerances (`forceTol`, `momentTol`).
  `F_ref_i` is the largest of the committed and predicted external force
  and the committed internal force. The optional **work criterion is not
  implemented**; its displacement criterion bounds the largest component
  of the last undamped correction. The arc-length-specific coupled
  correction criterion is `ArcLength::correction_tolerance`, relative to
  `s`. On the wire, `momentTol` defaults to `forceTol`, `relativeTol` to
  `1e-6`, and an arc-length stage that omits `convergence` gets
  `combined { forceTol: 1e-6, maxIter: 30 }`.
- **Shared Newton loop:** a zero iteration budget is now
  `AnalysisError::InvalidOption { field: "convergence.maxIter" }`, as the
  wire decoder already required. Two iteration counts change
  intentionally. A linear system converges on its one real solve, with no
  confirming second iteration, and an already-balanced predictor is
  accepted with no solve under the force tests.
- **Stop criteria:** `max_steps` is native-only, since a wire stage's
  `steps` is its cap. After a stop, `set_integrator` with the same path
  definition (scales, direction, seed, predictor) reopens the phase,
  keeping direction, radius and accumulated chord length and restarting
  the step count. Exact landing retakes the bracketing step with up to six
  regula-falsi radii.
- **Fixed thresholds:** sharp turn `cos < 0.866` (30 degrees),
  orientation `eps_orient = 1e-6`, landing tolerance `1e-8` of the
  bracketing step's change, Schur cancellation `1e-8`.

Tests: `core/tests/arclength_softening.rs` (§8: softening peak, load-control
baseline and rollback, series snap-back, reverse direction),
`core/tests/arclength_continuation.rs` (§9), `wasm-bridge/tests/
arclength_session.rs` and the arc-length section of `wasm-bridge/tests/
boundary-smoke.ts` (wire, session ordering, both profiles), and unit tests
in `bordered.rs`/`convergence.rs`.

Not yet covered by tests:

- continuation through an affine 3D rigid diaphragm (only `equal_dof`
  reduced coordinates are tested);
- a Pinching4 trial-isolation test (isolation is shown bit-for-bit with
  Hysteretic, after failed attempts);
- the follow-up cyclic-damage protocol;
- an asymmetric-bifurcation two-bar variant. The bifurcation warning is
  tested instead on a straight corotational cantilever passing its Euler
  load, since a pin-ended von Mises truss has no asymmetric bifurcation
  without added support flexibility.

## 1. Scope and recommended method

Add a static `Integrator::ArcLength` using a **scaled spherical constraint,
oriented predictor, full Newton correction of displacement and load factor,
and bounded adaptive step retries**. Solve the bordered system by block
elimination, with a direct bordered factorization as the fallback for a
singular `K` (§5). This avoids requiring an invertible structural tangent
at every point and keeps the method usable through ordinary load limits,
negative stiffness, snap-through, and snap-back.

This is equilibrium-path continuation, including unstable equilibria; it
does not predict the dynamic jump of a real structure. Target the existing
conservative nodal/element loads, material nonlinearity, and supported
geometric nonlinearity in both execution profiles. Cyclic loading still
requires a prescribed loading/reversal protocol: continuation direction alone
does not specify such a protocol. Bifurcation branch switching, follower
loads, nonhomogeneous prescribed displacement continuation, and a
dissipation-based constraint are later extensions. Abrupt complete failure
may destroy the continuous equilibrium path; report that failure rather
than regularizing stiffness silently.

The spherical constraint and augmented correction have established
precedents in [OpenSees's ArcLength class](https://opensees.berkeley.edu/OpenSees/html/classArcLength.html)
and the [FEniCS-arclength formulation](https://fenics-arclength.readthedocs.io/en/latest/math_prelim.html).
The scaling, state lifecycle, globalization, and API choices below are
proposals specific to Carapace. Do not assume matching parameter names imply
numerical equivalence with another implementation.

## 2. Existing code and required changes

| Location | Relevant behavior / planned change |
| --- | --- |
| `core/src/analysis/integrator.rs` | Currently a `Copy` configuration enum whose predictor returns only pseudo-time. Add arc-length configuration; keep vector/history state separately. Arc length must apply a displacement predictor too. |
| `core/src/analysis/state.rs` | `Analysis::step()` snapshots the domain and commits only on success. Extend the transaction to continuation history, radius, load factor, and retries. |
| `core/src/analysis/newton_loop.rs` | Currently solves a displacement system before the integrator corrector; cannot represent a singular-tangent bordered solve. Add a correction-system hook rather than forcing arc length through `Integrator::correct()`. |
| `core/src/analysis/convergence.rs` | Currently sees only `residual` and `du`; the ordinary check uses the residual assembled **before** the applied correction. Introduce an accepted-state convergence context and an all-required-tests criterion. |
| `core/src/analysis/solver.rs` | Sparse LU can handle indefinite/nonsymmetric operators. Build an `(n+1)` bordered sparse matrix; verify singular and nonfinite solve outcomes and backward error. |
| `core/src/model/domain.rs` | Reuse transformed assembly and displacement application. Add reduced-coordinate displacement access, scaling metadata, and inspection of active load-series kinds. |
| `core/src/analysis/builder.rs`, `error.rs`, `mod.rs` | Initialize/export continuation settings and history, validate combinations, and add structured failures. |
| `wasm-bridge/src/input_v1/sequence.rs`, `decode/stages.rs` | Add configuration and convergence wire variants and validation for both profiles (one shared stage compiler). |
| `wasm-bridge/src/input_v1/session.rs` | Accepted steps alone produce samples; expose useful retry/continuation diagnostics and preserve sample ordering when load factor decreases. |

Preserve the existing load/displacement integrator API. Remove the source
documentation assumption that all integrators only predict a fixed scalar.

## 3. Equations, loads, and units

Let `q` be independent DOFs after constraint transformation, `lambda` the
existing static pseudo-time/load parameter, and `h_n` committed constitutive
history. Within one attempted step, evaluate every trial from the same
`h_n`:

```text
P(lambda) = P_const + lambda p
r(q, lambda; h_n) = P(lambda) - F_int(q; h_n)
K = dF_int/dq
Delta q = q - q_n                 Delta lambda = lambda - lambda_n
g = Delta q^T W Delta q + beta^2 Delta lambda^2 - s^2
```

`p` is `assemble_reference_load_sensitivity()` for the combined active
linear patterns, including pattern scale factors and transformed nodal and
equivalent element loads. Constant and frozen patterns belong in `P_const`,
including gravity. First release accepts active `Linear`, `Constant`, and
frozen patterns only; reject any unfrozen `Path` series rather than assuming
its current slope remains valid across a reversal or breakpoint. Reject zero
combined reduced sensitivity, even if individual patterns are nonzero.
Always assemble actual `P(lambda)` for residuals and reactions.

Require equilibrium at entry to an arc-length phase. Check it using the
force criterion below; if not satisfied, return an initialization error
with the residual norm. Do not hide an initial fixed-load equilibration,
which could jump branches. A gravity phase must first converge and freeze
gravity. The new phase's initial load parameter follows the existing stage
reset convention; unfrozen preloads must be represented accordingly.

Use a dimensionless metric. For independent translations, `W_ii = 1/u_*^2`;
for rotations, `W_ii = 1/theta_*^2`; `beta = 1/lambda_*`, with all three
scales finite and strictly positive. Require explicit model scales; do not
mix lengths and radians in an unweighted Euclidean norm or guess units.
`s` is consequently dimensionless. Users can omit the rotation scale only
when there are no independent rotational DOFs. Allow per-equation overrides
later if needed for highly heterogeneous models.

A poorly chosen `beta` is the most common practical failure of spherical
continuation: a load term that dominates the metric makes the method behave
like load control near limit points. Offer an optional model-derived
auto-scale as an alternative to explicit scales, so users need not tune
them by hand. With `lambda_*` given, solve the first elastic tangent
`K v = p` and set `u_*` from the
translational part of `v` (for example `u_* = lambda_* * max_i |v_i|` over
independent translations, and `theta_*` likewise over rotations). This
derives scales from the model rather than guessing units, so it does not
violate the rule above. Record the derived scales in diagnostics, keep them
fixed for the phase, and fail initialization if the first tangent is
singular or the derived scale is zero or nonfinite.

Initially use the independent-coordinate metric: alias DOFs contribute once,
and a retained diaphragm rotation receives the rotational scale. Document
that this metric depends on the chosen independent coordinates and mesh;
it is not a physical mass norm. A later physical-node metric would require
`W = T^T W_full T`, including diaphragm lever arms and possibly off-diagonal
entries. Do not implement that accidentally as a diagonal approximation.

Changing reference load by a multiplier `c` and changing `lambda_*` by
`1/abs(c)` preserves the physical continuation metric; changing length units
requires changing `u_*` consistently. These are explicit invariance tests.

## 4. Predictor and branch orientation

Persist only accepted history: last total `(Delta q, Delta lambda)`, its
metric-normalized direction `t`, accumulated accepted chord length, and the
next proposed radius. At the first step:

```text
K v = p
t = sign * (v, 1) / sqrt(v^T W v + beta^2)
(Delta q_pred, Delta lambda_pred) = s t
```

The requested initial direction is increasing or decreasing `lambda`.
Check solve accuracy and finiteness. A singular first tangent with no
previous direction needs the `seedDirection` setting: a user-supplied
`(Delta q, Delta lambda)` direction in reduced coordinates. Normalize it in
the metric and use it as `t_prev` for the bordered tangent below, not as
the predictor itself. Without a seed, fail initialization clearly. Shrinking
the step repeatedly does not repair a singular first tangent.

For later steps, the **default predictor is the normalized last accepted
secant**: `t = (Delta q_n, Delta lambda_n) / s_n`. It needs no extra
factorization and is oriented by construction. At kinks in piecewise-linear
materials (the main first-release target) it is no worse than a tangent.
Optionally (`predictor: "tangent"`), and always for the step after a
cutback or a sharp turn, obtain a tangent from the bordered system at the
committed state, with `t_prev` the last accepted direction and
`M = diag(W, beta^2)`:

```text
[ K       -p                  ] [v_q     ] = [0]
[ (M t_prev)_q^T  (M t_prev)_l ] [v_lambda]   [1]
t = v / sqrt(v^T M v)
```

The normalization row gives `t^T M t_prev = 1/sqrt(v^T M v) > 0` by
construction. Assert this as a solve-quality check rather than flipping the
sign. This follows the path through a change in the sign of
`Delta lambda`; **do not choose the direction from determinant sign or
previous load-increment sign**. Apply both displacement and load predictor
before evaluating the first nonlinear corrector. If the tangent predictor
fails, fall back to the secant, then cut back if needed.

At acceptance, require the total increment `d` to have positive metric
projection onto the predictor direction, `d^T M t > 0`. Do not
additionally require a minimum turn angle. Kinks at material breakpoints
legitimately turn the path by large angles: 66.8 degrees for the softening
fixture and 84.4 degrees for the series snap-back fixture in §8. For the
step that crosses such a kink, the projection ratio does not depend on the
radius, so a tighter angle threshold would exhaust the cutbacks. Instead,
resolve ambiguity relative to the backward root. The accepted point must be
closer in the metric to the predictor point `x_n + s t` than to the
reflected point `x_n - s t`. That is the same as `d^T M t > 0` on the
sphere, so additionally require `d^T M t >= eps_orient * s^2` with a small
`eps_orient` (initially `1e-6`) that only guards against roundoff. A near-zero
projection triggers a smaller radius.
This protects tracking but does not resolve a true bifurcation or guarantee
the closest branch at a large radius. A discontinuous tangent with a metric
turn exceeding 90 degrees can defeat this orientation rule even with
cutbacks; report the ambiguity. Such a cusp needs an explicit seed/metric
change or an event-aware continuation extension.

## 5. Coupled Newton corrector and globalization

At every iteration evaluate the current trial, then solve:

```text
[ K       -p                   ] [delta q     ] = [ r ]
[ 2 Delta q^T W  2 beta^2 Delta lambda ] [delta lambda]   [-g ]
```

Signs follow the definition `r = P - F_int`: equilibrium linearization is
`K delta q - p delta lambda = r`. Scale unknowns by the displacement and
load scales and equilibrate rows before factorization. The bottom row/RHS
can be divided by `s^2` for the dimensionless constraint. Include structural
zeros in a stable bordered sparsity pattern so symbolic reuse works as
values and load-increment signs change. Keep the bordered cache separate
from an existing `n`-DOF structural factorization.

The border can remain nonsingular at a simple load limit even when `K` is
singular. Verify this directly in tests. An unconstrained mechanism or a
true branch point can still make the bordered system singular; shrinking
the radius does not guarantee a solution. The bottom row also degenerates
if an iterate returns to the committed point (`Delta q = 0`,
`Delta lambda = 0`). The applied predictor makes this unlikely, but detect
it and treat it as a failed attempt.

The constraint row `2 Delta q^T W` is dense over every DOF, and `p` may be
dense too. The implementation-step-2 measurement settled how to solve it.
Factoring the bordered matrix directly lets partial pivoting pull the dense
row in early, and fill grows far faster than `K`'s own. Measured with the
ignored `bordered_fill_and_factor_time_on_a_frame` test (release build,
corotational elastic frames, `n`-DOF LU = `SparseSolver`'s supernodal LU):

| Frame | `n` | `nnz(K)` | `n`-DOF LU | Block elimination: `nnz(L+U)`, time | Direct bordered LU: `nnz(L+U)`, time |
| --- | --- | --- | --- | --- | --- |
| 10 stories x 3 bays | 120 | 1,548 | 0.07 ms | 3,285, 0.14 ms | 8,477, 0.26 ms |
| 20 x 5 | 360 | 4,932 | 0.41 ms | 14,124, 0.63 ms | 70,385, 2.9 ms |
| 40 x 8 | 1,080 | 15,318 | 2.3 ms | 61,269, 2.3 ms | 604,298, 32 ms |

So **block elimination is the default**. Only `K` is factored, and the
border is eliminated through its Schur complement. For `[K c; r^T d]`:

```text
a = K^-1 f,  b = K^-1 c,  y = (h - r^T a) / (d - r^T b),  x = a - b y
```

For the corrector, `c = -p`, `r = 2 W Delta q` and `d = 2 beta^2 Delta lambda`.
The direct bordered LU is the **fallback**, used only when `K`'s
factorization breaks down (a zero pivot, as at an exactly singular `K` on
a smooth load limit, where the border keeps the system nonsingular) or when
the Schur complement cancels (`|d - r^T b| <= 1e-8 (|d| + sum |r_i b_i|)`).
On either path, every solve is checked against the *full* bordered system's
scaled backward error (target `1e-10`, with one refinement step). Both
paths use faer's low-level simplicial LU with cached COLAMD orderings. The
`U` diagonal and permutations give `sign(det A)` directly, and with block
elimination `det A = det K * schur`. The first release does not depend on
quadratic-root selection.

Use coupled backtracking with `0 < eta <= 1`, and **project every
candidate radially back onto the sphere**:

```text
d(eta) = Delta_base + eta * delta,     Delta_trial = s * d(eta) / ||d(eta)||_M
```

This replaces the original plan's "trial states may leave the sphere".
Without projection, the first implementation stalled at the softening
fixture's kink. There the full correction is tangent to the sphere and
comparable to `s`, so the quadratic constraint error after a full step is
`O(1)` relative to `s^2`. The merit function then rejected every step but
tiny damped ones, and the attempt stagnated at every radius. The scheme is
now an updated-normal-plane corrector with spherical projection. With
`g = 0` at the base, the constraint row makes `delta` tangent to the sphere,
so projection changes the step only at second order and `delta` stays a
descent direction for

```text
Phi = 0.5 * (||D_f^-1 r||_2^2 + (g / (s^2 arcTolerance))^2)
```

`D_f` holds the per-equation force thresholds, frozen for the attempt, so
`Phi <= ~1` means both tests pass, and `g` is roundoff after projection.
Start with full Newton. If needed, halve `eta` until a sufficient decrease
(Armijo constant `1e-4`) or until the minimum `eta = 1/128`. Every trial
satisfies the constraint by construction; acceptance still checks
`|g| / s^2 <= arcTolerance` independently. If a nonsmooth material gives no
decreasing trial, reject the attempt and reduce `s`. Each candidate is
evaluated from the same committed history `h_n`. If the merit hasn't
halved over five iterations, the attempt is rejected as stagnated.

Do not apply the existing work-based bisection/regula-falsi search to arc
length: it currently varies `q` while keeping the corrected scalar fixed.
Version one supports `Newton { tangent: Current, line_search: None }` only,
with the coupled backtracking owned by the continuation driver. Reject
`Linear`, stale-tangent strategies, `KrylovNewton`, and ordinary line-search
combinations explicitly. Coupled acceleration/reuse can be added after
they satisfy the constraint and singular-limit tests.

## 6. Convergence tests and numerical failure checks

### Accepted-state evaluation for every Newton analysis

Replace `check(residual, du)` with an iteration context containing the
**post-update** residual at the final accepted line-search state, actual
applied correction, scalar correction, total step increments, and optional
control residual. Transient callers must update velocity/acceleration
before residual evaluation. Reuse that assembly as the next iteration's
system where practical. Check the predictor before solving a correction:
an exact elastic predictor should pass without an unnecessary solve.
Count iteration budget as corrector solves, with an acceptance check after
the last permitted solve.

Keep the three existing absolute tests and their wire spellings. Fix the
state timing of `NormUnbalance`; document that this can change iteration
counts. `NormDispIncr` and `EnergyIncr` remain available for compatibility,
but neither alone certifies force equilibrium. Update `EnergyIncr`'s
documentation: work pairs force/translation and moment/rotation consistently,
but its value changes with work units and cancellation can make it small.

### New composite test

Add `ConvergenceTest::Combined` / wire `kind: "combined"`, with a required
force criterion and optional displacement/load correction and work criteria.
All configured criteria must pass. Use finite positive characteristic scales
and finite nonnegative relative tolerances with strictly positive absolute
tolerances. Provide these measures:

1. **Force equilibrium:** on each reduced equation require
   `abs(r_i) <= forceAbs_i + forceRel * F_ref_i`. Use separate absolute
   tolerances for force and moment equations. Freeze `F_ref_i` at attempt
   start using characteristic force/moment scales and committed/predicted
   external and internal forces, taking their maximum. Do not normalize by
   the initial residual (which can be zero), the current residual, or solely
   current load (which can decrease to zero). Report the maximum scaled
   component and optionally the L2 norm; avoid dilution by many quiet DOFs.
2. **Arc constraint, mandatory independently of test variant:**
   `abs(g)/s^2 <= arcRel`, initially `arcRel = 1e-8`. No force-only or
   work-only convergence setting can bypass this condition.
3. **Optional coupled correction:** with
   `c = sqrt(delta q^T W delta q + beta^2 delta lambda^2)`, require
   `c <= correctionAbs + correctionRel*s`. Include `delta lambda`, even
   when displacement correction is tiny. Use the last undamped correction
   as an additional stagnation diagnostic: a tiny accepted `eta` cannot
   make failed force/constraint tests pass. Exact predictor convergence
   treats the equilibrium correction as zero.
4. **Optional incremental work:** use the accepted correction and accepted
   residual, and a fixed characteristic work scale for abs/rel tolerance.
   Treat it as supplementary; do not use work as the sole acceptance gate.

For arc length, require force equilibrium even when an explicitly selected
legacy displacement/work test passes. A legacy force `tol` continues to
mean its absolute Euclidean residual tolerance; the arc configuration must
provide the extra force gate for legacy displacement/work selections, or
reject those selections with a clear validation message. Recommended first
release: reject those two legacy selections for arc length and require
`Combined` or `NormUnbalance`, always plus the arc gate. When convergence
is omitted in an arc-length wire stage, choose `Combined`, rather than the
ordinary static default. Publish concrete defaults for dimensionless
relative tolerances (`forceRel = 1e-6`, `maxIter = 30`); users supply force,
moment, and displacement scales/absolute tolerances appropriate to units.

Check all assembled entries, residuals, corrections, scalar values, norms,
and solve outputs for finiteness. Check linear-system backward error in
scaled coordinates; do not treat a finite LU result as proof of a usable
solve. Detect stagnation using force and constraint histories, not just
tiny increments. After five iterations without meaningful merit decrease,
cut back. Initially use a documented scaled backward-error target `1e-10`
and measure it against the achievable roundoff floor in tests; report poor
solve accuracy separately from nonlinear nonconvergence.

## 7. Adaptive radius and transactional state

Proposed settings: `initialRadius`, `minRadius`, `maxRadius`, `targetIterations`
(default 6), `maxRetries` (default 8), initial direction, optional
`seedDirection`, `predictor` (`"secant"` default, or `"tangent"`), metric
scales or auto-scale, `arcRel`, coupled-backtracking settings, and stop
criteria. Validate
`0 < minRadius <= initialRadius <= maxRadius`, finite positive scales,
positive iteration budgets, and bounded retry counts. A fixed-radius mode
sets all three radii equal and disables retries/growth.

For an accepted step with `N` corrector solves:

```text
s_next = clamp(s_used * clamp(sqrt(N_target/max(N,1)), 0.5, 1.5), s_min, s_max)
```

Suppress growth immediately after a retry, significant damping, or a sharp
turn. For a failed attempt restore the original step snapshot and halve
`s_used`. If the halved radius would be below `minRadius`, optionally try
`minRadius` once, then fail; never repeat indefinitely at the minimum.
Invalid configuration, unsupported loads, or failed initial equilibrium
are immediate errors, not retriable failures.

Only an accepted attempt commits materials/elements, increments step count,
updates `lambda`, stores direction/secant history, advances accumulated
chord length, and records a sample. Failed line-search candidates and retries
must never commit. Snapshot continuation state as well as the domain;
clear numeric factors after a failed attempt. On exhausted retries, restore
the entire pre-call state including the proposed next radius, so repeating
the same call has defined behavior. Diagnostics can survive rollback.

`set_integrator()` must clear continuation history when entering/leaving
arc length or changing the metric/load definition; retain it for unchanged
arc-length configuration. Changing radius bounds alone can preserve the
direction. Because `domain_mut()` permits model/load changes, invalidate
history and numerical caches conservatively when granting that mutable
access, and revalidate equilibrium/load sensitivity before continuation.

### Stop criteria

Because `lambda` can reverse, pseudo-time targets do not define when a
phase ends. Support these optional criteria. Any one that is met ends the
phase, after the step that meets it has been accepted:

- target value of one independent DOF (for example roof drift), in either
  direction;
- target `lambda`, and optionally its first zero crossing after the peak;
- maximum accumulated chord length;
- maximum accepted step count.

For displacement and `lambda` targets, optionally land on the target
exactly. When an accepted step brackets the target, roll back and retake
that step once with the radius found by interpolating the chord. Accept the
retaken step if it hits the target within tolerance; otherwise keep the
bracketing step and report the overshoot. Report which criterion ended the
phase. Reaching the step cap without meeting another criterion is not an
error.

### Bifurcation diagnostics

Do not use the determinant for orientation, but record its sign at each
accepted state. The sign is almost free: it comes from the product of the
signs of the LU's `U` diagonal and the parity of the row and column
permutations. Take it from the `n`-DOF factor, or from the bordered factor
corrected for the border. If the sign of `det K` changes while
`Delta lambda` does not reverse, the run has probably crossed a bifurcation
and stayed on the primary branch. Emit a structured warning with the step
index; it is not an error. The sign detects only odd numbers of eigenvalue
crossings. Counting negative eigenvalues (inertia) would need `LDL^T` and is
a later extension.

### Diagnostics

Return diagnostics: accepted radius, next radius, retries, corrector solves,
factorizations across all attempts including predictors, force measure, arc
error, accumulated chord length, sign of `det K`, any derived auto-scales,
and the criterion that ended the phase. Keep load factor as the physical path
parameter; use step/sample index for progress and storage ordering. Never
sort or deduplicate results by load factor: it can reverse and repeat.
Keep recorder tuple semantics compatible; add optional diagnostics rather
than substituting chord length for existing pseudo-time values. Keep retries
bounded in `advance()`; cancellation remains between accepted step calls
under the current API. Iteration-level yielding would require separate
session-state work if model sizes make this insufficient.

Suggested structured errors distinguish invalid arc configuration,
unsupported active load series, zero reference sensitivity, initial state
out of equilibrium, unavailable seed direction, nonfinite trial state,
inaccurate/singular bordered solve, branch-orientation failure, and exhausted
cutbacks. Include step, attempt count, radius, and last convergence measures
where applicable; update session error mapping for both profiles.

## 8. Degrading-strength integration test

Create `core/tests/arclength_softening.rs`. Use existing `ZeroLength` and
`Material::hysteretic`; no production material addition is required. Two
coincident nodes, one fully fixed, the other free only in x. Apply a linear
unit reference force to the free x DOF. Define the spring:

```rust
Material::hysteretic(
    10.0, 0.01, 6.0, 0.02, 2.0, 0.03,
   -10.0,-0.01,-6.0,-0.02,-2.0,-0.03,
    1.0, 1.0, 0.0, 0.0, 0.0,
)
```

Thus on monotonic positive deformation the independent oracle is:

```text
F(u) = 1000 u                         0 <= u <= 0.01
F(u) = 10 - 400 (u - 0.01)            0.01 < u <= 0.03
lambda = F(u)                        unit reference load
```

This is degrading **backbone strength**, not cycle-induced damage. It tests
negative tangent and continuation past peak strength without relying on
cyclic-damage formulas. Stay below `0.03` to avoid testing the material's
terminal plateau/tiny-tangent convention in this benchmark.

Proposed fixture configuration: `u_* = 0.01`, `lambda_* = 10`, fixed
dimensionless radius `0.05`, increasing initial direction, current Newton,
composite force tolerance `1e-9` absolute plus `1e-9` relative, and
`arcRel = 1e-9`, `maxIter = 30`. Step until `u >= 0.028`, with a hard
100-step cap. The configuration is illustrative until the new API exists;
implement the test with the final public builder API, not a private scalar
solver substitute.

At **every accepted step**, assert:

- finite states and convergence diagnostics;
- force equilibrium against the piecewise oracle and support reaction
  `reaction_fixed_x = -lambda`, using matching scale-aware tolerances;
- `u` strictly increases within numerical tolerance, including across the
  peak; `lambda` decreases on at least ten post-peak steps;
- `lambda` has an interior maximum with samples on both sides of `u=0.01`;
  for fixed `s=0.05`, peak undersampling is bounded by `1000*u_*s = 0.5`;
- the measured chord satisfies `(Delta u/0.01)^2 + (Delta lambda/10)^2
  = s_used^2` within `arcRel`, using **total step increments**;
- the final deformation lies in `[0.028, 0.029)` and force is below 3,
  demonstrating substantial strength loss without returning to unloading.

Also compute the oracle independently at `u=0.01`, `0.02`, `0.028`
(`F=10`, `6`, `2.8`) and compare the arc-length path to those landmarks
with interpolation or bracketing; do not assert it samples them exactly.
Run a separate increasing load-control baseline past `lambda=10`, verify
failure to converge and exact rollback. Success of arc length must not
depend on that baseline's iteration counts.

### Series-spring variant: actual snap-back

Add an elastic zero-length spring of stiffness `k_s=380` in series between
the degrading spring's free node (`x`) and a new loaded node (`y`). Fix all
other DOFs; the unit load acts at `y`. All nodes may be coincident.

```text
lambda = F(x)
y = x + F(x)/380
pre-peak: y = (69/19)x
post-peak: dy/dx = 1 + (-400)/380 = -1/19
```

Continue with `x` increasing through `0.01`: force **and loaded-node
displacement** must decrease afterward. At `x=0.01`, `(y,lambda)=(0.03631578947,10)`;
at `x=0.02`, `(y,lambda)=(0.03578947368,6)`; at `x=0.028`,
`(0.03536842105,2.8)` (use the exact series formula in assertions).
Use the metric on both free translations with the same scales and radius,
a 250-step cap, and stop when `x >= 0.028`. Assert both nodal equilibrium
equations, the series-force oracle, reactions, the two-DOF arc constraint,
and positive predictor projection. This catches direction selection based
only on load-factor sign or the loaded displacement. It is more demanding
than a single softening spring, which alone does not establish snap-back.

In the `(x/u_*, y/u_*, lambda/lambda_*)` metric, the path turns 84.4 degrees
at the kink: `cos = 0.097`. The step that crosses the kink therefore
exercises the orientation rule of §4 close to its 90-degree limit. Assert
that it succeeds without exhausting cutbacks. In this fixture the
displacement turning point coincides with the load peak at the material
kink. The smooth snap-back test in §9 covers a turning point where `K` is
nonsingular.

## 9. Further acceptance tests

| Test | Required evidence |
| --- | --- |
| Elastic closed form, 2D and 3D | Predictor radius and load factor match the analytic solution; no needless corrector; scales and negative initial direction work. |
| Smooth peak with exactly singular tangent | Test-only conservative spring/`ElementOps` fixture `F(u)=k*u*exp(-u/a)`. Start exactly at `u=a` in equilibrium: an initial displacement `a` plus a constant preload carrying the peak force, with a linear reference pattern to continue. A displacement-control phase can't reach the peak exactly, because its own corrector needs `K^-1` there. Then continue with `seed = (+1, 0)` in `(q, lambda)`. This proves that the bordered predictor and corrector do not require `K^-1`. Analytic peak is `k*a/e`. Also assert that omitting the seed fails initialization. |
| Smooth snap-back | Same exponential spring in series with an elastic spring of stiffness `k_s < k/e^2` (for example `k_s = 0.1k`). Then `dy/dx = 1 + F'(x)/k_s` changes sign smoothly at the two roots of `(x/a - 1) exp(-x/a) = k_s/k`, about `1.41a` and `3.0a` for `k_s = 0.1k`. These give a snap-back and a recovery turning point, both beyond the load peak at `x = a` and both with `K` nonsingular. Use the tangent and secant predictors. Assert the exact series oracle, a decreasing `y` between the turning points, and no cutbacks at a moderate radius. |
| Geometric snap-through | Shallow symmetric two-bar assembly using supported corotational geometry; intentional symmetry constraint or small declared imperfection. Trace descending and recovering load branches; compare equilibrium and a independently derived symmetric path/refined solution. Without the symmetry constraint, where an asymmetric bifurcation exists, assert the `det K` sign-change warning. |
| Stop criteria | Target DOF displacement, target and zero-crossing `lambda`, chord-length and step caps each end the phase and report the criterion. Exact landing hits the target within tolerance; on failure, the bracketing step is kept and the overshoot reported. Repeating `advance` after the stop is defined. |
| Auto-scale | Derived scales match the hand-computed values for the softening fixture and reproduce its explicit-scale path. A singular first tangent fails initialization. |
| Coupled solve accuracy | Indefinite and singular-`K` but nonsingular-border matrices; true singular border, near-singular border, and nonfinite outputs produce defined errors. |
| False convergence guards | Small `du` with large force residual; large `delta lambda` with tiny `du`; zero/cancelling work with large residual; balanced forces with violated arc constraint; last-budget correction newly creates a residual. None may commit. |
| Post-update shared-loop regression | Newton static and nonlinear Newmark converge at the final accepted state, including ordinary line search. Legacy API still decodes; updated iteration counts are intentional. |
| Cutback and rollback | Oversized radius and constrained iteration budget provoke retries; accepted path agrees with a small-radius run. Exhaustion preserves domain, material energy/damage/maxima, scalar, count, and history; records contain no failed samples. |
| Material-history isolation | Repeated trial/backtracking evaluations from the same committed Hysteretic/Pinching4 state are identical; forced retry plus successful smaller step agrees with a clean smaller-step run. |
| Frozen loads / phase transitions | Gravity stays constant; only lateral patterns contribute `p`; reaction balance includes both. Invalid initial equilibrium and active `Path` are rejected. |
| Constraints and scale changes | Equal DOF and affine 3D diaphragm use reduced coordinates consistently. Convert length/force units and reference-load amplitude with corresponding scales; recover the same physical path. |
| Refinement / orientation | Halved radii recover the same branch; record turn angles and load reversals. Repeat softening fixture in reverse initial direction using its negative backbone. |
| WASM session / storage ordering | Both profiles decode arc settings; generated types expose options; decreasing/repeated load factors survive recorder batches with monotonically increasing sample indices. `advance(1)` versus larger budgets produces identical accepted samples. |
| Invalid option combinations | Reject nonfinite/invalid radii and scales, zero sensitivity, unsupported algorithms/search, and incompatible legacy convergence tests in native and wire entry points. |

Smooth-peak and geometric tests complement the degrading-strength tests:
piecewise backbone peak crossing alone does not prove exact-singular-tangent
handling. A follow-up cyclic-damage test can use a prescribed reversal
protocol and Pinching4, but should distinguish physical history changes
from trial-state contamination.

## 10. Implementation sequence and completion gates

1. Introduce accepted-state convergence context, composite force tests, and
   nonfinite checks. Fix the shared loop timing; run static/transient
   regression tests before adding continuation.
2. Add metric/load validation, reduced-state access, sparse bordered
   assembly, scaling, and accurate-solve checks. Test singular-`K` cases.
   Measure bordered fill and factor time on a representative frame model
   now, and choose bordered LU or block elimination as the default.
3. Add transactional arc history, the secant/tangent predictors with
   `seedDirection`, coupled Newton and the arc acceptance gate. Deliver fixed-radius elastic, degrading-strength,
   series snap-back, and smooth-peak tests through the public API.
4. Add coupled backtracking, bounded cutbacks, adaptation, diagnostics,
   and rollback/history-isolation tests. Add geometric and refinement cases.
5. Add planar/spatial wire settings, decode errors, session diagnostics,
   and native boundary tests. Rebuild generated JS/TypeScript declarations
   through `wasm-pack`, then run Node/TypeScript boundary tests. Update
   `docs/algorithms.md` and README with the implemented API and restrictions.

Completion requires `cargo test --workspace` on native followed by the
repository's WASM build, Node smoke test, and generated TypeScript check.
Benchmark representative frame/truss systems to measure bordered fill and
factorization cost; an `(n+1)` sparse solve is not automatically cheap because
the border is globally coupled. Optimize only after the specified path,
convergence, and rollback tests pass. pysees UI/compiler exposure is a
separate integration change; the Carapace wire and session support belongs
in this implementation.
