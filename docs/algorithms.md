# Algorithm richness: line search, tangent strategy, Krylov acceleration, event-to-event stepping

Plan for M12 (implementation-plan §6). `Algorithm` today is `Linear` /
`NewtonRaphson` only (`core/src/analysis/algorithm.rs`) — full Newton,
current tangent re-formed and re-factored every single iteration, no
globalization (line search) and no way to trade iteration count for
factorization count. Every real nonlinear-analysis package (including
Xara/OpenSees, this project's own reference) treats those as separate,
composable axes, not a family of unrelated algorithm classes. This
document designs the same composability here, with each piece ported from
a specific, named Xara source file — not re-derived from a textbook
description, the same discipline `ForceBeamColumn`'s state determination
followed (see its module doc comment for why that matters: a hand
re-derivation of "the same algorithm" silently introduced a real
convergence bug that only surfaced under nonlinear load).

Scope: **static `Analysis` only.** `TransientAnalysis` has no Newton
corrector at all yet (a separately-flagged gap, implementation-plan §6's M6
note) — once it grows one, it can reuse whatever `TangentStrategy`/
`LineSearch` types come out of this milestone, but that's follow-on work,
not part of this one.

---

# Part I: Newton-family richness (tangent strategy, line search, Krylov acceleration)

## 1. What "richer algorithms" actually means — three independent axes

Xara's `EquiSolnAlgo` subclasses look like a large flat list
(`NewtonRaphson`, `ModifiedNewton`, `NewtonLineSearch`, `BFGS`, `Broyden`,
`AcceleratedNewton`, ...) but they're really three independent choices
multiplied together:

1. **Tangent strategy** — which stiffness matrix gets used at each
   iteration: re-form the *current* (trial-state) tangent every time (full
   Newton), form it once from the *initial* (undeformed/first-iteration)
   state and reuse it for the rest of the step (initial-stiffness
   iteration — quadratic convergence traded for no re-factorization cost),
   or form it once at the *start of the step* and reuse it (classic
   "Modified Newton" — cheaper than full Newton, better-conditioned than
   pure initial-stiffness for a step that's already far from the
   undeformed state).
2. **Globalization / line search** — after solving `K * du = r`, whether
   to accept `du` outright or scale it by some `eta` chosen to actually
   reduce the residual (guards against Newton overshooting badly for a
   very nonlinear or nearly-singular tangent).
3. **Acceleration** — whether a cheap subspace correction (built from the
   last few residual/solution differences, at the cost of a tiny dense
   least-squares solve, not a new factorization) is applied on top of a
   tangent-strategy-1/2/3 iterate to recover convergence speed lost by
   reusing a stale tangent.

Xara's class names are really `{tangent strategy} x {line search on/off}
x {accelerator on/off}` (see `NewtonRaphson.cpp`'s `prediction_tangent`/
`correction_tangent` constructor fields, `NewtonLineSearch.cpp`'s identical
fields plus a `LineSearch*`, and `AcceleratedNewton.cpp`'s
`Accelerator*` — the same iteration skeleton every time, just with
different pieces plugged in). Carapace should represent it as those three
orthogonal choices on `Algorithm`, not grow a combinatorial pile of enum
variants — closer to how `ConvergenceTest`'s three variants are one closed
enum already (§2.1), just with `Algorithm` carrying a couple of composable
fields instead of only a bare unit variant.

---

## 2. Xara reference map

| Concept | File(s) | What to port |
|---|---|---|
| Full Newton-Raphson (current behavior) | `SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.cpp/h` | Already have this; confirms the `prediction_tangent`/`correction_tangent` split is the right generalization axis. |
| Modified/initial-tangent Newton | `SRC/analysis/algorithm/equiSolnAlgo/ModifiedNewton.cpp/h` | Forms tangent **once**, before the iteration loop (not inside it) — the actual mechanical difference from `NewtonRaphson`, not just "which tangent." |
| Line search driver | `SRC/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.cpp/h` | Wraps the base Newton loop: after `update(dX)`, evaluates whether the new residual actually satisfies convergence; if not, hands off to a `LineSearch` strategy to rescale `dX`, re-tests, and only then commits the state. |
| Line search strategy base | `SRC/analysis/algorithm/equiSolnAlgo/search/LineSearch.h` | The `search(s0, s1, dU, G, Xs, integrator)` interface — `s0`/`s1` are `dU . residual` before/after the raw step (the 1-D "energy" the search drives toward zero). |
| Line search: bisection | `SRC/analysis/algorithm/equiSolnAlgo/search/BisectionLineSearch.cpp` | Bracket a sign change in `s(eta)` by geometric expansion (`eta *= 4`), then bisect. Most robust, the one to port first. |
| Line search: regula falsi | `SRC/analysis/algorithm/equiSolnAlgo/search/RegulaFalsiLineSearch.cpp` | Same bracketing, false-position instead of bisection for the refinement step — usually fewer iterations than bisection once bracketed. |
| Line search: secant | `SRC/analysis/algorithm/equiSolnAlgo/search/SecantLineSearch.cpp` | No bracketing — pure secant update on `s(eta)`; cheaper per iteration, less robust (can diverge without a bracket). |
| Line search: Crisfield interpolation | `SRC/analysis/algorithm/equiSolnAlgo/search/InitialInterpolatedLineSearch.cpp` | `eta *= s0/(s0 - s)`, from Crisfield's *Nonlinear FE Analysis*, cited in the file itself — cheapest, defer past the MVP two. |
| Krylov-accelerated modified Newton | `SRC/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.cpp/h` | The driver: form tangent once (like `ModifiedNewton`), then every iteration solve with the *stale* tangent and hand the raw correction to an `Accelerator` before applying it. |
| Krylov subspace accelerator | `SRC/analysis/algorithm/equiSolnAlgo/accelerator/KrylovAccelerator.cpp/h` | The actual math (§5 below) — a small dense least-squares solve (`dgels` in Xara) over the last `maxDimension` residual differences. `KrylovAccelerator2.cpp` is a later/alternate variant — not needed for a first port, noted only so nobody rediscovers it and wonders why it wasn't used. |
| Quasi-Newton (BFGS, Broyden) | `SRC/analysis/algorithm/equiSolnAlgo/BFGS.cpp/h`, `Broyden.cpp/h` | Explicitly **out of scope** for this milestone (§7) — a genuinely different update rule (secant condition on the *inverse* tangent, no re-factorization ever), not a combination of the three axes above. Revisit only if a model actually needs it. |

---

## 3. `TangentStrategy`: the first new type

```rust
/// Which tangent stiffness an `Algorithm` iteration solves against.
/// Closed enum (§2.1) — mirrors Xara's `TangentFlagType`
/// (`CURRENT_TANGENT` / `INITIAL_TANGENT` / `PREDICTOR_TANGENT`), which
/// every algorithm class above takes as a constructor argument rather
/// than hard-coding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TangentStrategy {
    /// Re-form and re-factor the tangent at the *current* trial state,
    /// every iteration. Full Newton-Raphson — today's only behavior.
    Current,
    /// Form the tangent once, at the state the step *started* from, and
    /// reuse it (still factored once) for every iteration of that step.
    /// Xara's `ModifiedNewton` with `CURRENT_TANGENT` passed at
    /// construction (misleadingly named from the caller's perspective —
    /// it means "current as of step start," not "current every
    /// iteration").
    ReuseAtStepStart,
    /// Form the tangent once, from the model's very first (undeformed)
    /// state, ever — reused for every iteration of every step. Xara's
    /// `ModifiedNewton` with `INITIAL_TANGENT`. Cheapest, worst
    /// convergence rate (linear, not quadratic); pairs naturally with the
    /// Krylov accelerator (§5) to claw the convergence rate back without
    /// paying for re-factorization.
    Initial,
}
```

`Current` is what `Algorithm::NewtonRaphson` already does — this is a
strict superset, not a rewrite of existing behavior.

### 3.1 The real prerequisite: `SparseSolver` needs factor/solve to split

`TangentStrategy::ReuseAtStepStart`/`Initial` are pointless as *performance*
features if `Analysis` still calls `Domain::form_tangent_and_residual`
(which assembles a fresh `K` every call) and `SparseSolver::solve` (which
re-factors from scratch every call — already flagged as a "future
optimization, not yet done" in `solver.rs`'s own doc comment) on every
iteration regardless of `TangentStrategy`. Implementing the enum without
this would still be *correct* (same iteration sequence as Xara) but
deceptive — it would look like it saves factorization cost and wouldn't.

This milestone is the natural point to finally do that split, since
`TangentStrategy` is the first real consumer:

```rust
impl SparseSolver {
    pub fn factor(&self, k: &SparseMatrix) -> Result<SparseFactorization, AnalysisError>;
}
pub struct SparseFactorization { /* faer's Lu, or equivalent */ }
impl SparseFactorization {
    pub fn solve(&self, r: &DVector<f64>) -> DVector<f64>;
}
```

`Analysis` holds `cached_factorization: Option<SparseFactorization>`,
invalidated (re-formed) exactly when `TangentStrategy` says to: every
iteration for `Current`, once per `step()` call for `ReuseAtStepStart`,
once ever (first `step()` call only) for `Initial`.

Getting the *residual* without re-assembling `K` also needs a new `Domain`
method — today `assemble_stiffness_triplets` computes stiffness and
resistance together in one element loop because every `Element` variant's
`form_tangent_and_resistance` returns both at once (no element currently
exposes a resistance-only path). Two honest options, in order of
preference:

1. **Do it properly**: add `Element::form_resistance(&self, node_i, node_j) -> SVector<f64, ELEMENT_DOF>` alongside `form_tangent_and_resistance`, implemented per element (usually just the second half of the existing method, e.g. `ElasticBeamColumn`'s `resistance_local` without ever building `k_local`; `ForceBeamColumn`'s `try_state_determination` already separates `q`/`k_basic` internally so this is mostly plumbing, not new math). Real assembly-cost savings on top of the factorization savings.
2. **Ship the factorization split alone first**: still call the existing combined assembly every iteration (so no assembly-cost savings yet), but only *factor* a fresh `K` when `TangentStrategy` says to — solve against the cached factorization otherwise. Still delivers the real, expensive half of the savings (factorization dominates assembly for anything past toy-sized models) with zero `Element` API surface change.

Recommendation: ship (2) first (stage A2 below), do (1) only if profiling on
a real `pysees` model shows assembly cost actually matters — same
"benchmark before optimizing further" discipline implementation-plan §9.1
already applies to `Material::evaluate`'s boxing cost.

---

## 4. `LineSearch`: a modifier, not a separate algorithm

```rust
/// Post-processes a raw Newton correction to actually reduce the
/// residual, rather than accepting it outright. Closed enum (§2.1).
/// Ported from `xara/SRC/analysis/algorithm/equiSolnAlgo/search/`.
#[derive(Debug, Clone, Copy)]
pub enum LineSearch {
    Bisection { tol: f64, max_iter: usize, max_eta: f64 },
    RegulaFalsi { tol: f64, max_iter: usize, max_eta: f64 },
}
```

Applied inside the Newton loop, after `du` is computed and *before* it's
committed to the domain — mirrors `NewtonLineSearch.cpp`'s structure
exactly:

1. Compute `s0 = du . residual` (the directional derivative of the
   1-D energy along `du`, at `eta = 0`).
2. Apply the raw `du` (`eta = 1`), re-form the residual, compute
   `s1 = du . residual_new`.
3. If the *ordinary* convergence test already passes at `eta = 1`, skip
   the search entirely (`NewtonLineSearch.cpp`'s `if (result < 1)` guard —
   don't pay for a line search step that wasn't needed).
4. Otherwise hand `(s0, s1, du)` to the strategy, which returns a rescaled
   `Xs` (`Bisection`: bracket `s(eta) = 0` by geometric expansion of `eta`
   then bisect, per `BisectionLineSearch.cpp`; `RegulaFalsi`: same bracket,
   false-position refinement, per `RegulaFalsiLineSearch.cpp`).
5. Re-test convergence against `(residual_at_Xs, Xs)`, same
   `ConvergenceTest` the outer loop already uses — no new test machinery
   needed, `ConvergenceTest::check`'s existing `(residual, du)` signature
   is already exactly what a line search needs to re-check.

Only `Bisection` and `RegulaFalsi` for the initial port (both bracket
first, so both are robust — `Secant`/`InitialInterpolated` can diverge
without a bracket and add little beyond what these two already cover for
the load cases `pysees` is likely to hit). Add the other two only if a
real model needs the extra speed.

---

## 5. Krylov-accelerated modified Newton

The actual algorithm (`KrylovAccelerator.cpp::accelerate`, distilled):
keep the last `k <= max_dimension` pairs of (raw correction `v_i`, the
change in residual it produced, `Av_i = r_{i-1} - r_i`). Given a new raw
correction `v_k` (from solving against the *stale* tangent) and its
residual `r`, solve the small `(numEqns x k)` least-squares problem

```text
minimize || r - Av * c ||   over c in R^k     (Av = [Av_0 .. Av_{k-1}])
```

and correct: `v_k += sum(c_i * v_i)`, i.e. use the history of how past
corrections actually moved the residual to extrapolate a better
correction than the stale tangent alone would give — a GMRES-flavored
idea, but built from already-computed Newton iterates instead of an
explicit Krylov-subspace matrix-vector product loop. When `k` exceeds
`max_dimension`, the accelerator signals the driver to re-form (and
re-factor) the tangent and reset the subspace — this is `is_finite`-style
bookkeeping, not new numerics.

`numEqns` is the full system size (can be large) but `k` is small (Xara
defaults `maxDimension = 3`), so this is a tall-skinny dense least-squares
problem — `nalgebra`'s QR decomposition (`Av.qr().solve(&r)`) covers this
exactly; no need for `faer` or a LAPACK binding the way `dgels` is used in
Xara (§2.6 already commits to one solver for the *global sparse* system;
this is an unrelated small-dense problem with its own natural solve, same
spirit as `FiberSection`'s hand-rolled 2x2 inversion needing no library).

```rust
/// Krylov-subspace acceleration of a stale-tangent Newton correction —
/// ported from `KrylovAccelerator.cpp`. A concrete struct, not a trait
/// object (§2.1/§2.6): exactly one accelerator implementation exists,
/// same reasoning as `SparseSolver` committing to one solver.
struct KrylovAccelerator {
    max_dimension: usize,
    v: Vec<DVector<f64>>,   // raw corrections
    av: Vec<DVector<f64>>,  // residual differences
}
impl KrylovAccelerator {
    fn new_step(&mut self); // reset subspace (called once per Analysis::step)
    fn accelerate(&mut self, v_raw: DVector<f64>, residual: &DVector<f64>) -> DVector<f64>;
    fn should_refactor(&self) -> bool; // subspace exceeded max_dimension
}
```

Exposed on `Algorithm` as its own variant (not a `LineSearch`-style
modifier field, since it fundamentally changes what "solve" means each
iteration, the same way Xara makes `AcceleratedNewton` its own algorithm
class rather than a flag on `ModifiedNewton`):

```rust
Algorithm::KrylovNewton { tangent: TangentStrategy, max_dimension: usize }
```

`tangent` is almost always `TangentStrategy::Initial` or
`ReuseAtStepStart` in practice (accelerating a `Current`-tangent Newton
that already re-factors every iteration defeats the point), but isn't
restricted to that — `AcceleratedNewton`'s constructor doesn't restrict it
either.

---

## 6. Revised `Algorithm` enum

```rust
#[derive(Debug, Clone, Copy)]
pub enum Algorithm {
    Linear,
    Newton {
        tangent: TangentStrategy,
        line_search: Option<LineSearch>,
    },
    KrylovNewton {
        tangent: TangentStrategy,
        max_dimension: usize,
    },
}
```

`NewtonRaphson` (today's bare unit variant) becomes
`Algorithm::Newton { tangent: TangentStrategy::Current, line_search: None }`
— a breaking rename, acceptable pre-M10 (no external API stability
promised yet, per §1's non-goals) but real: every existing
`Algorithm::NewtonRaphson` call site (tests, `wasm-bridge`) needs updating.
A `const NEWTON_RAPHSON: Algorithm = ...` or similar convenience isn't
worth it — the explicit struct literal is more honest about what's
actually configurable now, and every call site already spells out
`ConvergenceTest`'s fields the same way.

`Analysis::try_step`'s `match self.algorithm` grows the
`ReuseAtStepStart`/`Initial` tangent-caching logic and, for `Newton` with
`line_search: Some(_)`, the extra rescale-and-retest step between `update`
and the convergence check — both slot into the existing loop shape in
`state.rs` without needing a new control-flow structure, just more
branches inside it. `Algorithm::KrylovNewton` needs its own loop body
(closer to `AcceleratedNewton.cpp`'s shape than `NewtonRaphson`'s), since
"solve" now means "raw solve, then accelerate" rather than a single step.

---

## 7. Non-goals for this milestone

Same spirit as implementation-plan §1/§9 — named so a future session
doesn't rediscover and second-guess these:

- **No BFGS/Broyden.** Real quasi-Newton (secant condition on the inverse
  tangent operator itself, no factorization ever, not even an initial
  one) — a different family from the three axes here, not a
  generalization of them. `BFGS.cpp`/`Broyden.cpp` noted in §2 for when
  this becomes worth doing, not now.
- **No `KrylovAccelerator2`.** An alternate/improved variant of §5's
  accelerator exists in Xara; not adopted for the first port. Revisit only
  with a concrete case where `KrylovAccelerator` underperforms it.
- **No `Secant`/`InitialInterpolatedLineSearch`.** Bracketing strategies
  (`Bisection`/`RegulaFalsi`) cover the robustness case; these two trade
  robustness for a small speed gain, not needed until a specific model
  demands it.
- **No `TransientAnalysis` integration.** It has no Newton corrector at
  all yet (implementation-plan §6 M6 note) — out of scope here, follow-on
  work once it exists.
- **No reliability-analysis algorithms** (`FindDesignPointAlgorithm` and
  friends) — unrelated problem domain (§1 non-goals: no reliability
  analysis at all).

---

## 8. Staged implementation order

Verification for every stage: native `cargo test` then wasm32 build
(implementation-plan §7), same as every prior milestone. Each new
algorithm variant needs at least one test confirming it *converges to the
same answer* as `Algorithm::Newton { tangent: Current, line_search: None }`
on an already-covered nonlinear case (`m4_analysis.rs`'s or
`material_state.rs`'s `ZeroLength`+`ElasticPP` system, or
`m8_force_beam_column.rs`'s past-yield case) — same answer, different
(and ideally fewer, for the tangent-reuse/acceleration stages) iterations
or factorizations. That needs `StepResult` to grow `iterations: usize`
and, once §3.1's factorization split lands, `factorizations: usize` —
additive fields, not a breaking change to `StepResult` itself.

- **A1 — `TangentStrategy` enum + `Algorithm::Newton` restructure.**
  Rename `NewtonRaphson` → `Newton { tangent, line_search }`, add
  `TangentStrategy`, wire `ReuseAtStepStart`/`Initial` into
  `Analysis::try_step` *without* the factorization-caching optimization
  yet (§3.1 option 2's correctness-first version: still re-factor every
  iteration, just against whichever tangent the strategy says to form).
  Correct, not yet fast — matches this project's own precedent of landing
  correctness before the flagged optimization (e.g. M1's `SparseSolver`
  before its own re-factoring-across-iterations note in §5 decision #1).

- **A2 — `SparseSolver` factor/solve split + real caching.** `§3.1`'s
  `SparseFactorization` type; `Analysis` caches it per `TangentStrategy`.
  This is where `ReuseAtStepStart`/`Initial` actually start saving
  factorization cost — verify via the new `factorizations` counter on
  `StepResult`, not just "still converges."

- **A3 — `LineSearch` (`Bisection`, `RegulaFalsi`).** Needs a nonlinear
  test case where plain Newton without line search either fails to
  converge or takes visibly more iterations than with it — a coarse
  softening system (few fibers, like `ForceBeamColumn`'s own past-yield
  test had to be redesigned away from, ironically now useful *because*
  it's numerically harsh) is a natural candidate to construct
  deliberately for this test, not accidentally hit as a bug.

- **A4 — `Algorithm::KrylovNewton` + `KrylovAccelerator`.** The
  `nalgebra`-QR-based small least-squares solve (§5); verify against the
  same nonlinear case as A1-A3, checking `factorizations` stays near 1
  per step while `iterations` stays comparable to full Newton (the whole
  point: recover most of full Newton's convergence rate while paying
  `Initial`'s factorization cost).

---

# Part II: Event-to-event (EtE) stepping

A genuinely different algorithm class from everything in Part I — not a
Newton variant at all. OpenSees/Xara has no equivalent (confirmed: no
`EventToEvent`-anything anywhere under `xara/SRC/analysis/`), so unlike
every other milestone in this project there is no reference source to
port from. This section is designed from the published method and from
CSI's own (Perform-3D's vendor) public documentation of the same
technique in their other products, cited explicitly below — the same
"verify against a known-correct source, don't guess" discipline
implementation-plan §7 asks for, just aimed at literature instead of a
source tree.

## 9. What EtE is, and why it's worth having alongside Newton

For a structure built entirely from **piecewise-linear** force-deformation
relationships, the true response between any two "events" (points where
some component's tangent changes discontinuously — a fiber yielding, a
gap closing, a hysteretic backbone crossing a control point) is *exactly*
linear. That means the load factor (or displacement) increment to the
*next* event can be solved for directly, in one linear solve, with zero
iteration and zero convergence tolerance — not approximated by Newton and
then checked for convergence, but computed exactly. Advance to that exact
point, commit the triggered component(s), re-form the tangent (now
different, since something's regime changed), and repeat. For a model
that is genuinely piecewise-linear throughout, this produces the *exact*
solution to machine precision — a strictly stronger correctness property
than Newton-Raphson, which only guarantees convergence to within
`ConvergenceTest`'s tolerance.

This is the actual algorithmic basis of PERFORM-3D (per CSI's own
public documentation of the same technique in SAP2000/ETABS — Perform-3D
shares a numerical core with those products): PERFORM's component models
are deliberately *multi-linear by design* specifically so this strategy
applies exactly, not as an afterthought. It sidesteps needing
`ConvergenceTest`/iteration at all for the pure piecewise-linear case, and
it never overshoots a yield point the way a Newton step taken at a fixed
load increment can (a real robustness property, not just a performance
one — see the original paper's abstract, quoted in §10 below, on prior
event-to-event methods actually *violating* element force-deformation
limits).

## 10. The published method and CSI's production terminology

- **Karamchandani, A. and Cornell, C. A. (1992). "An Event-to-Event
  Strategy for Nonlinear Analysis of Truss Structures. I."** *Journal of
  Structural Engineering*, ASCE, 118(4), 895–909.
  [doi:10.1061/(ASCE)0733-9445(1992)118:4(895)](https://ascelibrary.org/doi/abs/10.1061/(ASCE)0733-9445(1992)118:4(895)).
  The originating paper for the *name* "event-to-event," restricted to
  truss structures (one uniaxial force-deformation relationship per
  element — the simplest possible case). Its abstract (confirmed via
  [TRID](https://trid.trb.org/view/365246)) is worth internalizing before
  implementing: it reviews *earlier* event-to-event approaches and shows
  they "lead to violations of element force-deformation characteristics
  and incorrect estimates of the response," then proposes a corrected
  method that "never violated" those characteristics while maintaining
  equilibrium, handles transient unloading before failure, and covers
  both proportional and nonproportional loading. **Read this paper (or
  its correction/caveats at minimum) before implementing §12's exact
  commit rule below** — this project doesn't have library access to the
  full text as of this writing, so the implementation stage must either
  obtain it or independently re-derive and stress-test the same
  correctness property (§14's test plan is designed to catch the failure
  mode the abstract describes: committing a component *past* its true
  event point).
- **CSI's own production terminology** (SAP2000/ETABS help, same vendor
  as Perform-3D):
  [Event-to-Event Stepping / Event-Lumping Tolerance](https://help.csiamerica.com/help/sap2000/26/26.0.0/SAP2000/WebHelp/Menus/Define/Load_Cases/Static/Nonlinear/Event_Lumping_Tolerance.htm).
  Confirms the exact terms this section adopts rather than inventing new
  ones: an **event** is "a change in stiffness of nonlinear elements";
  **Event-to-Event Stepping** subdivides load steps automatically at
  events; **Event-Lumping Tolerance (Relative)** "increase[s]" the load
  increment calculated to reach the first event "to include other events
  that would occur soon afterward" — larger tolerance means fewer event
  steps but more equilibrium unbalance; **Maximum Events per Step** caps
  how many events one step absorbs; **Minimum Event Step Size** and
  **Maximum Null Events per Step** guard against runaway tiny/zero-size
  event steps (their example: "interdependent hinges or catastrophic
  failure"); **"Event-to-Event Only"** is CSI's name for running this as
  the *sole* solution strategy (no Newton at all) versus using it only to
  subdivide steps within an otherwise-Newton scheme; **"Use Correction
  Step for Large Unbalance"** is CSI's own escape hatch back to iteration
  when pure event stepping accumulates too much unbalance — the same
  "optional final Newton cleanup" this plan proposes in §13.5, not an
  invention specific to this project.
- The user's own framing ("OpenSees has no equivalent... Perform3D does")
  is confirmed by the source search that grounded this section: no
  `EventToEvent`/`EtE` class exists anywhere in `xara/SRC/analysis/` (the
  same tree M12's Part I ported five algorithm classes from).

## 11. Scope: extending a truss-only method to fiber sections

Karamchandani & Cornell's method covers trusses — each element is a
single spring, so "the element's event" and "the structure's event" are
the same query repeated per element. Carapace's actual nonlinear elements
(`Truss`, `ZeroLength`, `DispBeamColumn`, `ForceBeamColumn`) are richer:
a `DispBeamColumn`/`ForceBeamColumn` has several `FiberSection`s, each
with several `Fiber`s, each carrying its own `Material` with its own
independent strain history. **Generalizing "the next event" from one
spring to a fiber-discretized section-integration element is this
project's own extension, not something the cited paper addresses** — call
this out explicitly in code comments when it's implemented, the same way
`ForceBeamColumn`'s doc comment is honest about which parts are a direct
port and which are this project's own derivation.

The generalization itself is mechanical, though, once stated precisely
(§12–13): a fiber's strain is a *linear function* of the section's
`(eps0, kappa)` deformation, which is in turn a *linear function* of the
element's nodal displacement, which is in turn a *linear function* of the
global displacement vector `K^-1 * reference_load` for a **fixed**
tangent `K`. Every one of those maps is linear between events by
construction (that's what "between events" *means* — nothing's regime
has changed yet), so "how far along the global unit response before some
fiber, somewhere, in some element, crosses a breakpoint" is a single
linear query, however many layers of indirection it passes through.

## 12. Per-material event query

New required method on the closed `Material` enum (§2.1 — every existing
variant gets a match arm, no default/fallback):

```rust
/// Given the material's *committed* strain and a trial strain-rate
/// direction (its sign is all that matters — magnitude is normalized
/// away), how much additional strain (always >= 0, in the direction of
/// `rate`'s sign) can be applied before this material's tangent changes
/// discontinuously? `None` means "no such point exists in this
/// direction" — either the material is genuinely always linear
/// (`Elastic`), or it's currently on an infinite plateau in that
/// direction (e.g. `ElasticPP` already yielded, still loading the same
/// way — zero tangent forever until a reversal, which is itself only
/// detected by `rate` changing sign on a *later* call, not by this one).
fn distance_to_event(&self, rate_sign: f64) -> Option<f64>;
```

Called with the material's own **current committed state already known
internally** (no `strain` parameter needed — unlike
`trial_stress_tangent(strain)`, this isn't evaluating a hypothetical
trial point, it's asking "starting from where I actually am, how far in
this direction"). Match arms by material family:

- **Exactly piecewise-linear leaves** (`Elastic` excepted — see below):
  `ElasticPP`, `Gap`, `Ent`, `Steel01`, `Hysteretic`, `Pinching4`. Each
  has a well-defined next breakpoint strain given its own committed state
  and history variables (`ElasticPP`: the yield-strain boundary at
  `ep ± eyp` relative to committed strain, `None` if already on the
  plateau in `rate_sign`'s direction; `Steel01`'s isotropic-hardening
  bookkeeping already tracks exactly the state this needs; `Hysteretic`/
  `Pinching4` are *defined* by explicit control points, so this is a
  direct lookup, arguably the easiest arms to write). Concretely: `Some`
  distance, always, except at an already-engaged infinite plateau.
- **`Elastic`**: `None` always — never a breakpoint, by definition.
- **Genuinely smooth curves**: `Steel02` (Menegotto-Pinto), `Concrete01`,
  `Concrete02`. `None` always, honestly — there is no exact finite-strain
  discontinuity to report. **This is the open design question — see §13.4.**
- **Composites** (`Parallel`, `Series`, `MinMax`): recurse, but not
  trivially — see below.
  - `Parallel`: every child sees the *same* strain (parallel springs
    share deformation), so `distance_to_event` is simply the minimum over
    children's own queries at the same `rate_sign`. Straightforward.
  - `Series`: every child sees the *same force*, strains sum — a
    child's *strain*-space event distance doesn't directly compose the
    way `Parallel`'s does, because the external deformation needed to
    reach a given internal force depends on that child's own
    (potentially different) current tangent. This needs the same
    flexibility-summation idea `ForceBeamColumn`'s section integration
    already uses (`f_total = sum(f_child)`, distance-to-event converted
    through the *series flexibility*, not compared strain-for-strain) —
    a real, self-contained sub-problem to design carefully at
    implementation time, not a two-line recursive call. Flagged as an
    open item, not resolved here.
  - `MinMax`: tracks which bound (min or max envelope) is currently
    active; needs both that active child's own event *and* the event of
    the *other* bound becoming newly binding (an envelope crossing) —
    another real sub-problem, not resolved here.

## 13. Aggregating events across a fiber section, element, and structure

### 13.1 `FiberSection`

```rust
/// Given the section's (eps0, kappa) deformation *rate* (a unit
/// direction, from the current linear predictor — see §13.3), the
/// minimum positive scalar `lambda` such that some fiber's strain
/// reaches a breakpoint, and which fiber it is (needed so `commit_event`
/// can commit exactly that fiber's material at exactly that point,
/// without waiting for a full section-level commit).
fn distance_to_event(&self, eps0_rate: f64, kappa_rate: f64) -> Option<(f64, usize)> {
    // Fiber i's strain rate is eps0_rate - y_i * kappa_rate (same linear
    // map `trial` already uses for the strain itself); ask that fiber's
    // Material for its own distance at that rate's sign, take the min.
}
```

### 13.2 `DispBeamColumn` / `ForceBeamColumn`

`DispBeamColumn`: direct — `strain_displacement`'s `(b_eps0, b_kappa)`
already gives each integration point's `(eps0, kappa)` as a linear
function of local nodal displacement, so the *rate* is just that same
linear map applied to the nodal displacement *rate* (which, in turn, is
`d_local_rate = T * d_global_rate`, `d_global_rate` being the relevant
slice of the structure-wide unit response — see §13.3). No new math, only
plumbing through machinery `DispBeamColumn::form_tangent_and_resistance`
already has.

`ForceBeamColumn`: **not** as direct — being force-based, its section
deformations aren't a simple linear function of nodal displacement the
way `DispBeamColumn`'s are; they come out of `state_determination`'s
Newton solve. Between events, though, `k_basic`/each section's `f_sec`
*are* fixed (that's what "between events" means), so the *rate* of each
section's `(eps0, kappa)` per unit change in the *basic deformation rate*
`v_rate` is `f_sec_i * (b(xi_i) * (k_basic * v_rate))` — a direct
consequence of the same linearization `try_state_determination` already
computes and discards every iteration (`f_secs[i]`, `k_basic` are already
right there in scope) — needs a dedicated method that reuses the last
converged `f_secs`/`k_basic` from `state_determination`'s final iteration
rather than redoing the Newton solve, which is a real, self-contained
piece of new plumbing (not a redesign of the element), flagged as its own
implementation task.

`Truss`/`ZeroLength`: single (or per-DOF) `Material`, so this is exactly
`Material::distance_to_event` called with that DOF's own strain rate sign
— the simplest case, matching the original truss-only paper directly.

`ElasticBeamColumn`: no `Material` at all (§3.1 of the implementation
plan) — always contributes `None`.

### 13.3 `Domain`

```rust
/// One linear solve for the whole structure's unit response to the
/// reference load pattern at the *current* (fixed, already-committed)
/// tangent, then the minimum event distance across every element's own
/// query (§13.2) converted to a *load factor* distance via that same
/// unit response. `None` if nothing in the structure has a pending event
/// in this direction (e.g. a fully-elastic model, or every remaining
/// component already on an infinite plateau) — the caller (`Analysis`)
/// then knows to fall back or take the full remaining step in one shot.
pub(crate) fn distance_to_next_event(&self, k: &SparseMatrix, du_hat: &DVector<f64>) -> Option<f64>;
```

### 13.4 Open decision: components with no exact event (smooth materials)

If *every* remaining active component reports `None` (an all-`Steel02`/
`Concrete01`/`Concrete02` model, or a piecewise-linear model where
everything left is already on an infinite plateau), pure EtE has nothing
to stop at before the step's full target — which is exactly wrong for a
smooth curve (it would take the *entire* remaining load increment in one
linear shot, silently ignoring the curve's curvature the same way
`Algorithm::Linear` would). Two honest options, not resolved here:

1. **Synthetic curvature-limited pseudo-events.** Give smooth materials a
   *non-`None`* `distance_to_event` too: estimate how far the strain can
   move before the *local tangent* changes by more than some tolerance
   (evaluate the tangent at the candidate endpoint, compare to the
   tangent at the start, bisect/scale to hit the tolerance), and report
   that as a pseudo-event distance. Keeps the *whole* structure inside
   one uniform EtE loop, at the cost of `Steel02`/`Concrete01/02` needing
   a genuinely new (and somewhat arbitrary, tolerance-parameterized)
   method that has no analogue in their existing `trial_stress_tangent`.
2. **Pre-tessellate smooth backbones into real piecewise-linear
   materials**, once, at construction time (a fixed or
   caller-specified resolution) — closer to what PERFORM-3D actually
   does by only ever offering multi-linear material definitions in the
   first place. Keeps `distance_to_event`'s contract uniformly exact for
   every material (no tolerance parameter smuggled into a "physics"
   method), at the cost of a real accuracy/resolution trade-off decided
   once at model-build time instead of adaptively during the solve, and
   doesn't help a model that specifically wants `Steel02`'s smooth
   Menegotto-Pinto shape rather than an approximation of it.

Recommendation: start with (1), scoped *only* to `Steel02`/`Concrete01`/
`Concrete02` (the three smooth leaves), with the tolerance as a field on
`Algorithm::EventToEvent` (§13.5) rather than per-material — keeps the
material catalog's existing contract (`trial_stress_tangent`/`commit`,
implementation-plan §3.2) otherwise untouched, confines the new concept
to the one place it's actually needed, and stays reversible: if (1) turns
out to behave badly in practice, (2) is still available without having
threaded a tolerance concept through `Material` itself.

### 13.5 `Algorithm::EventToEvent`

```rust
Algorithm::EventToEvent {
    /// CSI's "Event-Lumping Tolerance (Relative)": events within this
    /// fraction of the nearest one's distance are absorbed into the same
    /// step (§14) instead of triggering their own tangent re-formation.
    lumping_tolerance: f64,
    /// CSI's "Maximum Events per Step" — a lumped group larger than this
    /// forces a stop and re-formation regardless of `lumping_tolerance`.
    max_events_per_step: usize,
    /// §13.4 option 1's relative-tangent-change tolerance for smooth
    /// materials' pseudo-events. `None` disables it — only valid for a
    /// model verified to be exactly piecewise-linear throughout.
    smooth_material_tolerance: Option<f64>,
    /// CSI's "Use Correction Step for Large Unbalance": an optional
    /// final Newton pass (reusing `Algorithm::Newton`'s own loop, not a
    /// new implementation) once the step's target load factor is
    /// reached, closing residual drift from lumping and/or
    /// `smooth_material_tolerance`'s approximation. `None` is only safe
    /// for a verified exactly-piecewise-linear, unlumped model.
    correction: Option<ConvergenceTest>,
},
```

## 14. Event lumping, precisely

The key property worth stating explicitly (not obvious from CSI's
one-paragraph description, but follows directly from EtE's own math):
**lumping does not have to sacrifice each individual component's commit
exactness**, only the tangent's freshness in between. Concretely, within
one lumped group:

1. Solve once for `du_hat` (the unit response at the tangent formed at
   the *start* of the group).
2. Find `lambda_min`, the smallest individual event distance, and every
   other event whose own distance falls within
   `lambda_min * (1 + lumping_tolerance)` (capped at
   `max_events_per_step` events total).
3. Commit **each** of those triggered components at **its own** exact
   `lambda_i` (not the group's `lambda_min` or its max) — since `du_hat`
   was valid (nothing regime-changed) all the way out to the *largest*
   `lambda_i` in the group, every individual commit point is still
   exact, computed off the same single linear solve.
4. Advance the load factor by `max(lambda_i in the group)`, *not*
   `lambda_min` — components between `lambda_min` and that max were
   already committed correctly in step 3; what's approximate is only
   that the *tangent* stayed stale (reflecting the pre-group state, not
   the just-committed intermediate events) for that whole span, which is
   exactly the "may increase the equilibrium unbalance" caveat CSI's own
   docs give for a larger tolerance — an honest, bounded, named
   approximation, not silent drift.

This is precisely why §13.5's `correction` field matters for anything but
a verified-exact model: lumping's unbalance is real (if usually small)
and needs somewhere to go.

## 15. Testing / verification plan

Same "closed-form or oracle" discipline as implementation-plan §7, with
EtE's own correctness property giving a *stronger* oracle than usual:

- **Exact match, not tolerance match.** For a model built entirely from
  exactly-piecewise-linear materials (reuse `material_state.rs`'s
  `Truss` + `ZeroLength`/`ElasticPP` system, or a `Steel01`/`Hysteretic`/
  `Pinching4` variant of it), `Algorithm::EventToEvent` with
  `lumping_tolerance: 0.0` and `correction: None` must match a
  tightly-converged `Algorithm::Newton` run to a *much* tighter tolerance
  than Newton's own `ConvergenceTest` tolerance (e.g. `1e-13`, not
  `1e-9`) — proof it's actually solving exactly, not just "also
  converging."
- **Lumping's bounded error.** The same system, run with
  `lumping_tolerance: 0.0` vs. a nonzero value chosen to force at least
  two events into one group (needs a model with two components whose
  true event distances are close but not identical — a small,
  deliberately-constructed asymmetric two-fiber or two-spring case, not
  an accidental one) — the lumped run's final state should differ from
  the unlumped run's by an amount that shrinks as `lumping_tolerance`
  shrinks, and `correction: Some(_)` should close that gap. Needs a
  reformation/event counter exposed on `StepResult` (same idea as Part
  I's `factorizations` field) to also confirm lumping actually reduced
  the count, not just that the answer stayed close.
- **Smooth-material fallback.** A `Steel02`- or `Concrete01`-only model,
  `Algorithm::EventToEvent` with `smooth_material_tolerance: Some(_)` vs.
  a plain `Algorithm::Newton` run — expect ordinary `ConvergenceTest`-
  tolerance agreement here (not the exact-match bar above; §13.4 is an
  approximation by construction), and expect the answer to visibly
  degrade as `smooth_material_tolerance` is loosened, confirming the
  parameter actually does something rather than being ignored.
- **Never-violates-limits regression.** Directly targeting the original
  paper's own stated failure mode of *prior* event-to-event methods
  (§10): a test that pushes a component to (and numerically just past,
  by construction) its exact breakpoint and asserts the committed strain
  never exceeds the material's own valid domain for the regime it was
  committed into — e.g. `ElasticPP`'s committed strain never implies a
  stress magnitude exceeding `fy` by more than solver-precision.

## 16. Non-goals for this milestone

- **`DisplacementControl` + EtE.** The same unit-response scaling
  generalizes (stop at "target DOF displacement reached" instead of
  "target load factor reached," using identical event-detection
  machinery), but is real follow-on work, not core MVP — start with
  `LoadControl` only.
- **`Series`/`MinMax` composite event queries.** Flagged as open
  sub-problems in §12, not resolved here — a `Parallel`-only composite
  story is enough for a first, honestly-scoped implementation.
- **Nonproportional / multi-load-pattern loading.** Karamchandani &
  Cornell's Part I (and evidently a Part II — not tracked down here)
  address this; Carapace's own `Integrator`s are single-reference-load-
  pattern only today (implementation-plan §3.3), so this is naturally
  deferred until that's no longer true.
- **`TransientAnalysis` integration.** Out of scope here for the same
  reason as Part I's non-goals (no Newton corrector there at all yet) —
  worth flagging, though, that EtE is *arguably more* valuable for
  transient response than static (avoids any iteration inside a time
  step entirely for a piecewise-linear model), a strong motivation for
  revisiting once `TransientAnalysis` grows nonlinear support, not a
  reason to build it now.

## 17. Staged implementation order

- **E1 — `Material::distance_to_event` for the exactly-piecewise-linear
  leaves + `Parallel`.** `Elastic`/`ElasticPP`/`Gap`/`Ent`/`Steel01`/
  `Hysteretic`/`Pinching4`, plus `Parallel`'s straightforward recursion.
  `Steel02`/`Concrete01`/`Concrete02`/`Series`/`MinMax` all return `None`
  for now (correctly — they simply don't participate yet, not a
  half-finished implementation of an unresolved design). Unit-testable
  in complete isolation, same style as `FiberSection`'s own existing
  unit tests.
- **E2 — `FiberSection`/`Truss`/`ZeroLength`/`DispBeamColumn` event
  queries + `Domain::distance_to_next_event`.** `ForceBeamColumn`
  deliberately excluded from this stage (§13.2's harder sub-problem) —
  an EtE analysis restricted to `Truss`/`ZeroLength`/`DispBeamColumn`
  models is still a real, testable milestone on its own.
- **E3 — `Algorithm::EventToEvent` core loop**, unlumped
  (`lumping_tolerance: 0.0`), no smooth-material support
  (`smooth_material_tolerance: None`), `correction: None`. The exact-match
  test from §15 becomes possible here — this is the milestone's real
  correctness bar.
- **E4 — Event lumping** (§14) + the reformation/event counter on
  `StepResult`.
- **E5 — `ForceBeamColumn` event queries** (§13.2's harder case) +
  §13.4's smooth-material decision, whichever option §13.4 lands on, once
  it's actually been decided rather than deferred.
