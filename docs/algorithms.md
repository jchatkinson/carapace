# Solver algorithms and JavaScript configuration

Static and transient analysis support linear solves, Newton iteration with
current/reused/initial tangents, bisection or regula-falsi line search, and
Krylov-accelerated Newton. The same Rust iteration driver serves both analysis
types. These options are exposed through `CarapaceInputV1` and its generated
TypeScript definitions, in both the planar and spatial profiles.

Modal analysis is an eigenvalue solve; it does not use these iteration settings.
Event-to-event stepping is still unimplemented.

Static stages can also use arc-length continuation
(`integrator: { kind: "arcLength", ... }`, see [below](#arc-length-continuation)
and the [design record](arclength.md)) to trace equilibrium paths through
load limits, softening and snap-back.

## Input compatibility and defaults

The `algorithm` field accepts the original strings and configurable objects:

| Input | Behavior |
| --- | --- |
| `"linear"` or `{ kind: "linear" }` | One solve per step; no convergence iteration |
| `"newtonRaphson"` | Newton with the current tangent and no line search |
| `{ kind: "newton", ... }` | Newton with configurable tangent and optional line search |
| `{ kind: "krylovNewton", tangent, maxDimension }` | Newton with Krylov acceleration |

Static stages continue to require `algorithm`. Transient stages may omit it,
preserving the original linear Newmark solve. Both stage kinds accept an
optional `convergence` setting; the default is
`{ kind: "normUnbalance", tol: 1e-6, maxIter: 20 }`.
Linear solves ignore convergence settings. For nonlinear dynamic response,
select an iterative algorithm explicitly.

## Tangent strategy and line search

`tangent` is one of:

- `"current"`: rebuild and factor the effective tangent every iteration.
- `"reuseAtStepStart"`: factor once per step and reuse during that step.
- `"initial"`: reuse the initial tangent across steps within that stage.

Newton defaults to `"current"` when `tangent` is omitted. Krylov Newton
requires an explicit tangent and a positive integer `maxDimension`; using
`"reuseAtStepStart"` or `"initial"` trades factorization work for iteration.
Stage transitions construct a new analysis, so tangent caches reset there.

Newton's optional `lineSearch` is either
`{ kind: "bisection", tol, maxIter, maxEta }` or
`{ kind: "regulaFalsi", tol, maxIter, maxEta }`.
Krylov Newton has no line-search field, matching the Rust algorithm enum.
Line-search tolerance must be finite and positive, `maxIter` a positive
integer, and `maxEta` finite and at least 1.

Convergence kinds are `"normUnbalance"`, `"normDispIncr"`, and `"energyIncr"`;
each takes a finite positive `tol` and positive integer `maxIter`. Transient
`dt` must also be finite and positive. Counts in
the input format are `u32`. Absolute residual tolerances should account for
the force scale and floating-point cancellation in the effective dynamic
system. A tighter tolerance is not necessarily achievable for every model.

`{ kind: "combined", forceTol, momentTol?, relativeTol?, displacementTol?,
maxIter }` checks force equilibrium on every equation, and every configured
criterion must pass:
`|r_i| <= forceTol (momentTol on rotational equations) + relativeTol * F_ref_i`,
where `F_ref_i` is the largest of the committed external force, the
predicted external force, and the committed internal force on that
equation. `momentTol` defaults to `forceTol` and `relativeTol` to `1e-6`.
The optional `displacementTol` additionally bounds the largest component of
the last correction. Unlike a single norm, the per-equation check can't be
diluted by many quiet DOFs, and unlike `normDispIncr`/`energyIncr` it
certifies equilibrium.

Every test is evaluated at the **accepted state** of an iteration: the
residual is reassembled after the correction (and any line-search
rescaling) is applied. Earlier versions checked the residual the
correction was solved against. A linear system therefore converges on its
single solve, with no second confirming iteration, and a step whose
predictor is already balanced takes no solve at all under the force tests.
A zero `maxIter` is rejected.

## Examples

These objects go in `sequence.stages` or `sequence3.stages` of an otherwise
complete `CarapaceInputV1` payload:

```typescript
import type { StageSpec } from "./pkg/carapace_wasm.js";

const staticStage = {
  kind: "static",
  id: "pushover",
  steps: 100,
  integrator: { kind: "displacementControl", node: 1, dof: 0, increment: 0.01 },
  algorithm: {
    kind: "newton",
    tangent: "current",
    lineSearch: { kind: "bisection", tol: 1e-10, maxIter: 30, maxEta: 16 },
  },
  convergence: { kind: "normUnbalance", tol: 1e-6, maxIter: 100 },
  holdPatternsAfter: [],
} satisfies StageSpec;

const transientStage = {
  kind: "transient",
  id: "earthquake",
  steps: 1000,
  dt: 0.01,
  damping: { alphaM: 0.1, betaK: 0 },
  groundMotions: [], // supply acceleration time series here when needed
  algorithm: { kind: "krylovNewton", tangent: "reuseAtStepStart", maxDimension: 3 },
  convergence: { kind: "normUnbalance", tol: 1e-6, maxIter: 100 },
} satisfies StageSpec;
```

## Arc-length continuation

`{ kind: "arcLength", initialRadius, scales, ... }` follows the equilibrium
path with a scaled spherical constraint, correcting displacement and load
factor together, so load limits, negative stiffness, snap-through and
snap-back pass without special handling. A stage's `steps` is its step
cap; the optional `stop` criteria can end it earlier. The full design,
derivations and verification are in [arclength.md](arclength.md).

| Field | Meaning |
| --- | --- |
| `initialRadius`, `minRadius?`, `maxRadius?` | Dimensionless radius; omitted bounds default to `initialRadius` (fixed radius) |
| `scales` | `{ kind: "explicit", displacement, rotation?, load }` (rotation required only when rotational DOFs are active: a model of trusses and membranes has none, so it can be omitted) or `{ kind: "auto", load }` (translation/rotation scales from the first elastic tangent) |
| `direction?` | `"increasing"` (default) or `"decreasing"` initial load direction |
| `predictor?` | `"secant"` (default; the tangent is still used after a cutback or sharp turn) or `"tangent"` |
| `seed?` | `{ components: [{ node, dof, value }], load }` orienting a singular first tangent |
| `targetIterations?`, `maxRetries?` | Radius adaptation target (default 6) and cutbacks per step (default 8) |
| `arcTolerance?`, `correctionTolerance?` | Constraint tolerance `|g|/s^2` (default `1e-8`); optional bound on the last coupled correction relative to the radius |
| `backtracking?` | `{ armijo, minStep }` (defaults `1e-4`, `1/128`) |
| `stop?` | `{ displacement?: { node, dof, value, exact? }, loadFactor?: { value, exact? }, loadFactorZeroCrossing?, maxChordLength? }` |

Arc length requires `"newtonRaphson"` or `{ kind: "newton" }` with the
current tangent and no line search (its own coupled backtracking replaces
line search), and a `normUnbalance` or `combined` test. When `convergence`
is omitted it defaults to `{ kind: "combined", forceTol: 1e-6, maxIter: 30 }`.
The stage must start in equilibrium, and active load patterns must be
`linear`, `constant` or held constant (an active `path` series is
rejected). Each step reports `StepOutcome.continuation`: radius, retries,
corrector solves, factorizations, force and constraint measures,
accumulated chord length, the sign of `det K`, `bifurcationSuspected` when
that sign changed without a load reversal, and which stop criterion ended
the stage. The load factor can fall and repeat during continuation, so
recorder samples are ordered by their sample index, never by load factor.

```typescript
const pushover = {
  kind: "static",
  id: "softening pushover",
  steps: 500,
  integrator: {
    kind: "arcLength",
    initialRadius: 0.05, minRadius: 0.005, maxRadius: 0.2,
    scales: { kind: "explicit", displacement: 0.01, rotation: 0.01, load: 100 },
    stop: { displacement: { node: 12, dof: 0, value: 0.5, exact: true } },
  },
  algorithm: "newtonRaphson",
  holdPatternsAfter: [],
} satisfies StageSpec;
```

Runtime arc-length failures are structured too: `initialStateNotInEquilibrium`,
`unsupportedLoadSeries`, `zeroLoadSensitivity`, `missingSeedDirection`,
`cutbacksExhausted { step, attempts, radius, lastFailure }`.

## Errors

Malformed shapes throw during input deserialization. Invalid numeric solver
settings throw a structured `DecodeError` with
`{ kind: "invalidAnalysisOption", stage, field }`. Failure to converge while
stepping appears in `StepOutcome.error`; it stops the sequence and remains
set on subsequent calls. Already-produced recorder batches remain available.

## Verification and remaining gaps

Native tests in `core/tests/m_algorithm_richness.rs` and
`core/tests/m21_transient_corrector.rs` verify the solver implementations.
`wasm-bridge/tests/m10_carapace_input_v1.rs` verifies configurable algorithms
against nonlinear static and transient equilibrium. The generated-package
test `wasm-bridge/tests/boundary-smoke.ts` covers both profiles, legacy inputs,
object forms, validation, transient convergence failure, and arc length.
Arc-length continuation is verified by `core/tests/arclength_softening.rs`,
`core/tests/arclength_continuation.rs` and
`wasm-bridge/tests/arclength_session.rs`.

The JS interface still lacks initial nodal displacement/velocity input,
iteration/factorization diagnostics for the non-arc-length integrators,
and custom Series-material tolerance settings. The pysees compiler currently emits only the legacy static
algorithm strings; exposing these settings in its UI is separate work.

The [archived algorithm design](obsolete/algorithm-design.md) preserves the
original derivations, Xara references, section numbering used by source
comments, and implementation sequence. Its Part I is implemented. Part II
(§9–§17) describes the remaining event-to-event stepping proposal.
