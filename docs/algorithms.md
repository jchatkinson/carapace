# Solver algorithms and JavaScript configuration

Static and transient analysis support linear solves, Newton iteration with
current/reused/initial tangents, bisection or regula-falsi line search, and
Krylov-accelerated Newton. The same Rust iteration driver serves both analysis
types. These options are exposed through `CarapaceInputV1` and its generated
TypeScript definitions, in both the planar and spatial profiles.

Modal analysis is an eigenvalue solve; it does not use these iteration settings.
Event-to-event stepping is still unimplemented.

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
object forms, validation, and transient convergence failure.

The JS interface still lacks initial nodal displacement/velocity input,
iteration/factorization diagnostics, and custom Series-material tolerance
settings. The pysees compiler currently emits only the legacy static
algorithm strings; exposing these settings in its UI is separate work.

The [archived algorithm design](obsolete/algorithm-design.md) preserves the
original derivations, Xara references, section numbering used by source
comments, and implementation sequence. Its Part I is implemented. Part II
(§9–§17) describes the remaining event-to-event stepping proposal.
