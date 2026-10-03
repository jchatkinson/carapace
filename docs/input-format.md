# Model input format (`CarapaceInputV1`)

Reference for the JavaScript/TypeScript model API: the payload passed to
`decodeInput`, the session it returns, and the result and error shapes.

The Rust types in `wasm-bridge/src/input_v1/` are the source of truth. They
derive `serde` + `tsify::Tsify`, and this page summarizes the TypeScript
definitions generated from them. If this page and the types disagree, the types
win. To see the exact declarations:

```bash
wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg
# -> pkg/carapace_wasm.d.ts  (generated, gitignored, do not edit)
```

Field-level semantics live in the Rust doc comments, which tsify copies into the
`.d.ts`. `wasm-bridge/tests/boundary-smoke.ts` is a type-checked, working
example payload.

## Conventions

- **Naming.** All keys are `camelCase`. Tagged unions use a discriminator
  (`kind`, or `response` for recorders) with `camelCase` values.
- **Columnar tables.** Each table is a struct of parallel arrays, one entry per
  row (`nodeI[k]`, `nodeJ[k]`, `e[k]` all describe element `k`). Every array in
  a table must have the same length unless noted as *sparse* or *offset-indexed*.
- **Indices** are 0-based row indices into the named table. They are validated
  by the decoder and failures surface as a `DecodeError`.
- **Strided arrays.** Some flat arrays pack several values per row, for example
  `NodeTable.coords` (stride 2 in 2D, stride 3 in 3D).
- **Sparse tables** are lists of tuples keyed by a row index, such as
  `[row, dof, materialIndex]`, because most rows have no entry.
- **Offset-indexed tables** (fibers) flatten variable-length groups: group `k`
  occupies `sectionOffsets[k]..sectionOffsets[k+1]`, so `sectionOffsets` has
  `numSections + 1` entries.
- **2D and 3D profiles.** `header.space` selects 2D (`2`: 2 coordinates, 3 DOFs
  per node `ux, uy, rz`) or 3D (`3`: 3 coordinates, 6 DOFs
  `ux, uy, uz, rx, ry, rz`). Every table is always present, but the decoder
  reads only the set matching `header.space`: unsuffixed tables for 2D,
  `*3` tables for 3D. The other set is ignored (use empty tables).
  `materials`, `loadPatterns` and `nodalLoads` are shared; `nodalLoads3` reuses
  the `NodalLoadTable` type.
- **Plain JS values.** The boundary accepts ordinary objects and arrays. A
  transferable typed-array format is planned but not implemented.

## Entry points

```ts
decodeInput(value: CarapaceInputV1): WasmSession   // throws DecodeError
session.advance(stepBudget: number): StepOutcome
session.currentStageId(): string | undefined       // undefined once all stages done
session.free()
```

`advance` runs up to `stepBudget` steps and returns, so a caller can run it in a
loop and keep a worker responsive. Once `error` is set it is sticky: later
stages are not attempted.

The other exports (`axial_displacement`, `zero_length_ent_displacement`, ...)
are milestone verification functions, not part of the model API.

## `CarapaceInputV1`

| Field | Type | Profile |
|---|---|---|
| `header` | `Header` | both |
| `materials` | `MaterialSpec[]` | both (shared arena) |
| `loadPatterns` | `LoadPatternTable` | both |
| `nodalLoads` / `nodalLoads3` | `NodalLoadTable` | 2D / 3D |
| `nodes` / `nodes3` | `NodeTable` / `NodeTable3` | 2D / 3D |
| `fibers` / `fibers3` | `FiberTable` / `FiberTable3` | 2D / 3D |
| `trusses` / `trusses3` | `TrussTable` / `TrussTable3` | 2D / 3D |
| `elasticBeamColumns` / `…3` | `ElasticBeamColumnTable` / `…3` | 2D / 3D |
| `dispBeamColumns` / `…3` | `FiberBeamColumnTable` / `…3` | 2D / 3D |
| `forceBeamColumns` / `…3` | `FiberBeamColumnTable` / `…3` | 2D / 3D |
| `zeroLengths` / `…3` | `ZeroLengthTable` / `…3` | 2D / 3D |
| `zeroLengthSections` / `…3` | `ZeroLengthSectionTable` / `…3` | 2D / 3D |
| `equalDofs` / `…3` | `EqualDofTable` / `…3` | 2D / 3D |
| `rigidDiaphragms` / `…3` | `RigidDiaphragmTable` / `…3` | 2D / 3D |
| `elementLoads` / `elementLoads3` | `ElementLoadTable` / `…3` | 2D / 3D |
| `sequence` / `sequence3` | `SequenceSpec` / `SequenceSpec3` | 2D / 3D |

### `Header`

| Field | Type | Notes |
|---|---|---|
| `schemaVersion` | `number` | Wire format version, independent of the engine. |
| `space` | `number` | `2` or `3`; anything else is `unsupportedSpace`. |
| `engineVersion` | `string` | `carapace-core` version, recorded for provenance. |

## Nodes

`NodeTable` (2D) / `NodeTable3` (3D):

| Field | 2D | 3D |
|---|---|---|
| `coords` | stride 2 `(x, y)` | stride 3 `(x, y, z)` |
| `fixed` | one bitmask per node: bit 0 `ux`, 1 `uy`, 2 `rz` | bits 0..5 `ux, uy, uz, rx, ry, rz` |
| `massNodeIndex` | node indices that carry mass (sparse) | same |
| `mass` | stride 3 `(mx, my, mrz)`, parallel to `massNodeIndex` | stride 6, parallel to `massNodeIndex` |

## Materials

`materials: MaterialSpec[]` is one arena shared by all elements and fibers.
Elements and fibers refer to entries by arena index. Composite materials refer
to other arena entries; cycles are rejected (`cyclicMaterialReference`).

| `kind` | Parameters |
|---|---|
| `elastic` | `e` |
| `elasticPp` | `e`, `eyp` |
| `gap` | `e`, `gap` |
| `ent` | `e` (no tension) |
| `steel01` | `fy`, `e0`, `b`, `a1`..`a4` |
| `steel02` | `fy`, `e0`, `b`, `r0`, `cr1`, `cr2`, `a1`..`a4` |
| `concrete01` | `fpc`, `epsc0`, `fpcu`, `epscu` |
| `concrete02` | `fc`, `epsc0`, `fcu`, `epscu`, `rat`, `ft`, `ets` |
| `parallel` | `children: [materialIndex, factor][]` |
| `series` | `children: materialIndex[]` |
| `minMax` | `inner`, `minStrain`, `maxStrain` |
| `hysteretic` | `mom1p`/`rot1p`..`mom3p`/`rot3p`, `mom1n`/`rot1n`..`mom3n`/`rot3n`, `pinchX`, `pinchY`, `damfc1`, `damfc2`, `beta` |
| `pinching4` | `stress1p`/`strain1p`..`4p`, `…1n`..`4n`, `rDispP/N`, `rForceP/N`, `uForceP/N`, `gammaKParams`/`gammaDParams`/`gammaFParams` (4-tuples) with `gammaKLimit`/`gammaDLimit`/`gammaFLimit`, `gammaE`, `dmgCyc: "energyBased" \| "cycleBased"` |

## Elements

All element tables have `nodeI` / `nodeJ` (indices into the node table) and,
where present, a `density` array for mass.

| Table | Fields beyond `nodeI`, `nodeJ`, `density` |
|---|---|
| `TrussTable` / `TrussTable3` | `area`, `material` (arena index) |
| `ElasticBeamColumnTable` | `e`, `a`, `iz`, `transform: TransformSpec[]` |
| `ElasticBeamColumnTable3` | `e`, `g`, `a`, `j`, `iy`, `iz`, `transform: TransformSpec3[]` |
| `FiberBeamColumnTable` (disp and force) | `fiberSection` (index into `FiberTable.sectionOffsets`), `integration: IntegrationSpec[]`, `corotational: boolean[]` |
| `FiberBeamColumnTable3` (disp and force) | `g`, `j` (decoupled elastic torsion), `vecXz: [x,y,z][]`, `fiberSection`, `integration` |
| `ZeroLengthTable` | `materials: [row, dof, material][]` (sparse), `friction: [row, normalDof, shearDof, mu, k0, b][]` (at most one per row) |
| `ZeroLengthTable3` | same, with `friction: [row, normalDof, shearDof0, shearDof1, mu, k0, b][]` |
| `ZeroLengthSectionTable` | `fiberSection`, `materials: [row, dof, material][]` for DOFs the section does not drive (`uy`) |
| `ZeroLengthSectionTable3` | same; extra springs for `uy`, `uz`, `rx` |

Notes:

- `TransformSpec`: `"linear" | "pDelta" | "corotational"`.
  `TransformSpec3`: `{ kind: "linear3" | "pDelta3", vecXz: [x,y,z] }`. There is
  no 3D corotational transform. `vecXz` is a vector not parallel to the
  member axis; it fixes the local y/z orientation.
- `IntegrationSpec`: `{ kind: "legendre" | "lobatto", points }`.
- 3D fiber elements have no `corotational` flag and no `transform` field.

### Fiber sections

`FiberTable` (`y`, `area`, `material`) and `FiberTable3` (`y`, `z`, `area`,
`material`) hold fibers for all sections, flattened and offset-indexed by
`sectionOffsets` (see conventions). `material` is an arena index parallel to the
coordinate arrays. Fiber sections carry no torsion; 3D elements supply it
through their own `g`/`j`.

## Constraints

- `EqualDofTable` / `EqualDofTable3`: row `i` ties `constrained[i]` to
  `retained[i]` for each `[row, dof]` in the sparse `dofs`.
- `RigidDiaphragmTable`: `retained[]`, `constrained: [row, node][]`. Ties each
  listed node's `ux` to the retained node.
- `RigidDiaphragmTable3`: adds `normal: ("x" | "y" | "z")[]`, the diaphragm
  normal per row. Ties the two in-plane translations including the rotation
  lever-arm term.

## Loads

- `LoadPatternTable`: `series: TimeSeriesSpec[]`, `scaleFactor: number[]`.
  `TimeSeriesSpec` is `{ kind: "constant" }`, `{ kind: "linear", slope }`, or
  `{ kind: "path", times, factors }`.
- `NodalLoadTable`: `pattern`, `node`, `dof`, `value`, `stage`. `stage` is the
  index of the stage at whose start the load is registered, which matters when
  an earlier stage must not ramp this pattern.
- `ElementLoadTable` / `ElementLoadTable3`: `pattern`, `elementKind`,
  `elementIndex` (row in the table named by `elementKind`), `load`, `stage`.
  - `elementKind`: `"truss" | "elasticBeamColumn" | "dispBeamColumn" | "forceBeamColumn" | "zeroLength" | "zeroLengthSection"`.
  - `load`: `{ kind: "uniform", wx, wy }` (2D) or
    `{ kind: "uniform", wx, wy, wz }` (3D). Uniform force per length in the
    element's **local** axes: `wx` along the member axis, `wy`/`wz` transverse
    (3D local `y`/`z` come from the element's `vecXz`). All fields are
    required. The beam-column kinds (`elasticBeamColumn`, `dispBeamColumn`, `forceBeamColumn`) accept element loads today; any other
    `elementKind` fails decode with `unsupportedElementLoad`. Several rows
    for the same element and pattern sum.
  - `elementForce` recorders report member end forces including the fixed-end
    effect of element loads (the free end of a loaded cantilever reports zero).

## Analysis sequence

`SequenceSpec` / `SequenceSpec3`: `{ stages: StageSpec[], recorders: RecorderSpec[] }`.
Stages run in order.

### `StageSpec`

| `kind` | Fields |
|---|---|
| `static` | `id`, `steps`, `integrator`, `algorithm`, `convergence?`, `holdPatternsAfter: number[]` |
| `modal` | `id`, `modes` |
| `transient` | `id`, `steps`, `dt`, `damping: { alphaM, betaK }`, `groundMotions: GroundMotionSpec[]`, `algorithm?`, `convergence?` |

`GroundMotionSpec`: `{ direction, series: TimeSeriesSpec, scaleFactor }`.
`direction` is `0` for x, `1` for y, and `2` for z (3D only).

### Integrators

- `{ kind: "loadControl", increment }`
- `{ kind: "displacementControl", node, dof, increment }`
- `{ kind: "arcLength", ...ArcLengthSpec }`. See [arclength.md](arclength.md)
  for behavior. Fields: `initialRadius`, `minRadius?`, `maxRadius?` (both default
  to `initialRadius`), `targetIterations?`, `maxRetries?`,
  `direction?: "increasing" | "decreasing"`, `scales`, `predictor?: "secant" | "tangent"`,
  `seed?`, `arcTolerance?`, `correctionTolerance?`,
  `backtracking?: { armijo, minStep }`, and `stop?`.
  - `scales`: `{ kind: "explicit", displacement, rotation?, load }` or
    `{ kind: "auto", load }`.
  - `seed`: `{ components: { node, dof, value }[], load? }`.
  - `stop`: `{ displacement?: { node, dof, value, exact? }, loadFactor?: { value, exact? }, loadFactorZeroCrossing?, maxChordLength? }`.

### Algorithms and convergence

`algorithm` accepts the legacy bare strings `"linear"` / `"newtonRaphson"` or:

- `{ kind: "linear" }`
- `{ kind: "newton", tangent?, lineSearch? }`
- `{ kind: "krylovNewton", tangent, maxDimension }`

`tangent`: `"current" | "reuseAtStepStart" | "initial"`.
`lineSearch`: `{ kind: "bisection" | "regulaFalsi", tol, maxIter, maxEta }`.

`convergence`:

- `{ kind: "normUnbalance" | "normDispIncr" | "energyIncr", tol, maxIter }`
- `{ kind: "combined", forceTol, momentTol?, relativeTol?, displacementTol?, maxIter }`

### Recorders

Each recorder is one scalar channel, tagged by `response`. The 3D form
(`RecorderSpec3`) is identical except that `elementKind` names a 3D kind
and `component` indexes a width-12 local force vector (width-6 for 2D).

| `response` | Fields |
|---|---|
| `nodeDisp`, `nodeVel`, `nodeAccel`, `reaction` | `node`, `dof` |
| `elementForce` | `elementKind`, `elementIndex`, `component` |
| `modeShape` | `mode`, `node`, `dof` |
| `fiber` | `elementKind`, `elementIndex`, `point`, `fiber`, `quantity: "strain" \| "stress"` |

## Results: `StepOutcome`

Returned by `advance`:

| Field | Meaning |
|---|---|
| `done` | All stages finished (or stopped by `error`). |
| `stageComplete` | The current stage finished during this call. |
| `stepsTaken` | Steps taken in this call. |
| `loadFactor` | Static: integrator load factor. Transient: elapsed time. Modal: `0`. |
| `error` | `AnalysisErrorDetail \| undefined`; sticky once set. |
| `continuation` | Arc-length diagnostics for the last accepted step, or `undefined`. |
| `recorderBatches` | `RecorderBatch[]`: only the samples produced by this call. |

`RecorderBatch`: `{ recorderIndex, stageIndex, firstSample, samples: [pseudoTime, value][] }`.
`firstSample` plus `samples.length` gives the sample range, so a failed persist
can be retried without renumbering. With arc length the load factor can
decrease and repeat, so order samples by sample index, never by load factor.
The caller accumulates batches; there is no "history so far" accessor.

`ContinuationDetail`: `radius`, `nextRadius`, `retries`, `correctorSolves`,
`factorizations`, `forceMeasure`, `arcError`, `chordLength`, `detSign`,
`bifurcationSuspected`, `displacementScale`, `rotationScale`, `stop`
(`{ reason, landedExactly, overshoot }`, where `reason` is
`"displacementTarget" | "loadFactorTarget" | "loadFactorZeroCrossing" | "chordLength" | "stepCount"`).

### `AnalysisErrorDetail` (`kind`)

`failedToConverge {step}`, `singularSystem`, `invalidConstraint`,
`invalidModeCount {requested, freeDofs}`, `invalidOption {field}`,
`unsupportedLoadSeries`, `zeroLoadSensitivity`,
`initialStateNotInEquilibrium {measure}`, `missingSeedDirection`,
`cutbacksExhausted {step, attempts, radius, lastFailure}`, `continuationComplete`.

`lastFailure` is one of `"notConverged" | "stagnated" | "noDescent" | "nonFinite" | "singularSystem" | "inaccurateSolve" | "branchOrientation" | "degenerateConstraint"`.

## Errors: `DecodeError`

`decodeInput` throws a `{ kind, ... }` object (not an `Error` instance):

| `kind` | Extra fields |
|---|---|
| `invalidAnalysisOption` | `stage`, `field` |
| `unsupportedSpace` | `got` |
| `unknownNodeIndex`, `unknownMaterialIndex`, `unknownPatternIndex`, `unknownElementIndex`, `unknownStageIndex` | `table`, `row` |
| `cyclicMaterialReference` | `index` |
| `unknownFiberSectionIndex` | `row` |
| `unsupportedElementLoad` | `elementKind` |
| `invalidDof` | `table`, `row`, `dof` |
| `unknownConstraintRow` | `table`, `row` |

## Extending the format

Add or change the Rust types and decoder, then rebuild to regenerate the
TypeScript definitions, and update this page. A new binding is only needed for
a new operation exposed to JavaScript. A new material or solver behind the
existing session interface does not need one.
