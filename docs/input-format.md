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
  `NodeTable.coords` (stride `ndm`).
- **Sparse tables** are lists of tuples or small objects keyed by a row index,
  such as `[row, dof, materialIndex]`, because most rows have no entry.
- **Offset-indexed tables** (fibers, linear constraints) flatten variable-length
  groups: group `k` occupies `offsets[k]..offsets[k+1]`, so the offsets array has
  `numGroups + 1` entries.
- **2D and 3D.** One format serves both. `header.ndm` selects 2D (`2`: 2
  coordinates, 3 DOFs per node `ux, uy, rz`) or 3D (`3`: 3 coordinates, 6 DOFs
  `ux, uy, uz, rx, ry, rz`); anything else is `unsupportedNdm`. Entities whose
  data is the same in both (nodes, trusses, zero-lengths, fibers, constraints,
  loads, the sequence) have one table. Formulations that differ (beam-columns)
  have one table per formulation, suffixed `2d`/`3d`. The table of the other
  profile must be empty or omitted: a non-empty one is `tableNotInProfile`, and
  an element kind of the other profile in a load or recorder is
  `elementKindNotInProfile`. The decoder picks the profile once from `ndm`.
- **Every table is optional.** An omitted table is an empty one, so a model
  lists only what it uses.
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

Only `header` is required.

| Field | Type | Profile |
|---|---|---|
| `header` | `Header` | both |
| `nodes` | `NodeTable` | both |
| `materials` | `MaterialSpec[]` | both (shared arena) |
| `fibers` | `FiberTable` | both |
| `trusses` | `TrussTable` | both |
| `elasticBeamColumns2d` / `elasticBeamColumns3d` | `ElasticBeamColumn2dTable` / `ElasticBeamColumn3dTable` | 2D / 3D |
| `dispBeamColumns2d` / `dispBeamColumns3d` | `FiberBeamColumn2dTable` / `FiberBeamColumn3dTable` | 2D / 3D |
| `forceBeamColumns2d` / `forceBeamColumns3d` | `FiberBeamColumn2dTable` / `FiberBeamColumn3dTable` | 2D / 3D |
| `zeroLengths` | `ZeroLengthTable` | both |
| `zeroLengthSections` | `ZeroLengthSectionTable` | both |
| `equalDofs` | `EqualDofTable` | both |
| `rigidDiaphragms` | `RigidDiaphragmTable` | both |
| `rigidLinks` | `RigidLinkTable` | both |
| `linearConstraints` | `LinearConstraintTable` | both |
| `loadPatterns` | `LoadPatternTable` | both |
| `nodalLoads` | `NodalLoadTable` | both |
| `elementLoads` | `ElementLoadTable` | both |
| `sequence` | `SequenceSpec` | both |

### `Header`

| Field | Type | Notes |
|---|---|---|
| `schemaVersion` | `number` | Wire format version, independent of the engine. |
| `ndm` | `number` | `2` or `3`; anything else is `unsupportedNdm`. |
| `engineVersion` | `string` | `carapace-core` version, recorded for provenance. |
| `recordInitial` | `boolean?` | See the analysis sequence. |

## Nodes

`NodeTable`:

| Field | Meaning |
|---|---|
| `coords` | stride `ndm`: `(x, y)` in 2D, `(x, y, z)` in 3D |
| `fixed` | one bitmask per node, bit `k` = DOF `k`: `ux, uy, rz` in 2D; `ux, uy, uz, rx, ry, rz` in 3D |
| `massNodeIndex` | node indices that carry mass (sparse) |
| `mass` | stride `ndf` (3 in 2D, 6 in 3D): one value per DOF, parallel to `massNodeIndex` |

A DOF is an equation only if an element stiffens it, a constraint uses it as a master, it carries
a nonzero mass, or `fixed` pins it. DOFs nothing uses (a truss's rotations) need not be fixed.
A nodal load on such a DOF is rejected: `advance` returns `error.kind = "invalidModel"` with
`error.error = { kind: "loadOnInactiveDof", node, dof }`, and `decodeInput` rejects element
geometry problems the same way (`DecodeError` `invalidModel`).

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
| `TrussTable` (2D and 3D) | `area`, `material` (arena index) |
| `ElasticBeamColumn2dTable` | `e`, `a`, `iz`, `transform: TransformSpec[]` |
| `ElasticBeamColumn3dTable` | `e`, `g`, `a`, `j`, `iy`, `iz`, `transform: TransformSpec3[]` |
| `FiberBeamColumn2dTable` (disp and force) | `fiberSection` (index into `FiberTable.sectionOffsets`), `integration: IntegrationSpec[]`, `corotational: boolean[]` |
| `FiberBeamColumn3dTable` (disp and force) | `g`, `j` (decoupled elastic torsion), `vecXz: [x,y,z][]`, `fiberSection`, `integration` |
| `ZeroLengthTable` (2D and 3D) | `materials: [row, dof, material][]` (sparse), `friction?: FrictionRow[]`, `orient?: OrientRow[]` |
| `ZeroLengthSectionTable` (2D and 3D) | `fiberSection`, `materials: [row, dof, material][]` for DOFs the section does not drive (`uy` in 2D; `uy`, `uz`, `rx` in 3D), `orient?: OrientRow[]` |

Notes:

- `FrictionRow`: `{ row, normalDof, shearDofs, mu, k0, b }`, at most one per row.
  `shearDofs` has one entry in 2D and two in 3D (both coupled to the same normal
  force); any other length is `invalidRow`.
- `OrientRow`: `{ row, x: [x1, x2, x3], yp?: [yp1, yp2, yp3] }`, at most one per
  row, following OpenSees' `-orient`. In 2D, `yp` must be absent, `x3` must be 0,
  local x is `(x1, x2)` and local y is local x turned 90° counter-clockwise, so
  `rz` is unchanged. In 3D, `yp` is required: local z is `x × yp`, local y is
  `z × x`. Every per-DOF material, friction coupling and the section's
  axial/flexural directions are evaluated on the relative displacement and
  rotation projected onto those axes (translations and rotations share the
  frame). Rows without an entry use the global axes, and the field may be
  omitted. A zero vector (or, in 3D, parallel vectors; in 2D, `x3 != 0`) decodes
  to `invalidOrientation`; a `yp` in 2D, or none in 3D, to `invalidRow`.
  Element forces recorded for an oriented element are still in global DOFs.
- `TransformSpec` (2D): `"linear" | "pDelta" | "corotational"`.
  `TransformSpec3`: `{ kind: "linear3" | "pDelta3", vecXz: [x,y,z] }`. There is
  no 3D corotational transform. `vecXz` is a vector not parallel to the
  member axis; it fixes the local y/z orientation.
- `IntegrationSpec`: `{ kind: "legendre" | "lobatto", points }`.
- 3D fiber elements have no `corotational` flag and no `transform` field.

### Fiber sections

`FiberTable` (`y`, `z?`, `area`, `material`) holds fibers for all sections,
flattened and offset-indexed by `sectionOffsets` (see conventions). `material`
is an arena index parallel to the coordinate arrays. `z` is omitted (or empty)
for a 2D model and parallel to `y` in a 3D model (biaxial bending); a mismatch
is `invalidRow`. Fiber sections carry no torsion; 3D elements supply it through
their own `g`/`j`.

## Constraints

All constraints are linear multi-point constraints resolved by the transformation
handler (the session chooses it when the model has any).

- `EqualDofTable`: row `i` ties `constrained[i]` to `retained[i]` for each
  `[row, dof]` in the sparse `dofs`.
- `RigidDiaphragmTable`: `retained[]`, `constrained: [row, node][]`, and in 3D
  only `normal?: ("x" | "y" | "z")[]` (one per row; omitted means `"y"` for
  every row; a 2D model must omit it). In 2D it ties each listed node's `ux` to
  the retained node. In 3D it ties the two in-plane translations (perpendicular
  to the normal) including the rotation lever-arm term.
- `RigidLinkTable`: `master[]`, `slave[]` node pairs. The slave moves with the
  master as a rigid body: its rotations equal the master's and its translations
  equal the master's plus the small-rotation lever arm, taken from the node
  coordinates (`core::Domain::rigid_link`).
- `LinearConstraintTable`: general `u[slaveNode, slaveDof] = Σ coeff · u[node,
  dof]`, as flattened sparse rows: `slaveNode[]`, `slaveDof[]`, `termOffsets[]`
  (`rows + 1` entries, ending at the term count; may be omitted for an empty
  table), and the parallel term arrays `termNode[]`, `termDof[]`, `termCoeff[]`.
  Row `i` uses terms `termOffsets[i]..termOffsets[i + 1]`. Structural problems
  (a slave defined twice, a fixed slave, a dependency cycle, a prescribed
  master) are reported as `invalidModel`.

## Loads

- `LoadPatternTable`: `series: TimeSeriesSpec[]`, `scaleFactor: number[]`.
  `TimeSeriesSpec` is `{ kind: "constant" }`, `{ kind: "linear", slope }`, or
  `{ kind: "path", times, factors }`.
- `NodalLoadTable`: `pattern`, `node`, `dof`, `value`, `stage`. `stage` is the
  index of the stage at whose start the load is registered, which matters when
  an earlier stage must not ramp this pattern.
- `ElementLoadTable`: `pattern`, `elementKind`, `elementIndex` (row in the table
  named by `elementKind`), `load`, `stage`.
  - `elementKind` (shared with recorders): `"truss" | "elasticBeamColumn2d" |
    "elasticBeamColumn3d" | "dispBeamColumn2d" | "dispBeamColumn3d" |
    "forceBeamColumn2d" | "forceBeamColumn3d" | "zeroLength" | "zeroLengthSection"`.
    A kind of the other profile is `elementKindNotInProfile`.
  - `load`: `{ kind: "uniform", wx, wy, wz? }`. Uniform force per length in the
    element's **local** axes: `wx` along the member axis, `wy`/`wz` transverse
    (3D local `y`/`z` come from the element's `vecXz`). `wz` is 3D only (omitted
    means 0; a nonzero value in 2D is `invalidRow`). The beam-column kinds accept
    element loads today; any other `elementKind` fails decode with
    `unsupportedElementLoad`. Several rows for the same element and pattern sum.
  - `elementForce` recorders report member end forces including the fixed-end
    effect of element loads (the free end of a loaded cantilever reports zero).

## Analysis sequence

`SequenceSpec`: `{ stages: StageSpec[], recorders: RecorderSpec[] }`.
Stages run in order.

### `StageSpec`

| `kind` | Fields |
|---|---|
| `static` | `id`, `steps`, `integrator`, `algorithm`, `convergence?`, `holdPatternsAfter: number[]` |
| `modal` | `id`, `modes` |
| `reset` | `id` |
| `transient` | `id`, `steps`, `dt`, `damping: { alphaM, betaK }`, `groundMotions: GroundMotionSpec[]`, `algorithm?`, `convergence?` |

`GroundMotionSpec`: `{ direction, series: TimeSeriesSpec, scaleFactor }`.
`direction` is `0` for x, `1` for y, and `2` for z (3D only).

`reset` reverts the model to its as-built state (OpenSees `reset`): displacements,
velocities and all element/material history return to their starting values and time to
`0`. Like OpenSees it does not undo `holdPatternsAfter` freezes: a held pattern stays
applied at its frozen factor, so the next stage sees that load in full from its first
step (patterns that were not held ramp again from `0`). It takes no steps and records
nothing.

`header.recordInitial` (optional, default `false`): record one sample of every supported
recorder at the start of each `static`/`transient` stage, before its first step, at
load factor/time `0`: the stage's initial conditions, i.e. step 0 of that analysis.

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

Each recorder is one scalar channel, tagged by `response`. `elementKind` is the
shared kind enum above. `dof` is bounded by the profile's DOFs per node
(`invalidDof`). `component` of an `elementForce` recorder indexes the element's
local force vector (width 6 in 2D, 12 in 3D) and of an `elementLoad` recorder
the load (`wx, wy` in 2D; `wx, wy, wz` in 3D); one past the width is
`invalidRecorderComponent`.

| `response` | Fields |
|---|---|
| `nodeDisp`, `nodeVel`, `nodeAccel`, `reaction` | `node`, `dof` |
| `elementForce`, `elementLoad` | `elementKind`, `elementIndex`, `component` |
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

### Modal results: `modalResults()`

A `modal` stage records nothing in the step timeline. After it finishes,
`session.modalResults()` returns `{ stages: ModalStageResult[] }` (every modal
stage finished so far, in order; empty until one completes):

`ModalStageResult`: `{ stageIndex, stageId, ndf, modes: ModeResult[], totalMass: number[] }`.
`totalMass[d]` is the total mass `rᵀMr` in global translation direction `d`.

`ModeResult`: `{ frequency, shape, participation, massRatio }`, ascending
frequency.
- `frequency`: circular frequency ω in rad/time (period `T = 2π/ω`).
- `shape`: node-major over every node in input node-table order,
  `shape[node * ndf + dof]`; fixed DOFs are `0`. Mass-normalized (`φᵀMφ = 1`)
  and sign-fixed so the largest-magnitude entry is positive.
- `participation[d]`: `φᵀ M r_d`. `massRatio[d]`: `participation[d]² / totalMass[d]`.
  Models with rigid diaphragms get approximate ratios (a slave node's
  translation is not exactly 1 on its retained equation).

Massless DOFs (e.g. rotations with no rotational inertia) are statically
condensed; at least `modes` DOFs must carry mass.

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
| `unsupportedNdm` | `got` |
| `unknownNodeIndex`, `unknownMaterialIndex`, `unknownPatternIndex`, `unknownElementIndex`, `unknownStageIndex` | `table`, `row` |
| `cyclicMaterialReference` | `index` |
| `unknownFiberSectionIndex` | `row` |
| `unsupportedElementLoad` | `elementKind` |
| `invalidDof` | `table`, `row`, `dof` |
| `unknownConstraintRow` | `table`, `row` |
| `tableNotInProfile` | `table` (a 2D/3D formulation table of the other profile is not empty) |
| `elementKindNotInProfile` | `table` (an element kind of the other profile in a load or recorder) |
| `invalidRecorderComponent` | `recorder` (index), `component`, `width` |
| `invalidRow` | `table`, `row`, `reason` (a row shape that is wrong for the profile or its sparse encoding) |
| `invalidOrientation` | `table`, `row` |
| `invalidModel` | `error` (a `core` model error, e.g. `{ kind: "loadOnInactiveDof", node, dof }`) |

## Extending the format

Add or change the Rust types and decoder, then rebuild to regenerate the
TypeScript definitions, and update this page. A new binding is only needed for
a new operation exposed to JavaScript. A new material or solver behind the
existing session interface does not need one.

Rules for new tables:

- Every core capability that a model can use gets a wire table in the same
  change that adds it; nothing is reachable from core only.
- One table per element *formulation*, not per profile. If the 2D and 3D
  versions have the same fields (truss, zero-length), they share a table and the
  decoder branches on `ndm` only where the element constructors differ. If the
  fields differ, add `…2d`/`…3d` tables, a matching `ElementKind` variant for
  each, and a `tableNotInProfile` check in the other profile's decoder.
- Add the table to `CarapaceInputV1` with `#[serde(default)]` and
  `#[tsify(optional)]`, add its decoding to `decode/shared.rs` if it is
  dimension-agnostic (or `decode/elements_2d.rs` / `elements_3d.rs`), and a
  `DecodeError` for any new way a row can be wrong.
