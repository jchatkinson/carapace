# Model input format (`CarapaceInputV1`)

Reference for the JavaScript/TypeScript model API: the payload passed to
`decodeInput`, the session it returns, and the result and error shapes.

The Rust types in `wasm-bridge/src/input_v1/` are the source of truth. They
derive `serde` + `tsify::Tsify`, and this page summarizes the TypeScript
definitions generated from them. If this page and the types disagree, the types
win. To see the exact declarations:

```bash
scripts/build-wasm.sh
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
- **Plain JS values in, typed arrays for bulk out.** The boundary accepts ordinary objects and arrays
  as input (a typed-array input format is not implemented). The large per-step output, the `bulk`
  recorders' `progress` and `values`, comes back as transferable `Float64Array`s.

## Entry points

```ts
decodeInput(value: CarapaceInputV1): WasmSession   // throws DecodeError
session.advance(stepBudget: number): StepOutcome
session.currentStageId(): string | undefined       // undefined once all stages done
session.modelOutput(): ModelOutput                 // resolved embedded-bar segments, routing, sr
session.free()
```

`advance` runs up to `stepBudget` steps and returns, so a caller can run it in a
loop and keep a worker responsive. Once `error` is set it is sticky: later
stages are not attempted.


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
| `planeMaterials` | `PlaneMaterialSpec[]` | 2D (arena for `triangles`/`quads`) |
| `triangles` | `TriangleTable` | 2D |
| `quads` | `QuadTable` | 2D |
| `embeddedBars` | `EmbeddedBarSpec[]` | 2D (rebar polylines embedded in `triangles`/`quads`) |
| `shellSections` | `ShellSectionSpec[]` | 3D (arena for `shell4s`) |
| `shell3s` | `Shell3Table` | 3D |
| `shell4s` | `Shell4Table` | 3D |
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
| `csfmRebar` | `es`, `fy`, `ft`, `epsU`, `fctm`, `ecm`, `diameter`, `rhoEff`, `lambda`, `tensionStiffening`. The CSFM embedded-rebar law; see [CSFM](#csfm-concrete-and-embedded-rebar). |
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
| `TriangleTable` (2D) | `nodeIds` (stride 3, counter-clockwise), `thickness`, `material` (index into `planeMaterials`), `density` |
| `QuadTable` (2D) | `nodeIds` (stride 4, counter-clockwise), `thickness`, `material`, `density`, `formulation: ("full" \| "enhanced")[]` |
| `Shell3Table` (3D) | `nodeIds` (stride 3; the node order sets the local normal by the right-hand rule), `section` (index into `shellSections`) |
| `Shell4Table` (3D) | `nodeIds` (stride 4; the node order sets the local normal by the right-hand rule), `section` (index into `shellSections`) |
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
- **Continuum elements** (`triangles` = `Tri3`, constant strain, one Gauss point; `quads` = `Quad4`,
  bilinear, 2x2 Gauss points) carry only translations (`ux`, `uy`), so their nodes need no rotational
  fixity. Nodes must be counter-clockwise; geometry is validated (coincident nodes, a clockwise
  or degenerate triangle, a quad with a non-positive Jacobian at any corner) and reported as
  `invalidModel`. `density` gives an exact lumped mass per translational DOF (a quad's masses are the
  row sums of the consistent mass). `formulation: "enhanced"` selects Wilson-Taylor incompatible
  modes, which relieve bending locking and need linear materials; it is not a mixed formulation and
  makes no volumetric-locking guarantee. Mass-bearing continuum elements cannot touch a rigid-link or
  diaphragm slave node (`invalidModel`: `massOnConstrainedDof`).
- `PlaneMaterialSpec` (the `planeMaterials` arena; strain is `[eps_x, eps_y, gamma_xy]` with
  engineering shear): `{ kind: "isotropic", e, nu, state: "planeStress" | "planeStrain" }`,
  `{ kind: "orthotropic", ex, ey, nuXy, gXy, angle }` (plane stress; `angle` in radians
  counter-clockwise from global x to material axis 1), `{ kind: "elasticMatrix", d: [d11, d12, d13,
  d22, d23, d33] }` (packed symmetric `D`). Invalid constants are `invalidPlaneMaterial`.
- **Shells** (`shell4s` = `Shell4`, the MITC4 flat shell, OpenSees `ShellMITC4`; `shell3s` = `Shell3`, the DKT/Allman flat triangle, OpenSees `ShellDKGT`; 3D only) carry all
  six DOFs per node and use 2x2 Gauss points. The local frame follows OpenSees: `e1` from the edge
  midpoints, `e3 = e1 x e2` (the node order sets the sign), and a quad is treated as flat in that
  frame (warp is tolerated as OpenSees does). `ShellSectionSpec` (the `shellSections` arena):
  `{ kind: "elasticMembranePlate", e, nu, h, rho }` (invalid constants are `invalidShellSection`).
  Mass is a row-sum lumped translational mass `rho h` (OpenSees uses a consistent mass), so modal
  frequencies agree with OpenSees only approximately. A `shell3s` triangle has no drilling penalty and no transverse shear (its `gamma_xz`, `gamma_yz` and `Qx`, `Qy` are zero) and its local x axis is the edge from node 0 to node 1. A triangle has 4 Gauss points `(1/3,1/3,1/3)` with weight `-9/16` and `(1/5,3/5,1/5)` permutations with `25/48`; a Gauss-point recorder on a shell has 4 points
  and 8 components per quantity, in OpenSees' section order and sign: strain `[eps_x, eps_y,
  gamma_xy, kappa_x, kappa_y, 2 kappa_xy, gamma_xz, gamma_yz]`, `stress` the resultants `[Nx, Ny, Nxy,
  Mx, My, Mxy, Qx, Qy]` in the local axes (OpenSees' bending resultants have the opposite sign to
  `M = D kappa`). An `elementForce` recorder reports the 24 global nodal forces.
- `PlaneMaterialSpec` also has `{ kind: "csfmConcrete", fck, ecm?, law, softening, tensileRegularization? }`
  (see [CSFM](#csfm-concrete-and-embedded-rebar)).
- 3D fiber elements have no `corotational` flag and no `transform` field.

### CSFM concrete and embedded rebar

The Compatible Stress Field Method pieces: a rotating principal-strain concrete, the CSFM rebar
law, and rebar polylines embedded in a 2D continuum mesh. Units are N, mm, MPa throughout.

**`csfmConcrete`** (a `planeMaterials` entry; 2D only):

| Field | Meaning |
|---|---|
| `fck` | Characteristic strength, `0 < fck <= 50` MPa. |
| `ecm?` | Replaces the derived `Ecm = 22000 ((fck + 8) / 10)^0.3`. Finite and positive. |
| `law` | `"parabolaRectangle"` (parabola to `eps_c2 = -0.002`, then a plateau) or `"elasticIdeallyPlastic"`. |
| `softening` | Kaufmann compression softening `kc2 = min(1, 1 / (0.8 + 140 eps_1))` on or off. |
| `tensileRegularization?` | `gamma`, `0 < gamma < 1`. Omitted gives the exact zero-tension law. |

It is a stateless rotating principal-strain law with `nu = 0`: stress depends on the current
principal strains only, and unloading retraces loading. The tangent is the analytic consistent
tangent, which is nonsymmetric when softening is on (the solver uses LU for it). Concrete has no
tensile strength. `tensileRegularization` is a **numerical regularization**, not a material
property: it adds a small tension stiffness `gamma * fctm (1 - exp(-eps / eps_reg))` so a cracked
element does not make the stiffness singular, and results (especially bar stress between cracks)
depend on it slightly. Choose it deliberately and report it with the results. `fck` outside the
range, a bad `ecm` or a bad `gamma` is `invalidPlaneMaterial`. The enhanced quad accepts linear
materials only, so `csfmConcrete` with `formulation: "enhanced"` is `invalidModel`
(`invalidElement`).

**`csfmRebar`** (a `materials` entry):

| Field | Meaning |
|---|---|
| `es`, `fy`, `ft`, `epsU` | Bare-steel bilinear properties; `ft > fy` and `epsU > fy / es`. |
| `fctm`, `ecm` | The host concrete's tensile strength and modulus, given explicitly (not referenced). |
| `diameter` | Bar diameter (mm). |
| `rhoEff` | Effective reinforcement ratio, `0 < rhoEff < 1`. |
| `lambda` | Crack-spacing calibration factor. |
| `tensionStiffening` | `true`: tension follows the TCM/POM average-strain relation. `false`: symmetric bare bilinear steel. Compression is always bare steel. |

With tension stiffening, `rhoEff` at or above `rho_crit = fctm / (fy - (es / ecm - 1) fctm)` routes to
the Tension Chord Model (TCM, crack spacing `sr`), below it to the Pull-Out Model (POM). The law is
stateless and nonlinear-elastic: unloading retraces the loading curve. Invalid steel, ratio or
bond inputs, and a TCM startup bridge with `2 sigma_min > fy`, are `invalidMaterial`
(`index`, `reason`). A bar is rupture-flagged strictly past `eps_m(ft)` (tension stiffening) or
`|eps| > epsU` (bare steel or compression).

**`embeddedBars`**: `{ polyline: [x0, y0, x1, y1, ...], area, material, hosts? }[]`, one row per
bar; `material` indexes `materials` (a `csfmRebar` for CSFM, but any uniaxial material works) and
`area` is the bar cross-section (mm^2). The polyline is split at every mesh-edge crossing so each
piece lies in one `triangles`/`quads` host, and each piece becomes an `embeddedBar` element whose
strain is the host displacement field projected on the piece (each piece is bound to the host
containing its midpoint). `hosts` is `"all"` (default) or `{ kind: "tri3" | "quad4", index }[]`,
which restricts the elements that may host the bar. Errors: fewer than two vertices, a bad area
or a non-finite coordinate is `invalidRow` (`embedded_bars`); an unknown material is
`unknownMaterialIndex`; an unknown host element is `unknownElementIndex`; a piece whose midpoint
is in no allowed host (outside the mesh, or across a hole) is `invalidModel` with
`embeddedBarOutsideMesh`, one in several overlapping hosts is `embeddedBarAmbiguousHost`, and a
path with no length is `embeddedBarDegeneratePath`. The table in a 3D model is `tableNotInProfile`.

**`modelOutput()`** is available right after `decodeInput`:
`{ embeddedBars: { firstSegment, segments, routing?, sr? }[] }`, one group per `embeddedBars` row.
`segments[s]` is `{ start, end, length, hostKind, hostNodes }` (undeformed end points in drawing
order; the host's nodes as `nodes` rows). Results are indexed `(bar, segment)`; the flat
`embeddedBar` element index of segment `s` of bar `b` is `groups[b].firstSegment + s`. For a
`csfmRebar` bar, `routing` (`"tcm"` or `"pom"`) and `sr` (crack spacing, mm; `0` for POM) report
what the material resolved to.

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
  `{ kind: "path", times, factors }`, or `{ kind: "boundedPath", times, factors, useLast }`.
  `path` retains its endpoint-holding behavior. `boundedPath` is zero before its first sample
  and after its last sample (unless `useLast` is true), for OpenSees-compatible accelerograms.
  Both require nonempty, equally sized finite arrays and strictly increasing times.
- `NodalLoadTable`: `pattern`, `node`, `dof`, `value`, `stage`. `stage` is the
  index of the stage at whose start the load is registered, which matters when
  an earlier stage must not ramp this pattern.
- `ElementLoadTable`: `pattern`, `elementKind`, `elementIndex` (row in the table
  named by `elementKind`), `load`, `stage`.
  - `elementKind` (shared with recorders): `"truss" | "elasticBeamColumn2d" |
    "elasticBeamColumn3d" | "dispBeamColumn2d" | "dispBeamColumn3d" |
    "forceBeamColumn2d" | "forceBeamColumn3d" | "zeroLength" | "zeroLengthSection" | "tri3" | "quad4"`.
    A kind of the other profile is `elementKindNotInProfile`.
  - `load`: `{ kind: "uniform", wx, wy, wz? }`. Uniform force per length in the
    element's **local** axes: `wx` along the member axis, `wy`/`wz` transverse
    (3D local `y`/`z` come from the element's `vecXz`). `wz` is 3D only (omitted
    means 0; a nonzero value in 2D is `invalidRow`). The beam-column kinds accept
    element loads today; any other `elementKind` fails decode with
    `unsupportedElementLoad`. Several rows for the same element and pattern sum.
  - Continuum loads (2D `tri3`/`quad4` only; a beam, truss or spring takes none of them, and a
    continuum element takes no `uniform`, both `unsupportedElementLoad`):
    `{ kind: "body", bx, by }` (force per unit volume, global axes),
    `{ kind: "edgeTraction", edge, tx, ty }` (force per unit edge length, global axes) and
    `{ kind: "edgePressure", edge, pressure }` (positive into the element along the inward
    normal). Edge `k` joins local node `k` to node `k + 1` (a triangle has edges 0..3, a quad
    0..4; past that is `invalidRow`). Body forces are integrated with the same weights as the
    lumped mass, edges with a 2-point Gauss rule. Rows add up per element and pattern.
  - Shell loads (`shell3` and `shell4` only; a shell takes no `uniform`): `{ kind: "shellPressure", pressure }`
    (per unit area, positive along the local normal `e3`) and `{ kind: "shellBody", bx, by, bz }`
    (a body acceleration in global axes, force per area `rho h b`; a gravity vector points down).
    Both are consistent nodal forces `integral(N_i) * value`, translations only. OpenSees'
    `ShellMITC4` has no pressure load and applies `-selfWeight` with the opposite sign, and `ShellDKGT` applies twice
    the self-weight (a bug in its integration weight), so exporters expand pressure, and a triangle's self-weight, into
    nodal loads and negate the quad's body vector.
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
| `transient` | `id`, `steps`, `dt`, `damping: { alphaM, betaK, modalAnchors? }`, `groundMotions: GroundMotionSpec[]`, `algorithm?`, `convergence?` |

Optional `damping.modalAnchors: { mode1, mode2, ratio }` uses positive, 1-based
mode numbers and a damping ratio in `[0, 1]`. At transient-stage entry, the engine
solves the requested modes from the current domain (including preceding gravity)
and overrides `alphaM` and `betaK` to match that ratio at both modal frequencies.
No preceding modal stage is required. Unavailable modes fail the run; invalid
mode numbers or ratios fail decoding. Omitting `modalAnchors` preserves direct
coefficient behavior.

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
- `{ kind: "loadControlCutbacks", increment, minIncrement }`: load control with adaptive
  substepping, `0 < minIncrement <= increment`. A step that fails to converge is reverted and
  retried at half its size, down to `minIncrement`, after which the stage fails with
  `loadControlExhausted`. The stage's `steps` counts **nominal** increments: it ends when
  `steps * increment` of load factor has been accepted, however many substeps that took, and the
  uncovered remainder of a cut-back step is requested again by the next substep. Every accepted
  substep counts against an `advance` budget and records one sample of every recorder. This is a
  Carapace extension, not OpenSees' `LoadControl dLambda numIter minLambda maxLambda`, which
  rescales the next step by iteration count without ever reverting a failed one; a plain
  `loadControl` stays exactly OpenSees' fixed-increment semantics.
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
  - `stop`: `{ displacement?: { node, dof, value, exact? }, loadFactor?: { value, exact? }, loadFactorZeroCrossing?, maxChordLength?, materialLimit? }`.
  - `stop.materialLimit`: `{ exact?, concreteCrushingLimit?, concreteTensileLimit?, rebarStrainLimit? }`
    stops when the controlling CSFM material-limit utilization reaches `1` (crushing at
    `|eps_cu2| = 0.0035` and rebar rupture are always checked; the optional positive strain limits
    replace the crushing limit or add a concrete tensile-strain / user rebar-strain limit). With
    `exact` the bracketing step is retaken so the stage ends on the limit
    (`continuation.stop.reason = "materialLimit"`). A non-positive or non-finite limit is
    `invalidAnalysisOption` with field `integrator.stop.materialLimit`.

### Algorithms and convergence

`algorithm` is a tagged object:

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
the load (3D: `wx, wy, wz`; 2D: 16 components, `[wx, wy, bx, by, t0x, t0y, p0, ..., t3x, t3y, p3]` for beam uniform, body force, then traction and pressure of each edge); one past the width is
`invalidRecorderComponent`.

| `response` | Fields |
|---|---|
| `nodeDisp`, `nodeVel`, `nodeAccel`, `reaction` | `node`, `dof` |
| `elementForce`, `elementLoad` | `elementKind`, `elementIndex`, `component` |
| `modeShape` | `mode`, `node`, `dof` |
| `fiber` | `elementKind`, `elementIndex`, `point`, `fiber`, `quantity: "strain" \| "stress"` |
| `gaussPoint` | `elementKind` (`tri3`/`quad4`), `elementIndex`, `point` (`tri3`: 0; `quad4`: 0..4 in node order of the 2x2 rule), `quantity: "strain" \| "stress" \| "principalStrain" \| "principalStress" \| "kc2"`, `component` (0..3 = x, y, xy; engineering shear for strain). The enhanced quad reports the recovered values including the internal modes. Static and transient stages only. For `principalStrain`/`principalStress`, components `0, 1, 2` are `eps_1 >= eps_2` (`sigma_1 >= sigma_2`) and `theta`, the counter-clockwise angle in radians from global x to the first principal direction; a shell records nothing for them. `kc2` has one component (`0`) and needs a `csfmConcrete` point (`recorderNotApplicable` otherwise). `elementKind: "embeddedBar"` addresses a segment by its flat index and reports axial strain/stress in component 0. |
| `embeddedBar` | `bar`, `segment`, `quantity: "axialForce" \| "strain" \| "stress" \| "ruptured"` (`axialForce = stress * area`, tension positive; `ruptured` is 1 past the rupture strain, else 0). Static and transient stages only. |
| `bulk` | `table: "triangles" \| "quads" \| "embeddedBars" \| "nodes"`, `quantity: "strain" \| "stress" \| "principalStrain" \| "principalStress" \| "kc2" \| "axialForce" \| "ruptured" \| "displacement" \| "reaction"`. Records every element and Gauss point of the table at every accepted step as flat arrays, in `StepOutcome.bulkBatches` instead of `recorderBatches` (this recorder's slot in `recorderBatches` stays empty). `axialForce`/`ruptured` belong to `embeddedBars` (`strain`/`stress` work on all three element tables, one value per segment); principal values and `kc2` belong to the continuum tables, and `kc2` needs every point to be `csfmConcrete`. `nodes` takes `displacement` or `reaction` only (and they take only `nodes`): every node in input node-table order, one point, `ndf` components (`ux, uy, rz` in 2D; the six DOFs in 3D), a fixed DOF reporting its prescribed displacement (`0` at a support), and a free DOF's reaction being only the convergence residual. A reaction frame costs one pass over the model, not one per DOF; in a `transient` stage it includes the ground-motion inertial load, as the scalar `reaction` recorder does, and costs a sweep per DOF. Static and transient stages only. |

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
| `bulkBatches` | `BulkBatch[]`: this call's `bulk` recorder output. |
| `stepDiagnostics` | `StepDiagnostic[]`: the solver work of each `static` step this call took (see below). |
| `limitHistory` | `LimitSample[]`: the controlling CSFM material-limit utilization at each sample (see below). |

`RecorderBatch`: `{ recorderIndex, stageIndex, firstSample, samples: [pseudoTime, value][] }`.
`firstSample` plus `samples.length` gives the sample range, so a failed persist
can be retried without renumbering. With arc length the load factor can
decrease and repeat, so order samples by sample index, never by load factor.
The caller accumulates batches; there is no "history so far" accessor.

`BulkBatch`: `{ recorderIndex, stageIndex, firstSample, elements, points, components, progress, values }`.
`progress` and `values` are `Float64Array`s, not `number[]`: each is a fresh buffer of its own (not a view of
wasm memory), so a worker can list `values.buffer` in `postMessage`'s transfer list and hand the frames to the
main thread without a copy. `values` holds `progress.length` consecutive frames of `elements * points * components` numbers,
in element, then point, then component order (element = table row, or the flat segment index for
`embeddedBars`); `progress[i]` is frame `i`'s load factor or time. `firstSample` continues across
`advance` calls like `RecorderBatch`'s.

`LimitSample`: `{ progress, utilization, kind, elementKind, elementIndex }` for each sample a
`static` stage recorded in the call (every accepted step, plus the initial state with
`header.recordInitial`). `utilization >= 1` means a limit is reached. `kind` is
`"concreteCrushing" | "concreteTensileStrain" | "rebarRupture" | "rebarStrainLimit"` (the default
report checks crushing at `|eps_cu2| = 0.0035` and rupture only; the other two appear only through
`stop.materialLimit`). The location is element-level: `elementIndex` is the row in the table of
`elementKind` (the flat segment index for `embeddedBar`, mapped to `(bar, segment)` by
`modelOutput()`). The list is empty for a model with no `csfmConcrete` or `csfmRebar`, and for
`transient` and `modal` stages.

`StepDiagnostic`: `{ progress, iterations, factorizations, cutback? }` for each step a `static` stage took in
the call, in order, so entry `k` belongs with frame `k` of every `bulkBatches` entry and with `limitHistory[k]`
(`header.recordInitial`'s initial sample has a frame and a limit entry but no step, so it shifts them by one).
`progress` is the load factor the step ended at; `iterations` the Newton iterations (corrector solves under arc
length) of the accepted attempt; `factorizations` the tangent factorizations over the step, failed attempts
included. Under `loadControlCutbacks`, `cutback` is `{ nominalIncrement, acceptedIncrement, retries }`: the step
asked for `nominalIncrement`, halved `retries` times, and landed `acceptedIncrement`, the rest being asked again
by the next step. `cutback` is absent for any other integrator. The list is empty for `transient` and `modal`
stages.

`ContinuationDetail`: `radius, `nextRadius`, `retries`, `correctorSolves`,
`factorizations`, `forceMeasure`, `arcError`, `chordLength`, `detSign`,
`bifurcationSuspected`, `displacementScale`, `rotationScale`, `stop`
(`{ reason, landedExactly, overshoot }`, where `reason` is
`"displacementTarget" | "loadFactorTarget" | "loadFactorZeroCrossing" | "chordLength" | "stepCount" | "materialLimit"`).

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
`cutbacksExhausted {step, attempts, radius, lastFailure}`, `continuationComplete`,
`loadControlExhausted {step, acceptedLoadFactor}` (a `loadControlCutbacks` stage halved down to
`minIncrement` and still failed; `acceptedLoadFactor` is the last committed load factor).

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
| `invalidGaussPoint` | `recorder`, `point`, `count` (the element's Gauss-point count, 0 for a kind without any) |
| `invalidPlaneMaterial` | `index`, `reason` |
| `invalidMaterial` | `index`, `reason` (a `materials` entry with unusable constants, today a `csfmRebar`) |
| `invalidEmbeddedBarSegment` | `recorder`, `bar`, `segment`, `count` (an `embeddedBar` recorder past the bar's piece count) |
| `recorderNotApplicable` | `recorder`, `reason` (a quantity the target cannot produce: `kc2` off a `csfmConcrete` point, an empty `bulk` table, a `bulk` quantity the table does not have) |
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
