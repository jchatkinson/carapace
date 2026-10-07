# CSFM integration feasibility

Status: investigation only. No code has been written, and no performance has
been measured; the performance notes below are estimates from reading both
codebases.

Question: can the Compatible Stress Field Method (CSFM) work in the sibling
`CalcsApp` repo (`lib/fem`, TypeScript) be ported into Carapace and run on
the Rust solver?

Answer: yes. The material and rebar mechanics port easily, and Carapace's
solution algorithms are already better than the TS ones. The cost is mostly in
the Rust `Domain`, which has no continuum elements today.

Source of the CSFM work: `CalcsApp/lib/fem`, with design history in
`CalcsApp/docs/concrete/csfm/` (start with `CSFM_CONTEXT.md` and
`phase8-solver-robustness.md`).

## What exists in CalcsApp

CSFM there is a 2-D plane-stress concrete membrane analysis: a constrained
Triangle mesh of Q4/T3 elements, point loads, supports, nodal springs, and
embedded rebar, solved by incremental Newton-Raphson.

| Piece | Size | Port? |
|---|---|---|
| Q4/T3 plane-stress element (`shape`, `bmatrix`, `quadrature`) | ~100 lines | Yes. Standard. |
| CSFM concrete law (`nonlinear-material.ts`, `makeConcreteContinuumLaw`) | ~250 lines | Yes. Stateless rotating principal-strain law, parabola-rectangle compression, Kaufmann kc2 softening, optional tensile regularization. |
| Rebar laws (tension stiffening POM/TCM, bare bilinear steel) | ~200 lines | Yes. Pure functions. |
| Embedded-rebar discretization (`rod.ts`) | ~140 lines | Yes. Preprocessing: split each bar at element edges, locate the host element. |
| Newton engine (`nonlinear.ts`: line search, secant-Newton, retry bisection) | 685 lines | No. Carapace already covers it (below). |
| Sparse LU and SuperLU WASM glue | n/a | No. Carapace uses faer. |
| Triangle meshing, Valtio store, worker | n/a | Stays in JS. The mesh crosses the boundary as arrays. |

## What Carapace already provides

- Arc-length continuation with softening and snap-back. CalcsApp phase 8 lists
  its absence as the remaining limitation (force control stalls at the limit
  point).
- Displacement control, Krylov-Newton, regula-falsi/bisection line search,
  tangent reuse with a cached factorization, `Combined` per-equation
  convergence, staged analysis, modal and transient analysis.
- Crisfield secant-Newton does not need a port. CalcsApp's own phase 8
  measurements show full Newton with line search is the best configuration,
  and Krylov-Newton covers the same ground.
- Springs, supports, prescribed displacements. A `ZeroLength` to a fixed
  ground node replaces `spring.ts`; fixed DOFs are exact rather than
  penalized.

## Gaps in Carapace

1. **Element framework: resolved.** Elements now have a variable node list and push their own
   contribution through sinks; Quad4/Tri3 exist (see [architecture.md](architecture.md) and
   [domain-refactor-plan.md](domain-refactor-plan.md)). No separate `membranes` store is needed.
2. **Plane material seam exists; the CSFM law does not.** `PlaneMaterial` provides the 3-component
   trial/commit seam (per-Gauss-point copies), with elastic variants only. A
   `PlaneStressMaterial` for the CSFM law is still needed. The CSFM law is stateless, so
   Carapace's trial/commit model fits trivially, and the catalog can grow
   history later.
3. **Constraints: resolved for the general case, still not needed for rebar.** General linear
   multi-point constraints and `rigid_link` now exist ([architecture.md](architecture.md)). Embedded
   bars should still contribute directly to the host membrane's stiffness (a bar segment's axial strain
   is linear in the host element's nodal DOFs), with no extra DOFs or constraints.
4. **Nodes carry `rz` (resolved).** DOF activation leaves `rz` inactive on membrane-only nodes and lets a beam and a membrane share a node. Original note: Planar nodes are `(ux, uy, rz)`. Membrane nodes would
   fix `rz`, which costs nothing, but a frame element cannot share a membrane
   node (there is no drilling DOF). Open question: keep planar nodes with `rz`
   fixed, or add a 2-DOF node profile.
5. **Wasm input and results: partly resolved.** `planeMaterials`, `triangles`, `quads` and
   `gaussPoint` recorders (strain, stress) exist. Still missing: embedded bars and CSFM-specific
   results (principal strains/stresses, kc2, bar stress). Mechanical work; the fiber
   response recorder (`element_fiber_responses`) is the template.

## Risks and notes

- **Tangent correctness.** This is where the TS lib went wrong: a secant used
  as a tangent (`phase7b`), then a sparse-LU pivoting bug misread as a
  material problem (`phase8`). The concrete tangent there is a regime-aware
  finite difference. In Rust, derive it analytically and keep the finite
  difference as a test oracle; build a generic tangent-vs-finite-difference
  check into the tests from the start.
- **Tensile regularization (default gamma 0.1).** phase 7b concluded the LU
  bug, not the law, caused the stalls. faer's LU is robust, so the
  regularization may be reducible. Re-test; do not assume.
- **Nonsymmetric tangent.** kc2 depends on the maximum principal strain, so
  the tangent is nonsymmetric. faer's general LU handles that.
- **Stop criteria.** `findCriticalLoad` bisects on material stop criteria.
  Carapace's arc-length stop criteria with exact landing is the likely
  replacement.
- **Validation.** Dump golden vectors from the TS tests (material points, rod
  laws, the span-beam and deep-beam models) and assert Rust agrees before
  changing any behavior. This is the same pattern as the OpenSees comparison
  in `comparison/`.
- **Units.** CalcsApp runs mm/N/MPa. Carapace is unit-agnostic, so nothing
  changes.

## Performance (estimate, unmeasured)

- Assembly and material evaluation should improve a lot. The TS lib allocates
  `number[][]` per Gauss point and makes about 13 `evaluateConcrete` calls per
  Gauss point for the finite-difference tangent. Rust with stack-allocated 3x3
  matrices and an analytic tangent should be roughly an order of magnitude
  faster there.
- The linear solve is already SuperLU in WASM, and faer is not guaranteed to
  beat it. For meshes up to a few thousand DOF, assembly dominates and Rust
  should win clearly. For very large meshes the gain shrinks.
- Packaging: `pkg/carapace_wasm_bg.wasm` is about 2 MB, and CalcsApp does not
  depend on Carapace yet, so it needs a build/consume step.
- The spike should report a real number (span beam and deep beam) before the
  rest is committed to.

## Suggested order

1. **Spike.** Standalone module in `core`: T3/Q4, the CSFM law, embedded bars,
   a static `Domain` path, golden tests against the TS lib, one timing
   comparison.
2. **Domain integration.** Membrane store, commit and recorders, then
   arc-length on a softening CSFM model.
3. **Wire protocol.** New input tables, decode, results.
4. **CalcsApp cutover.** Replace `solveNonlinear` in `femWorker.js` behind a
   flag; retire `nonlinear.ts` and the TS sparse LU once parity holds.

Steps 1 and 2 are the real work. Step 3 is mechanical and step 4 is mostly
packaging.
