# CSFM in Carapace: implementation and test plan

Status: plan. Nothing in this document is implemented yet.

Goal: move the Compatible Stress Field Method (CSFM) analysis from CalcsApp (`../calcs.app/lib/fem`, TypeScript)
into Carapace, so that Carapace becomes the solver CalcsApp ships. This plan covers only the CSFM procedure and its
tests in Carapace. Moving CalcsApp's other ad-hoc FEM calculations to Carapace is a later, separate step.

This plan replaces the "Suggested order" section of [csfm-integration.md](csfm-integration.md). That document's
feasibility findings still apply. Since it was written, the element framework, `PlaneMaterial`, Tri3/Quad4, DOF
activation and Gauss-point recorders have been built. What is still missing is the CSFM-specific work below.

Sources in CalcsApp:

| What | Where |
|---|---|
| Overview and checklist | `docs/concrete/csfm/CSFM_CONTEXT.md` |
| Concrete law | `lib/fem/nonlinear-material.ts` (`evaluateConcrete`, `resolveConcreteProperties`) |
| Continuum adapter, rebar laws, model resolution | `lib/fem/csfm-nonlinear.ts` (`makeConcreteContinuumLaw`, `makeRodLaw`, `resolveModel`, `findCriticalLoad`) |
| Embedded rebar | `lib/fem/rod.ts` (`resolveRebarPaths`, `findHostBinding`), assembly in `nonlinear.ts` (`rodContribution`) |
| Newton driver | `lib/fem/nonlinear.ts` (`solveLoadPath`, `attemptIncrement`) |
| Worker integration | `app/concrete/csfm/femWorker.js` (`solveNonlinear`) |
| Tests | `test/fem/csfm-nonlinear.test.ts`, `rod.test.ts`, `spring.test.ts`, `sparse.test.ts` |
| Design history | `docs/concrete/csfm/phase5…phase8*.md`. Read `phase7a`, `phase7b` and `phase8` before changing a law or the solver. |

## 1. What has to be reproduced

The CalcsApp procedure, which is the parity target:

1. **Model.** A 2D plane-stress membrane with thickness `t`, meshed as T3/Q4 by the caller in mm/N/MPa. It has fixed
   DOFs, nodal point loads, nodal springs (`kx`, `ky`), and embedded rebar polylines.
2. **Concrete** (`evaluateConcrete` + `makeConcreteContinuumLaw`). This is a stateless rotating principal-strain law
   with `nu = 0`.
   - Properties come from `fck` (0 < fck ≤ 50 MPa): `Ecm = 22000 ((fck+8)/10)^0.3` (an override is allowed),
     `fctm = 0.3 fck^(2/3)`, `eps_c2 = -0.002`, `eps_cu2 = -0.0035`, `eta_fc = min(1, (30/fck)^(1/3))`.
   - Compression is either a parabola-rectangle (peak `eta_fc fck` at `eps_c2`) or elastic-ideally-plastic (`Ecm`,
     capped at `eta_fc fck`).
   - Kaufmann softening applies to the minor principal direction when `eps_1 > 0`:
     `kc2 = min(1, 1 / (0.8 + 140 eps_1))`.
   - Tension is zero. Optional numerical regularization adds `sigma_t = gamma fctm (1 - exp(-eps/eps_reg))`, with
     `eps_reg = gamma fctm / Ecm` derived rather than configured, so the tension tangent at the origin is `Ecm`. An
     optional Hermite crack-transition band joins compression and tension. Absent regularization must give the exact
     zero-tension law.
   - At the exact origin (all strains < 1e-9) the tangent is seeded as isotropic
     `diag(E0, E0, E0/2)`, where `E0` is the compression origin tangent.
   - The Cartesian tangent in TS is a regime-aware finite difference. It is not rotated secants (phase 7a §2).
3. **Rebar** (`resolveModel`, `makeRodLaw`, `makeBareRodLaw`). This is stateless and depends on the path.
   - In compression, or for every strain when tension stiffening is off: bilinear bare steel (`Es`, `fy`, `ft`,
     `eps_u`) with rupture.
   - In tension with tension stiffening on: the inverse of the Tension Chord Model (TCM) average strain, used when
     `rho_eff ≥ rho_crit`, or the Pull-Out Model (POM) otherwise. Inputs are `tau_b0 = 2 fctm`, `tau_b1 = fctm`,
     crack spacing `sr = lambda d fctm (1 - rho)/(2 tau_b0 rho)` (default `lambda = 0.67`), and the branch boundaries
     computed at resolve time.
   - TCM has a C1 cubic startup bridge below `sigma_min / Es`. POM has an origin guard.
   - The tangent is the analytic reciprocal slope (phase 7b §1), bounded to `[Es/1e4, 10 Es]`.
   - Past `eps_m(ft)` the bar is ruptured: stress `ft`, with a small positive secant tangent.
   - Validation: `0 < rho_eff < 1` and `ft > fy`.
4. **Embedding.** Each polyline is split at every mesh-edge crossing. Each piece lies in one host element. Its two
   end points are bound to host shape-function values, and its axial strain is linear in the host DOFs. No extra
   DOFs or constraints are added.
5. **Procedure.** Load control is applied in `n` equal increments of the full load (default 20). Each increment
   runs full Newton with the Crisfield slope line search.
   - Convergence requires both a relative energy criterion (`energyTol = 1e-8`) and a free-DOF residual below
     `1e-3 · max(||F||, 1)`.
   - A failed increment is halved, down to `minLoadFraction = 1/128`. The unconverged remainder is put back in the
     schedule.
   - An aborted run is still a useful result: it reports the last accepted fraction and the retry history.
6. **Outputs used by the UI.** Displacements, Gauss-point stress and strain (used for contours, section cuts and
   principal directions), bar force and stress per rod segment, reactions, spring forces, per-path TCM/POM routing,
   and per-increment progress.

The following stay in CalcsApp: meshing (Triangle), quad recombination, unit conversion (in/kip/ksi to mm/N/MPa),
contour shaping, nodal stress recovery and section-cut resultants. Carapace receives a mesh and returns Gauss-point
and bar results.

## 2. Design decisions

### 2.1 Concrete: `PlaneMaterial::CsfmConcrete`

Add a new `PlaneMaterial` variant, boxed if needed to keep the enum small, in a new
`core/src/model/materials/csfm_concrete.rs`.

- **Constructor.** It takes `fck`, an optional `ecm`, `law`, `softening`, and an optional
  `TensileRegularization { gamma, transition_strain }`. It resolves the derived properties once. It validates eagerly
  and returns `PlaneMaterialError` variants for an out-of-range `fck`, `0 < gamma < 1` violated, or a non-positive
  transition strain.
- **Stateless.** `commit` returns `self`, and `is_linear()` is false, so `Quad4Formulation::Enhanced` keeps
  rejecting it.
- **Analytic consistent tangent.** It replaces the TS finite difference. Let `(eps_1, eps_2, theta)` be the
  principal strains and direction, and `T(theta)` the stress rotation. Then

  ```
  D = T^T [[dσ1/dε1, 0,       0 ],
           [dσ2/dε1, dσ2/dε2, 0 ],
           [0,       0,       G*]] T_eps,      G* = (σ1 - σ2) / (2 (ε1 - ε2))
  ```

  - `dσ2/dε1` is the kc2 coupling. It is nonzero only when softening is on, `eps_1 > 0` and `eps_2 < 0`. Because of
    it, `D` is **nonsymmetric**.
  - When `|ε1 - ε2|` is near zero, `G*` takes its limit `dσ/dε / 2` for the shared branch.
  - The branch at the kink `eps_1 = 0` follows the TS regime rule: `eps_1 ≤ 0` counts as compression.
  - The exact origin keeps the TS isotropic seed.
  - Write `1 - exp(-x)` as `-expm1(-x)`.
- **Solver fit.** The faer solver already falls back from Cholesky to LU for a nonsymmetric matrix. No solver change
  is needed. A test still pins the behavior (§4, H).
- **Limit-state queries** for §2.5: principal strains, `kc2` and a compression utilization `-eps_2 / eps_cu2`.

### 2.2 Rebar: `Material::CsfmRebar`

Add a new uniaxial `Material` variant in `core/src/model/materials/csfm_rebar.rs`. Making it a uniaxial `Material`
means the existing fiber, truss and probe machinery can use it, and PySees' Material Preview gets it free through
`material_probe.rs`.

- **Constructor.** It takes the steel properties, the bond and concrete inputs (`fctm`, `Ecm`, `d`, `rho_eff`,
  `lambda`) and `tension_stiffening: bool`. It resolves routing, `sr`, `sigma_min` and the branch boundaries once,
  exactly as `resolveModel` does.
- **Behavior.** It is stateless and nonlinear-elastic: unloading retraces the loading curve, as in TS. This is
  monotonic-only by construction. Document that, and do not imitate hysteresis.
- **Rupture.** Keep the TS response (stress `±ft`, small positive tangent) so that parity holds. Also expose a
  `ruptured(strain)` query so that rupture is reported (§2.5) and not just absorbed.
- **Reuse the TS logic for the inverse.** The inverse is a closed-form quadratic per branch (`chooseQuadratic`).
  Keep it closed-form rather than switching to a Newton inverse, so that it stays bit-close to TS. Errors become
  `MaterialError`, never panics.

### 2.3 Embedded bars: `EmbeddedBar` element (2D only)

Add a new `Element` variant in `core/src/model/elements/embedded_bar.rs`.

- **Nodes.** Its nodes are the **host element's** 3 or 4 nodes, so it fits `MAX_ELEMENT_NODES = 4` and the existing
  sink path (`sink.add::<6>` or `sink.add::<8>`). It stores, for each end, the shape-function values
  `N_a[i], N_b[i]` in the host, plus the unit direction, length, area and a `Material`.
- **Strain.** `eps = d · (Σ N_b,i u_i - Σ N_a,i u_i) / L`. The tangent is `(E_t A / L) Aᵀ d dᵀ A`.
- **`dof_mask`.** `ux` and `uy` of the host nodes only. A bar therefore never activates `rz`.
- **Binding.** **Bind both ends to the host that contains the segment midpoint.** TS binds each end point on its own
  and can pick a neighbor for an end that sits on a shared edge. This changes nothing: shape functions of nodes off
  that edge vanish on it. But it keeps every segment within one element's ≤4 nodes. Verify this in the tests (§4, D)
  instead of assuming it.
- **Responses.** `local_force` returns axial force. Use the existing fiber/Gauss pattern for strain and stress, and
  add a `ruptured` flag.

Discretization (polyline to segments) goes in a **core** model-building helper,
`model/continuum/embed.rs::embed_polyline(&Domain, polyline, …) -> Result<Vec<EmbeddedBar>, ModelError>`, and not
in JS. Reasons:

- Carapace then owns and tests it.
- PySees can use it too.
- The wire format takes polylines, not precomputed shape functions.

The helper needs:

- edge-crossing splits;
- a T3 closed-form inverse and a Q4 Newton inverse of the isoparametric map (port `naturalCoordinatesAt`);
- tolerance scaled to element size;
- an error for a point outside the mesh (`ModelError::EmbeddedBarOutsideMesh`) and for an empty or degenerate path.

It handles a bar through a mesh node and a bar along a mesh edge. Start with brute force and a bounding-box filter,
as TS does. Add a spatial index only if the benchmarks call for it.

### 2.4 Supports, springs and loads

- **Supports.** Exact fixed DOFs, not penalty. Golden comparisons therefore check displacements, Gauss states and bar
  forces. They do not compare penalty-inflated reactions in TS: recompute those as `internal - external`.
- **Springs.** Each spring becomes a `ZeroLength` with `Elastic` materials to a fixed ground node, one per sprung
  node. A spring-only support (no fixed DOFs) must still be accepted.
- **Loads.** Nodal loads go in one load pattern, scaled by the load factor.

### 2.5 Procedure: static analysis, cutbacks, limit states, arc-length

1. **Load-control cutbacks (a generic core feature that is still missing).** Today only arc-length has
   transactional cutbacks. Add optional adaptive substepping to `Integrator::LoadControl`: on non-convergence,
   revert, halve, and retry down to `min_increment`; then restore the remainder.
   - This is the TS `solveIncrementWithSession` and remainder logic, and it is useful beyond CSFM.
   - Report each cutback in the step result, as the source of CalcsApp's `retryHistory`, and fail with
     `CutbacksExhausted` while keeping the last accepted state.
   - Optional growth back toward the nominal step (Crisfield eq. 20) can come later.
2. **Algorithm and convergence.** Use the existing Newton, `TangentStrategy::Current`, and Carapace's `LineSearch`
   (`RegulaFalsi`/`Bisection` with `max_eta > 1`).
   - Do **not** port the Crisfield line search or secant-Newton up front. TS phase 8 found that full NR plus line
     search is the best configuration, and the line search was worth about 8% in iterations.
   - Port the Crisfield line search only if the fixtures show Carapace's line search converging worse (§4, F).
   - Use `ConvergenceTest::Combined` with a relative force tolerance. Its per-equation reference makes it better
     than the TS norm-with-`max(·,1)` floor, which phase 8 §6 calls a dimensional assumption.
   - Parity tests use tight tolerances on both sides, so that the different convergence tests cannot account for
     any difference.
3. **Material limit states.** CSFM's answer is a load at which a material limit is reached, not "converged at 100%".
   - Add an element query `limit_utilization() -> Option<(f64, LimitKind)>`, where values ≥ 1 mean violated. Kinds
     are concrete crushing (`-eps_2 / eps_cu2`, optionally against a user strain limit), rebar rupture
     (`eps / eps_m(ft)`), and an optional rebar strain limit. A `Domain` reduction returns the controlling element
     and point.
   - Add a stop criterion `material_limit: Option<{ exact: bool }>` to `ArcStopCriteria`, and the same check after
     each accepted load-control step.
   - Exact landing reuses arc-length's existing bracketing and retake. The utilization is a scalar target, exactly
     like `DisplacementTarget`.
   - This replaces `findCriticalLoad`'s full re-solve bisection.
4. **Arc-length past the peak.** Load control cannot pass the flexural limit point. TS stops at 97.5% on the span
   beam at 50 kip per point (phase 8 §5).
   - Run the CSFM fixtures under `Integrator::ArcLength` and check that it passes the peak.
   - Crisfield warns that arc-length can fail under steep softening, and this law has steep softening (kc2).
     Treat arc-length as validated only once §4, G passes.

### 2.6 Results and recorders

- The existing `gaussPoint` strain and stress recorder stays.
- Add `quantity: "principalStrain" | "principalStress"` (components `1, 2, theta`) and `"kc2"` for `CsfmConcrete`
  points.
- Add an `embeddedBar` recorder: axial force, strain, stress, ruptured.
- Add a bulk form that records every element and point of a table each step as flat arrays. CalcsApp plots every
  accepted increment, and per-point recorders would make the batches too large.
- The model output reports the resolved TCM/POM routing and `sr` per embedded bar group, as `resolveModel`'s
  `resolvedRodMaterials` does now.
- Each step reports the material-limit utilization and its location.

### 2.7 Wire format (`input_v1`)

- `planeMaterials` gains `{ kind: "csfmConcrete", fck, ecm?, law, softening, tensileRegularization? }`. A regularized
  material must be labeled as numerical regularization in docs, not described as a material property (phase 7a §3).
- `materials` gains `{ kind: "csfmRebar", es, fy, ft, epsU, fctm, ecm, diameter, rhoEff, lambda, tensionStiffening }`.
  An alternative is to resolve `fctm` and `ecm` from a referenced concrete. **Decision needed (§6).**
- A new table `embeddedBars`: `{ polyline: [x, y, …], area, material, hosts?: "all" | element index list }`.
  It is decoded through `embed_polyline`. Results are indexed by `(bar, segment)`, and the segment table is returned
  so the UI can draw it.
- Load-control `cutbacks`, `stop.materialLimit`, and the recorders in §2.6.
- New `DecodeError` variants for every new failure. Mirror the change in PySees (`carapaceInputV1.ts`,
  `compileInputV1.ts`) as `AGENTS.md` requires, even if PySees has no CSFM UI: the type must compile.
- Document all of it in `input-format.md`.

## 3. Phases

Each phase ends with `cargo test`, `cargo fmt` and `cargo clippy` clean. Work goes test-first: the fixtures from
phase 0 define "correct" before any law is written.

| Phase | Deliverable | Exit criterion |
|---|---|---|
| 0. Fixtures and harness | Golden export script in CalcsApp. A generic tangent-vs-finite-difference checker in Carapace. | Fixtures committed. The checker passes on existing materials and elements. |
| 1. Materials | `CsfmConcrete`, `CsfmRebar` | §4 A–C pass. |
| 2. Embedding | `EmbeddedBar`, `embed_polyline` | §4 D passes. |
| 3. Static procedure | Load-control cutbacks. Fixtures run natively. | §4 E–F and H pass. Parity with TS on every fixture. |
| 4. Limit states and arc-length | `limit_utilization`, `material_limit` stop, arc-length runs | §4 G passes. |
| 5. Wire and results | Tables, recorders, decode errors, docs, PySees type mirror | §4 I passes. `boundary-smoke.ts` case added. |
| 6. Benchmark | CSFM benchmarks against CalcsApp TS | Numbers committed to `docs-site/benchmarks/`. |
| 7. CalcsApp cutover (outside this plan's scope) | `femWorker.js` calls Carapace behind a flag | See §5. |

**Phase 0 detail.**

- **Export script.** `../calcs.app/test/fem/export-csfm-golden.ts`, run with CalcsApp's existing test runner. It
  imports `lib/fem` and writes JSON to `core/tests/csfm/fixtures/`. Each file records the CalcsApp commit it was
  generated from. Never hand-edit a fixture; regenerate it.
- **Dependencies.** Add `serde` and `serde_json` as **dev-dependencies** of `carapace-core` to read the fixtures.
  This does not touch the wasm build.
- **Fixture contents:**
  - material points (strain → stress, plus the TS finite-difference tangent) on a grid covering every regime,
    rotated states, and both laws with and without softening, regularization and the transition band;
  - rebar (strain → stress and tangent) on every TCM/POM branch, at the startup bridge, at rupture and in
    compression;
  - embedding (polyline → segments, host and `N`);
  - model runs: per-increment displacements, Gauss strains and stresses, bar forces, accepted fractions.
- **The checker.** A `core/tests` helper that, for any `PlaneMaterial` or `Material` at a strain, or any element at a
  displacement state, compares the returned tangent with a central difference of the returned stress or resistance.
  The step should scale with magnitude and should not straddle a regime kink. Every new law and element must pass
  it. This is the check TS lacked until phase 7b (`csfm-integration.md`, "Tangent correctness").

## 4. Tests

Test layout follows `AGENTS.md`:

- Material laws: `core/tests/materials/csfm_concrete.rs` and `csfm_rebar.rs`.
- Element: `core/tests/continuum/embedded_bar.rs`.
- A new area binary `core/tests/csfm/` (`main.rs`, one module per fixture family) for model-level runs.
- `wasm-bridge/tests/wire/csfm.rs`.

Each file gets a `//!` header with what it verifies and against which reference. Each test gets a `///` comment with
its reference and tolerance, which the docs-site test report shows.

There are two kinds of reference, and they are kept apart:

- **Parity** with the TS implementation (the golden fixtures). This shows the port is faithful.
- **Independent** references (closed form, hand checks, published examples). This shows the method is right.
  TS is not an independent reference.

**A. Concrete law (parity and closed form)**

- Stress matches the fixtures to 1e-12 relative at every grid point.
- Tangent matches the TS finite difference to ≤1e-5 relative, which is TS's finite-difference noise.
- Tangent matches the Rust central-difference checker to ≤1e-7 relative.
- Closed form:
  - parabola peak `-eta_fc fck` at `eps_c2`;
  - elastic-ideally-plastic yield;
  - `kc2` at chosen `eps_1`;
  - `Ecm`, `fctm` and `eta_fc` for fck = 20, 30, 50;
  - regularization origin tangent equal to `Ecm` (phase 7a correction) and the cap `gamma fctm`.
- Regularization absent gives exactly zero tension.
- The transition band is C1 at `±width`.
- The origin seed is `diag(E0, E0, E0/2)`.
- A rotated state couples normal and shear terms.
- The kc2 coupling term makes `D` nonsymmetric when softening is on, and symmetric when it is off.
- Invalid `fck` and `gamma` are errors, not panics.

**B. Rebar law (parity and closed form)**

- Routing: `rho_crit = fctm / (fy - (Es/Ecm - 1) fctm)`, checked on both sides of it.
- Every branch round-trips: TS forward `eps_m(sigma)` → Rust inverse → sigma, to 1e-9.
- The TCM startup is C1 at `sigma_min / Es`, and the stress at the startup strain is `2 sigma_min`. The stress just
  above zero carries no `sigma_min` jump.
- The POM origin guard holds.
- The tangent equals `1/slope` and passes the finite-difference checker on every branch, but is clamped to
  `[Es/1e4, 10 Es]`.
- Rupture gives stress `ft` and a positive tangent, and `ruptured()` is true.
- Compression is bilinear. Bare mode is bilinear in both directions.
- `rho_eff` out of `(0, 1)` and `ft ≤ fy` are errors.
- Parity with the TS law at every fixture strain to 1e-12.

**C. Independent material-level validation**

These check the method, not the port. Sources are in `../calcs.app/docs/concrete/csfm/`.

- A single-element reinforced panel under pure shear and under biaxial loads, compared with the Cracked Membrane
  Model (Kaufmann and Marti 1998) solution, using the closed-form CMM equations for an orthogonally reinforced
  panel. The tolerance is set from the model gap, not tuned to pass.
- A uniaxial tension chord (bar in a concrete prism) compared with the TCM average strain–stress relation: the
  member-level response must reproduce the material law it embeds.

**D. Embedded bar element and embedding**

- Patch tests on T3, Q4 and a distorted Q4 under uniform strain: bar strain equals `dᵀ ε d` of the host for any
  position and angle, to 1e-12.
- Rigid translation gives zero force.
- A bar coincident with a mesh line matches a `Truss` on shared nodes exactly.
- A bar crossing many elements gives a force independent of how it was split (cut it at extra points).
- Mask conformance (added to the test in `elements/mod.rs`). The finite-difference checker passes on a nonlinear
  `CsfmRebar`.
- Midpoint-host binding gives the same strain as TS's per-end binding (fixture), including ends on shared edges and
  at nodes.
- `embed_polyline` against the fixtures: segment count, points and hosts. Ported cases from `rod.test.ts`:
  - 36 segments on the deep beam and 30 on the span beam;
  - a bar through a node;
  - a bar along an edge;
  - zero-length vertices dropped;
  - outside the mesh is an error.

**E. Linear sanity on the CSFM fixtures**

- Linear deep beam with embedded bars (`expectLinearDeepBeamToBalance`): one Newton correction balances it to
  1e-10.
- Linear-elastic displacements match TS `solveLinearElastic` on the deep beam, with penalty effects excluded by
  comparing free DOFs at a matched penalty-free reference.

**F. Nonlinear model runs (parity with `csfm-nonlinear.test.ts`)**

Each fixture must reach the TS outcome, and must match TS's converged state at each increment: displacements to 1e-6
relative and bar stress to 1e-5. The cases:

- Q4 first increment from zero displacement (origin seed).
- 2x2 in patch, first 5%.
- T3 2x2 in patch at 50 kip/node, first 5%.
- 8 in patch 1x1 and 2x2, each with the six option combinations (plain or transverse bar × regularization × tension
  stiffening × softening). All six reach 100%.
- Pinned deep beam (12x6 Q4, 4 bars, 7.5 kip, gamma 0.01): reaches 100%.
- Span beam (15x7 Q4, 2+2 #5, 10 kip, gamma 1e-3): reaches 100%. The maximum bar stress is within 0.5–1.5x the
  67 MPa hand check (`M ≈ 130 kip·in`, `jd ≈ 22 in`). That is an independent check: a low value means the
  regularization is carrying the tie.
- Gamma sweep on the span beam, gamma = 1e-1, 1e-2, 1e-3, 1e-4. All converge. Bar stress is gamma-independent to 2%
  for gamma ≤ 1e-3 (TS: 50.5 vs 51.3 MPa), and gamma = 0.1 visibly under-reports it (TS: 6.3 MPa). This test
  records the default question in §6.
- Load-control cutbacks: a step forced to fail retries at half size, restores the remainder, and reaches the same
  final state as a fine schedule. Exhausted cutbacks keep the last accepted state.
- Solver comparison, recorded and not asserted: iterations per increment for Carapace's line searches against the TS
  counts in phase 8 §4. If Carapace needs over 1.5x TS's iterations on the span beam at 30 kip, port the Crisfield
  slope line search.

**G. Limit states and arc-length**

- Span beam at 50 kip/point under load control stops near 97.5% (TS). The material limit reports the bar at yield
  plateau.
- Under arc-length the same model passes the peak. The peak load factor matches the load-control limit to 2%, and
  the descending branch continues to a stop criterion.
- `material_limit` with `exact`: lands on utilization 1 within the landing tolerance, for rebar rupture (a lightly
  reinforced tie) and for concrete crushing (an over-reinforced strut).
- The bracketing load agrees with TS `findCriticalLoad` on the same fixture to its `loadFractionTolerance`.

**H. Solver**

- Dump a cracked CSFM tangent from the span beam as a fixture. faer's LU solves it with relative residual < 1e-10.
  This is the phase 8 lesson: assert the linear-solve residual literally.
- A nonsymmetric CSFM tangent goes to LU and not Cholesky.

**I. Wire**

- Native vs wire, identical results for the span beam, including recorders.
- Each new `DecodeError` case: invalid `fck` or `gamma`, invalid rebar, a bar outside the mesh, a bar on a 3D
  profile, and `csfmConcrete` with the enhanced quad.
- Add a case to `boundary-smoke.ts`.

**Benchmarks (phase 6).** These are not part of CI.

- A comparison script times CalcsApp's TS `solveLoadPath` with SuperLU against Carapace native and Carapace wasm on
  the span beam, the deep beam, and a refined deep beam (`createPinnedDeepBeamFixture(48, 24)`).
- Report wall time, increments and iterations, and commit the JSON next to the existing benchmarks.
- `csfm-integration.md` predicts about 10x on assembly. This is where that prediction is measured.

## 5. CalcsApp cutover (summary only)

This is outside this plan, but it constrains its choices.

- **Packaging.** Build the wasm with `scripts/build-wasm.sh` and vendor `pkg/` into CalcsApp's `public/wasm/carapace/`,
  the way SuperLU is vendored today.
- **Worker.** `femWorker.js` builds a `CarapaceInputV1` from the same SI mesh, supports, springs, loads and
  `RebarPath`s. Triangle output must be counter-clockwise, since Carapace rejects clockwise elements. It maps the
  results back to the current result shape.
- **Rollout.** Run Carapace behind a `model.fem.solver` flag beside `solveLoadPath` until the UI regressions agree.
  Then retire `nonlinear.ts`, `csfm-nonlinear.ts`, `rod.ts` and the TS sparse LU/SuperLU path.

## 6. Open decisions

1. **Regularization default.** The CalcsApp UI ships `gamma: 0.1`, phase 7a recommends 0.01, and phase 7b shows
   gamma = 0.1 lets the concrete carry the span-beam tie (bar 6.3 vs about 67 MPa). Carapace should have **no**
   default and require an explicit value. The CalcsApp default should drop to ≤ 1e-3 at cutover. This needs the
   owner's decision.
2. **Rebar concrete inputs.** `csfmRebar` can either take `fctm` and `Ecm` explicitly, or reference a `csfmConcrete`.
   Explicit is simpler and keeps `Material` self-contained. A reference prevents mismatches. Recommendation:
   explicit in core, with an optional reference resolved by the wire decoder.
3. **Fracture-energy tension softening** (phase 7b, "Optional future change"). Out of scope. It has a negative
   tangent and would need arc-length, which phase 4 provides.
4. **Crisfield line search and secant-Newton.** Port them only on evidence from §4, F.
5. **Compression failure criterion.** TS has none in the worker path; `findCriticalLoad` exists but is unused by the
   UI. Default concrete limit: `eps_cu2` for parabola-rectangle. Confirm whether elastic-ideally-plastic should use
   the same limit.
6. **Bond-slip.** Out of scope. `rod.ts` keeps a `RodCoupling` seam for it. Here, `EmbeddedBar` would gain slip DOFs,
   which needs nodes it owns and is a later element.
