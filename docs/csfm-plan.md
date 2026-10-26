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
     `eps_reg = gamma fctm / Ecm` derived rather than configured, so the tension tangent at the origin is `Ecm`.
     Absent regularization must give the exact zero-tension law.
   - CalcsApp used to have an optional Hermite crack-transition band (`crackTransitionStrain`). It was applied to
     `sigma_max` only, so at equal principal strains the two directions used different laws, and the stress
     depended on the arbitrary principal direction. It was unused by the worker and has been **removed**: passing
     it is now an error (§6). Both principal directions now use one scalar law.
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
   - TCM has a cubic startup bridge below `sigma_min / Es`, which is C1 only when `2 sigma_min ≤ fy`. Past that the
     bridge overshoots yield and the stress jumps; for example, `lambda = 2` gave about 573.5 MPa just below the
     startup strain and 512.2 MPa at it. `resolveModel` now **rejects** `2 sigma_min > fy` (§6).
     - With TCM routing, `2 sigma_min ≤ lambda (fy - (Es/Ecm - 1) fctm)`, so this rejects only `lambda > 1`, which
       is outside the 0.5–1.0 range the UI offers.
   - POM has an origin guard: stress 0 and tangent `Es` up to a strain of 1e-9.
   - The tangent is the analytic reciprocal slope (phase 7b §1), bounded to `[Es/1e4, 10 Es]`.
   - Rupture, with a strict `>` comparison:
     - tension with tension stiffening: past `eps_m(ft)`, stress `ft`;
     - compression, and bare steel in either direction: `|eps| > eps_u`, stress `±ft`.

     In both cases the tangent is a small positive secant.
   - Validation: `0 < rho_eff < 1` and `ft > fy`.
4. **Embedding.** Each polyline is split at every mesh-edge crossing. Each piece lies in one host element. Its two
   end points are bound to host shape-function values, and its axial strain is linear in the host DOFs. No extra
   DOFs or constraints are added.
5. **Procedure.** Load control is applied in `n` equal increments of the full load (default 20). Each increment
   runs full Newton with the Crisfield slope line search.
   - Convergence: an iteration is accepted **as soon as** the free-DOF residual is at most
     `residualForceTol · max(||F||, 1)` (default `1e-3`), with no energy check ([nonlinear.ts:500](../../calcs.app/lib/fem/nonlinear.ts)).
     After the update, the iteration is accepted when the relative energy (`energyTol = 1e-8`, floor `1e-18`) and
     the residual both pass. An iteration diverges when the residual exceeds `residualForceCapRatio · ||F||`
     (default 100).
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
  `TensileRegularization { gamma }`. It resolves the derived properties once. It validates eagerly and returns
  `PlaneMaterialError` variants for an out-of-range `fck` or `0 < gamma < 1` violated. There is no transition band
  (§1).
- **Stateless.** `commit` returns `self`, and `is_linear()` is false, so `Quad4Formulation::Enhanced` keeps
  rejecting it.
- **Analytic consistent tangent.** It replaces the TS finite difference. Let `(eps_1, eps_2, theta)` be the
  principal strains and direction. Let `T_eps(theta)` map Cartesian engineering strain `[eps_x, eps_y, gamma_xy]` to
  the principal frame. This is the same matrix and convention as the `Orthotropic` arm in
  [plane.rs](../core/src/model/materials/plane.rs), which computes `D = Tᵀ D' T`, so stress transforms by its
  inverse transpose. Then

  ```
  D = T_eps^T C_p T_eps,     C_p = [[dσ1/dε1, 0,       0 ],
                                    [dσ2/dε1, dσ2/dε2, 0 ],
                                    [0,       0,       G*]],     G* = (σ1 - σ2) / (2 (ε1 - ε2))
  ```

  - `G*` is the engineering-shear rotation term for distinct principal strains.
  - `dσ2/dε1` is the kc2 coupling. It is nonzero only when softening is on, `eps_1 > 0`, `eps_2 < 0`, the `kc2 = 1`
    cap is inactive (`eps_1 > 0.2/140`), and the compression law is not on a flat branch: past `eps_c2` for the
    parabola-rectangle, or capped for elastic-ideally-plastic. Because of it, `D` is **nonsymmetric**.
  - The scalar branch derivatives respect the same caps: zero on the plateau, and `Ecm` or zero for
    elastic-ideally-plastic depending on whether the strength cap is active.
  - At `ε1 = ε2` both directions use the same scalar law (without the transition band, this always holds), so
    `G* → (dσ/dε)/2` and `D` does not depend on `theta`. Below a relative gap of about 1e-8, use that limit.
  - The branch at the kink `eps_1 = 0` follows the TS regime rule: `eps_1 ≤ 0` counts as compression.
  - The exact origin keeps the TS isotropic seed. That seed is a numerical policy, not the derivative of the law
    there (§4, A).
  - Write `1 - exp(-x)` as `-expm1(-x)`.
- **Solver fit.** The faer solver already falls back from Cholesky to LU for a nonsymmetric matrix. No solver change
  is needed. A test still pins the behavior (§4, H).
- **Limit-state queries** for §2.5: principal strains, `kc2`, and the utilizations in §2.5.

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
  `ruptured(strain)` query, with TS's sign- and mode-dependent thresholds (§1), so that rupture is reported (§2.5)
  and not just absorbed.
- **Errors happen at construction only.** `Material::trial_stress_tangent` returns `(f64, f64)`, and element
  assembly cannot fail ([materials/mod.rs](../core/src/model/materials/mod.rs), [ops.rs](../core/src/model/elements/ops.rs)).
  TS's `chooseQuadratic` instead throws during evaluation. The constructor therefore validates everything that
  makes each branch inverse well posed: each branch is monotone between the precomputed boundaries, so a root
  exists. Evaluation is then total, with no error channel added to the API.
  - Construction also **rejects** `2 sigma_min > fy` as a new `MaterialError`, rather than reproducing the
    discontinuous TS bridge (§1, §6).
  - Compute the quadratic roots in the cancellation-free form (`q = -(b + sign(b) sqrt(disc))/2`, roots `q/a` and
    `c/q`) and pick the root by branch. Do not keep TS's textbook formula just to stay bit-close.
  - Tests check the inverse's residual on every branch, and a `debug_assert` guards the branch selection.

### 2.3 Embedded bars: `EmbeddedBar` element (2D only)

Add a new `Element` variant in `core/src/model/elements/embedded_bar.rs`.

- **Nodes.** Its nodes are the **host element's** 3 or 4 nodes, so it fits `MAX_ELEMENT_NODES = 4` and the existing
  sink path (`sink.add::<6>` or `sink.add::<8>`). It stores, for each end, the shape-function values
  `N_a[i], N_b[i]` in the host, plus the unit direction, length, area and a `Material`.
- **Strain.** `eps = d · (Σ N_b,i u_i - Σ N_a,i u_i) / L`. The tangent is `(E_t A / L) Aᵀ d dᵀ A`.
- **`dof_mask`.** `ux` and `uy` of the host nodes only. A bar therefore never activates `rz`.
- **Binding.** **Bind both ends to the host that contains the segment midpoint.** TS binds each end point on its own
  and can pick a neighbor for an end that sits on a shared edge. This changes nothing on a **conforming** mesh,
  where neighbors share the edge's nodes, because shape functions of nodes off that edge vanish on it. But it keeps
  every segment within one element's ≤4 nodes.
  - Document the conforming-mesh precondition. Reject a midpoint that falls inside more than one host (overlapping
    elements) as ambiguous.
  - Verify the equivalence in the tests (§4, D) instead of assuming it.
- **Kinematics.** A segment measures the *average* axial strain between its two ends. Under a non-affine host field
  (a Q4 in bending), different splits of one bar give different strain distributions. Splitting is exact only for
  affine fields, which is the basis of the §4, D tests. TS has the same formulation.
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
- a T3 closed-form inverse and a Q4 Newton inverse of the isoparametric map.
  - Do **not** port `naturalCoordinatesAt` as it is. It returns after 10 iterations without checking that it
    converged or that the coordinates are finite ([postprocess.ts](../../calcs.app/lib/fem/postprocess.ts)).
  - Require finite coordinates and a physical residual `|x(ξ,η) - x|` below the tolerance.
  - Report non-convergence as an error.
- separate tolerances: a physical one scaled to element size, and a natural-coordinate one for the inside test;
- the segment midpoint is checked as well as the end points. This catches a segment whose ends are inside the mesh
  but which crosses a hole, a notch, or the gap between disconnected regions.
- an error for a point outside the mesh (`ModelError::EmbeddedBarOutsideMesh`) and for an empty or degenerate path.

It handles a bar through a mesh node and a bar along a mesh edge. Start with brute force and a bounding-box filter,
as TS does. Add a spatial index only if the benchmarks call for it.

### 2.4 Supports, springs and loads

- **Supports.** Exact fixed DOFs, not penalty.
  - TS imposes supports with a penalty (`constraintPenalty`), so its free-DOF solution carries a small
    penalty-dependent error. Ignoring the constrained entries does not remove that error.
  - Phase 0 therefore measures it: export each fixture at two penalty scales (10x apart). Parity tolerances are set
    above the observed difference, or the export adds a penalty-free mode that eliminates fixed DOFs.
  - TS reactions are already physical, `internal - external` ([femWorker.js](../../calcs.app/app/concrete/csfm/femWorker.js)),
    so Carapace's reactions compare to them directly.
- **Springs.** Each spring becomes a `ZeroLength` with `Elastic` materials to a fixed ground node, one per sprung
  node. A spring-only support (no fixed DOFs) must still be accepted.
- **Loads.** Nodal loads go in one load pattern, scaled by the load factor.

### 2.5 Procedure: static analysis, cutbacks, limit states, arc-length

1. **Load-control cutbacks (a generic core feature that is still missing).** Today only arc-length has
   transactional cutbacks. Add optional adaptive substepping to `Integrator::LoadControl`: on non-convergence,
   revert, halve, and retry down to `min_increment`; then restore the remainder.
   - This is the TS `solveIncrementWithSession` and remainder logic, and it is useful beyond CSFM.
   - **Scheduling contract.** One `Analysis::step()` accepts **one** (sub)increment, so every accepted state is
     visible to recorders and progress reports.
     - The integrator owns the pending remainder. The stage ends when the accepted load factor reaches the stage
       target, **not** after a count of `step()` calls.
     - Today the static runner decrements `steps_remaining` once per call
       ([session.rs](../wasm-bridge/src/input_v1/session.rs) `step_once`). With cutbacks it must instead finish when
       the integrator reports the target reached. `steps` stays the *nominal* increment count, and each step result
       reports nominal and accepted increments separately.
     - The same applies to the `advance(stepBudget)` budget, which counts accepted substeps.
   - Report each cutback in the step result, as the source of CalcsApp's `retryHistory`.
   - Failure after partial progress fails with `CutbacksExhausted`. The last accepted state stays committed and
     recorded, and the error carries the accepted load factor, as TS's aborted runs do.
   - Optional growth back toward the nominal step (Crisfield eq. 20) can come later.
2. **Algorithm and convergence.** Use the existing Newton, `TangentStrategy::Current`, and Carapace's `LineSearch`
   (`RegulaFalsi`/`Bisection` with `max_eta > 1`).
   - Do **not** port the Crisfield line search or secant-Newton up front. TS phase 8 found that full NR plus line
     search is the best configuration, and the line search was worth about 8% in iterations.
   - Port the Crisfield line search only if the fixtures show Carapace's line search converging worse (§4, F).
   - Production uses `ConvergenceTest::Combined`. Its frozen per-equation reference makes it better than the TS
     norm-with-`max(·,1)` floor, which phase 8 §6 calls a dimensional assumption.
     - Defaults (in N and mm): `relative_tol = 1e-6`, `force_tol = 1e-6 · max|F|` (absolute, positive as `Combined`
       requires), and `displacement_tol = None`. Tune these on the fixtures in phase 3 and record the result.
   - Parity runs use the same increment schedule on both sides and tolerances well below the parity tolerance:
     - TS: `residualForceTol = 1e-10`, `energyTol = 1e-14`.
     - Carapace: `Combined` with `relative_tol = 1e-10` and `force_tol = 1e-10 · max|F|`.

     Phase 0 checks that this is enough by solving one TS fixture at both 1e-10 and 1e-12 and confirming that the
     change is well below the §4 parity tolerance. Without that check, "tight" is not evidence.
3. **Material limit states.** CSFM's answer is a load at which a material limit is reached, not "converged at 100%".
   - Add an element query `limit_utilization() -> Option<(f64, LimitKind)>`, where values ≥ 1 mean the limit is
     reached. `eps_cu2 = -0.0035`, so signs are written out. The kinds:
     - **Concrete crushing:** `max(0, -eps_2) / limit_c`, where `limit_c = |eps_cu2|` or a user compressive strain
       limit.
     - **Concrete tension strain** (optional; part of `findCriticalLoad`'s `StopCriteria`): `max(0, eps_1) / limit_t`.
     - **Rebar rupture:**
       - tension with tension stiffening: `eps / eps_m(ft)`;
       - compression: `-eps / eps_u`;
       - bare steel: `|eps| / eps_u`.
     - **Rebar strain limit** (optional, user value).

     A `Domain` reduction returns the controlling element and point.
   - Utilization 1 is "landed on the limit". TS's `ruptured` status uses a strict `>`. The query reports both, and
     tests use states just below and just above the limit, not exactly at it.
   - Add a stop criterion `material_limit: Option<{ exact: bool }>` to `ArcStopCriteria`, and the same check after
     each accepted load-control step.
   - Exact landing reuses arc-length's existing bracketing and retake. The utilization is a scalar target, exactly
     like `DisplacementTarget`.
   - This replaces `findCriticalLoad`'s full re-solve bisection.
4. **Arc-length past the peak.** Load control cannot pass the flexural limit point. TS stops at about 97.5% on the
   span beam at 50 kip per point (phase 8 §5).
   - That figure is a historical solver outcome: the last converged step under 20 increments with halving. It is
     neither an independently located peak nor a material limit; the bars there are at about 473 MPa, between `fy`
     and `ft`, with positive hardening.
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
| 1. Materials | `CsfmConcrete`, `CsfmRebar` | §4 A–B pass. |
| 2. Embedding | `EmbeddedBar`, `embed_polyline` | §4 D passes. |
| 3. Static procedure | Load-control cutbacks with the scheduling contract. Fixtures run natively. | §4 C, E, F and H pass. Parity with TS on every fixture. |
| 4. Limit states and arc-length | `limit_utilization`, `material_limit` stop, arc-length runs | §4 G passes. |
| 5. Wire and results | Tables, recorders, decode errors, docs, PySees type mirror | §4 I passes. `boundary-smoke.ts` case added. The built wasm agrees numerically with native. PySees rebuilt and its regression checks pass. |
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
    rotated states, repeated principal strains, and both laws with and without softening and regularization;
  - rebar (strain → stress and tangent) on every TCM/POM branch, at the startup bridge, at rupture and in
    compression;
  - embedding (polyline → segments, host and `N`);
  - model runs: per-increment displacements, Gauss strains and stresses, bar forces, accepted fractions.
- **Solver-side checks** (§2.4 and §2.5):
  - the penalty sensitivity of each model fixture;
  - the convergence-tolerance sensitivity of one fixture.

  Record both in the fixture README, because they set the parity tolerances in §4.
- **Reference calculations** for §4, C, written **before** any Carapace result exists, with the tolerance fixed in
  advance. Each one is a short script plus a note under `core/tests/csfm/reference/`:
  - the CMM panel equations;
  - the panel's parameters and boundary conditions;
  - the expected quantities;
  - the tolerance, with its justification.
- **The checker.** A `core/tests` helper that, for any `PlaneMaterial` or `Material` at a strain, or any element at a
  displacement state, compares the returned tangent with a central difference of the returned stress or resistance.
  - The step scales with magnitude.
  - The checker refuses a stencil that straddles a regime or yield kink.
  - It is meant for **smooth branches**. Deliberate numerical policies (§4, A and B) are tested explicitly instead.
  - Every new law and element must pass it on its smooth branches. This is the check TS lacked until phase 7b
    (`csfm-integration.md`, "Tangent correctness").

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

**Tolerances are absolute plus relative.** Every comparison is componentwise:
`|a - b| ≤ abs + rel · max(|a|, |b|)`. A pure relative tolerance is undefined at zero stresses and displacements.
Defaults per quantity:

| Quantity | `abs` | `rel` |
|---|---|---|
| Stress | 1e-9 MPa | as stated per test |
| Strain | 1e-15 | as stated per test |
| Displacement | 1e-12 · max\|u\| | as stated per test |

A test that departs from these says why in its `///` comment.

**A. Concrete law (parity and closed form)**

- Stress matches the fixtures to `rel` 1e-12 at every grid point.
- **Smooth-branch tangent.** The Rust tangent matches the Rust central-difference checker to `rel` 1e-7. It matches
  the TS finite difference to `rel` 1e-5, which is TS's noise, **only at points whose TS stencil stays on one
  branch**.
  - The fixture export records which points those are, and excludes points within a stencil width of a regime or
    yield kink.
  - TS's stencil only protects the tension/compression regime and can still average across the `eps_c2` or
    elastic-ideally-plastic yield kink.
- **Numerical policies, tested on their own.** The origin seed is exactly `diag(E0, E0, E0/2)` at strains below
  1e-9, and it is not checked against a finite difference.
- **Repeated principal strains and objectivity.** At `ε1 = ε2` (compression, regularized tension, and zero
  tension), and just off it, the stress and `D` do not depend on the principal direction chosen. Rotating the
  strain state rotates the stress and `D` consistently.
- Closed form:
  - parabola peak `-eta_fc fck` at `eps_c2`;
  - elastic-ideally-plastic yield;
  - `kc2` at chosen `eps_1`;
  - `Ecm`, `fctm` and `eta_fc` for fck = 20, 30, 50;
  - regularization origin tangent equal to `Ecm` (phase 7a correction) and the cap `gamma fctm`.
- Regularization absent gives exactly zero tension.
- A rotated state couples normal and shear terms.
- The kc2 coupling term makes `D` nonsymmetric when softening is on and the cap is inactive. `D` is symmetric when
  softening is off, when `kc2 = 1` (`0 < eps_1 ≤ 0.2/140`), and on the compression plateau.
- Invalid `fck` and `gamma` are errors, not panics.

**B. Rebar law (parity and closed form)**

- Routing: `rho_crit = fctm / (fy - (Es/Ecm - 1) fctm)`, checked on both sides of it.
- **Round trips.** TS forward `eps_m(sigma)` → Rust inverse → sigma, to `rel` 1e-9, on each inverse branch.
  - This **excludes** TCM stresses in `[0, sigma_min]`. Those all map forward to strain 0, so the inverse is not
    unique there.
  - It also excludes the startup bridge, which deliberately replaces the inverse below `sigma_min / Es`.
- **The inverse's own residual** is below `rel` 1e-12 on every branch. This checks the stable quadratic form.
- **TCM startup.** Where it is admitted (`2 sigma_min ≤ fy`), it is C1 at `sigma_min / Es`, the stress at the
  startup strain is `2 sigma_min`, and the stress just above zero carries no `sigma_min` jump.
  `2 sigma_min > fy` is a construction error. The test uses TS's discontinuous case: `fck = 30`, `fy = 500`,
  `ft = 550`, `rho = 0.01`, `lambda = 2`.
- **Smooth-branch tangent.** On each branch, away from the clamp, the tangent equals `1/slope` and passes the
  finite-difference checker.
- **Numerical policies, tested on their own:**
  - the clamp to `[Es/1e4, 10 Es]`, where it differs from the derivative of the unchanged stress curve;
  - the POM origin guard (stress 0, tangent `Es`);
  - rupture: stress `±ft`, positive secant tangent, and `ruptured()` true. Rupture is checked just past each
    threshold in §1, including compression and bare tension at `eps_u`.
- Compression is bilinear. Bare mode is bilinear in both directions.
- `rho_eff` out of `(0, 1)` and `ft ≤ fy` are errors.
- Parity with the TS law at every fixture strain to `rel` 1e-12. The cancellation-free roots may differ from TS by
  rounding; any larger difference is a bug.

**C. Independent validation (phase 3: needs embedding)**

These check the method, not the port. Sources are in `../calcs.app/docs/concrete/csfm/`. Each case runs against the
reference calculation written in phase 0. Its equations, parameters, boundary conditions, expected quantities and
tolerance are fixed before any Carapace run.

- **Cracked Membrane Model panel** (Kaufmann and Marti 1998). An orthogonally reinforced panel, uniform and
  single-element, under pure shear and under one biaxial case.
  - Displacements are prescribed so the strain is uniform. The comparison is stress and bar stress against the CMM
    equilibrium and compatibility solution at the same strain.
  - The CMM and this CSFM law differ by design: CSFM has no concrete tension, and its crack spacing and kc2 are
    simplified. So the reference note states which quantities must match closely (equilibrium identities, to
    rounding) and which only within a stated band (the stresses, with the band derived from those differences).
- **Tension chord.** A bar in a uniform concrete prism, pulled axially. The average bar stress against the average
  strain must reproduce the TCM relation it embeds, to `rel` 1e-8, for both TCM and POM routing.

**D. Embedded bar element and embedding**

- Patch tests on T3, Q4 and a distorted Q4 under uniform strain: bar strain equals `dᵀ ε d` of the host for any
  position and angle, to 1e-12.
- Rigid translation gives zero force.
- A bar coincident with a mesh line matches a `Truss` on shared nodes exactly.
- **Split invariance, affine fields only.** Under an affine host field, a bar crossing many elements gives the same
  forces however it is split (cut at extra points). Under a non-affine Q4 field (`ux = xy`), splitting legitimately
  changes the segment strains: a diagonal at 0.5 becomes 0.25 and 0.75. That case is a convergence-under-subdivision
  test instead.
- Embedding errors:
  - a segment crossing a hole or notch, with both ends inside;
  - a non-converging Q4 inverse on a badly distorted element;
  - overlapping hosts.
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
- Linear-elastic displacements match TS `solveLinearElastic` on the deep beam, within the penalty sensitivity
  measured in phase 0 (or exactly, against the penalty-free export).

**F. Nonlinear model runs (parity with `csfm-nonlinear.test.ts`)**

Each fixture runs with the parity settings of §2.5 and the same increment schedule as TS. It must reach the TS
outcome, and must match TS's converged state at each increment:

- displacements to `rel` 1e-6;
- bar stress to `rel` 1e-5.

Both tolerances are confirmed in phase 0 to exceed the penalty and convergence-tolerance sensitivity. Near a limit
point (a nearly singular tangent), parity applies only to increments at or below 90% of the last TS-accepted
fraction. The cases:

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
- **Load-control cutbacks.** A step forced to fail retries at half size, restores the remainder, and reaches the
  same final state as a fine schedule. Exhausted cutbacks keep the last accepted state committed and report its load
  factor.
- Solver comparison, recorded and not asserted: iterations per increment for Carapace's line searches against the TS
  counts in phase 8 §4. If Carapace needs over 1.5x TS's iterations on the span beam at 30 kip, port the Crisfield
  slope line search.

**G. Limit states and arc-length**

These are three separate tests. Passing a peak is not reaching a material limit.

- **The peak.** The span beam at 50 kip/point runs under arc-length. The test asserts that `lambda` rises, reaches a
  maximum and then falls, so the peak is passed. The peak load is located by refining the radius around the maximum
  until it changes by less than 0.5% between refinements, and that refined value is recorded.
  - The only assertion against TS: the load-control run stops at or below the peak, consistent with TS's last
    accepted step (about 97.5%).
  - The 97.5% figure is not treated as the peak (§2.5).
- **Material limits with exact landing.** `material_limit` with `exact` lands on utilization 1 within the landing
  tolerance. One model per limit kind, each built so that kind controls:
  - a lightly reinforced tie (tensile rupture);
  - an over-reinforced strut (crushing);
  - a bare-steel case (`|eps| = eps_u`).

  Utilization is evaluated just below and just above each landing, to distinguish landing from TS's strict-`>`
  rupture.
- **`findCriticalLoad` parity.** The bracketing load agrees with TS on the same fixture and stop criteria, including
  the concrete tensile-strain criterion, to its `loadFractionTolerance`.

**H. Solver**

- Dump a cracked CSFM tangent from the span beam as a fixture. faer's LU solves it with relative residual < 1e-10.
  This is the phase 8 lesson: assert the linear-solve residual literally.
- A nonsymmetric CSFM tangent goes to LU and not Cholesky.

**I. Wire**

- Native vs wire, identical results for the span beam, including recorders.
- **Cutback scheduling over the wire.** A nominal 20-step stage that needs cutbacks reaches its target load factor
  and records every accepted substep. The `advance` budget counts accepted substeps.
- **Native vs built wasm.** The span beam runs through the real `pkg/` build under Node (`comparison/load_wasm.mjs`)
  and matches native to rounding.
- **PySees.** After the wire change, rebuild and copy the bundle (`../pysees/scripts/build-carapace.sh`). PySees'
  type check and `carapaceVsOpenSees.test.ts` must pass.
- Each new `DecodeError` case:
  - invalid `fck` or `gamma`;
  - invalid rebar, including `2 sigma_min > fy`;
  - a bar outside the mesh or crossing a hole;
  - a bar on a 3D profile;
  - `csfmConcrete` with the enhanced quad.
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

## 6. Decisions

Items 1, 2 and 5 were decided by the owner on 2026-10-25. The rest are standing scope notes.

1. **Regularization default (decided).** Carapace has **no** default: `tensileRegularization`, when present, must
   give `gamma`. Phase 7b showed that gamma = 0.1 lets the concrete carry the span-beam tie (bar 6.3 vs about
   67 MPa), so the CalcsApp UI default is now `gamma = 0.001`.
2. **Rebar concrete inputs (decided).** `csfmRebar` takes `fctm` and `Ecm` explicitly in core, which keeps `Material`
   self-contained. The wire decoder may later resolve them from a referenced `csfmConcrete`, as a convenience.
3. **Fracture-energy tension softening** (phase 7b, "Optional future change"). Out of scope. It has a negative
   tangent and would need arc-length, which phase 4 provides.
4. **Crisfield line search and secant-Newton.** Port them only on evidence from §4, F.
5. **Compression failure criterion (decided).** TS has none in the worker path; `findCriticalLoad` exists but is
   unused by the UI. The default crushing limit is `eps_cu2 = -0.0035` for **both** compression laws, and a user
   limit can override it.
6. **TS defects found while planning (fixed in CalcsApp).** Carapace matches the fixed behavior, and the phase 0
   export runs after these fixes.
   - **Crack-transition band:** removed from CalcsApp, and passing it is an error (§1). Carapace does not
     implement it. If it is ever wanted, it must apply to both principal directions.
   - **TCM startup bridge:** both reject `2 sigma_min > fy` (§1). CalcsApp's `csfm-nonlinear.test.ts` covers three
     cases:
     - the rejected `lambda = 2` case;
     - the UI's worst case, `lambda = 1` just above `rho_crit`, which stays continuous;
     - equal principal strains giving equal stresses in every regime.
7. **Bond-slip.** Out of scope. `rod.ts` keeps a `RodCoupling` seam for it. Here, `EmbeddedBar` would gain slip DOFs,
   which needs nodes it owns and is a later element.
