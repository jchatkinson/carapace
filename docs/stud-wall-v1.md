# Stud wall v1: implementation spec

Status: implementation spec. It supersedes the parts of `docs/stud-wall-panel.md` it contradicts. That
document remains the background (prior work, rationale, RAM SS lessons) and the list of deferred features.
How to use the element (layout, laws, supports, analysis, results, assumptions, limitations) is in the user
guide, [`docs/stud-wall.md`](stud-wall.md).

v1 is a 2D macro-element whose frame remains linear elastic and whose per-pier sheathing can be linear or
nonlinear. Its purpose is design:

- gravity takedowns that match hand calculations;
- lateral analysis that matches segmented shear-wall design (SDPWS 2021 Eq. 4.3-1; CSA O86:19 11.7.1.2 and
  Annex A.11.7.1 with the `O86` law).

Hold-downs and anchorage are ordinary exterior elements, and the framing is linear elastic; only the per-pier
sheathing law can be nonlinear.

## 1. What v1 includes and defers

| In v1 | Deferred |
|---|---|
| Studs and stud packs (one row each), top chord `beam` or `direct`, bottom chord, openings with header, jacks and cripples (windows and doors) | Sills and cripples below windows, force transfer around openings (segmented design only) |
| Stations (exterior nodes) on either chord | Automatic station generation (a PySees/decoder concern) |
| Linear sheathing (`ga`), CSA O86:19 nonlinear-elastic shear and nail-slip law, and uniaxial material laws, each one constant-shear field per full-height pier | Hysteretic member framing, staged installation |
| Top point and top line loads in any number of wall load cases (selected by patterns) | Self-weight, mass, out-of-plane output (PySees post-processing) |
| Core API and core tests | Wire format, CLI, wasm and PySees mirror (a later phase) |

## 2. Frame and exterior nodes

- **Corners.** Four corner nodes, counter-clockwise from bottom-left: `c0` (local `(0, 0)`), `c1` (`(L, 0)`),
  `c2` (`(L, H)`) and `c3` (`(0, H)`).
  - Local `x` runs along `c0 → c1`, and local `y` is `x` rotated +90°.
  - The rectangle may have any orientation in the global plane.
  - `validate` checks that the corner coordinates form that rectangle, with lengths `L` and `H` matching the
    layout to a relative 1e-9.
- **Stations.** Extra exterior nodes, each given as `(chord: Top | Bottom, x, node)`, whose node sits at that
  chord point.
- **Exterior DOFs are translations only.** `dof_mask` is `ux, uy` for every exterior node. The element never
  connects a node's rotation. Plate and rim splices are not moment connections, so chords are pinned to the
  outside world; chord rotations at corners and stations stay internal. Exterior node order is corners `c0..c3`,
  then stations in input order.

## 3. Layout (input)

All quantities are in consistent units, and the layout is defined in the local frame.

```text
StudWallLayout {
  length L, height H,
  studs: [{ x, ea }],                 // sorted by x; first at 0, last at L; a pack is one row (summed EA)
  top: { ea, ei, distribution: Beam | Direct },   // ei ignored for Direct
  bottom: { ea, ei },
  openings: [{ x0, x1, head_height, header_ei, jack_ea, door: bool }],
  sheathing: Option<{ law }>,         // Linear { ga }, O86 { bv, nail_spacing, nail_diameter }, or Material(Material)
  foundation: bool,                   // default false; internally fix interior bottom-chord u and v
}
```

**Validation.** Each of these errors names the offending item:
- `L`, `H` and every `ea` and `ei` are positive and finite (`top.ei` only with `Beam`).
- Studs are strictly increasing in `x`, start at 0 and end at `L`.
- For each opening: `0 < x0 < x1 < L` and `0 < head_height < H`. Openings don't overlap, and at least one
  full-height stud lies in `[x1 of the previous opening, x0]`.
- A stud row exactly at an opening's `x0` or `x1` is the king stud. It is a full-height stud that shares the
  jack's bottom-chord node (one end post), so the adjacent pier ends at the opening edge.
- With `Direct`, a top station's `x` must coincide with a stud or cripple. A hinged chord node needs a vertical
  member below it.
- A bottom station may not lie strictly inside a door.
- With `foundation: true`, every interior bottom-chord position (stud or jack) is internally fixed in `u` (sill
  anchorage) and `v` (bearing), except where a bottom station is placed: the station replaces the internal
  support (for a hold-down at a pier end). Bottom corners remain exterior nodes.
- A singular internal stiffness is reported as an error at construction, never regularized. It can't occur
  for a valid layout; it guards against layout bugs.

## 4. Internal model

### 4.1 Nodes and DOFs

All DOFs are in the local frame:

- **Top chord nodes** sit at every full-height stud, every cripple and every top station.
  - With `Beam`, the nodes also include each opening's `x0` and `x1`.
  - DOFs are `u, v, θ` with `Beam`, and `u, v` with `Direct`.
- **Bottom chord nodes** sit at every full-height stud, at each opening's `x0` and `x1`, and at every bottom
  station. DOFs are `u, v, θ`.
- **Header nodes** sit at `x0`, at each cripple `x` and at `x1`, at height `head_height`. DOFs are `v, θ` only;
  the header has no axial DOF, since v1 has no lateral load path through openings.

### 4.2 Members

| Member | Model |
|---|---|
| Chord segment | Between consecutive chord nodes. `Beam`: Euler–Bernoulli (`ea` on `u`, `ei` on `v, θ`). `Direct` top chord: an axial bar on `u` only, hinged at every node. |
| Bottom chord at a door | No segment between `x0` and `x1`. |
| Chord segments over an opening | Beam segments lying within `[x0, x1]` (top `Beam` chord, bottom chord under a window) are hinged at `x0` and `x1`: the segment end gets its own rotation DOF, and loads on the segment put their fixed-end moment on it. Piers on either side are tied axially only. |
| Full-height stud | Vertical bar between the bottom and top chord nodes at its `x`: `ea/H` on the two `v`s. |
| Jack | Vertical bar from the bottom chord node at `x0` (or `x1`) to the header end node: `jack_ea/head_height`. |
| Cripple | A stud row strictly inside `(x0, x1)` becomes a vertical bar from the header node to the top chord node: `ea/(H − head_height)`. |
| Header | Euler–Bernoulli beam on `v, θ` between consecutive header nodes. Header ends are simply supported on the jacks, bilaterally. |

Below an opening there are no members; window sills are deferred.

**Stability.** Every vertical bar couples only `v`. Chord `u` is stiffened by chord axial stiffness, chained to
the corners. `v` and `θ` are stiffened by studs, jacks, cripples and bending. So `K_ii` is non-singular, and an
unsheathed wall has one exterior zero-energy mode, racking. The model must restrain it, for example with `ux`
fixes or a lateral system. There is no `rackingStiffness`.

### 4.3 Sheathing law: one shear field per pier

**Piers.** A pier is a maximal run of full-height studs not interrupted by an opening. Its extent is
`[xa, xb]`, from its first to its last full-height stud, so `b = xb − xa`. A pier with fewer than two studs
has no field. Each pier has one constant shear flow `q_p` over `[xa, xb] × [0, H]`.

**Compatible strain.** `u_top(x)` and `u_bot(x)` are interpolated linearly between chord nodes, and
`v̄(x) = (v_top(x) + v_bot(x))/2` is taken at a boundary stud:

```text
γ_p = [∫_xa^xb u_top dx − ∫_xa^xb u_bot dx] / (b H) + (v̄(xb) − v̄(xa)) / b  =  β_pᵀ u
```

The integrals are trapezoid sums over every chord segment in `[xa, xb]`, including those created by stations.

**Lag compliance.** Each edge member's force varies linearly at slope `q`. Writing `ℓ_k` for each chord
segment's length inside the pier:

```text
c_p = H³/12 · (1/EA(xa) + 1/EA(xb)) + Σ_top ℓ_k³/(12 EA_top) + Σ_bottom ℓ_k³/(12 EA_bot)
```

Interior studs of a pier get equal and opposite traction from the two sides, so they have no lag term.

**Law and series lag.** With `A_p = b H`, `q_p` is shear flow and `γ_sh,p` is the sheathing strain. The
frame's lag compliance is in series with the sheathing law:

```text
β_pᵀu = γ_sh,p + c_p q_p/A_p
```

For `Linear { ga }`, `q = ga γ_sh`, giving the existing closed form
`q = A βᵀu/(A/ga + c)` and rank-one stiffness `A² β βᵀ/(A/ga+c)`. This path retains the original
one-time condensation. `O86` uses `γ_sh = q/Bv + a q |q|`, with
`a = 0.0025·(0.013·nail_spacing/nail_diameter²)²`; its inverse is
`q = 2γ_sh/(1/Bv + sqrt(1/Bv² + 4a|γ_sh|))`. `Material` evaluates `q = σ(γ_sh)` from each pier's
committed uniaxial material state and commits only after a converged analysis step. A singular local Jacobian
or failed local solve has no element-failure channel; the implementation propagates NaNs into the element
contribution.

**Nonlinear laws: local pier solve.** The frame is condensed without sheathing (`K_ii` is non-singular
without it, see 4.2). Write `B_i` for the internal rows of the pier strain vectors `β_p` (free internal DOFs
only) and `A = diag(A_p)`. Eliminating the internal DOFs leaves a P×P system in the pier sheathing strains:

```text
γ_sh + S q(γ_sh) = y,     S = B_iᵀ K_ii⁻¹ B_i A + diag(c_p/A_p),     y = G u_e + Σ_k λ_k y_k
```

- `y` is the pier strain of the unsheathed frame: `G = B_eᵀ − B_iᵀ K_ii⁻¹ K_ie` (the exterior map) and
  `y_k = B_iᵀ K_ii⁻¹ p_i,k` (load case `k`; `λ_k` are the effective case factors, the sum over active patterns).
- A local Newton solve (fixed-size buffers, hand-written LU with partial pivoting, residual `≤ 1e-14 (1 + max|y|)`,
  at most 40 iterations) starts from `γ_sh = y` and evaluates each law and its tangent `D = diag(dq/dγ_sh)`.
  `Material` laws are evaluated as trials from their committed state.
- Exterior resisting force: `K_c,frame u_e + C q`, with `C = G ᵀ A`, the condensed force map of the pier
  fields. The case loads are applied by `assemble_load` as in the linear path.
- Consistent tangent: `K_e = K_c,frame + C D (I + S D)⁻¹ G` (implicit derivative of the local system; `S`
  is not symmetric, so the order matters).
- Internal recovery: `u_i = Σ_k λ_k w_k − R u_e − K_ii⁻¹ B_i A q`.

Nonlinear walls are capped at 16 piers and, like linear walls, 256 exterior nodes. They need an iterative
algorithm such as Newton; `Algorithm::Linear` performs one tangent solve and does not iterate to equilibrium.
The dense nonlinear tangent is built in a per-thread buffer that is reused between calls (it allocates only
when it first grows), plus fixed local solve buffers for a 16×16 Jacobian (2 KiB) and its vectors. `prepare`
caches rotated frame and pier maps, so prepared assembly does not otherwise allocate. Preparation and unprepared calls may allocate, as in the linear path.

**Member force recovery.** Each member's mean force comes from its own stiffness, and the field adds a linear
part. Tension is positive; `q_L`, `q_R` are the fields left and right of a stud, with 0 where there is none:

- Stud: `N(y) = N̄ + (q_L − q_R)(y − H/2)`.
- Top chord segment, midpoint `x_m`: `N(x) = N̄ + q (x − x_m)`.
- Bottom chord segment: `N(x) = N̄ − q (x − x_m)`.

**SDPWS check.** For one pier with rigid chords, end posts `EA` and rigid anchorage, the top deflection
under `V` is exactly `δ = 2vH³/(3EA b) + vH/ga` with `v = V/b`. These are the first two terms of SDPWS 2021
Eq. 4.3-1 in consistent units; hold-down springs outside the element give the third.

### 4.4 Condensation

For linear sheathing, condensation is built once, in the constructor, in the local frame. It depends only on the layout and the station list;
orientation comes from the corner coordinates at call time. The linear path retains its 256-exterior-node limit.

- `K_c = K_ee − K_ei K_ii⁻¹ K_ie` (dense is fine; walls have at most a few hundred internal DOFs).
- `R = K_ii⁻¹ K_ie`, stored for recovery.
- Per load case `k`: `f_k = p_e,k − K_ei K_ii⁻¹ p_i,k` and `w_k = K_ii⁻¹ p_i,k`.
- Recovery: `u_i = Σ_k λ_k w_k − R u_e`, where `λ_k` are the effective case factors.

Every element method gives the same answer whether or not `prepare` has run.

For nonlinear sheathing, the constructor instead condenses the bare frame and stores the pier coupling and
compliance maps. The pier Newton solve and consistent correction run at each tangent evaluation.

## 5. Loads: pattern-owned cases

- **Case field.** `ElementLoad` gains `wall: Option<WallCaseLoad { case, factor }>`, scalable like the other
  fields (the same case adds factors; different cases on one element in one pattern are `WallCaseLoad::MIXED`, which
  `Domain::validate` rejects). Other elements reject a `wall`; the stud wall rejects every other field and a case
  index it does not have. Across patterns the wall receives the list of active `(case, factor)` pairs.
- **Case contents.** The stud wall holds, per case, a list of `WallLoad`:
  - `TopPoint { x, fx, fy }`;
  - `TopLine { x0, x1, wx0, wx1, wy0, wy1 }`, linear along `x`, in local axes.

  A pattern's loads on a wall go in one case. The pattern then carries `wall = Some(case)`, and its factor
  scales the case.
- **Discretization** into internal and exterior consistent loads, per top chord segment:
  - `Beam`: fixed-end forces of the beam segment for points and trapezoid pieces.
  - `Direct`: simple-span (lever-rule) reactions to the two end nodes, for both `fx` and `fy`.
  - A load exactly at a node goes to that node.
- **Assembly.**
  - `assemble_load` emits `Σ_k wall[k] · f_k` on the exterior DOFs.
  - `assemble_tangent` emits `K_c` and resisting force `K_c u_e`. It never subtracts the load.
  - `local_load_force` returns `Σ_k wall[k] · f_k` in local axes. For nonlinear walls, `local_force` reads the
    committed pier flows cached by `commit` and returns `K_c u_e + C q`; `element_end_force` subtracts the slot
    load once, including its influence on the committed pier state.

## 6. Responses (core API)

`StudWall::response(&self, nodes, load: Option<&ElementLoad>) -> StudWallResponse` uses the current nodal
state (`Domain::element_stud_wall_response(id, pseudo_time)` passes the effective slot factors). For nonlinear
laws it runs the local pier solve from the committed law state, so after a converged, committed step it gives
the committed response. It reports, in local axes with tension positive:

- **Studs** (full-height): axial force at the bottom and top.
- **Jacks and cripples:** axial force.
- **Header:** end reactions (the jack forces), plus moment and shear at each header node.
- **Chords:** axial force, moment and shear at each chord node, per chord.
- **Piers:** `xa`, `xb`, `q`, `γ`.
- **Pier overturning:** `moment_base` and `moment_top` are `Σ (x − x_mid) · N` over the full-height studs in the
  pier, using each stud's axial force (tension positive) at its bottom and top end, CCW positive. By statics
  this is the external overturning moment of everything above that section: a +x lateral load at the top gives
  negative values, and `moment_top = moment_base + q·b·H`. `jd = xb − xa`; `chords = [−moment_base/jd,
  +moment_base/jd]` (tension positive, so +x gives tension in the left post) is the overturning-only CSA O86 /
  SDPWS M/jd design force. It is exact however the element shares overturning among studs; the gravity share
  of a post comes from its stud force or a gravity analysis.
- **Base** (bottom-chord exterior nodes: `c0`, `c1` and bottom stations): sums of the element's resisting
  forces, `Fx` and `Fy`, and the moment about `(L/2, 0)`, `Σ (x_i − L/2) · Fy_i`.

## 7. Engine infrastructure (prerequisite)

- **More nodes per element:** more than `MAX_ELEMENT_NODES = 4` nodes.
- **Run-time-sized matrices:** element matrices whose size is not a compile-time constant. `TangentSink` and
  `VectorSink` get slice-based paths.
- **Constraints:**
  - No behaviour change for existing elements, and bit-identical results.
  - No heap allocation in the assembly loop for existing elements.
  - Serial and parallel assembly stay bit-identical, scattering in element order and in a fixed order
    within each element.

## 8. Verification (v1)

| Test | Reference | Tolerance |
|---|---|---|
| Mask conformance. Rigid-body translations and rotation give zero force. Rank with corners only (8 DOFs): 5 sheathed, 4 unsheathed. Rotated wall equals unrotated. | Kinematics | 1e-10 relative |
| Unsheathed, `Beam` and `Direct`, with openings, stations and doors. Condensed element vs the same frame built explicitly from `Truss` and `ElasticBeamColumn` (hinges as coincident nodes plus `equalDOF`). Displacements, reactions, member forces. | Explicit Carapace model | 1e-9 relative |
| Loaded condensation. Points between nodes, partial trapezoids, three patterns with different factors in different slots. Exterior forces, reactions and recovered member forces. | Explicit model | 1e-9 |
| One pier, rigid chords (`ea × 1e6`), end posts `EA`: top deflection and post forces. | `2vH³/(3EAb) + vH/ga`; post force `±qH` at the base, 0 at the top | 1e-6 |
| One pier with O86 nonlinear-elastic law, solved by Newton. | CSA O86:19 11.7.1.2 shear and nail-slip terms | 1e-6 |
| Elastic `Material` versus `Linear` on window and door walls with stations and slot loads. | Linear-law condensation | 1e-8 (the two condensation paths differ by roundoff, about 1e-9) |
| Slot-loaded nonlinear `element_end_force` after commit. | Response exterior resistance | 1e-10 of applied load |
| Two-pier door wall with slot load, O86 and hysteretic material. | Central-difference tangent | 1e-6 relative (2e-6 absolute floor for near-zero terms) |
| Cyclic Steel01 sheathing under displacement control and trial reversion. | Standalone Steel01 at matching shear strains | 1e-5 relative (rigid stand-in conditioning) |
| Two piers with a door and different `ga`; an unloaded station inserted inside a pier leaves the response unchanged. | Consistency | 1e-10 |
| Foundation supports versus fixed bottom stations, including window and door walls; per-position reactions and member responses. | Equivalent restrained frame | 1e-10 relative |
| Pier base/top stud-force moments and M/jd chord forces, 2- and 5-storey walls in every variant. | Section statics; O86 A.11.7.1 hold-down force | 1e-6 relative (rigid stand-in conditioning) |
| O86 law per storey against the matching secant walls, 2- and 5-storey, every variant (Newton, displacement-increment test). | Secant-at-design-load identity; code idealized also against the O86 A.11.7.1 hand calculation | 1e-8 with real plates; 1e-6 with rigid stand-ins (`EA × 1e6` plates, `RIGID` links); 1e-5 hand calculation |
| Hand-calc takedowns (`stud_wall_takedown.rs`): uniform, point, packs, window and door headers, factored combinations, stacked storeys, wall on a beam. | Hand calculations | 1e-9 (1e-6 for the continuous-beam top plate) |
| Window wall with king studs and hold-downs at every pier end (`stud_wall_openings.rs`), `Direct` and rim-stiff `Beam` top chords. | Two separate pier walls joined by an axial drag strut (segmented design) | 1e-10 |
| User-guide example (`stud_wall_guide_example.rs`): window wall with king studs on a foundation, O86 law, hold-downs at every pier end, gravity plus wind with Newton. | Statics: base reactions, pier shear sum, M/jd | 1e-8 |
