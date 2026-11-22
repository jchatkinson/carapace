# Stud wall v1: implementation spec

Status: implementation target. It supersedes the parts of `docs/stud-wall-panel.md` it contradicts. That
document remains the background (prior work, rationale, RAM SS lessons) and the list of deferred features.

v1 is a **linear** 2D macro-element whose purpose is design:

- gravity takedowns that match hand calculations;
- lateral analysis that matches segmented shear-wall design (SDPWS 2021 Eq. 4.3-1, three-term form).

Everything nonlinear (hold-downs, bearing, slip) lives in ordinary exterior elements (`zeroLength`).

## 1. What v1 includes and defers

| In v1 | Deferred |
|---|---|
| Studs and stud packs (one row each), top chord `beam` or `direct`, bottom chord, openings with header, jacks and cripples (windows and doors) | Sills and cripples below windows, `chordOverOpening: released`, force transfer around openings |
| Stations (exterior nodes) on either chord | Automatic station generation (a PySees/decoder concern) |
| Linear sheathing: one constant-shear field per full-height pier, stiffness `ga` | Code nonlinear-elastic laws (`sdpws4`, `o86`), hysteretic sheathing, staged installation |
| Top point and top line loads in up to 8 pattern slots | Self-weight, mass, out-of-plane output (PySees post-processing) |
| Condensation once at construction; internal recovery for responses | Internal Newton, failure channel, scratch arena |
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
  sheathing: Option<{ ga }>,          // shear flow q = ga * gamma (force/length), both faces summed
  foundation: bool,                   // default false; internally fix interior bottom-chord u and v
}
```

**Validation.** Each of these errors names the offending item:
- `L`, `H` and every `ea` and `ei` are positive and finite (`top.ei` only with `Beam`).
- Studs are strictly increasing in `x`, start at 0 and end at `L`.
- For each opening: `0 < x0 < x1 < L` and `0 < head_height < H`. Openings don't overlap, and at least one
  full-height stud lies between consecutive openings.
- No stud row sits exactly at an opening's `x0` or `x1`, since the jacks occupy those positions.
- With `Direct`, a top station's `x` must coincide with a stud or cripple. A hinged chord node needs a vertical
  member below it.
- A bottom station may not lie strictly inside a door.
- With `foundation: true`, bottom stations are rejected; every interior bottom-chord position (stud or jack)
  is internally fixed in `u` (sill anchorage) and `v` (bearing). Bottom corners remain exterior nodes.
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
| Full-height stud | Vertical bar between the bottom and top chord nodes at its `x`: `ea/H` on the two `v`s. |
| Jack | Vertical bar from the bottom chord node at `x0` (or `x1`) to the header end node: `jack_ea/head_height`. |
| Cripple | A stud row strictly inside `(x0, x1)` becomes a vertical bar from the header node to the top chord node: `ea/(H − head_height)`. |
| Header | Euler–Bernoulli beam on `v, θ` between consecutive header nodes. Header ends are simply supported on the jacks, bilaterally. |

Below an opening there are no members; window sills are deferred.

**Stability.** Every vertical bar couples only `v`. Chord `u` is stiffened by chord axial stiffness, chained to
the corners. `v` and `θ` are stiffened by studs, jacks, cripples and bending. So `K_ii` is non-singular, and an
unsheathed wall has one exterior zero-energy mode, racking. The model must restrain it, for example with `ux`
fixes or a lateral system. There is no `rackingStiffness`.

### 4.3 Sheathing: one shear field per pier

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

**Stiffness.** With `A_p = b H`, the term is stationary in `q_p`:

```text
q_p = A_p β_pᵀ u / (A_p/ga + c_p),     K += A_p² β_p β_pᵀ / (A_p/ga + c_p)
```

**Member force recovery.** Each member's mean force comes from its own stiffness, and the field adds a linear
part. Tension is positive; `q_L`, `q_R` are the fields left and right of a stud, with 0 where there is none:

- Stud: `N(y) = N̄ + (q_L − q_R)(y − H/2)`.
- Top chord segment, midpoint `x_m`: `N(x) = N̄ + q (x − x_m)`.
- Bottom chord segment: `N(x) = N̄ − q (x − x_m)`.

**SDPWS check.** For one pier with rigid chords, end posts `EA` and rigid anchorage, the top deflection
under `V` is exactly `δ = 2vH³/(3EA b) + vH/ga` with `v = V/b`. These are the first two terms of SDPWS 2021
Eq. 4.3-1 in consistent units; hold-down springs outside the element give the third.

### 4.4 Condensation

Built once, in the constructor, in the local frame. It depends only on the layout and the station list;
orientation comes from the corner coordinates at call time.

- `K_c = K_ee − K_ei K_ii⁻¹ K_ie` (dense is fine; walls have at most a few hundred internal DOFs).
- `R = K_ii⁻¹ K_ie`, stored for recovery.
- Per load slot `k`: `f_k = p_e,k − K_ei K_ii⁻¹ p_i,k` and `w_k = K_ii⁻¹ p_i,k`.
- Recovery: `u_i = Σ_k λ_k w_k − R u_e`, where `λ_k` are the effective slot factors.

Every element method gives the same answer whether or not `prepare` has run.

## 5. Loads: pattern slots

- **Slot field.** `ElementLoad` gains `wall: [f64; 8]` (`WALL_SLOTS = 8`). It is additive and scalable like the
  other fields. Other elements reject a nonzero `wall`; the stud wall rejects every other field.
- **Slot contents.** The stud wall holds, per slot, a list of `WallLoad`:
  - `TopPoint { x, fx, fy }`;
  - `TopLine { x0, x1, wx0, wx1, wy0, wy1 }`, linear along `x`, in local axes.

  A pattern's loads on a wall go in one slot. The pattern then carries `wall = e_slot`, and its factor scales
  the slot.
- **Discretization** into internal and exterior consistent loads, per top chord segment:
  - `Beam`: fixed-end forces of the beam segment for points and trapezoid pieces.
  - `Direct`: simple-span (lever-rule) reactions to the two end nodes, for both `fx` and `fy`.
  - A load exactly at a node goes to that node.
- **Assembly.**
  - `assemble_load` emits `Σ_k wall[k] · f_k` on the exterior DOFs.
  - `assemble_tangent` emits `K_c` and resisting force `K_c u_e`. It never subtracts the load.
  - `local_load_force` returns `Σ_k wall[k] · f_k` in local axes, and `local_force` returns `K_c u_e` in local
    axes, so `element_end_force` behaves as for a loaded beam.

## 6. Responses (core API)

`StudWall::response(&self, nodes, load: Option<&ElementLoad>) -> StudWallResponse` uses the committed nodal
state. It reports, in local axes with tension positive:

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
| Two piers with a door and different `ga`; an unloaded station inserted inside a pier leaves the response unchanged. | Consistency | 1e-10 |
| Foundation supports versus fixed bottom stations, including window and door walls; per-position reactions and member responses. | Equivalent restrained frame | 1e-10 relative |
| Pier base/top stud-force moments and M/jd chord forces, 2- and 5-storey walls in every variant. | Section statics; O86 A.11.7.1 hold-down force | 1e-6 relative (rigid stand-in conditioning) |
| Hand-calc takedowns (separate test file, written after the element). | Hand calculations | see file |
