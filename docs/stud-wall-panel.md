# Stud-wall macro-element

Status: background and long-term design. **What is implemented is specified in `docs/stud-wall-v1.md`**
(linear elastic framing, per-pier sheathing laws: linear, CSA O86 nail slip, or any uniaxial material), which
supersedes this document wherever they differ; how to use it is in the user guide, `docs/stud-wall.md`. This
document keeps the rationale and the features still deferred (nonlinear framing,
out-of-plane output, wire format). Revised after review (Codex; Bentley's RAM Structural System walls FAQ,
KB0112985).

A 2D element for light-frame walls (wood or cold-formed steel) that covers one wall segment for one storey and
is defined by its framing layout: studs and stud packs, top and bottom chords (plates, rim), openings with
headers, and optional sheathing. Internally it is a small frame of studs, chord beams, headers and sheathing
shear fields, condensed onto a few exterior nodes, so users model walls, not studs.

The same element serves:

- **Gravity takedown**, without relying on sheathing: non-uniform loads, stud packs, rim vs blocking, headers.
- **Lateral (in-plane) analysis**, with sheathing: reproduces the bending and shear terms of the SDPWS and CSA
  O86 deflection equations in the elastic case, with any uniaxial hysteretic material for the sheathing.
- **Out-of-plane wind** as per-stud output, with no extra DOFs.

## 1. Scope

In v1:

- 2D (ndm=2), rectangular wall segments, small displacements, no in-plane second-order effects.
- Studs are pinned at the chords; their in-plane bending is ignored.
- **Sheathed walls have elastic framing.** In an element with sheathing (installed from the start or at a
  stage), studs, jacks, cripples, chords and headers must be linear elastic. Wall nonlinearity comes from the
  sheathing material and from exterior hold-down and anchorage springs. Nonlinear framing materials are
  allowed only in elements without sheathing (Section 5.4 explains why).
- Openings carry gravity through headers and jacks. Laterally, only full-height sheathed bays get shear fields
  (Section 5.5).

Out of scope for v1, listed in Section 14: force transfer around openings, 3D, non-rectangular walls,
perpendicular-to-grain bearing, shrinkage, sills, and CFS stud and track stability.

## 2. Prior work

| Work | What it does | Relation |
|---|---|---|
| Kasal, Leichti & Itani (1994) | Condensed wall "superelements" in whole-house models | Same idea: a condensed wall element. This design adds gravity load paths and staged sheathing. |
| Garvey (1951); NASTRAN `CSHEAR` + `CROD` | Shear panels with constant shear flow; edge members whose force varies linearly | The per-bay sheathing field (Section 5.4). |
| Pian (1964), assumed-stress hybrid elements | Stiffness from an equilibrium stress field | How sheathing and framing are coupled. |
| Folz & Filiatrault (2001, 2004), CASHEW / SAWS | Nail-level wall model condensed to one hysteretic spring | Source of the sheathing backbone. |
| Judd & Fonseca (2005); Pang & Rosowsky; Christovasilis & Filiatrault (Timber3D) | Nail-spring, beam-spring and condensed-panel wood models | Higher-fidelity calibration sources. |
| Leng, Schafer & Buonopane (2013); Buonopane et al. (2015); Kechidi & Bourahla (2016) | CFS-NEES models, fastener-based CFS walls, `CFSWSWP` material | CFS backbones and the frame + spring idealization. |
| Bentley RAM Structural System (RAM SS) walls FAQ (KB0112985) | Shell-element walls in a commercial program | A catalogue of pitfalls this design addresses: corner flange effects, lintel coupling, pilaster load sharing, gravity-induced wall shear, walls on beams (Section 10). |

(Academic references are from memory; check details before quoting them externally.)

## 3. Why a macro-element

Two smoother options were rejected:

- **p-version edge modes:** polynomial displacement along the chords.
- **Smeared studs:** studs as an elastic foundation under a top-chord beam.

Both fail on what drives light-frame design: stud packs and beam bearings (point stiffnesses and point
loads), openings (discontinuities), and per-stud nonlinearity. The basis that represents these exactly is the
framing layout itself.

A side benefit is no mesh dependence. RAM's FAQ shows wall force splits changing with shell mesh size. Here,
the internal model is fixed by the framing, and only the station choice (Section 7) is a modelling decision.

## 4. Input

| Item | Content |
|---|---|
| Geometry | Four corner nodes (counter-clockwise from bottom-left). Defines `L`, `H` and the local frame: `x` along the bottom chord, `y` up. |
| Studs | `[x, area, material, eiOop?]` rows, sorted by `x`. A stud pack is one row; this assumes the plies are fastened to deform together (Section 12). The first and last rows are the end posts (`x = 0`, `x = L`). `eiOop` is used only for out-of-plane output. |
| Top chord | `ea`, `ei`: rim plus double top plate (rim edgewise, plates flatwise, non-composite by default), or plates only where there is only blocking. Elastic. `distribution: "beam" \| "direct"` (Section 6.2). |
| Bottom chord | `ea`, `ei` of the sole plate. Elastic. |
| Sheathing | `material` (uniaxial: sheathing shear strain `γ_sh` → shear flow `q`, force per length) and `install: "always" \| "stage"`, or absent. |
| Framing racking | Optional `rackingStiffness` (force per length per unit wall racking strain; Section 6.3). |
| Openings | `[x0, x1, headerHeight, header: {ea, ei}, jack: {area, material}, chordOverOpening: "continuous" \| "released"]`. Studs with `x0 < x < x1` become cripples above the header. |
| Stations | Extra exterior nodes on the top or bottom chord at given `x` (Section 7). |
| Weight and mass | Optional weight and mass per wall area, on the net area (openings excluded). Applied and lumped at the chord nodes by tributary length. |

Many walls share a layout, so the wire format keeps layouts in an arena and the element rows reference them
(Section 11).

## 5. Internal model

### 5.1 Internal nodes and DOFs

- **Chord nodes** on the top (`y = H`) and bottom (`y = 0`) chords at every stud `x`, every opening edge and
  every station. DOFs `u, v, θ` in the local frame (`u, v` only with `distribution: "direct"`).
- **Header nodes** at each opening's `x0`, `x1` and its cripple positions. DOFs `v, θ`. Header `u` is slaved
  to the chords (Section 5.3).

Corners and stations are exterior; every other node is internal and condensed (Section 5.6).

### 5.2 Members

| Member | Model |
|---|---|
| Chord segment between adjacent chord nodes | Euler–Bernoulli beam (`ea`, `ei`). Loads between nodes enter as fixed-end forces, which is exact. With `chordOverOpening: "released"`, the top chord has a moment release at both opening edges. |
| Full-height stud | Axial bar from the bottom chord node to the top chord node at its `x`, pinned at both ends. |
| Jack stud | Axial bar from the bottom chord node at `x0` (or `x1`) to the header end node. |
| Cripple | Axial bar from a header node to the top chord node at the same `x`. |
| Header | Euler–Bernoulli beam between header nodes, pinned in bearing on the jacks. Bearing is bilateral in v1, so header uplift is resisted. |
| Sheathed bay | Hybrid shear field between two adjacent full-height studs (Section 5.4). |

King studs are ordinary full-height studs listed next to the opening.

### 5.3 Stability of the internal frame

With positive elastic framing stiffness, every internal DOF is stiffened without sheathing:

- Chord `u` by the chord's axial stiffness.
- Chord `v` and `θ` by chord bending and the studs.
- Header `v` and `θ` by header bending, jacks and cripples.

Header `u` has no stiffness of its own, so it is interpolated linearly in `y` from the chord `u` values at the
same `x`. Openings have no shear fields in v1, so this constraint carries no force.

So `K_ii` is non-singular, and the only mechanism is exterior: the wall racking when it has no sheathing
(Section 6.3). That holds for positive elastic stiffness only. In unsheathed walls with nonlinear studs, a
zero-tangent branch (gap, no-tension, rupture) can make `K_ii` singular. That is reported as an element
failure (Section 5.6), never regularized.

### 5.4 Sheathing: per-bay hybrid shear fields

**Field.** Each sheathed bay `j` (width `b_j`, between full-height studs `j` and `j+1`, area `A_j = b_j H`) has
a constant shear flow `q_j`. Its edge members carry the axial forces that keep that field in equilibrium:

- A stud between bays `j−1` and `j` varies linearly over its height with slope `q_j − q_{j−1}`. An end post
  varies with `±q_j`.
- Each chord segment inside a bay varies along its length with slope `±q_j`.

**Bay strain from the actual chord segmentation.** The work of `q_j` on the boundary displacements is
`q_j A_j γ_j`, with

```
γ_j = (1/(b_j H)) [ ∫_bay u_top dx − ∫_bay u_bot dx ] + (v̄_{j+1} − v̄_j) / b_j
```

- The integrals use the chords' piecewise-linear axial interpolation, so they are trapezoid sums over every
  chord segment in the bay, including segments created by stations.
- `v̄_i` is the mean of stud `i`'s end `v` values.

This replaces the earlier four-corner average, which was wrong whenever a station split a bay.

**Lag flexibility.** Each member's force is its mean (carried by the bar or beam) plus a linear part driven by
`q`. The cross terms vanish because `∫ s ds = 0` over each member. The linear part adds complementary energy
`½ qᵀ C q`, with:

```
C = Σ_studs H³/(12 EA_i) d_i d_iᵀ + Σ_chord segments k in bay j ℓ_k³/(12 EA_chord) e_j e_jᵀ
```

- `d_i = e_j − e_{j−1}` for an interior stud, `±e_j` for an end post.
- `ℓ_k` is each chord segment's own length, so stations reduce the chord terms consistently rather than
  counting the same flexibility twice.

**Why framing must be elastic.** The split into a mean force (bar) plus an independent quadratic lag term is
exact only for linear elastic members. A nonlinear stud with force `N̄ + Δq·s` yields or unloads at different
heights, and the split no longer holds. Using initial tangents would quietly produce member forces that
don't satisfy the stud's own material law. Hence the v1 rule in Section 1. Lifting it means integrating the
varying member force at several points with their own material states, and deriving the consistent
flexibility.

**Local solve, using the forward material law.** The unknowns are the bays' sheathing strains `γ_sh`, not the
shear flows, so pinching, plateaus and descending branches stay well defined. With `D = diag(A_j)`,
`q = σ(γ_sh)` from each bay's material (its normal strain-to-stress evaluation from the committed state), and
`γ0` the installation offsets (Section 6.1):

```
r(γ_sh) = D (γ − γ0 − γ_sh) − C σ(γ_sh) = 0
```

- Jacobian: `−(D + C K_t)`, with `K_t = diag(k_t,j)` the material tangents.
- Bay tangent seen by the displacements: `dq/dγ = K_t (D + C K_t)⁻¹ D`.
- Materials are committed at `γ_sh`, not at the total bay strain.
- A singular or non-converging local system (for example, a steep descending branch with stiff coupling) is
  an element failure (Section 5.6). Tangents are never floored.

**What this does and doesn't claim.** The field is in equilibrium and the coupling between neighbouring bays
is exact for the assumed field. That is not the same as the assumed field being accurate for real fasteners.

- In the idealized case (equal bays, equal sheathing, rigid chords) racking gives equal `q_j`, so overturning
  goes to the end posts as the codes assume.
- With unequal bays, materials or chord flexibility, the `q_j` differ, and the element reports them as they
  are.
- Under a point load on a pack, the pack force decays linearly within the assumed field. Whether real walls
  do the same is a calibration question.

### 5.5 Openings, laterally

v1 gives shear fields only to full-height bays; bays over and under openings are unsheathed. That defines the
resisting regions but is **not by itself segmented design**, so the earlier "conservative" claim is
withdrawn. Segmented behaviour also needs:

- **Anchorage at each pier boundary.** Bottom-chord stations at both sides of every opening are generated
  automatically, so hold-down and bearing springs can attach at the pier end posts (Section 7).
- **A defined coupling between piers.** The continuous top chord and the header couple piers by frame
  action, which RAM's FAQ discusses for lintel beams (beam to wall corners vs rigid-link extension).
  - `chordOverOpening: "released"` puts moment releases in the top chord at the opening edges, giving
    pin-connected piers, close to segmented design.
  - `"continuous"` keeps the coupling.
  - Either way, the element reports chord moments at the opening edges so the coupling is visible.

Force transfer around openings (sheathing above and below, strap forces in the header and sill) is a
follow-up that needs header and sill `u` DOFs and partial-height bay fields.

### 5.6 Condensation, state determination and failure

Split DOFs into exterior `e` (corners, stations) and internal `i`. Each evaluation starts from the committed
internal state, as `ForceBeamColumn` does, so there is no mutable trial cache and assembly stays `&self`:

1. Solve internal equilibrium `R_i(u_e, u_i) = f_int,i(u_e, u_i) − p_i = 0` by Newton. Inside each iteration,
   the bay system of Section 5.4 is solved for `γ_sh`. Retry with step subdivision.
2. Condense: `K = K_ee − K_ei K_ii⁻¹ K_ie`. The exterior resisting force `f_e = f_int,e(u_e, u_i)` includes the
   effect of the internal loads, so `p_i` is never condensed into a separate exterior load (Section 8).
3. On `commit`, store `u_i`, `γ_sh` and the material states.

**Failure channel.** Exhausted iterations, a singular `K_ii` or bay Jacobian, or an invalid material state
must fail the step so the analysis can cut back, leaving committed internal state, installation offsets and
material histories unchanged.

The element interface has no way to report this today: `assemble_tangent` and `commit` return nothing
(`core/src/model/elements/ops.rs`). `ForceBeamColumn` currently returns its last state once it has
subdivided 256 times, even if unconverged (`force_beam_column.rs:291`). A separate task is adding an
element-failure path; this element depends on it.

**Workspace and cost.**
- Assembly is `&self` on `Send + Sync` elements, so the element can't own mutable scratch. Scratch comes from
  an assembly context: one preallocated arena per worker, sized at `prepare` from the largest element's
  needs. This is an infrastructure addition shared with Section 11.
- Material trial evaluation currently allocates (`Pinching4` and `Hysteretic` box their next state). A
  separate task is removing that; until it lands, this element inherits the allocation.
- Eliminating the `γ_sh` unknowns couples all bays densely. Instead, keep them in a bordered local system:
  banded framing unknowns (two chord chains coupled stud by stud, bandwidth about 9) plus a small dense block
  of `n_bays` sheathing unknowns. The condensed exterior matrix is small and dense either way. Benchmark a
  12 m wall with 30 studs and 3 openings before committing to the approach.

## 6. Sheathing, gravity and racking

### 6.1 Not relying on the sheathing

> Implemented (2026-11-26) as one stage action for every wall rather than a per-wall `install` property: a sequence with an
> `installSheathing` stage builds every sheathed or floor-zone wall uninstalled, and the floor-zone shear bays are staged with
> the sheathing. Uninstalled fields keep 1e-6 of their stiffness instead of a separate `rackingStiffness`. See
> [stud-wall.md §6](stud-wall.md#staged-sheathing-gravity-by-the-framing-lateral-by-the-sheathing).

Sheathing spreads gravity load between studs. That effect is real but unquantified: sheathing is often
installed after most of the dead-load deformation has happened through direct bearing.

| Mode | Behaviour |
|---|---|
| No sheathing | No shear fields. Gravity goes down through studs, chords, headers and jacks only (a stick takedown with beam distribution). Nonlinear framing is allowed. |
| `install: "always"` | Shear fields active from the start, with `γ0 = 0`. Framing must be elastic. |
| `install: "stage"` | Inactive until an `installSheathing` stage action, which records each bay's `γ0` from the **committed** state and activates the fields stress-free. Dead load applied before installation stays in the framing. Framing must be elastic. |

Installation rules: each element is installed once. Installing twice, installing an element without
sheathing or with `install: "always"`, or naming an unknown element is a `DecodeError`.

### 6.2 Point loads on packs: direct bearing

RAM's FAQ notes that engineers want pilasters designed for at least the full load from the column above, even
when a wall analysis shares it out. Here the analogue is a stud pack under a beam. Even without sheathing, a
stiff rim spreads the pack's point load into neighbouring studs.

`distribution: "direct"` on the top chord gives RAM-gravity-wall behaviour ("not fanned out"):

- The chord has no bending stiffness, and its internal `θ` DOFs and corner `rz` are dropped (the mask
  follows).
- A load between studs goes to the two adjacent studs by the lever rule.
- A load over a stud or pack goes straight down it.

PySees can run the takedown both ways and envelope the two for pack and stud design.

### 6.3 Racking without sheathing

Without sheathing, the wall has one exterior racking mode, which the model must close explicitly. Either:

- `rackingStiffness` adds an elastic spring on the wall racking strain
  `γ_w = (ū_top − ū_bot)/H + (v̄_R − v̄_L)/L`, where `ū` is the chord's mean `u` and `v̄_L`, `v̄_R` are the
  end posts' mean `v`. This strain is invariant under rigid rotation. The spring stays active after
  installation, and its shear is included in the reported wall shear. Or:
- the model restrains the mode (for example `ux` fixes in a gravity-only model).

Otherwise the analysis stops with the solver's singular-matrix error. Vertical loads on pinned vertical
studs don't excite the racking mode, so a small `rackingStiffness` regularizes the model without changing
gravity results.

### 6.4 Analysis history for design forces

Design forces are differences along one defined analysis history. They are not general superposition, since
nonlinear sheathing and anchorage are history-dependent. Recommended history (PySees default):

1. Dead load, sheathing not installed.
2. Remaining gravity (superimposed dead, live): the **takedown state**. Stud, pack, jack and header gravity
   design forces are read here, with `distribution` run both ways if requested.
3. `installSheathing`.
4. Lateral stages. Lateral design forces are each stage's committed end state minus its committed start
   state.

Eccentric gravity creates wall shear and overturning even without lateral load (RAM's FAQ, "Shear and Moment
in Walls due to Gravity Loads"). The stage increments above keep that out of lateral design forces, and
the outputs label gravity-stage shears so users aren't surprised.

## 7. Exterior nodes and stations

**Corners** carry `ux`, `uy` and `rz`. `rz` is the chord rotation; it is dropped with `distribution: "direct"`.

**Stations** are chord nodes promoted to exterior nodes. They are needed wherever something outside the
element attaches along a chord:

- **Stacked storeys.** The upper wall's sole plate bears on the floor line. Each floor line gets a station at
  every stud position of the walls above and below (PySees computes the union; this is the default).
- **Walls on beams.** A wall bearing on a global beam needs a station at every stud on its bottom chord,
  shared with beam nodes, or the stud loads can't reach the beam. RAM's FAQ describes the reverse problem:
  a wall stiff enough to span over its supporting beam as a deep beam. Takedown without sheathing avoids
  this, and the stations make the load path explicit.
- **Pier anchorage.** Bottom-chord stations at both sides of every opening, generated automatically
  (Section 5.5), plus any the user requests for hold-downs.
- **Global beams and posts** in the wall plane that frame into the wall.

**Sharing a node is a connection idealization.** It enforces full compatibility, including tension transfer
and rotation continuity. Adding the two chords' bending stiffnesses is valid for that idealization only. It
does not model unilateral bearing, fastener slip, rim offset or floor depth. Where those matter, use separate
coincident nodes joined by `zeroLength` springs (compression-only bearing, slip).

The same applies between adjacent segments in a wall line. Sharing a corner node makes the chord
moment-continuous, but plate splices aren't moment connections. For a hinge, use separate corner nodes tied
in translation only (`equalDOF ux uy`).

**Posts at wall ends.** A global column placed at a wall end shares the corner node with the wall's end
post. Both then carry load, and self-weight is counted twice if both include it (RAM's FAQ, "Combining Columns
and Walls"). Model a post either as the end-post stud row or as a global column with that row's area set to
zero, never both. PySees validates this; the engine can't tell the two apart.

## 8. Loads

### 8.1 Kinds

| Load | Use |
|---|---|
| `wallTopPoint { x, fx?, fy }` | Beam reaction or girder bearing at any `x`. |
| `wallTopLine { x0, x1, wy0, wy1 }` | Floor/roof joist loads (trapezoidal). |
| `wallSelfWeight { w }` | Weight per net wall area, applied at the chord nodes by tributary length. |
| `wallOopPressure { p }` | Out-of-plane pressure. No in-plane load; feeds the out-of-plane output (Section 9.2). |

### 8.2 Representation: pattern slots

The catalog's `Load` type must be `Copy`, additive and scalable (`ops.rs`). Loads at different positions
can't be combined by adding their coordinates, so the positions stay out of the `Load` value:

- At decode, each pattern's loads on a wall are discretized into a consistent load vector over the wall's
  internal and exterior DOFs (fixed-end forces for loads between nodes). The vector is stored in the element
  in **pattern slot** `k`.
- The element's `Load` is `StudWallLoad { factors: [f64; K], oop: f64 }`. Adding and scaling act on the
  factors, so the effective load `Σ_k factors[k] · p_k` is exact, `Copy` and allocation-free.
- `K` is a fixed capacity (8 to start). A wall loaded by more than `K` patterns is a `DecodeError` naming the
  wall. Out-of-plane pressure is a scalar, so it adds directly.

### 8.3 Residual contract

- **Internal parts** `p_i` enter only the internal equilibrium (Section 5.6). Their effect reaches the global
  system through `f_e`.
- **Exterior parts** `p_e` (loads exactly at corners or stations) go through `assemble_load` as ordinary
  nodal loads.
- Nothing is subtracted twice: `assemble_tangent` returns `f_e` without `p_e`, and `local_load_force` returns
  only the exterior part.

Recorders use the same split. Exterior forces at corners and stations are `local_force`. Stud, bay, chord
and header responses use a dedicated response interface (Section 9), not the nodal-force vector.

## 9. Output

### 9.1 In-plane

Variable width, sized at decode from the layout and bounds-checked like other recorders.

- **Wall, at the base section** (the bottom chord line):
  - shear `V = Σ q_j b_j` plus the framing racking spring's shear;
  - axial force;
  - overturning moment about the mid-length of the base, `x = L/2`;
  - racking strain `γ_w`.

  The sign convention is documented with the response, following the element's local frame. PySees sums
  segments into wall lines.
- **Studs, jacks, cripples:** force at the top and the bottom (`N̄ ± Δq·H/2` in sheathed elastic walls).
- **Bays:** `q_j`, `γ_j`, `γ_sh,j`.
- **Chords and headers:** moment and shear at the internal nodes, including at opening edges; header end
  reactions.

### 9.2 Out-of-plane, per stud

Requires a `wallOopPressure` and `eiOop`. Each member spans its own length: full-height studs `H`, jacks to
the header, cripples from the header to the top chord.

- **Tributary width:** half the bay on each side. v1 assumes that jambs collect half of each adjacent
  opening's width. The real path through the header, sill and jamb is a follow-up.
- **First-order results** (always reported): simple-span moment `w·trib·ℓ²/8`, end reactions, deflection
  `5w·trib·ℓ⁴/(384EI)`.
- **Amplified moment** uses the beam-column secant formula with a **constant** compression `P`, taken as
  the member's maximum compression along its length. This is conservative, and it is an approximation
  whenever the force varies along the member.
  - With `k = √(P/EI)`: `M = (w·trib·EI/P)(sec(kℓ/2) − 1)`.
  - Evaluated with a series expansion as `P → 0`.
  - A member in tension reports first-order values.
  - At `P ≥ P_cr = π²EI/ℓ²` it reports an instability flag instead of a number.

Pairing out-of-plane pressures with gravity cases, and code checks, belong to PySees.

## 10. Lessons from RAM Structural System (KB0112985)

| RAM SS issue | Handling here |
|---|---|
| Walls sharing corner nodes act as a 3D flange system; overturning splits between parallel and perpendicular walls, and depends on mesh size. | 2D in v1. The 3D follow-up must make corner coupling an explicit choice per corner (shared post vs separate), because code design assumes independent segments with their own hold-downs. |
| Lintel coupling depends on whether the beam meets wall corners only or extends with rigid links. | `chordOverOpening` and explicit hinges between segments (Sections 5.5, 7). |
| Pilasters attract less load than the column above delivers. | `distribution: "direct"` and an envelope (Section 6.2). |
| Column and wall overlapping counts self-weight twice. | One owner per post (Section 7). |
| Eccentric gravity creates wall shear; diaphragms make it internal. | Defined analysis history and labelled gravity-stage shears (Section 6.4). |
| A wall over a beam spans as a deep beam and unloads the beam. | Takedown without sheathing, plus bottom stations shared with beam nodes (Sections 6.1, 7). |
| Reported forces are defined at the storey's base section about the centroid. | Same definition (Section 9.1). |

## 11. Carapace integration

| Item | Plan |
|---|---|
| Element | `core/src/model/elements/stud_wall/`, `Element::StudWall`, 2D only. A generic condensed-substructure core (bars, beams, hybrid bay fields, bordered local solve) plus a stud-wall generator, so the substructure can be tested directly. |
| Variable connectivity | Shared infrastructure, not an escape hatch for one element. `NodeList` (an `ArrayVec` capped at 4, `ops.rs`) becomes a borrowed `&[NId]` that elements store at construction. Tangent and vector sinks gain slice-based `add` paths. DOF activation, `dof_mask` conformance, constraints, reactions, mass and recorders iterate the slice. Serial and parallel scatter keep element order and a fixed order within each element, so results stay bit-identical. Closed enums and generic dispatch are unchanged. |
| Assembly scratch | Per-worker arena in the assembly context (Section 5.6). |
| Failure | Depends on the element-failure path being added separately (Section 5.6). |
| `dof_mask` | `ux`, `uy`, `rz` at corners and stations; `ux`, `uy` with `distribution: "direct"`. |
| `validate` | **Geometry:** rectangle within tolerance and counter-clockwise; studs sorted, within `[0, L]`, with posts at both ends; openings inside `(0, L)`, not overlapping, header below `H`, at least one stud between openings; stations on a chord within `[0, L]`. **Properties:** positive areas, `ea` and `ei` (except `ei` with `"direct"`). **Materials:** elastic framing whenever sheathing is present. Errors are `ModelError::InvalidElement` with the reason. |
| `prepare` | Generate the substructure; factor its linear parts; record the scratch size. |
| Loads | Pattern slots (Section 8.2). |
| Stage action | `installSheathing { elements: "all" \| index[] }` (Section 6.1). |
| Wire | `studWallLayouts: StudWallLayoutSpec[]`, `StudWallTable` rows (`nodeIds` stride 4, `layout`, `stations: [chord, x, node][]`, `weight?`, `areaMass?`), and the load kinds above. Document in `docs/input-format.md`; mirror in PySees (`carapaceInputV1.ts`, `compileInputV1.ts`), including the station union, the post-ownership check and the default analysis history. |
| Materials | `pinching4` and `hysteretic` cover sheathing backbones now. Add `saws` (Folz–Filiatrault) separately, cross-checked against OpenSeesPy. |

## 12. Assumptions and calibration (stated, not modelled)

- **Stud packs:** the plies are fastened to share load, so their areas add. Loose plies are a design check.
- **Sheathing calibration:** each bay's `q–γ_sh` law is a wall-level property, from full-height wall tests
  or code `G_a`. Its initial tangent is SDPWS's apparent shear stiffness `G_a` (panel shear plus fastener
  slip) in consistent units. It already includes panel seams, blocking, the fastening schedule and two-sided
  sheathing (sum the two sides' laws). Aspect-ratio reductions (SDPWS 4.3.3.2) apply per pier: calibrate on
  the pier's full-height width, not the stud bay.
- **CFS:** track flexibility, stud-to-track connections, local and distortional buckling, and effective
  section properties are not modelled. Use effective axial stiffness and check stability in design.
- **Wood:** creep, moisture, shrinkage and perpendicular-to-grain bearing are not modelled.
- **Connections:** header bearing is bilateral, and plate and rim splices are continuous unless the user
  models a hinge (Section 7).

## 13. Verification

| Test | Reference | Tolerance |
|---|---|---|
| Mask conformance; rigid-body modes give zero force; rank with corners only (12 DOFs, 3 rigid-body modes): 9 with sheathing, 8 without; elastic framing | Kinematics | 1e-12 |
| One- and two-bay boundary work and `C`, derived independently: unequal widths and `EA`, end-post signs, symmetry, energy, force and moment balance | Hand derivation | 1e-12 |
| An unloaded station inserted inside a bay leaves the response unchanged | Consistency of the bay strain and `C` | 1e-10 |
| Without sheathing, linear: condensed element vs the same frame from `Truss` and `ElasticBeamColumn` | Carapace explicit model and OpenSeesPy | 1e-10 |
| Without sheathing, nonlinear studs with reversals | Explicit model | 1e-8 |
| Loaded condensation: point loads, partial trapezoids, three patterns with different factors; exterior forces, reactions and internal recovery; finite-difference tangent at fixed load; no double subtraction | Explicit model; finite differences | 1e-8 |
| Pack point load, `"beam"` vs `"direct"` distribution; blocking vs rim | Beam on discrete springs; lever rule | 1e-10 |
| Header over an opening, `"continuous"` vs `"released"` chord | Simply supported beam with cripple loads; frame closed form | 1e-10 |
| SDPWS 2021 Eq. 4.3-1, bending and shear terms: one storey; end posts `EA`; interior studs present; chord `EA` → 1e6× and chord `EI` → 1e-6× (controlled limits, never zero); rigid base with hold-downs at the posts; `rackingStiffness = 0`; top load split equally to both corners | `VH/(L G_a) + 2VH³/(3EA L²)` in consistent units | 1e-6 |
| CSA O86:19 A.11.7.1 a), three storeys with different storey loads, same limits | Per storey: `θ_i = M_i H_i/(EI) + V_i H_i²/(2EI)`, `Δ_b,i = M_i H_i²/(2EI) + V_i H_i³/(3EI)`, storey drift `Δ_b,i + H_i Σ_{j<i} θ_j` plus shear, with `(EI) = EA L²/2` and `M_i` the moment at the top of storey `i`; compare every storey, not just the roof | 1e-6 |
| Nonlinear sheathing (`pinching4`, then `hysteretic` with a descending branch) with **flexible** elastic framing, cyclic | Explicit reference: studs split into `k` bars, each sub-bay carrying a `zeroLength` shear spring; convergence in the end-post force norm as `k` grows | trend, documented |
| Stiff-framing limit of the previous test | OpenSeesPy `zeroLength` with the same material scaled by `L` and `H` | 1e-9 |
| Failure and rollback: forced local non-convergence, singular `K_ii`, failed global step and retry | Retried run bit-identical to a clean run; committed state unchanged after failure | exact |
| Installation lifecycle: deformation before installation, more gravity after, unloading; repeat installation and invalid selections rejected; offsets taken from committed state | `q = 0` at installation; superposition in the linear case | 1e-10 |
| Out-of-plane: closed forms, `P → 0` continuity, tension, the `P_cr` flag; members with different spans | Closed form | 1e-12 |
| Infrastructure: constraints at stations, rotated walls, mass conservation, recorder bounds, 3D rejected, native/wire parity, wasm smoke test, no allocation in assembly, serial/parallel bit-identical above the parallel threshold | — | exact |

## 14. Follow-ups and open decisions

**Follow-ups:**
- Force transfer around openings.
- Sills and out-of-plane load paths through the header, sill and jamb.
- Raked walls (Garvey's non-uniform shear flow).
- Nonlinear framing in sheathed walls (distributed member states, Section 5.4).
- Perpendicular-to-grain bearing springs.
- 3D with explicit corner coupling (Section 10).
- Stud in-plane bending.

**Open decisions:**
- Default `rackingStiffness` that PySees emits for unsheathed stages, or `ux` restraints instead.
- Pattern-slot capacity `K`.
- Whether PySees envelopes `"beam"` and `"direct"` distribution by default for takedown.
- Exact (union) stations as the only mode, or a reduced-station option for large models.
