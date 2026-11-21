# Stud-wall macro-element

Status: design only, nothing implemented.

A 2D element for light-frame walls (wood or cold-formed steel) that covers one wall segment for one storey and
is defined by its framing layout: studs and stud packs, top and bottom chords (plates, rim), openings with
headers, and optional sheathing. Internally it is a small frame of studs, chord beams, headers and sheathing
shear fields. All of it is condensed onto a few exterior nodes, so users model walls, not studs.

The same element serves:

- **Gravity takedown**, without relying on sheathing: non-uniform loads, stud packs, rim vs blocking, headers.
- **Lateral (in-plane) analysis**, with sheathing: matches the SDPWS and CSA O86 wall deflection equations in
  the elastic case, and takes any uniaxial hysteretic material for the sheathing.
- **Out-of-plane wind** as per-stud output, with no extra DOFs.

## 1. Scope

In v1:

- 2D (ndm=2), rectangular wall segments, small displacements.
- Studs are pinned at the chords; their in-plane bending is ignored.
- Openings carry gravity through headers and jacks. Laterally they use the segmented-wall assumption: only
  full-height sheathed bays resist shear (Section 5.5).

Out of scope for v1, listed in Section 11:

- Force transfer around openings
- 3D
- Non-rectangular walls
- Perpendicular-to-grain bearing and shrinkage
- Sills

## 2. Prior work

| Work | What it does | Relation |
|---|---|---|
| Kasal, Leichti & Itani (1994) | Condensed wall "superelements" in whole-house models | Same idea: a condensed wall element. This design adds gravity load paths and staged sheathing. |
| Garvey (1951); NASTRAN `CSHEAR` + `CROD` | Shear panels with constant shear flow; edge members whose force varies linearly | The per-bay sheathing field (Section 5.4). |
| Pian (1964), assumed-stress hybrid elements | Stiffness from an equilibrium stress field | How the sheathing–stud coupling is made exact. |
| Folz & Filiatrault (2001, 2004), CASHEW / SAWS | Nail-level wall model condensed to one hysteretic spring | Source of the sheathing backbone. |
| Judd & Fonseca (2005); Pang & Rosowsky; Christovasilis & Filiatrault (Timber3D) | Nail-spring, beam-spring and condensed-panel wood models | Higher-fidelity calibration sources. |
| Leng, Schafer & Buonopane (2013); Buonopane et al. (2015); Kechidi & Bourahla (2016) | CFS-NEES models, fastener-based CFS walls, `CFSWSWP` material | CFS backbones and the frame + spring idealization. |

OpenSees users build walls from trusses, rigid links and a shear spring by hand, or model every stud.
Neither gives gravity load paths and lateral response from one wall definition. (References are from
memory; check details before quoting them externally.)

## 3. Why a macro-element

Two smoother options were considered and rejected:

- **p-version edge modes:** polynomial displacement along the chords.
- **Smeared studs:** studs as an elastic foundation under a top-chord beam, with exact beam-on-foundation
  shape functions.

Both fail on the features that drive light-frame design: stud packs and beam bearings (point stiffnesses and
point loads), openings (discontinuities), and per-stud nonlinearity. The basis that represents these exactly
is the framing layout itself. So the element generates the stud-level frame internally and condenses it.
Users enter one element per wall segment per storey; the solver sees only corner and station DOFs.

## 4. Input

| Item | Content |
|---|---|
| Geometry | Four corner nodes (counter-clockwise from bottom-left). Defines `L`, `H` and the local frame: `x` along the bottom chord, `y` up. |
| Studs | `[x, area, material, eiOop?]` rows, sorted by `x` (`eiOop` is only for out-of-plane output). A stud pack is one row with the plies' combined area. The first and last rows are the end posts (`x = 0`, `x = L`). |
| Top chord | `ea`, `ei`: rim plus double top plate (rim edgewise, plates flatwise, non-composite by default), or plates only where there is only blocking. Elastic. |
| Bottom chord | `ea`, `ei` of the sole plate. Elastic. |
| Sheathing | `material` (uniaxial: shear strain `γ` → shear flow `q`, force per length) and `install: "always" \| "stage"`, or absent. |
| Framing racking | Optional `rackingStiffness` (force per length per unit `γ`): racking stiffness of the unsheathed frame. See Section 6. |
| Openings | `[x0, x1, headerHeight, header: {ea, ei}, jack: {area, material}]`. Studs with `x0 < x < x1` become cripples above the header. |
| Stations | Extra exterior nodes on the top or bottom chord at given `x` (Section 7). |
| Mass | Optional mass per wall area, lumped to the exterior nodes by tributary length. |

Many walls share a layout, so the wire format keeps layouts in an arena and the element rows reference them
(Section 9).

## 5. Internal model

### 5.1 Internal nodes and DOFs

- **Chord nodes** on the top (`y = H`) and bottom (`y = 0`) chords at every stud `x`, every opening edge and
  every station. DOFs `u, v, θ` in the local frame.
- **Header nodes** at each opening's `x0`, `x1` and its cripple positions. DOFs `v, θ`. Header `u` is slaved
  to the chords (Section 5.3).

Corners and stations are exterior. Every other node is internal and condensed (Section 5.6).

### 5.2 Members

| Member | Model |
|---|---|
| Chord segment between adjacent chord nodes | Euler–Bernoulli beam (`ea`, `ei`); loads between nodes enter as fixed-end forces, which is exact. |
| Full-height stud | Axial bar from the bottom chord node to the top chord node at its `x`, pinned at both ends; any uniaxial material. |
| Jack stud | Axial bar from the bottom chord node at `x0` (or `x1`) to the header end node. |
| Cripple | Axial bar from a header node to the top chord node at the same `x`. |
| Header | Euler–Bernoulli beam between header nodes, pinned where it bears on the jacks. |
| Sheathed bay | Hybrid shear field between two adjacent full-height studs (Section 5.4). |

King studs are ordinary full-height studs that the user lists next to the opening.

### 5.3 Stability of the internal frame

Every internal DOF is stiffened without sheathing:

- Chord `u` by the chord's axial stiffness.
- Chord `v` and `θ` by chord bending and the studs.
- Header `v` and `θ` by header bending, jacks and cripples.

Header `u` has no internal stiffness of its own, so it is interpolated linearly in `y` from the chord `u`
values at the same `x`. Because v1 treats openings as unsheathed laterally (Section 5.5), this constraint
carries no force and changes nothing. So `K_ii` is non-singular with or without sheathing. The only mechanism
is exterior: the whole wall racking when it has no sheathing (Section 6).

### 5.4 Sheathing: per-bay hybrid shear fields

Each sheathed bay `j` (width `b_j`, between full-height studs `j` and `j+1`) has a constant shear flow `q_j`.
Its generalized strain, from the four bay corners, is:

```
γ_j = (ū_top,j − ū_bot,j) / H + (v̄_{j+1} − v̄_j) / b_j
```

- `ū` is the mean chord `u` at the bay's two stud positions.
- `v̄_i` is the mean of stud `i`'s end `v` values.

**Equilibrium.** Following Garvey, a stud between bays `j−1` and `j` carries a force that varies linearly
along its height with slope `q_j − q_{j−1}`. Each chord segment carries a force that varies with `q_j` along
the bay. So the edge members' complementary energy has a mean-force part, carried by the bars and beams above,
plus

```
U*_lag = ½ qᵀ C q,   C = Σ_studs H³/(12 EA_i) d_i d_iᵀ + Σ_bays (b_j³/(12 EA_top) + b_j³/(12 EA_bot)) e_j e_jᵀ
```

- `d_i = e_j − e_{j−1}` for an interior stud.
- `d_i = ±e_j` for an end post.
- `C` uses the stud materials' initial tangents.
- The cross terms with the mean forces vanish (`∫ s ds = 0`).

**Conjugacy.** Conjugate to `q_j` is `A_j γ_j`, with `A_j = b_j H`. The sheathing state for a given set of
bay strains solves

```
A_j (γ_j − γ0_j) = A_j γ_sh,j(q_j) + (C q)_j          for all sheathed bays j
```

- `γ_sh,j` is the inverse of the sheathing material.
- `γ0_j` is the strain recorded at installation (Section 6).

This is a small system of size `n_bays` with a tridiagonal Jacobian `diag(A_j / k_t,j) + C`, solved by
Newton. With elastic sheathing it is linear. The bays' contribution to the internal stiffness is
`Bᵀ diag(A) (diag(A/k_t) + C)⁻¹ diag(A) B`, where `B` holds the rows of `γ_j`.

**What this gets right:**
- Under uniform racking all `q_j` are equal. Interior studs then see no net shear, so overturning goes to the
  end posts, as both codes assume.
- Under a point load on a pack, the neighbouring bays have opposite `q`, and the pack's force decays down the
  height exactly. That is the `(q₁ − q₂)` coupling the earlier 4-node design could only approximate.

### 5.5 Openings, laterally

v1 uses the segmented-wall assumption: bays over and under openings are unsheathed. Only full-height bays get
shear fields. This is conservative and matches segmented design. Force transfer around openings (sheathing
above and below, header and sill straps) is a follow-up that needs header and sill `u` DOFs and partial-height
bay fields.

### 5.6 Condensation and state determination

Split DOFs into exterior `e` (corners, stations) and internal `i`. Given exterior displacements and the
effective element load `p` (the trait already passes this to `assemble_tangent` and `commit`; see
`core/src/model/elements/ops.rs`):

1. Solve internal equilibrium `R_i(u_e, u_i) = f_int,i − p_i = 0` for `u_i` by Newton. Start from the
   committed internal state, as `ForceBeamColumn` does. If it doesn't converge, subdivide the step.
2. Condense: `K = K_ee − K_ei K_ii⁻¹ K_ie`, and `f = f_e − p_e`.
3. On `commit`, store `u_i`, the bay shear flows and the material states.

With linear materials, step 1 is a single solve and the condensation is exact.

**Cost:**
- Internal DOFs are about `6 × n_studs + 2 × n_header_nodes`. The chords and studs form a banded structure:
  two chains coupled stud by stud, so a banded factorization costs `O(n b²)` with `b ≈ 9`.
- Workspace is allocated when the element is constructed, never in the assembly loop.
- The sheathing system adds a dense `n_bays²` block, which is small (a 6 m wall has about 15 bays).

## 6. Sheathing and gravity: not relying on the sheathing

Sheathing spreads gravity load between studs. That effect is real but unquantified: sheathing is often
installed after most of the dead-load deformation has happened through direct bearing. The element therefore
lets the sheathing be absent, or present only from a chosen stage onward.

| Mode | Behaviour |
|---|---|
| No sheathing | No shear fields. Gravity goes to the foundation through studs, chords, headers and jacks only, giving a classic stick takedown with exact beam distribution. |
| `install: "always"` | Shear fields active from the start, with `γ0 = 0`. |
| `install: "stage"` | Inactive until an `installSheathing` stage action, which records each bay's `γ0` and activates the fields stress-free. Dead load applied before installation stays in the framing. |

Without sheathing the wall has one exterior racking mode. The model must close it explicitly. Either:
- the element's `rackingStiffness` adds a small elastic racking spring, or
- the model restrains it (for example `ux` fixes in a gravity-only model).

Otherwise the analysis stops with the solver's singular-matrix error rather than degrading silently. Vertical
loads on pinned vertical studs don't excite the racking mode, so a small `rackingStiffness` regularizes the
model without changing gravity results.

**Recommended practice (PySees default):**
- Design takedown stages run with the sheathing not installed. Stud, pack, jack and header design forces come
  from these stages.
- Sheathing is installed before the lateral stages.
- Stud design forces combine the takedown forces with the lateral-stage increments (the overturning).
- Live load applied after installation is shared with the sheathing in the lateral model. Nothing in the
  gravity design reads that model.

## 7. Exterior nodes and stations

**Corners** carry `ux`, `uy` and `rz`. `rz` is the chord rotation, which gives bending continuity to
neighbouring walls in the same line and to beams.

**Stations** are chord nodes promoted to exterior nodes. They are needed wherever something outside the
element attaches along a chord:

- **Stacked storeys.** The upper wall's sole plate bears continuously on the floor line. For exact load paths,
  each floor line gets a station at every stud position of the walls above and below. Both walls' chords then
  meet at each station, so the floor line acts as both chords together, with their bending stiffnesses added.
  PySees computes this union (the default). With fewer stations, upper studs between stations reach the
  floor line only through chord bending. That approximation is acceptable for aligned (stack) framing.
- **Global beams and posts** in the wall plane that frame into the wall at a point.

Loads don't need stations; they are element loads (Section 8). In the exact default, global DOFs are about
3 per stud per floor line. Studs' internals, sheathing, headers and cripples are all condensed.

## 8. Loads

| Load | Use |
|---|---|
| `wallTopPoint { x, fx?, fy }` | Beam reaction or girder bearing at any `x`. |
| `wallTopLine { x0, x1, wy0, wy1 }` | Floor/roof joist loads (trapezoidal). |
| `wallSelfWeight { w }` | Weight per wall area, applied at the chord nodes by tributary length. |
| `wallOopPressure { p }` | Out-of-plane pressure. Produces no in-plane load; feeds only the out-of-plane output (Section 10). |

Loads land on internal chord nodes, or as fixed-end forces between them. They enter internal equilibrium
through the effective element load (Section 5.6), so nonlinear stages stay consistent.

## 9. Carapace integration

| Item | Plan |
|---|---|
| Element | `core/src/model/elements/stud_wall.rs`, `Element::StudWall`, 2D only. A generic condensed-substructure core (bars, beams, hybrid bay fields, internal Newton) plus a stud-wall generator, so the substructure can be tested directly. |
| Node count | Needs `4 + n_stations` nodes. `MAX_ELEMENT_NODES` is 4 (`ops.rs`), and `TangentSink::add` takes a const `N`. Add a dynamic-size path (`add_dyn` with slices), used only by this element, rather than raising the cap for every element. |
| `dof_mask` | `ux`, `uy`, `rz` at corners and stations. All three are stiffened: chord axial, studs, chord bending. |
| `validate` | Rectangle within tolerance and counter-clockwise; studs sorted, within `[0, L]`, with posts at both ends; openings inside `(0, L)`, not overlapping, header below `H`, at least one stud between openings; stations on a chord within `[0, L]`; positive areas, `ea` and `ei`. Errors are `ModelError::InvalidElement` with the reason. |
| `prepare` | Generate the substructure; factor its linear parts; size the workspace. |
| `assemble_tangent`, `commit` | Section 5.6, reading the effective element load. |
| Stage action | `installSheathing { elements: "all" \| index[] }`: records `γ0` and activates the bay fields. Rejected on elements without sheathing or with `install: "always"`. |
| Wire | `studWallLayouts: StudWallLayoutSpec[]` (studs, chords, sheathing, openings), `StudWallTable` rows (`nodeIds` stride 4, `layout`, `stations: [chord, x, node][]`, `areaMass?`), and the element loads above. Document in `docs/input-format.md`; mirror in PySees (`carapaceInputV1.ts`, `compileInputV1.ts`), including the station union and the takedown/lateral stage defaults. |
| Materials | `pinching4` and `hysteretic` cover sheathing backbones now. Add `saws` (Folz–Filiatrault) as a separate change, cross-checked against OpenSeesPy. |

## 10. Output

Variable width, sized at decode from the layout:

- **Wall:** base shear `V = Σ q_j b_j`, overturning moment, racking `γ`.
- **Studs and jacks:** force at the top and the bottom. The hybrid field gives both: `N̄ ± Δq·H/2`.
- **Bays:** `q_j`, `γ_j`.
- **Chords and headers:** moment and shear at the internal nodes; header end reactions.
- **Out-of-plane, per stud** (when a `wallOopPressure` is present):
  - tributary width (jambs take half the opening);
  - simple-span moment `w·trib·H²/8` and the reactions;
  - deflection `5wH⁴/(384EI)`;
  - moment amplified by the stud's compression using the exact beam-column formula
    `M = (w·trib·EI/P)(sec(kH/2) − 1)`, `k = √(P/EI)`.

  Out-of-plane stud EI comes from the layout (`eiOop` per stud row). Pairing out-of-plane pressures with
  gravity cases, and code checks, belong to PySees.

## 11. Verification

| Test | Reference | Tolerance |
|---|---|---|
| Mask conformance; rigid-body modes give zero force; rank with corners only (12 DOFs, 3 rigid-body modes): 9 with sheathing, 8 without (the racking mode) | Kinematics | 1e-12 |
| Without sheathing, linear: condensed element vs the same frame built from `Truss` and `ElasticBeamColumn` | Carapace explicit model, and OpenSeesPy | 1e-10 |
| Same, with nonlinear studs (`steel01` packs in compression) | Explicit model | 1e-8 |
| Point load on a pack over blocking-only vs rim: stud forces | Beam on discrete springs, closed form | 1e-10 |
| Header over an opening: moment and jack reactions | Simply supported beam with cripple point loads | 1e-10 |
| One-storey cantilever, sheathing on, stiff interior studs, `ei → 0`: drift `VH/(L G_a) + 2VH³/(3EA L²)` | SDPWS 2021 Eq. 4.3-1 (bending and shear terms), consistent units | 1e-8 |
| Three storeys with stations, different storey loads | CSA O86:19 A.11.7.1 a) bending plus shear | 1e-8 |
| Pack under a point load with sheathing: stud force decay | Explicit model with studs split into `k` segments, converging as `1/(4k²)` to the hybrid result | trend |
| `installSheathing`: dead load before installation leaves `q = 0`; lateral after installation is unchanged by the dead load (linear case) | Superposition | 1e-10 |
| Cyclic `pinching4` sheathing, stiff framing: `V` vs drift | OpenSeesPy `zeroLength` with the same material scaled by `L` and `H` | 1e-9 |
| Out-of-plane output | Closed forms above | 1e-12 |
| Wire: decode round trip; each validation error | Expected `DecodeError` / `ModelError` | exact |

## 12. Limitations, follow-ups and open decisions

**Follow-ups:**
- Force transfer around openings: header and sill straps, partial-height sheathed bays.
- Sills and out-of-plane window loads on jambs and sills.
- Non-rectangular (raked) walls, which need Garvey's non-uniform shear flow.
- Perpendicular-to-grain bearing springs under packs, and shrinkage.
- 3D: the wall in its own plane, with out-of-plane DOFs supported by diaphragms.
- Stud in-plane bending, for very slender piers.

**Open decisions:**
- Default `rackingStiffness` that PySees emits for takedown stages, or `ux` restraints instead.
- Whether the lag matrix `C` should track nonlinear stud tangents. v1 uses initial tangents; studs in sheathed
  walls stay nearly elastic.
- Exact (union) stations as the only mode, or a reduced-station option for large models.
