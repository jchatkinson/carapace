# Stud wall bearing base (design record)

Status: implemented with this document. User guide sections: [stud-wall.md §2, §4](stud-wall.md).

## Problem

A stud wall connects to what is below it only at nodes. With `foundation: true` every bottom-plate stud position
bears on rigid ground. Without it, the bottom plate touches nothing between its stations, and the modeller must:

1. put a bottom station at every stud,
2. find a node of the wall below (or a beam) at the same plan point, which needs a top station there, and
3. join each pair with a vertical `ZeroLength` (anchorage or bearing) and an `equal_dof`.

That only works when the studs above line up with the studs below. In real platform framing they often don't: a
second-storey wall set back over a cantilever, a party wall split differently at each storey, a wall on joists.
The load still gets through, because the studs sit on a **sandwich**: one or two bottom plates on the subfloor,
on the rim and joists (or blocking and joists), on the double top plate of the wall below. That sandwich spreads
each stud's load to the studs below on either side.

Found on the BC Rowhouse template (PySees): the grid 3 second-storey party wall (8.9 m, no openings) was held
only at its two pier ends. Its bottom plate spanned 8.9 m between them and its top sagged 178 mm under gravity.
The main-floor walls, without `foundation`, sagged 18 mm between their ground anchors.

The aim: **a wall carries its own connection to whatever it stands on.** A modeller (or PySees) names the nodes it
stands on and never builds per-stud links.

## Decisions

### D1. The bearing is inside the element

The sandwich is modelled inside the stud-wall element, as springs from each bottom-plate position to the
surface below, condensed with the rest of the internal frame. The alternatives were rejected:

- *External links per stud* (today's recipe, automated by a preprocessor): needs aligned nodes, so extra studs or
  stations in the other wall; the links go stale when a wall is edited; and they clutter the model.
- *A separate interface element* between two walls: it would need the other walls' internal plate DOFs, which are
  condensed away.

### D2. The surface below is a set of base nodes, at any positions

`layout.bearing = Some(WallBearing { kv, kh, hold_down })` and stations with the new chord `WallChord::Base`.
A base station is an exterior node on the wall's bottom line (local `y = 0`), at any `x`, including outside
`0..L`, so a wall can bracket its ends with nodes beyond them. Base nodes are what the wall stands on: the
exterior top row of the wall(s) below (corners and top stations, at every stud with a floor zone), beam-column
nodes, or fixed ground nodes. Their number and spacing are unrelated to this wall's studs.

Several base nodes may share one `x` (the corners of two walls that meet end to end, a beam node on a wall
corner). They form one base position, and a load there is split equally between them.

### D3. Linear interpolation between base positions

The surface between consecutive base positions moves linearly: `v_b(x) = N_i v_i + N_j v_j`, likewise `u_b`.
A stud between two studs below therefore splits its load by the lever rule. This is the `Direct` hand takedown,
and it is consistent with the wall below, which carries loads on its top row between stations by the same
rule. The lower top plate's continuity over one stud space is left out, as in the hand method.

Every bearing position must lie within the base positions' range (`x_first <= x <= x_last`), otherwise the
element is rejected with the position named. There is no extrapolation and no snapping.

### D4. A spring at every bottom-plate bearing position

Bearing positions are the bottom-plate nodes under a vertical member: every stud row not inside an opening and
both jacks of every opening (the king stud and jack share the node). Under a door there is none (the plate is
cut); under a window the plate spans between the jacks. The bottom corners `c0`, `c1` are bearing positions too
(they are the end posts); their springs act on the exterior DOFs directly.

At each position, between the plate node and the interpolated base point:

- **vertical** `kv` (force/length, per position): the series crushing stiffness of the sandwich under one stud
  (plates, subfloor, rim or blocking, top plate below, perpendicular to grain). Bilateral, as the element already
  assumes for interior studs.
- **horizontal** `kh` (force/length, per position): sill nailing and friction. This replaces both the
  foundation's internal `u` fix and the `equal_dof` plan tie between storeys. It is a penalty: give it a large
  value for the conventional rigid tie.

Bottom stations stay exterior nodes for other connections; with `bearing` they get the springs too, unlike
`foundation`, where a station replaces the internal support. The bottom plate under a window has no bearing
position between the jacks (nothing stands on it).

### D5. Hold-downs are inside the element too

`hold_down: Option<Material>` (an `Anchorage` material, or any uniaxial material) replaces `kv` at every pier end
(`layout.pier_ends()`), between the end-post bottom node and the interpolated base point. Deformation is the
upper post's rise relative to the surface (`v_p − v_b`, uplift positive), the convention of the external
anchorage link, so tension is positive. With `None` the pier ends bear like any other stud (`kv`).

A hold-down is nonlinear (tension and bearing differ). It becomes a component of the element's existing
nonlinear internal solve, which so far held only the piers' shear fields: each component has a strain vector
`β` over the element DOFs, a law `f(e)` and its own state, and the solve finds `e + C f(e) = y` with
`C = βᵢᵀ K_ii⁻¹ β` (the same algebra, with area 1 and no own compliance). Linear sheathing walls keep their
piers in the condensed stiffness, and only the hold-downs are components. Consequence: a wall with hold-downs
needs an iterative algorithm (Newton), as the rowhouse already uses.

`MAX_WALL_PIERS` becomes the limit on components (piers plus hold-downs) and rises from 16 to 32.

Rejected: keeping hold-downs as external `ZeroLength`s at bottom stations. The `ZeroLength` needs a node below at
the same `x`, which brings back the alignment problem for the posts. A hold-down rod or strap actually anchors
to whatever is below the post, which the interpolated point represents.

### D6. Responses

- `base_nodes`: one entry per bearing position (`x`, `fx`, `fy`: force on the wall, gravity positive), from the
  spring forces (and the hold-down force at pier ends). `base` sums them. Corner and station external forces are
  not base reactions in this mode.
- `hold_downs: Vec<WallHoldDownResponse { x, deformation, force }>`: tension positive.
- The base nodes receive the equal and opposite forces through the condensed stiffness, so the wall below is
  loaded exactly.

### D7. Foundation stays

`foundation: true` keeps its meaning (rigid internal supports, external hold-downs). `bearing` and `foundation`
are exclusive. Bearing on two fixed ground nodes at `x = 0` and `x = L` is rigid ground with internal
hold-downs: it reproduces `foundation: true` plus pier-end anchorage links (verified, see below), so a
preprocessor can use `bearing` for every storey.

### D8. Wire format

`WallLayoutSpec.bearing: Option<WallBearingSpec { kv, kh, holdDown: Option<AnchorageSpec> }>` (the building's
anchorage spec, inline) and `WallChordSpec::Base`. The building expansion rejects a layout's own `bearing`, as it
rejects its own `floor`; it keeps generating links.

## PySees (consumer)

- A stud wall stores its bearing (`kv`, `kh`) and hold-down law in its args; `drawStudWall` no longer writes
  ground nodes, `zeroLength` links or `equalDOF`s.
- `compileInputV1` finds each wall's base nodes when it compiles: the nodes at the wall's bottom elevation on its
  line (other walls' top rows, beam-column nodes), keeping those within `0..L` and the nearest one beyond each
  end. A wall at the lowest stud-wall elevation stands on two generated fixed nodes. A wall end with no base
  node beyond it is a compile error naming the wall. A gap between base positions longer than three of the wall's
  stud spaces is a warning (the interpolation treats the sandwich as rigid across it).
- Walls saved with explicit links (older projects) keep compiling as before; bearing applies to walls whose
  args carry it.

## Verification

| Check | Tolerance |
|---|---|
| Aligned stack (wall above on the floor-zone top row of the wall below, studs aligned): bearing vs bottom stations at every bearing position joined by `ZeroLength`s with the same springs, gravity and lateral | 1e-8 |
| Bearing on two fixed nodes (`kv = kh = 1e12`) vs `foundation: true` with pier-end anchorage links | 1e-6 |
| Offset studs: base reactions equal the lever-rule split of each stud force | 1e-9 |
| Coincident base nodes share equally | 1e-12 |
| Equilibrium: base reactions sum to the applied top loads | 1e-9 |
| Hold-down in tension: force `tension·uplift` | 1e-9 |
| Position outside the base range, bearing with foundation, base station without bearing: named errors | exact |

Tests: `core/tests/elements/stud_wall_bearing.rs`.
