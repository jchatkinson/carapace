# Floors for stud-wall buildings: plan

Status: Phases 1 (floor band), 2 (one-way gravity transfer), 3 (pier-end anchorage material) and 4a (walls and
bands in 3D) are implemented in the core Rust API; 4b (diaphragms) is implemented; the builder remains. User guide: [floors.md](floors.md). Per-position bearing and
anchorage links of the band are
external `ZeroLength` elements (the initial plan's internal vertical links are deliberately not part of the
band). Companion to `stud-wall.md` (user guide) and `stud-wall-v1.md` (wall spec).

## 1. What a floor has to do

In a light-frame building the floor does three jobs, and the stud wall needs all three:

1. **Gravity, one way.** Joists span one way between bearing lines (walls, beams, headers). Each support gets
   its tributary share, which then runs down through the wall framework: top plate → studs, or top plate →
   cripples → header → jacks over an opening. Designers expect hand-calculation reactions here.
2. **Floor depth.** The joist and rim depth `d` sits between the storeys. It belongs in the overturning lever and
   the storey height (O86 Fig. 11.7: `h` = wall height + joist depth), and it is where bearing, crushing,
   shrinkage and the hold-down connection between storeys live. Without it, hold-down forces and drift are low
   (about 10% for a 254 mm floor on a 2.44 m wall).
3. **Diaphragm.** The sheathing collects the lateral load and delivers it to the wall lines: flexible (by
   tributary area), rigid (by stiffness, with torsion) or semi-rigid (by an in-plane shear stiffness). The
   FPInnovations guide reports up to 60% difference in wall-line loads between these idealizations, so
   designers bracket them. The rim joists act as chords and collectors.

These pieces depend on each other in one direction: gravity transfer needs only plan geometry; the floor depth
lives inside a wall line; the diaphragm needs more than one wall line at once.

## 2. The key decision: wall lines or a building

The `StudWall` is a 2D element: one wall line in elevation. A 2D model can represent jobs 1 and 2 completely.
Job 3 needs plan geometry: several wall lines in different vertical planes connected by the floor.

| | 2D wall lines + distribution step | 3D building |
| --- | --- | --- |
| How the diaphragm works | A separate step splits each storey's shear among the wall lines (tributary, or by line stiffness plus torsion), then each line is solved in 2D | Walls sit in vertical planes of a `Domain3`; the floor ties their tops; distribution, torsion and interaction come out of one solve |
| Rigid diaphragm with nonlinear O86 walls | Iterate the line stiffnesses until the split converges | Automatic (`Domain3::rigid_diaphragm` exists) |
| Gravity | Same one-way transfer (Section 3) | Same |
| Engine work | Small | `StudWall` made generic over 2D/3D (Section 6) |
| Matches | WoodWorks Shearwalls | ETABS/SAP-style building models |

**Recommendation:** build it in phases so that every phase is useful in 2D and the parts carry into 3D
unchanged:

- **Phase 1:** the floor band, which fixes floor depth (Section 4).
- **Phase 2:** one-way gravity transfer (Section 3).
- **Phase 3:** a pier-end anchorage law: crushing, shrinkage, tension-only hold-downs (Section 5).
- **Phase 4:** `StudWall` in 3D and the diaphragm (Section 6).

Phases 1–3 close the WoodWorks gaps A and B for single wall lines. Phase 4 adds the building.

## 3. Gravity: one-way joist transfer (Phase 2)

### Approach: transfer loads, don't model joist stiffness

Joists are treated as simple spans (optionally continuous over supports that don't deflect), and their
reactions are applied as loads on the supporting members. This is how designers do it, and it keeps takedowns
equal to hand calculations: the floor's share of load reaching a header is what a designer would compute,
whatever the header's stiffness.

The alternative, joists as beams on flexible supports, sends load toward stiff supports and away from
deflecting headers. Its results stop matching hand calculations, so it is not the first choice.

Because the transfer is statically determinate, it needs no element stiffness. `carapace-core` provides a pure
transfer function in `model::floor`; a future builder/input layer can resolve its reactions into existing loads:

- onto walls as `WallLoad::TopLine` (trapezoidal, partial) and `WallLoad::TopPoint`, in the wall's own load
  slot;
- onto beams and girders as their existing element loads.

From there, the wall and header framework already does the rest: a `TopLine` over an opening goes top plate →
cripples → header → jacks.

### `Floor` input

- A plan outline (rectangles first, then polygons).
- The span direction.
- Supports: wall-line segments (wall id and local `x` range), beams, and free edges (cantilevers).
- Area loads per slot (D, L, S…), plus line and point loads on the floor (partitions, a post from above).

### Algorithm

1. Cut the floor into strips parallel to the span. Each strip's edges are where some support starts or ends,
   so the loads come out piecewise linear.
2. Find the supports each strip crosses, in order.
3. For `Simple` continuity, each span shares its load to its two ends. A cantilever past the last support
   loads it and relieves the back span, which may lift the back support (reported, not hidden). For
   `Continuous`, solve each strip as a continuous beam on rigid supports, a tiny 1D solve.
4. Emit a `TopLine` or beam load for each support and strip.
5. Report every support's tributary width and reaction per slot, so the takedown can be checked by hand.

### Girders

A girder (a beam element supporting joists) receives line loads from the transfer. Its end reactions reach the
wall through the analysis, at a top-plate station where it bears: point loads go to the station, a post below
it, or a header.

### Platform framing

Floor joists bear on the top plate of the wall below. The wall above bears on the subfloor, and its stud loads
pass through the floor band (Section 4) to the wall below.

### Tests

Each compares against a hand calculation, at 1e-9:

- two-span joists, simple and continuous;
- a cantilever with uplift on the back span;
- joists over a window header (header reaction and jack loads);
- a girder with a post below;
- a two-storey takedown through the band.

### Phase 2 implementation and deviations

The current core API supports rectangular floor bounds, area patches, cross-joist line loads, joist-aligned
line loads and point loads, with eight wall load slots. Supports are clipped to the floor width. Reactions can
be mapped to `WallLoad`s or returned as beam segments. The transfer is not on the wire format because the
stud-wall builder and its identifiers are not on the wire format yet. Polygon outlines, automatic support
discovery remain future builder work. The end-to-end test maps D/L floor reactions to a windowed stud wall,
applies a factored combination, and checks the header jack and base forces against the hand takedown. The API
is usable before a builder can discover supports automatically.

## 4. Floor depth: the floor band (Phase 1)

### What it is

`FloorBand` is a 2D element between the top plate of the wall below and the bottom plate of the wall above.
It is a strip of depth `d` along the wall line.

- **Nodes.**
  - a bottom row (the lower wall's top corners and stations);
  - a top row (the upper wall's bottom corners and stations);
  - both rows at the same positions `x_i` (the union of both walls' stations).

  Exterior DOFs are `ux, uy` only, like the wall.
- **Shear across the depth.** A constant shear field over each bay between posts, with stiffness from the rim
  joist and blocking (`GA_band`, rigid by default). The lateral load enters at the top row (the sheathing
  level), so the overturning moment at the lower wall's base comes out as `V·(H + d)`: the O86 lever, with no
  post-processing.
- **Vertical link at each position.** This is bearing: the rim, plates and squash blocks perpendicular to grain.
  It is linear in Phase 1 and gets the anchorage law in Phase 3.
- **Rim along the line.** The rim's `EA` acts as the chord or collector along the floor edge.

### Reuse, not a new frame

The band is the stud wall's own frame:

- the vertical links play the studs;
- the rim plays the plates;
- the shear field is the sheathing field.

The implemented `FloorBand` uses direct constant-shear bay stiffness and axial rim bars. There are no internal
vertical springs or condensed DOFs: users add bearing and hold-down `ZeroLength` elements between corresponding
bottom and top nodes. Omitting those links leaves a mechanism. This keeps anchorage laws available to evolve
independently in Phase 3.

### Results

- The base moment and M/jd at the band level, for hold-downs between storeys.
- Storey drift measured over `H + d`.
- Bearing forces at each post.

A small helper reports drift per storey from the node displacements, so it matches the code's definition of
storey height.

### Hold-downs between storeys

A hold-down between storeys runs from the upper wall's pier-end post, through the band, to the lower wall's
post. In Phase 1 that is a `ZeroLength` between the band's top and bottom nodes at that position, in parallel
with the bearing link. In Phase 3 the two merge into one anchorage law.

### Tests

- One storey on a band: hold-down force `V(H + d)/jd`, exact.
- Two and five storeys against O86 A.11.7.1 with `h = H + d`, at the existing tolerances.
- A band with a rigid shear field equals a hand statics lever.
- A takedown through the band, with bearing at each stud, equals a hand calculation.

## 5. Pier-end anchorage law (Phase 3)

Implemented as `Material::Anchorage`, one uniaxial law used wherever a pier end meets what is below it: the
foundation, a band, or the storey below. It replaces separate plain `ZeroLength` hold-down and bearing springs.

- **Tension:** hold-down elongation, either the manufacturer's `(T/T_cap)·Δ_max` or a stiffness, plus rod
  stretch `TL/(A·E)` and fastener slip.
- **Compression:** bottom-plate crushing, a nonlinear function of `C/Q_r`, in the same form WoodWorks uses.
- **Shrinkage:** a constant slack, `0.002·ΔMC·L_shrink`, which closes before the hold-down engages. It is
  reported separately because it is not amplified by `RdRo/Ie` for seismic drift.

The material is available for explicit `ZeroLength` links now. Automatic pier-end placement (`anchorage` on a
wall layout generating stations and links at `pier_ends()`) moves to the Phase 4 builder. The material has no
OpenSees counterpart and is not added to the Python API or `input_v1` wire format in this phase; wire exposure
can be considered when the wall/builder input format is designed.

Tests:

- tension only, with load reversal;
- crushing and shrinkage each checked against the O86 11.7.1.2 `d_a` hand values in the WoodWorks sample
  reports.

## 6. The building in 3D (Phase 4)

### `StudWall` in a vertical plane of a `Domain3`

The wall's frame stays 2D. Only the mapping changes: local `(x, y)` → global `(along, up)`. In 3D the up axis
is `Y`, as for `rigid_diaphragm`, and `along` is a unit vector in the `X–Z` plane.

- `dof_mask` declares exactly the global translations the wall stiffens: `uy`, plus `ux` and/or `uz` for its
  direction.
- Out-of-plane stiffness is zero, as in design. A node used only by one wall line needs out-of-plane support
  from a crossing wall or the diaphragm.

This follows the rule "share assembly, split leaf kinematics": one frame builder, two thin node mappings.

Test: walls along `X`, along `Z` and at 30° each equal the 2D wall exactly.

**Phase 4a status: done.** One implementation for 2D and 3D: `StudWall<K>`/`FloorBand<K>` with a shared
`WallPlane` map (`core/src/model/elements/wall_plane.rs`), 3D wall-slot loads (`ElementLoad3::wall`) and `Domain3`
response accessors. `core/tests/elements/stud_wall_3d.rs` checks 3D against 2D: walls along X, Z and 30°, the
window wall with pier-end hold-downs, a gravity takedown, the O86 nonlinear law, the floor-band stack (1e-7,
stand-in limited) and an L corner sharing its post between planes, plus DOF numbering (no `uz` for an X wall).

### Diaphragm options

- **Rigid:** `Domain3::rigid_diaphragm` on the band-top nodes of a level. It gives distribution by stiffness,
  with torsion. Accidental torsion is a load case (an offset storey force).
- **Flexible:** tributary distribution. This is the gravity transfer's algorithm turned sideways. It applies to
  each wall line's band top as loads, with no stiffness.
- **Semi-rigid:** a `Diaphragm` element. It is a stud-wall frame lying flat, with these pieces:
  - the sheathing field over each bay;
  - the rim joists as chords and collectors;
  - the joists as the transverse members.

  It uses the same O86 nail-slip law, and is checked against the O86 diaphragm deflection equation (to be read
  from the code PDF when this is built). It needs a 2D-in-3D mapping like the wall, in the horizontal plane.

**Phase 4b status: done.** `Domain3` now has a 3D semi-rigid `Diaphragm` element, rigid-diaphragm documentation and flexible load-transfer helpers. The element is axis-aligned in global X-Z at a specified Y elevation.

### A `Storey`/`Building` builder in the input layer

Users won't wire stations by hand. Given walls per storey, floors and bands, the builder should:

- place the stations at every stud and pier end;
- share or link the nodes;
- place the bands and anchorage;
- run the gravity transfer;
- attach the diaphragm.

Tests:

- a rigid diaphragm over two unequal wall lines, with an eccentric load, against a hand stiffness-and-torsion
  split;
- a flexible diaphragm against tributary;
- a two-storey box building's takedowns against hand calculations.

## 7. Wire format and PySees

Each wire-visible phase extends `input_v1` and must be mirrored in PySees (`carapaceInputV1.ts`, `compileInputV1.ts`),
`docs/input-format.md`, the CLI schema and the wasm build:

- Phase 1: the band;
- Phase 2: `Floor`;
- Phase 3: the anchorage law (core only; no OpenSees material-name or wire counterpart);
- Phase 4: 3D walls, the diaphragm and the builder.

The stud wall itself is not on the wire yet, so that comes first, or together with Phase 1.

## 8. Open decisions

1. 2D wall lines with a distribution step, or a 3D building for the diaphragm. Recommended: phases 1–3 in 2D,
   then 3D.
2. Gravity by load transfer (hand-calculation reactions) or by joist stiffness. Recommended: transfer, with
   `Simple` and `Continuous` options.
3. Which diaphragm idealizations, and in what order. Recommended: rigid and flexible first, because they bracket
   the result; semi-rigid later.
4. Whether `FloorBand` is a separate element or an option on the wall above (a `floor_depth` below its bottom
   plate). Recommended: a separate element, because the band is shared by two walls and later by the
   diaphragm.
