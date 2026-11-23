# Stacked stud-wall buildings

`carapace-input::building` expands a native Rust `BuildingSpec` (plan X-Z with Y up by default, N/mm; camelCase serde
form, which `carapace expand` and `expandBuilding` read) into an ordinary explicit model: the same tables a hand-written
`CarapaceInputV1` has (nodes, materials, `studWalls`, `zeroLengths`, `equalDofs`, `rigidDiaphragms`, `linearConstraints`,
`loadPatterns`, `nodalLoads`, `elementLoads`, `floorLoads`) plus a map of where things went. The wire decoder has no
building table and no building recorders, so a client can show, edit and re-export the generated model like any other.

```rust,ignore
use carapace_input::building::{expand, ExpandOptions, UpAxis};
let Expansion { mut input, map } = expand(&spec, &ExpandOptions { up: UpAxis::Y })?;
input.sequence = my_sequence; // empty on return; the load rows are at stage 0
let mut session = carapace_input::input_v1::decode(input)?;
```

`build(&spec)` is `expand` (Y up) plus a decode, returning the native `Building`: the `Domain3` with every load
registered and handles (`walls`, `levels`, `links`, `floor_reactions`) for `wall_response`, `storey_drift`,
`holddown_forces` and `apply_combination`. The same code makes both, so the native tests below exercise the expansion.

## Expansion

`expand(spec, ExpandOptions { up }) -> Result<Expansion, BuildError>`; `Expansion { input, map }` (camelCase, tsify types in
wasm). From the command line `carapace expand spec.json [--up y|z] [--sequence seq.json] [--model-only]` prints it;
in wasm `expandBuilding(spec, up)` returns it (errors reject with `{ kind: "invalidBuilding", reason }`).

`input` has `header.ndm = 3` and these tables, all explicit:

| Table | Content |
| --- | --- |
| `nodes` | Corners, top stations (every stud and pier end), bottom stations, one fixed ground node (all three translations) per foundation pier end, the diaphragm masters; coincident nodes are not merged |
| `materials` | `elastic` 1e12 for the ground links' horizontal DOFs, `anchorage` laws (hold-down at pier ends, bearing elsewhere), deduplicated |
| `studWalls` | One row per wall, storey by storey in spec order; `layout.foundation` on storey 0, `layout.floor` from the storey's `floorDepth`, `bandGa`, `rimEa` |
| `zeroLengths` | Ground-to-pier-end anchorage links (storey 0), then the vertical link between each wall and the wall below (anchorage at pier ends, compression-only bearing elsewhere) |
| `equalDofs` | Each upper bottom node to the lower top node in the wall-plane translations |
| `rigidDiaphragms` | Rigid levels: one retained master at the centroid of the level nodes, with `normal` = the vertical axis |
| `linearConstraints`, `nodes.fixed` | Flexible levels: the plan-normal translation of every top node is fixed (wall along a plan axis) or tied to the in-plane one (skew wall) |
| `loadPatterns` | One `linear` pattern per used load slot (ascending), `scaleFactor = combination[slot]` |
| `nodalLoads` | Rigid levels: the lateral force and its moment on the master; flexible: lumped line reactions on the top stations (all at stage 0) |
| `elementLoads` | The spec walls' `loads` as `wallLoad` rows of their slot's pattern |
| `floorLoads` | Each floor as one row per slot it loads (`pattern` = that slot's pattern), supports given as `studWalls` rows, so the decoder does the transfer |

`input.sequence` is empty. Every load row has `stage: 0`, so append a stage (and recorders) before decoding; without one
the decoder reports `unknownStageIndex` for the first load row.

`map: ExpansionMap { storeys, slotPatterns }`, all indices into the tables above:

- `storeys[s].walls[i]`: `name`, `element` (`studWalls` row), `corners` (bottom start, bottom end, top end, top start),
  `topNodes` / `bottomNodes` (`{ x, node }`, wall-local x) and `pierLinks` (`{ x, element }`, the `zeroLengths` row of the
  anchorage at that pier end: to the ground on storey 0, to the wall below above it);
- `storeys[s].linksBelow`: the `zeroLengths` rows joining the storey's walls to the storey below (empty on storey 0);
- `storeys[s].master`: the rigid level's retained node (`null` for flexible levels); `storeys[s].nodes`: the level nodes;
- `slotPatterns[k]`: the `loadPatterns` row of load slot `k`, or `null` if no load uses it.

Results come from the generic recorders on those rows: `studWallPier` and `studWallBase` with `elementKind: "studWall"`
and `elementIndex = element`; storey drift from `nodeDifference` with `node = storeys[s].master` and
`otherNode = storeys[s-1].master` (omitted on storey 0), the DOFs below; hold-down forces from `elementForce`
(component 1, the vertical) on the `pierLinks` rows. Flexible levels have no master: difference the wall's top-row and
bottom-row `nodeDisp` instead.

### Up axis

`UpAxis::Y` (the default, and what `build` and the native helpers use): plan `(x, z)` maps to model `(x, z)`,
elevations go along `+Y`, the diaphragm normal is `Y`, a vertical link acts on DOF 1 and a level's drift is DOFs
`ux, uz, ry` = 0, 2, 4. `UpAxis::Z`: plan `(x, z)` maps to model `(x, y)` (the spec's second plan coordinate is the
model's Y), elevations go along `+Z`, the normal is `Z`, vertical links act on DOF 2 and the drift DOFs are `ux, uy, rz` =
0, 1, 5. Node numbering, rows and the map are identical in both; only coordinates and DOF indices change.

Handedness: the plan frame `(x, z)` with `y` up is the mirror of `(x, y)` with `z` up (Y-up: `x` to `z` turns clockwise seen
from above; Z-up: `x` to `y` counter-clockwise). A rotation of the level about the vertical axis therefore changes sign:
for the same physical motion `rz = -ry`, and the generated moment of a lateral force about the master is
`M_y = f_x (z - z_m) - f_z (x - x_m)` about `+Y` and `M_z = -M_y` about `+Z`. Translations, shears, pier forces and drift
magnitudes are identical. Tested (`input/tests/building/expansion.rs`): the same building expanded with `Y` and `Z` agrees
to 1e-10, with `rz = -ry`.

## Generated model

Each `WallSpec` defines a vertical wall plane by a plan start, unit direction and local `WallLayoutSpec`.
Stations are generated at stud rows and pier ends. Storey zero has `foundation: true`; interior pier ends are
connected to fixed ground nodes with the selected anchorage law. At an upper storey, matching named walls must
have identical plan start, direction and length. A positive `floorDepth` gives the storey's walls a
[floor zone](stud-wall.md#2a-floor-zone) (`depth = floorDepth`, `ga = bandGa`, `rimEa`), including the roof storey:
the wall spans storey grid to storey grid, its exterior top row (a station at every stud) lies at `H + floorDepth`,
and the wall above starts on new nodes coincident with it. Each pair of coincident nodes is tied in the wall plane
(`equal_dof`) and joined by a vertical `ZeroLength` link: the anchorage law at pier ends, linear bearing
(`bearing`, compression only) elsewhere. A zero depth links the walls the same way with no zone. A roof zone
needs no links. A wall layout's own `floor` is rejected: the zone comes from the storey. The wall above has bottom
stations at every stud of the wall below; its pier ends must be among them (`BuildError` otherwise).

Rigid diaphragms retain one node at the centroid of the level nodes. A lateral load's force is applied at that
node, and its plan offset becomes the moment about the vertical (Y-up: about `+Y`, `M_y = f_x (z - z_m) - f_z (x - x_m)`,
`r × F`), the convention of `Domain3::rigid_diaphragm` (a rotation `θ` moves a point by `θ dz` in x and `-θ dx` in z).
Flexible diaphragms span between the wall lines that resist the load (lines parallel to it), distribute it by
the floor transfer statics, and lump each reaction onto the wall's exterior top stations. Nodal loads are summed per
node and DOF before they are added to the domain (`add_nodal_load` replaces an existing value for the same node
and DOF). Unsupported plan-normal motion is restrained at the level nodes. Semi-rigid diaphragms are not
selected through this builder; add a core `Diaphragm` explicitly.

`FloorSpec` supports rectangular plan bounds `[x0, x1, z0, z1]`, X or Z span (load coordinates `s` run along the span and `t` across it, so for a Z span `s` is z and `t` is x), named bearing walls perpendicular to the span (joists bear
on the top plate of this storey's walls: platform framing), area,
line and point loads, and simple or continuous transfer. Reactions remain available as `Building::floor_reactions`
and are also mapped into the supporting wall's load cases. The expansion makes one constant load pattern per slot
that carries a load (wall, floor or lateral), and each wall has one load case per slot that loads it (the decoder
groups a pattern's `wallLoad` rows and floor reactions on a wall); a slot's pattern selects that case on every wall.
The spec keeps `slot` on its loads and `combination: [f64; 8]` (default all ones), which becomes the patterns' scale
factors. `Building::apply_combination([f64; 8])` resets pattern `k`'s scale factor to `combination[k]` (slots without a
load are ignored) before constructing a native analysis, scaling D/L/W or other combinations independently.

## Assumptions and limits

- Walls in successive storeys stack at the same line and length. Walls bearing on a floor between supports are
  rejected; there is no automatic collector, transfer beam or post path.
- Floor bounds are rectangles, supports are full wall lines, and a floor has one axis-aligned span direction.
- A storey has one diaphragm idealization, shared by all its levels in this phase. Flexible distribution is a
  tributary statics model; use a core semi-rigid diaphragm element when floor stiffness matters.
- Pier ends receive anchorage links. Other stations receive linear bearing (compression only).
- `build` and the native helpers (`storey_drift`, `holddown_forces`) are Y-up; for Z-up use `expand` and the map.
- Plan-normal translation at flexible diaphragm stations is restrained. Crossing wall lines are not merged into
  common nodes by this builder.
- Loads and material parameters must be finite; invalid names, slots, wall alignment and support geometry return
  `BuildError`.

## Verification

`input/tests/building/builder.rs`, against hand statics:

| Case | Reference | Tolerance |
| --- | --- | ---: |
| One storey, one X wall + two Z returns, rigid: X wall carries the whole storey shear | statics | 1e-9 |
| Symmetric two-storey box with a 254 mm floor zone: each X wall takes half of each storey; lower base moment `5000(2H + d) + 5000(H + d)` | statics, O86 Fig. 11.7 levers | 1e-8 |
| Eccentric X force, rigid, four separate equal walls: shears `F/2 ∓ M d/J`, `+Y` rotation | rigid-body stiffness + torsion split | 1e-8 |
| Flexible, three parallel lines under uniform load: 3/6/3 kN | tributary reactions | 1e-9 |
| Floor D + L at 1.25D + 1.5L to a window wall: base 32e6 N, header 4.8e6 N per end from the cripple tributaries | hand takedown | 1e-9 |
| Misaligned stacking, unknown/non-perpendicular supports, bad direction | named `BuildError` | exact |

`input/tests/building/expansion.rs` runs the expansion through the wire decoder:

| Case | Reference | Tolerance |
| --- | --- | ---: |
| Two-storey eccentric rigid box with floor zone and gravity floor, and a flexible three-line building: expand, serialize, parse, decode, solve: wall pier q, base moment and reactions, storey drift | the native `build` solution | 1e-12 |
| Eccentric box, torsion case, flexible case with `up = Z` against `Y`: reactions, pier forces, drift; `rz = -ry` | the Y-up expansion | 1e-10 |
| `studWallPier`, `studWallBase`, hold-down `elementForce` and master `nodeDifference` recorders addressed through the map | `Building::wall_response`, `holddown_forces`, `storey_drift` | 1e-12 |
| Floor load rows: front/back base axial `w Lx Lz / 2` | statics | 1e-9 |
| A `building` field in a model fails to parse; tables, pattern scale factors, normals and ground fixities of the expansion; bad slots and factors | exact | exact |

(The wire round trip needs `serde_json`'s `float_roundtrip`: its default float parser can move a coordinate by an ulp, which
the 1e12 link stiffnesses amplify to 1e-11.) `cli/tests/cli.rs` checks `carapace expand` on the example, bad specs, and that the
committed `examples/models/two_storey_box_building.json` is what `expand` generates; `wasm-bridge/tests/boundary-smoke.ts`
expands for both axes and runs the results.

## Regenerating the example model

`examples/models/two_storey_box_building.spec.json` is the building spec, `two_storey_box_building.sequence.json` the stage and
recorders (rows 4 = the upper front wall, nodes 53/52 = the two masters, from the map), and
`two_storey_box_building.json` the runnable model they expand to:

```bash
cargo run -q -p carapace-cli --release -- expand examples/models/two_storey_box_building.spec.json \
  --sequence examples/models/two_storey_box_building.sequence.json --model-only --pretty \
  -o examples/models/two_storey_box_building.json
```

Run it after changing the spec, the sequence or the generator (the CLI test fails on a stale file).
