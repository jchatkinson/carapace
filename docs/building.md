# Stacked stud-wall buildings

`carapace-input::building` builds a 3D `Domain3` from a native Rust `BuildingSpec`. The spec uses global Y up,
plan X-Z, and N/mm. Its camelCase serde representation is also accepted as the optional `building` field in
`CarapaceInputV1` for the 3D profile.

```rust,ignore
let building = carapace_input::building::build(&spec)?;
let mut domain = building.domain;
```

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
node, and its plan offset becomes the moment about `+Y`, `M_y = f_x (z - z_m) - f_z (x - x_m)` (`r × F`), the
convention of `Domain3::rigid_diaphragm` (a rotation `θ` moves a point by `θ dz` in x and `-θ dx` in z).
Flexible diaphragms span between the wall lines that resist the load (lines parallel to it), distribute it by
the floor transfer statics, and lump each reaction onto the wall's exterior top stations. Nodal loads are summed per
node and DOF before they are added to the domain (`add_nodal_load` replaces an existing value for the same node
and DOF). Unsupported plan-normal motion is restrained at the level nodes. Semi-rigid diaphragms are not
selected through this builder; add a core `Diaphragm` explicitly.

`FloorSpec` supports rectangular plan bounds `[x0, x1, z0, z1]`, X or Z span (load coordinates `s` run along the span and `t` across it, so for a Z span `s` is z and `t` is x), named bearing walls perpendicular to the span (joists bear
on the top plate of this storey's walls: platform framing), area,
line and point loads, and simple or continuous transfer. Reactions remain available as `Building::floor_reactions`
and are also mapped into the supporting wall's load cases. The builder makes one constant load pattern per slot
that carries a load (wall, floor or lateral), and each wall has one load case per slot that loads it; a slot's
pattern selects that case on every wall. The wire format keeps `slot` on the building's loads and
`combination: [f64; 8]`. Call `Building::apply_combination([f64; 8])` before constructing an analysis: it sets
pattern `k`'s scale factor to `combination[k]` (slots without a load are ignored), scaling D/L/W or other
combinations independently.

## Assumptions and limits

- Walls in successive storeys stack at the same line and length. Walls bearing on a floor between supports are
  rejected; there is no automatic collector, transfer beam or post path.
- Floor bounds are rectangles, supports are full wall lines, and a floor has one axis-aligned span direction.
- A storey has one diaphragm idealization, shared by all its levels in this phase. Flexible distribution is a
  tributary statics model; use a core semi-rigid diaphragm element when floor stiffness matters.
- Pier ends receive anchorage links. Other stations receive linear bearing (compression only).
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

The wire decoder builds this same model, applies `combination`, and maps `buildingWallPier`, `buildingWallBase` and rigid-level
`buildingStoreyDrift` recorders. `buildingStoreyDrift` is rejected for flexible levels because they have no master
node. PySees mirrors only the tables its compiler emits; it has no stud-wall API and needs no change for these
tables.
