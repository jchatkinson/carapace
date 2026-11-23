# Plan: stud walls as ordinary elements (Carapace and PySees)

Status: C1–C4 complete; P1–P5 proposed, 2026-11-23. Supersedes the PySees part of `floor-system-plan.md` and the `building`-as-model
approach.

## Why

The first PySees hookup stored a stud-wall building as one opaque `building` entity that Carapace expanded at
decode time. Nothing in it could be inspected or edited except through the template form, and the model explorer,
deformed shape and node/element results had nothing to work with. A stud wall must be modelled like any other
element: nodes on the storey grid, elements between them, supports, constraints, materials and load patterns, all
of which the user can edit after a template has generated them.

## Decisions

1. **The floor band folds into the wall.** A wall spans storey grid to storey grid. Its layout gains a floor
   zone of depth `d` at the top: sheathing, studs and chords stay over the framed height `H`; the zone carries the
   storey force down to the framing as the band does now (couple into the chords, O86 Fig. 11.7 levers `H + d`),
   and the rim's EA adds to the top chord (one drag strut over openings, which also settles the rim + top plate
   coupling question). The zone keeps the band's shear stiffness (`floor_ga`) and its vertical load path is rigid posts at the top stations; crushing and hold-downs stay in the links between the walls. The `FloorBand` element and its tables
   are removed.
2. **Geometry defines the wall plane.** A wall's local x runs `c0 → c1` and local y `c0 → c3`, in 2D and 3D, with no
   up axis assumed, so Z-up PySees models and Y-up models both work. The `plane` field goes. The semi-rigid
   diaphragm likewise takes its plane from its nodes.
3. **Wall loads belong to load patterns**, like any element load. No slot cap.
4. **One-way floor loads are a pattern load** that names its support walls and is expanded at decode time (as
   shell pressure is in PySees).
5. **Building generation is an expansion, not a model.** The builder outputs an explicit wire model (nodes,
   walls, links, materials, constraints, patterns), exposed natively, in the CLI and in wasm. The decoder no longer
   accepts a `building` table.
6. **OpenSees export** writes everything else and a warning comment in place of each Carapace-only element or
   material (stud wall, diaphragm, anchorage), as zeroLength friction does.
7. **PySees** drops the building entity, the Building card and the `building` project directive (never
   committed); it keeps what carries over (wire types, the wall drawing, the results tables, recorder plumbing).

## On the wire tables

Each element formulation has its own table (`trusses`, `elasticBeamColumns2d`, `shell4s`, ...). That is the
right shape: decoding is a typed, closed set (it maps onto the closed `Element`/`Element3` enums), each table
validates its own columns, the TypeScript types are exact, and columnar tables are compact for thousands of simple
elements. The stud-wall tables follow it, with two deviations and one real mistake:

- **Row objects instead of columns** (`studWalls*`, like `embeddedBars`): justified. A wall's layout is nested
  and variable-length (studs, openings, stations); columns would need offset arrays like `fibers`, which is worse
  to write, read and validate.
- **Twin 2D/3D tables**: not justified here. Beams split because their columns differ (3D adds `g`, `j`, `iy`,
  `vecXz`). A wall row is identical in both profiles apart from `plane`, which decision 2 removes. Merge into one
  `studWalls` table whose profile comes from `header.ndm`, and one element kind `studWall`.
- **Generator tables in the decoder** (`building`, and `floorLoads` as an unowned table): the mistake. They mix
  describing a model with generating one, hide the generated model from the client (the PySees problem), and needed
  an exclusivity rule (`buildingWithTables`), which is the symptom. Generation moves to an explicit expansion
  (decision 5); floor loads become a pattern-owned load (decision 4), which is description, not generation.

## Carapace phases (one commit each)

**C1. Plane from geometry — done.** Local x comes from `c0 → c1`, local y from `c0 → c3`
(orthonormal within 1e-9 relative, else an error). The DOF mask takes the global translation slots with nonzero
components in either axis. `StudWall::new` and `new3` lose their plane argument; the wire `plane` field and the
2D/3D twin tables merge into `studWalls` and kind `studWall` (recorders and element loads follow). `Diaphragm` takes
its normal from its grid nodes (axis-aligned grids in any of the three planes). The builder and every existing test
keep passing; new tests: the same wall in a Y-up and a Z-up model, and a skew wall, agree to 1e-10; a Z-up diaphragm
equals the Y-up one.

**C2. Floor zone in the wall; FloorBand removed — done.** `StudWallLayout { floor: Option<FloorZone { depth, ga, rim_ea }>, .. }`: exterior top
corners and stations sit at `H + d`; the zone is condensed internally. Verification: equal to today's wall +
FloorBand with a stiff band (1e-9); the 2- and 5-storey stacks and the O86 Fig. 11.7 hand checks
(`floor_band.rs`) pass on the new form; drift and hold-down results of the guide example unchanged. Remove
`FloorBand`, its tables, `floorBandBay`, docs.

*Outcome.* Implemented as briefed (`FloorZone`, `floor_zone.rs` bay math shared by the wall, rigid posts by
sharing the exterior node's `v` DOF with the framing node, builder on coincident nodes + `ZeroLength3` links +
`equal_dof`). Top slot loads: the horizontal part goes on the exterior row by the lever rule, the vertical part
stays on the top plate (a vertical force has the same line of action either way, so gravity hand values are
unchanged). The equivalence to the old band model holds to 1e-8 (stacks 1e-6): the recorded references carry
that model's own conditioning noise (its statics residual is 3e-8).

**C3. Pattern-owned wall loads (done).** The wall keeps an internal list of load cases (each a prepared load vector);
`ElementLoad` carries `wall: Option<(case, factor)>` (Copy, no cap) instead of `[f64; 8]`. Wire: an element load
`{ kind: "wallLoad", load: WallLoadSpec }` in any pattern; the decoder gives each pattern's loads on a wall their
own case. Floor loads: `floorLoads` rows gain `pattern` and their reactions become that pattern's wall cases.
`WALL_SLOTS`, `wallSlots` and slot fields go. Verification: the takedown and window-header tests in pattern form,
two patterns with different scale factors equal one pattern with combined loads (1e-12).

As built: `StudWall::new(.., cases: Vec<Vec<WallLoad>>)` (no `add_load_case`: the decoder gathers a wall's cases before
building it); `ElementLoad.wall: Option<WallCaseLoad { case, factor }>`. A single `ElementLoad` cannot add two
patterns' different cases, so the effective load of a wall is a list (`effective_load` in `load_pattern.rs`) that
`ElementOps::assemble_tangent_cases`, `commit_cases` and `stud_wall_response` take; `Add` of two different cases on
one element in one pattern gives `WallCaseLoad::MIXED`, which `Domain::validate` rejects. The decoder (`wall_load_plan`)
groups each (pattern, wall) pair's `wallLoad` rows and floor reactions into one case and emits one `ElementLoad` per
pair, so it never produces `MIXED`. `FloorLoad::slot` stays a generic case index (`wall_loads` returns a map). The
building builder keeps its wire format and makes one pattern per used slot.

**C4. Building as an expansion (done).** `building::expand(&BuildingSpec) -> CarapaceInputV1` (explicit tables, walls on
the storey grid with floor zones, anchorage zeroLengths, rigid diaphragm rows or flexible-diaphragm loads, one
pattern per load case, plus an `up` axis option so PySees gets Z-up). Native API, `carapace expand` in the CLI, and
`expandBuilding` in wasm. The decoder drops the `building` table, `buildingWithTables` and the building recorders
(`buildingWallPier` etc. become `studWallPier` on the expanded indices; the expansion returns a name → element
index map). Verification: the builder tests run on the expansion; a round trip expand → decode → solve equals the
current builder results.

## PySees phases

**P1. Element and material.** `studWall` and `diaphragm` element types and an `Anchorage` uniaxial material, all
flagged Carapace-only. Wall form: corners, stations, and a layout editor (studs, openings, plates, sheathing law,
floor zone) as a custom form block. The explorer lists them like any element; the viewport draws walls with studs,
openings and the floor zone (reusing the current drawing code).

**P2. Compiler and loads.** `studWalls`, `diaphragms`, anchorage materials, 3D `rigidDiaphragm` with its normal,
`wallLoad` eleLoads and floor loads in patterns.

**P3. Results.** Real nodes give the deformed shape, reactions and node plots for free. Stud-wall element results
(piers, base reactions, chords) in the element-force table and the plot channel picker.

**P4. Export and project files.** Warning comments for Carapace-only items; project files save the new element
and material (round trip test); importing a script does not see them.

**P5. Template.** The stud-wall building template calls `expandBuilding` and converts the tables to ordinary
entities (nodes, walls, zeroLengths, materials, fixes, constraints, patterns, analysis blocks). After that it is the
user's model. Remove the building entity, card, directive and building-only results code.

Each phase keeps `npm test` green, including `carapaceVsOpenSees.test.ts` (unaffected: stud walls have no OpenSees
counterpart; the template gets its own statics tests as now).
