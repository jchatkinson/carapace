# Floors

Floor components for stud-wall models (core Rust API and `CarapaceInputV1` tables, 2D and 3D). The plan and its rationale are in
[floor-system-plan.md](floor-system-plan.md); the walls themselves are in [stud-wall.md](stud-wall.md).
For automatic stacked building generation, see [building.md](building.md).

## 1. Floor depth: the floor zone of the wall

The floor or roof depth `d` between storeys is part of the stud wall that carries it: `StudWallLayout::floor =
Some(FloorZone { depth, ga, rim_ea })`. The old `FloorBand` element, its `floorBands` table and the
`floorBandBay` recorder are gone. The model, the O86 Fig. 11.7 levers, the verification and the wire form are in
[stud-wall.md](stud-wall.md#2a-floor-zone); in short:

- the wall's exterior top corners and top stations sit at `H + d`, framing (studs, plates, sheathing) over `0..H`
  as before; the zone is condensed inside the wall (bays of constant shear field `ga`, a rim bar `rim_ea`, rigid
  vertical posts to the framing), so a wall spans storey grid to storey grid;
- the zone has **no vertical stiffness** of its own beyond the rigid posts: bearing and hold-downs are
  `ZeroLength` links *between* walls, from the lower wall's exterior top node to the upper wall's coincident
  bottom node (tie the two horizontally with `equal_dof`). Without those links the stacked system is a mechanism;
- apply lateral diaphragm loads as nodal loads on the exterior top row (the sheathing level), or as top slot
  loads, which act there when a zone exists.

**Where the floor-depth couple goes.** In each bay the shear field takes the horizontal force from the exterior
row to the framing top, half to each node, and puts its vertical couple `q·a·d` on the framing top: the rim bears
on the wall below. So the wall carrying the zone gets the lever `V · (H + d)` at its base, and a link across the
floor (the hold-down between storeys) carries the moment at the top of the zone, the base of the wall above:
`R21 = F21 (H + h2)/L1` in Fig. 11.7, where `h2` is the depth of the floor or roof above that wall. Give the roof
storey a zone too if the roof framing depth belongs in its lever.

**Stacking.** Wall below (with its zone) → vertical links at each top station → wall above. The 2- and 5-storey
tests match an exact hand calculation that walks up the stack with the A.11.7.1 storey terms over `H`, a rigid
rotation over `d`, and true lever arms (1e-5). Compared with the code's own A.11.7.1 with `H_i = H + d` (treating
the floor depth as flexible wall height), the rigid zone gives 4–6% less drift; leaving out the floor depth
gives a further 4–6% less.

**Openings.** A zone bay over an opening of its own wall has no shear field, so each pier takes its own shear.
The rim still runs over the opening as the floor-level drag strut, and the wall's own top plate also crosses it,
so the drag runs in two struts a depth `d` apart whose couple ties the piers together slightly. Under a diaphragm
load spread along the exterior top with a real rim, each pier's moment is within 1% of `V_p (H + d)` (0.5% in the
test). With a very stiff rim and the whole load entering at one end it reaches about 5%.

Choosing `ga`: the zone's shear stiffness is the rim joist and blocking with their nailing; it is usually much
stiffer than the walls, and a large value (a rigid stand-in such as `1e10` N/mm) is the normal choice. Very large
stand-ins limit the solution to about 1e-9 relative. `StudWallResponse::floor_bays` / `floor_rim` report each
bay's `gamma` and shear flow `q` (zero over an opening) and the rim force per bay (tension positive).

## 2. Gravity: one-way joist transfer

`model::transfer_floor` sends one-way joist-floor gravity loads to bearing lines using hand-calculation
reactions. It creates no floor element and does not use support stiffness. Define a rectangular `FloorBounds`
in span-aligned `(s, t)` coordinates, list the bearing lines as `FloorSupport`s, then provide `FloorLoad`s by
slot. `FloorContinuity::Simple` shares each span by statics; `Continuous` solves a constant-EI beam on rigid
supports. EI cancels from the reaction calculation. Reactions preserve support input order.

```rust
let floor = FloorBounds { s0: 0.0, s1: 3600.0, t0: 0.0, t1: 2400.0 };
let supports = [
    FloorSupport { s: 0.0, t0: 0.0, t1: 2400.0 },
    FloorSupport { s: 3600.0, t0: 0.0, t1: 2400.0 },
];
let reactions = transfer_floor(
    floor,
    &supports,
    &[FloorLoad::Area { slot: 0, w: 2.0, region: None }],
    FloorContinuity::Simple,
)?;
let wall_top_loads = reactions.per_support[0].wall_loads(0.0, false);
```

For this 3.6 m by 2.4 m floor with pressure 2 force/area, each wall receives 3.6 force/length over its
2.4 m support length: 8.64 force total per wall. A wall with a window receives that `TopLine` on its top
plate; the existing stud-wall takedown sends the opening tributary to the header and jacks. Put dead and live
loads in separate slots, then select/factor them with `ElementLoad::wall`, as with any other wall load.

Floor joists in platform framing bear on the top plate of the wall below, so map their reactions to that
wall. `wall_loads` converts downward-positive reactions to the wall's upward-positive local axes. Set
`reversed` when the floor's increasing `t` direction runs opposite the wall's increasing local `x`.
`beam_segments(slot)` returns uniform line-load pieces for a beam; split the beam at each piece endpoint before
applying `ElementLoad::beam_uniform` or `beam_taper`.

This phase supports rectangles only. Supports must lie within the span bounds and overlap the floor width;
their `t` extents are clipped to the floor. Loaded strips or concentrated joist lines need at least two
supports. Non-finite values, zero/negative extents, loads outside the floor and slots outside `0..WALL_SLOTS`
return `FloorError`. Free floor edges are allowed and become overhangs when they extend beyond an end support;
negative support reactions are returned as uplift without clamping. The idealizations are simple spans or a
continuous beam on rigid supports, constant EI, and no joist/support stiffness interaction.

## 3. Pier-end anchorage

Use one `Material::anchorage(tension, slack, bearing)` on a vertical `ZeroLength` at every pier end. The same
link supplies tension-only hold-down response after slack closes and compression bearing/crushing below the
post. Positive relative `uy` is uplift; negative is compression. For a between-storey link, order its nodes so
positive `uy` means the upper post rises relative to the lower one. The constructor validates all parameters.

`tension` is the effective series stiffness in N/mm: combine device stiffness (`T_cap / delta_max`), rod
stiffness (`A E / L`) and fastener slip as series springs. Use zero for bearing-only support. `slack` is mm of
framing shrinkage and any physical gap; an estimate is `0.002 * (MC_fab - MC_service) * L_shrink`, with
moisture contents in percent. `AnchorageBearing::Linear { k }` uses N/mm. `AnchorageBearing::O86Crush { qr }`
uses `qr = K_D K_S K_Z f_cp A` in N and crush in mm. The crush curve is the WoodWorks Shearwalls 11 design-report
legend “Crush – Deformation of bottom plate at compression end of wall segment.” CSA O86:19 11.7.1.2 defines
`d_a` as total anchorage displacement; WoodWorks composes elongation, slip, shrinkage and crush. One link at
each pier end lets the wall produce the rotation `(uplift + crush) / jd` directly. `anchorage_split` reports the
current slack and remaining strain portions for drift reporting.

The slack portion is constant shrinkage drift and is not amplified by `RdRo/Ie` in NBC seismic drift. For a
linear wall after all hold-downs are engaged, it may also be isolated by two runs: drift with the specified
slack minus drift with `slack = 0`.

Analysis: the law is nonlinear (slack, bilinear tension, crush curve), so use Newton with load steps. With
rigid stand-ins elsewhere, use an absolute displacement-increment test no tighter than about 1e-7 mm. Verified
in `core/tests/elements/stud_wall_anchorage.rs`: a single-storey wall's drift equals the O86 11.7.1.2 hand total
with `d_a = slack + T/k + crush(C)` from statics, with and without uplift, both directions (1e-5, observed 2e-8),
and the slack adds exactly `H·slack/L`. The window wall with anchorage at all four pier ends equals two separate
pier walls (1e-10).

## 4. Diaphragms

For the rigid and flexible approaches below, this guide uses the conventional Y-up building model with floors in
the X-Z plane. The semi-rigid `Diaphragm` element can use any axis-aligned global plane. Choose one of three
idealizations and bracket wall-line forces when they matter: the FPInnovations guide reports differences up to 60%
between idealizations.

### Rigid

`Domain3::rigid_diaphragm(master, slaves)` constrains each slave's `ux` and `uz` to the retained node's
translation and in-plane rotation about global `Y`. Apply the storey force and any accidental torsion as
loads on the retained node. This distributes forces by wall-line stiffness and torsion. The retained node
must have rotational `ry` active through the rigid constraint.

Verified (`rigid_diaphragm_matches_stiffness_and_torsion_split`): four 3D walls (two X lines of unequal length, two
Z lines) under an eccentric X force give wall-line shears equal to the hand rigid-body split
`V_i = k_i r_i · K⁻¹ F`, `K = Σ k_i r_iᵀ r_i`, with each line's stiffness `k_i` from a run of that wall alone (1e-9).

```rust,ignore
domain.rigid_diaphragm(master, &level_nodes);
```

### Flexible

`flexible_diaphragm_loads` wraps the simple-span statics in [`transfer_floor`](../core/src/model/floor.rs).
Supports are wall lines perpendicular to the lateral load, and area loads represent lateral force per floor
area. Returned line reactions are distributed to those supports without using wall stiffness. Use
`lump_line_reaction` to distribute each uniform reaction piece to ordered exterior top nodes by midpoint
tributary length; the nodal forces sum exactly to the line reaction.

```rust,ignore
let reactions = flexible_diaphragm_loads(bounds, &wall_lines, &wind_area_loads)?;
let nodal = lump_line_reaction(&top_node_coordinates, reaction_line)?;
```

### Semi-rigid

`Diaphragm` is a 3D-only finite element on a rectangular grid whose nodes share one global coordinate (X, Y or Z).
That shared coordinate defines the normal and `level`. The in-plane axes cycle after the normal:
`(u, v) = ((normal + 1) % 3, (normal + 2) % 3)`, so a Y-normal grid uses `(Z, X)` and a Z-normal grid uses
`(X, Y)`. Nodes are ordered `i + j * nu`, with `i` indexing `u` and `j` indexing `v`; each node activates only
the two in-plane translations. Each bay has an independent `WallShearLaw` or `None` for an opening. Engineering
shear is the symmetric bilinear average `gamma = du/dv + dv/du`; the bay energy is `area * integral(q dgamma)`.
Every grid segment has an axial bar with its line's `EA`: `ea_u` for constant-`v` lines and `ea_v` for
constant-`u` lines. Interior grid lines need positive axial stiffness. There are no element loads: lump floor-level
lateral loads to grid nodes, or use the flexible transfer helper when a distribution by statics is wanted.

```rust,ignore
let layout = DiaphragmLayout {
    u: vec![0., 3000., 6000.], v: vec![0., 2400.], level: 3000.,
    bays: vec![Some(WallShearLaw::Linear { ga: 2.0e5 }); 2],
    ea_u: vec![1.0e8; 2], ea_v: vec![1.0e8; 3],
};
let id = domain.add_element(Element3::Diaphragm(Diaphragm::new(grid, layout)?));
let response = domain.element_diaphragm_response(id).unwrap();
```

The response reports bay `gamma` and shear flow `q`, axial force on every segment, boundary shear flow, and
chord force. For `Linear { ga: Bv }`, compare mid-span deflection with CSA O86:19 11.7.2 (p. 156):
`Delta_d = 5 v L^3/(96 E A L_D) + v L/(4 B_v) + 0.00061 L e_n + sum(Delta_c x)/(2 L_D)`.
Verified (`simply_supported_diaphragm_converges_to_o86_11_7_2`): a 12 m × 6 m diaphragm, one bay deep, with
`ga = B_v` converges to the first two terms at second order in the number of bays (1.4e-4 at 64 bays). The shear
term is exact for nodal load lumping; the chord term carries the discretization error.

Nail slip: with the `O86` law the model integrates the slip of the actual shear distribution. For `e_n ∝ v²` and a
shear that falls linearly to zero at mid-span, that is `0.0025 · e_n · L/6 ≈ 0.000417 L e_n`, i.e. 0.68 of the
code's `0.00061 L e_n` (the test reproduces 0.683). The code term is the more conservative of the two.
