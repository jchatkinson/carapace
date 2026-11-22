# Floors

Floor components for stud-wall models (core Rust API, 2D). The plan and its rationale are in
[floor-system-plan.md](floor-system-plan.md); the walls themselves are in [stud-wall.md](stud-wall.md).

## 1. Floor depth: the floor band

`FloorBand` represents the floor depth `d` between the top plate of the wall below and the bottom plate of the
wall above. Its two node rows share station positions. The local frame runs from the first to last bottom node;
it may have any orientation in the 2D plane. Each node contributes `ux, uy` only, so the band never connects
wall rotations.

`FloorBandLayout::between(&below, &above, d, ga, rim_ea)` forms stations from both walls' stud and pier-end
positions. It marks bays lying inside an opening below as `None`, because segmented design has no shear transfer
there. Supply the bottom and top node rows to `FloorBand::new` in station order. Geometry is checked when the
analysis is built.

The band has a constant shear field in each active bay and an axial rim bar along each bay. It has **no vertical
stiffness**: add a `ZeroLength` between corresponding bottom and top nodes for bearing, and hold-down springs at
pier ends. Without vertical links the stacked system is a mechanism. Apply lateral diaphragm loads as nodal loads
on the band's top row (the sheathing level), spread along it the way the diaphragm delivers it.

**Where the floor-depth couple goes.** In each bay the shear field takes the horizontal force from the top row
to the bottom row, half to each node, and puts its vertical couple `q·a·d` on the **bottom row only**: the rim
bears on the wall below. So:

- the wall below gets the O86 Fig. 11.7 lever `V · (H + d)` at its base;
- a link across the band (the hold-down between storeys) carries the moment at the **top** of the band, the
  base of the wall above: `R21 = F21 (H + h2)/L1` in Fig. 11.7, where `h2` is the depth of the floor or roof
  above that wall.

For the roof storey, put a band on top of the top wall too if the roof framing depth should be in its lever,
and apply the roof force at that band's top row.

**Stacking.** Wall below → band (bottom row = the wall's top corners and top stations; top row = the wall
above's bottom corners and bottom stations) → wall above. The 2- and 5-storey tests match an exact hand
calculation that walks up the stack with the A.11.7.1 storey terms over `H`, a rigid rotation over `d`, and
true lever arms (1e-5). Compared with the code's own A.11.7.1 with `H_i = H + d` (treating the floor depth as
flexible wall height), the rigid band gives 4–6% less drift; leaving out the floor depth gives a further 4–6% less
than the band.

**Openings in the wall below.** `between` gives the bays over an opening no shear field, so each pier takes its
own shear. The rim still runs over the opening as the floor-level drag strut. The wall's own top plate also
crosses it (axially), so the drag runs in two struts a depth `d` apart, and their couple ties the piers
together slightly. Under a diaphragm load spread along the band with a real rim, each pier's moment is within
1% of `V_p (H + d)` (0.5% in the test). With a very stiff rim and the whole load entering at one end, it reaches
about 5%.

```rust
let layout = FloorBandLayout::between(&below_layout, &above_layout, 254.0, 1e10, 1e12)?;
let band = FloorBand::new(bottom_nodes, top_nodes, layout)?;
let id = domain.add_element(Element::FloorBand(band));
// Add a vertical bearing ZeroLength at each station and hold-down springs at pier ends.
// Apply diaphragm shear as nodal loads on the top-node row.
let response = domain.element_floor_band_response(id).unwrap();
```

Choosing `ga`: the band's shear stiffness is the rim joist and blocking with their nailing; it is usually much
stiffer than the walls, and a large value (a rigid stand-in such as `1e10` N/mm) is the normal choice. Bearing
links: perpendicular-to-grain stiffness of the plates and rim at each stud, or a rigid stand-in. Very large
stand-ins limit the solution to about 1e-9 relative.

`FloorBandResponse` reports each bay's `gamma` and shear flow `q` (zero for inactive bays), plus rim axial force
per segment (tension positive). The band is linear and stateless; it has no element loads.

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
