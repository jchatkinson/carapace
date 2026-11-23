# Stud wall element: user guide

`StudWall` is a 2D macro-element for one storey of a light-frame wall (wood or cold-formed steel). It is
defined by its framing: studs and stud packs, plates, openings with headers, jacks and cripples, and an
optional sheathing law. The studs, plates and headers form an internal frame that is condensed onto a few
exterior nodes. You model walls; the element models the studs.

It is built for design:

- **gravity takedowns** that match hand calculations (tributary or continuous top plate, stud packs, headers,
  jacks and cripples);
- **lateral analysis** that matches segmented shear-wall design: the SDPWS 2021 Eq. 4.3-1 and CSA O86:19
  11.7.1.2 deflection terms for one storey, and the O86 Annex A.11.7.1 multi-storey method for stacked walls;
- **design forces** per pier: shear flow `q` and the M/jd chord (hold-down) force.

Status: core Rust API only (`carapace-core`, 2D `Domain` and vertical-plane `Domain3`). `FloorBand` supports stacking storeys with floor depth.
Neither element is yet in the JSON wire format, the CLI, the Python package or PySees. Implementation details
are in [`stud-wall-v1.md`](stud-wall-v1.md); the long-term
design and its rationale in [`stud-wall-panel.md`](stud-wall-panel.md).

Contents: [1 Geometry](#1-geometry-and-nodes) · [2 Layout](#2-layout) · [3 Sheathing laws](#3-sheathing-laws) ·
[4 Supports](#4-supports-anchorage-and-hold-downs) · [5 Loads](#5-loads) · [6 Analysis](#6-analysis-setup) ·
[7 Results](#7-results) · [8 Example](#8-example) · [9 Assumptions](#9-built-in-assumptions) ·
[10 Limitations](#10-limitations) · [11 Errors and limits](#11-errors-and-limits) ·
[12 Accuracy](#12-accuracy-and-verification)

## 1. Geometry and nodes

- **Corners.** `StudWall::new([c0, c1, c2, c3], stations, layout, loads)`. The corners go counter-clockwise from
  bottom-left: `c0` at local `(0, 0)`, `c1` at `(L, 0)`, `c2` at `(L, H)`, `c3` at `(0, H)`. Local `x` runs
  along the bottom plate from `c0` to `c1`; local `y` is up the studs. The rectangle may have any orientation
  in the global plane; `validate` checks it against `L` and `H` (relative 1e-9) when the analysis is built.
- **Stations** (`WallStation { chord: Top | Bottom, x, node }`) are extra exterior nodes on a plate, at local
  `x`. Use them where something else connects to the plate between the corners: a floor or the wall above
  (to transfer bearing at each stud), a drag strut, a beam bearing on the top plate, a sill anchor you want to
  model as a spring. Every station's node must lie at its plate point.
- **Exterior DOFs are `ux` and `uy` only.** The element never connects nodal rotations: plate splices and
  plate-to-floor connections are not moment connections. A node used only by stud walls has no rotational
  equation (Carapace activates only DOFs some element stiffens), so you need not fix `rz`.
- Corners and stations must be distinct nodes. Two walls side by side or stacked can share corner nodes.

### 3D models

Use `StudWall::new3(..., plane)` for a wall in a vertical plane of `Domain3`. `WallPlane::vertical(ex, ez)` takes
the unit plan direction from `c0` to `c1`; global `Y` is up. Corners and stations must lie in that plane at the
same local `(x, y)` positions as in a 2D model. For example, the point at local `(x, y)` is
`p0 + x * (ex, 0, ez) + y * (0, 1, 0)`.

The wall stiffens global `uy` and the plan translation components present in its direction (`ux` when `ex != 0`,
`uz` when `ez != 0`). It never stiffens rotations. A node on one wall line has no stiffness normal to that wall;
provide that support with a crossing wall, a diaphragm constraint, or a fixity. Wall response quantities and
`local_force` remain in wall-local axes. `ElementLoad3::wall` carries the same eight top-plate slot factors as
`ElementLoad::wall`.

## 2. Layout

All values are in consistent units (N and mm throughout this guide), in the wall's local frame.

| Field | Meaning |
|---|---|
| `length`, `height` | `L` (bottom plate, `c0 → c1`) and `H` (stud length, plate centreline to centreline). |
| `studs: Vec<WallStud { x, ea }>` | Every stud row, strictly increasing in `x`, the first at `0` and the last at `L` (the end posts). A pack of plies is one row with the summed `EA`. Studs are axial bars, pinned at both plates. |
| `top_ea` | Top plate axial stiffness `EA` (all plies). |
| `top_distribution` | `Direct`: the top plate is hinged at every stud, so gravity load goes to studs by the lever rule (simple-span tributary distribution, the usual hand calculation). `Beam`: the top plate is a continuous beam over the studs (continuity reactions). |
| `top_ei: Option<f64>` | Top plate `EI`; `Some` with `Beam`, ignored with `Direct`. |
| `bottom_ea`, `bottom_ei` | Bottom plate `EA` and `EI`. The bottom plate is always a beam. |
| `openings: Vec<WallOpening>` | `x0`, `x1` (clear opening edges, where the jacks are), `head_height` (underside of header), `header_ei`, `jack_ea` (all plies of one side), `door` (cuts the bottom plate between `x0` and `x1`). |
| `sheathing: Option<WallSheathing { law }>` | One shear law for the whole wall, applied to each pier ([Section 3](#3-sheathing-laws)). `None` for a gravity-only wall. |
| `foundation: bool` | The wall bears on a rigid foundation: every interior bottom-plate position (each stud and jack) is fixed internally in `u` (sill anchorage) and `v` (bearing). The bottom corners remain exterior nodes, and a bottom station replaces the internal support at its position, for hold-downs at interior pier ends ([Section 4](#4-supports-anchorage-and-hold-downs)). |

**Openings.**
- Each opening gets two jacks (at `x0` and `x1`, from the bottom plate to the header) and a header from `x0` to
  `x1`, simply supported on the jacks. Stud rows strictly inside `(x0, x1)` become cripples from the header to
  the top plate. Below a window there are no members (sills and lower cripples are not modelled).
- **Put a stud row exactly at `x0` and at `x1`: these are the king studs.** The king stud and the jack share
  one bottom-plate node and form the pier's end post at the opening, which is where its hold-down goes. Give
  the row the king stud's `EA`; the jack's is `jack_ea`.
- The top and bottom plates are hinged at `x0` and `x1` for the segments over the opening. Across an opening
  the plates tie the piers axially (drag) only, and the header transfers no lateral force
  ([Section 9](#9-built-in-assumptions)).
- At least one full-height stud must lie between consecutive openings (king studs count).

**Piers.** A pier is a run of full-height studs not interrupted by an opening, from its first to its last
full-height stud (`xa`, `xb`, width `b = xb − xa`). A pier needs at least two studs. Jacks are not full-height
studs, so **without a king stud at the opening edge a pier stops at the last regular stud**, and the strip next
to the opening is unsheathed in the model. `layout.pier_ends()` lists the pier end positions. Only piers carry
shear: sheathing above and below openings is ignored, which is segmented design
([Section 9](#9-built-in-assumptions)).

## 3. Sheathing laws

The law relates a pier's shear flow `q` (force per length, N/mm) to its sheathing shear strain `γ_sh`
(radians). It is calibrated per wall, for the whole sheathing assembly. One law applies to every pier of the
wall, each pier with its own state.

| Law | Behaviour | Use for |
|---|---|---|
| `Linear { ga }` | `q = ga · γ_sh`. Linear: load cases superpose, and the linear algorithm is exact. | SDPWS apparent shear stiffness `Ga` (convert to N/mm); a code secant stiffness at a known design shear; gravity-plus-lateral models where superposition matters. |
| `O86 { bv, nail_spacing, nail_diameter }` | `γ_sh = q/Bv + 0.0025 · e_n` with `e_n = (0.013 · q · s / d_F²)²` (odd in `q`): the shear and nail-slip terms of CSA O86:19 11.7.1.2 divided by `H`, with `e_n` from Annex A.11.7. Nonlinear elastic (no hysteresis). | Code deflection at any load level, with each pier at its own shear. |
| `Material(m)` | `q = σ(γ_sh)` from any uniaxial `Material` (`Steel01`, `Hysteretic`, `Pinching4`, ...), with history. Initial tangent must be positive. | Pushover and cyclic analysis with a calibrated backbone. |

**O86 law: units and scope.** The constant 0.013 makes the law unit-specific: `bv` in N/mm (Tables 9.1 to 9.3;
7000 N/mm for gypsum wallboard), `nail_spacing` (panel-edge spacing) and `nail_diameter` in mm, so the model
must be in N and mm. A.11.7 applies to nails in wood-based sheathing on **dry lumber**. Adapting it:

- **Green lumber** (A.11.7 doubles `e_n`): `e_n` scales with `s²`, so use `nail_spacing × √2`.
- **Two identical sides:** each side carries `q/2`. Use `bv × 2` and `nail_spacing / 2`.
- **Gypsum wallboard** (`e_n = 0.76 mm` constant) and **unblocked walls** (O86 11.7.1.3 divides the whole
  deflection by `J_us`) cannot be expressed exactly. Use a `Linear` secant stiffness at the design shear, or a
  `Material` backbone.
- O86 evaluates `v` as the maximum specified shear at the top of the wall. The law uses each pier's own `q`, so
  piers of one wall can sit at different points of the curve. That is the intent, but it is a modelling
  choice, not a clause.

**Calibration notes for any law.**
- The law is the whole wall's sheathing response: panel shear plus fastener slip, panel seams, blocking and
  both sides. Plate and stud axial flexibility are added by the element; don't include them in `ga`.
- Aspect-ratio adjustments (for example SDPWS 4.3.3.2) are per pier: compute them on the pier width `b`.
- The hold-down (anchorage) term is not in the law. It comes from your hold-down springs
  ([Section 4](#4-supports-anchorage-and-hold-downs)).
- For SDPWS `Ga` in kips/in, convert to N/mm (1 kip/in = 175.13 N/mm).

## 4. Supports, anchorage and hold-downs

**Racking must be restrained.** An unsheathed wall has a zero-stiffness racking mode. Gravity-only walls need a
lateral restraint (an `ux` fix at a top or bottom corner, or the lateral system).

**Anchor every pier end.** Segmented design puts overturning anchorage at both ends of every pier, including
the end posts beside each opening. `layout.pier_ends()` gives their positions: `0` and `L` are the corners
`c0`/`c1`, and each other position needs a bottom station with its own hold-down. Without one, that post is
anchored by whatever holds the bottom plate there. On a foundation that support is rigid, which leaves out
that hold-down's elongation (`d_a`) and makes the wall too stiff: 23% too stiff for a 3.6 m wall with one
window in a probe. The M/jd forces are right either way.

Use `Material::anchorage(tension, slack, bearing)` on a vertical `ZeroLength` at every pier end, at the
foundation and between storeys. It combines tension-only hold-down action and compression bearing/crushing in
one path-independent law. Positive relative `uy` is uplift. In a between-storey link, order nodes so positive
`uy` means the upper post rises relative to the lower one. Use the same link for each pier end of a `FloorBand`.

`tension` is the effective hold-down stiffness in N/mm: combine device stiffness (`T_cap / delta_max`), rod
stiffness (`A E / L`) and fastener slip in series. Set it to zero for bearing only. `slack` is shrinkage plus
physical gaps in mm; estimate shrinkage as `0.002 * (MC_fab - MC_service) * L_shrink` (moisture contents in
percent). For compression use `AnchorageBearing::O86Crush { qr }` with `qr = K_D K_S K_Z f_cp A` in N, or
`Linear { k }` with `k` in N/mm. The O86 crush curve is the WoodWorks Shearwalls 11 design-report legend
“Crush – Deformation of bottom plate at compression end of wall segment.” O86:19 11.7.1.2's `d_a` is total
anchorage-system displacement; WoodWorks adds hold-down elongation, slip, shrinkage and crush. A single link at
each pier end gives the wall's rotation `(uplift + crush) / jd` directly.

For seismic drift, separate the slack contribution: shrinkage drift is not amplified by `RdRo/Ie`. The
`anchorage_split` helper reports slack and remaining deformation. For a linear wall after all hold-downs engage,
the two-run difference (specified slack minus zero slack) gives the same constant shrinkage drift.

**On a foundation** (`foundation: true`), interior studs and jacks bear and are anchored internally, and their
reactions are reported in `base_nodes`. Put hold-downs at the corners and at bottom stations on the other
pier ends: a bottom station replaces the internal support at its position. For load reversal, every pier end
needs a hold-down (both directions).

**On a floor or another wall**, leave `foundation: false` and connect the bottom plate with stations: a
station at each stud and each opening edge, joined to the plate or wall below (shared nodes, or `ZeroLength`
springs for bearing and slip), with hold-downs at the pier ends. For stacked storeys with floor depth, use a
`FloorBand` between the walls ([floors.md](floors.md)); its bearing and hold-down links are external.

**Anchor the bottom plate along its length**, not only at its corners. If it is anchored only at the corners,
it acts as a drag strut, and the wall comes out 1–5% more flexible than the code equations. `foundation: true`
does this; between storeys, connect the bottom-plate stations sideways as well.

If a separate tension-only idealization is needed, `Material::Ent` on `ZeroLength::new(c0, ground)` resists
uplift when the wall node is first; add separate compression bearing. That legacy recipe is an alternative to
`Anchorage`. The wall's reported pier M/jd force remains the design chord force ([Section 7](#7-results)).

## 5. Loads

Wall loads are defined once, at construction, in up to `WALL_SLOTS = 8` slots. A load pattern switches a slot
on with a scale factor.

- `WallLoad::TopPoint { x, fx, fy }`: point load on the top plate at `x`.
- `WallLoad::TopLine { x0, x1, wx0, wx1, wy0, wy1 }`: line load on the top plate from `x0` to `x1`, linear
  from `(wx0, wy0)` to `(wx1, wy1)` per unit length.
- All loads are in the wall's local axes: `fx`/`wx` along the wall (shear), `fy`/`wy` along the studs
  (negative is down).
- Put each load case in its own slot (dead, live, wind ...). In a pattern, add an
  `ElementLoad { wall: e_slot, ..Default::default() }` for the wall (`wall[slot] = 1.0`); the pattern's series
  and factor then scale that slot. Several patterns can scale the same slot.
- Loads at exterior nodes (corners and stations) can also be ordinary nodal loads.
- The wall rejects other element-load kinds (beam, body, edge), and other elements reject a nonzero `wall`.
- There is no self-weight or mass: include wall weight as a top line load, or as nodal loads.

### Floor loads

Floor gravity loads reach a wall as top-plate loads computed by the one-way floor transfer; see
[floors.md](floors.md#2-gravity-one-way-joist-transfer).

## 6. Analysis setup

| Sheathing | Algorithm | Notes |
|---|---|---|
| none or `Linear` | `Algorithm::Linear` (or Newton) | One step. Results from separate runs superpose. |
| `O86` or `Material` | `Algorithm::Newton` | `Algorithm::Linear` takes one tangent step and does **not** converge a nonlinear wall; the result would be wrong without an error. Results do not superpose: analyse each load combination as one run. |

- Use load control with a few steps (for example 10 × 0.1) for nonlinear walls. The O86 law is smooth and
  monotonic, so Newton converges quickly.
- `ConvergenceTest::NormDispIncr` is the robust choice when the model contains near-rigid stand-in springs:
  an unbalance norm becomes ill-conditioned there.
- For cyclic `Material` laws, use displacement control or small load steps and commit each step: the pier
  material states advance only at `commit`.
- Static analysis only: the element has no mass, so it contributes nothing to modal or transient analysis.

## 7. Results

`domain.element_stud_wall_response(id, pseudo_time) -> Option<StudWallResponse>` recovers the internal
response at the current state. Its argument is the pseudo-time at which slot factors are evaluated (use the
analysis time). All quantities are in local axes, with **axial tension positive**.

| Field | Contents |
|---|---|
| `studs` | Full-height studs: `x`, axial force at the `bottom` and `top`. They differ inside a sheathed pier by the shear lag `(q_L − q_R) · H`. |
| `jacks`, `cripples` | Axial force per member. |
| `header_reactions`, `header` | Header end reactions (the jack forces), and moment and shear at each header node. |
| `top_chord`, `bottom_chord` | Plate axial force, moment and shear at each plate node. Moments are positive sagging (as for the header). |
| `piers` | Per pier: `xa`, `xb`, `q` (N/mm), `gamma`; `moment_base`, `moment_top`; `jd = xb − xa`; `chords`. |
| `base_nodes` | Force on the wall at each supported bottom position: `c0`, `c1`, bottom stations, and with `foundation` every internal bearing position, sorted by `x`. Gravity gives positive `fy`. |
| `base` | `[Fx, Fy, M]`: the sums over `base_nodes`, with `M` about `(L/2, 0)`. |

**Design use.**
- **Shear:** pier `q` is the design shear flow for sheathing and nailing.
- **Hold-down and chord tension and compression:** use `chords = [−moment_base/jd, +moment_base/jd]`, the
  M/jd force of CSA O86 and SDPWS. It holds exactly by statics, however the element shares overturning among
  studs. It is the overturning share only: add the post's gravity force from a gravity analysis (or from
  `studs`) per your load combination. A +x (rightward) lateral load gives tension in the left post.
  - `moment_base` and `moment_top` are the moment of the stud axial forces about the pier mid-width, CCW
    positive. They equal the overturning from everything above that section, so a +x load gives negative
    values, and `moment_top = moment_base + q·b·H`.
- **Do not size hold-downs from the end-stud force** when interior studs bear (`foundation: true` or bearing
  stations). Interior studs then carry part of the overturning, and the end-stud force is 4–11% below M/jd in
  the 2- and 5-storey tests.
- **Gravity:** `studs`, `jacks`, `cripples` and the header results are the takedown.
- `domain.element_end_force(id, time)` gives the exterior nodal forces (local, `[fx, fy]` per exterior node:
  corners, then stations).

## 8. Example

A 3.6 m wall with a window (king studs at its edges), on a foundation, with the O86 law, a hold-down at every
pier end, and roof dead load plus wind run together. The same code runs as a test in
`core/tests/elements/stud_wall_guide_example.rs`, so it stays current.

```rust
let (l, h) = (3600.0, 2440.0); // mm
let ea_ply = 9500.0 * 38.0 * 140.0; // N, one 38x140 SPF ply
let mut studs: Vec<WallStud> = (0..=9)
    .map(|i| WallStud {
        x: i as f64 * 400.0,
        ea: ea_ply,
    })
    .collect();
studs[0].ea = 3.0 * ea_ply; // 3-ply end posts
studs[9].ea = 3.0 * ea_ply;
// The rows at the window edges (x0 = 1600, x1 = 2400) are the king studs: with the jacks they form the
// pier end posts at the opening.
let layout = StudWallLayout {
    length: l,
    height: h,
    studs,
    top_ea: 2.0 * ea_ply, // double top plate
    top_ei: None,
    top_distribution: ChordDistribution::Direct,
    bottom_ea: ea_ply,
    bottom_ei: 9500.0 * 140.0 * 38f64.powi(3) / 12.0,
    openings: vec![WallOpening {
        x0: 1600.0,
        x1: 2400.0,
        head_height: 2000.0,
        header_ei: 9500.0 * 2.0 * 38.0 * 235f64.powi(3) / 12.0,
        jack_ea: ea_ply,
        door: false,
    }],
    sheathing: Some(WallSheathing {
        law: WallShearLaw::O86 {
            bv: 11_000.0,
            nail_spacing: 100.0,
            nail_diameter: 3.25,
        },
    }),
    foundation: true,
};

let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
let ground = d.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
// A hold-down at every pier end (segmented design): the corners, plus a bottom station at each interior
// pier end, which replaces the foundation support there. Linear springs, so they also act in compression.
let mut post = |x: f64| {
    let n = d.add_node(Node::new([x, 0.0]).fix(0));
    d.add_element(Element::ZeroLength(
        ZeroLength::new(ground, n).with_material(1, Material::Elastic { e: 20_000.0 }),
    ));
    n
};
let (c0, c1) = (post(0.0), post(l));
let stations: Vec<WallStation> = layout
    .pier_ends()
    .into_iter()
    .filter(|&x| x > 0.0 && x < l) // 1600 and 2400
    .map(|x| WallStation {
        chord: WallChord::Bottom,
        x,
        node: post(x),
    })
    .collect();
let c2 = d.add_node(Node::new([l, h]));
let c3 = d.add_node(Node::new([0.0, h]));

// Slot 0: roof dead load. Slot 1: wind shear along the top plate.
let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
loads[0].push(WallLoad::TopLine {
    x0: 0.0,
    x1: l,
    wx0: 0.0,
    wx1: 0.0,
    wy0: -5.0,
    wy1: -5.0,
});
loads[1].push(WallLoad::TopLine {
    x0: 0.0,
    x1: l,
    wx0: 4.0,
    wx1: 4.0,
    wy0: 0.0,
    wy1: 0.0,
});
let wall = StudWall::new([c0, c1, c2, c3], stations, layout, loads).expect("valid layout");
let id = d.add_element(Element::StudWall(wall));

for (slot, factor) in [(0, 1.0), (1, 1.0)] {
    let pattern = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1.0 }, factor);
    let mut wall = [0.0; WALL_SLOTS];
    wall[slot] = 1.0;
    d.add_element_load(
        pattern,
        id,
        ElementLoad {
            wall,
            ..Default::default()
        },
    );
}

// Nonlinear sheathing law: Newton, in load steps (load cases do not superpose).
let mut a = AnalysisBuilder::new()
    .constraint_handler(ConstraintHandler::Transformation)
    .integrator(Integrator::LoadControl { increment: 0.1 })
    .algorithm(Algorithm::Newton {
        tangent: TangentStrategy::Current,
        line_search: None,
    })
    .test(ConvergenceTest::NormDispIncr {
        tol: 1e-10,
        max_iter: 25,
    })
    .build(d);
for _ in 0..10 {
    a.step().expect("converged step");
}
let r = a.domain().element_stud_wall_response(id, 1.0).unwrap();
for pier in &r.piers {
    // q (N/mm), M/jd chord forces (N, tension positive) for the hold-down and post design.
    println!(
        "pier {}..{}: q = {:.2}, chords = {:?}",
        pier.xa, pier.xb, pier.q, pier.chords
    );
}
```

Output: the piers carry `q` = 5.21 and 5.06 N/mm (they share the 14.4 kN by stiffness), and each pier
reports its M/jd chord forces for its two hold-downs.

## 9. Built-in assumptions

**Framing.**
- Small displacements; no in-plane second-order (P-Δ) effect.
- Linear elastic framing: studs, plates, jacks, cripples and the header never yield or buckle. Check stud and
  plate capacity, stability and bearing perpendicular to grain in design.
- Studs, jacks and cripples are axial bars pinned at both ends; their bending is ignored.
- Stud packs share load fully (one row, summed `EA`).
- Plates are beams (`Beam` top plate, bottom plate) or hinged bars (`Direct` top plate), continuous along each
  pier (splices are not modelled), and hinged at opening edges for the segments over an opening.
- Interior studs bear bilaterally on what is below them (no uplift separation), and the header is bilaterally
  supported on the jacks.

**Sheathing.**
- One constant shear flow per full-height pier: the field is in equilibrium with the plates and end posts,
  and overturning goes to the posts as the codes assume.
- Sheathing over and under openings is ignored (segmented design): there is no force transfer around
  openings. Laterally, each pier is its own shear wall; the opening framing (header, jacks, cripples) carries
  gravity only.
- Piers are tied across an opening axially only (drag through the top plate). Both plates are hinged at the
  opening edges, and the header has no axial DOF and sits on the jacks, so neither plate bending nor the
  header couples the piers. With hold-downs at every pier end, a wall with openings gives exactly the same
  drift, pier shears and M/jd forces as separate pier walls joined by a drag strut, for `Direct` and `Beam`
  top plates alike (`stud_wall_openings.rs`).
- An unsheathed wall gets no in-plane shear stiffness from the opening framing either. A gravity-only wall
  needs bottom supports on both sides of a door (bearing under the studs or a foundation); otherwise a pier
  pinned at one corner can rotate about it.
- Shear lag: the plates' and end posts' axial flexibility under the linearly varying edge forces is included
  in series with the law.
- Interior studs inside a pier get equal and opposite sheathing traction, so they carry no lag force.
- `FloorBand` has no vertical bearing or anchorage stiffness. Those connections must be modelled externally; omit
  them and the stacked model has a mechanism.

**Code alignment.**
- With rigid plates, end posts `EA` and rigid anchorage, one pier reproduces the SDPWS/O86 bending and shear
  terms exactly. The hold-down term comes from your springs.
- Real interior studs that bear make the wall about 2% stiffer in bending than the code equation, which
  ignores them. This is accepted: it is a more complete model of the same wall.
- Creep, moisture, shrinkage, temperature and the CFS track and connection effects are not modelled.

## 10. Limitations

Not modelled, or not yet available:
- JSON wire format, CLI, Python and PySees (core Rust API only).
- Out-of-plane loading and per-stud out-of-plane output (a 3D wall has only its in-plane translations).
- Mass and self-weight (no dynamics; add wall weight as a load).
- Sills and cripples below windows; force transfer around openings and the perforated shear-wall method
  (segmented design only).
- Nonlinear or hysteretic framing members; staged sheathing installation.
- Internal hold-downs, tension-only bearing and slip inside the element: use exterior `ZeroLength`s.
- A built-in floor band is core Rust only; wire format and automatic storey builder are not implemented.
- Unilateral interior-stud bearing (uplift of interior studs).
- Automatic station generation.
- O86 gypsum `e_n`, unblocked `J_us` and the no-hold-down `d_a` equation (11.7.1.2); SDPWS's nail-slip form.

**Element failure.** If a nonlinear wall's internal pier solve fails to converge, the element returns NaN
forces and the global step fails to converge. There is no wall-specific error message: check for piers with
very soft or descending `Material` branches.

## 11. Errors and limits

`StudWall::new` returns `StudWallError` (with the offending item named) for:
- non-positive or non-finite `L`, `H`, `EA` or `EI` values, or an invalid law parameter;
- studs not strictly increasing, or not starting at `0` and ending at `L`;
- openings out of range (`0 < x0 < x1 < L`, `0 < head_height < H`), overlapping, or without a full-height
  stud between them;
- with `Direct`, a top station not at a stud or cripple; a bottom station inside a door;
- corners or stations reusing a node; a singular internal stiffness (a layout bug, never regularized);
- too many exterior nodes: 256 for linear walls; for `O86`/`Material` walls, 32 exterior nodes and 16 piers.

At analysis build, `validate` checks the corner rectangle and that each station node lies at its plate point.

## 12. Accuracy and verification

Each check below is a test in `core/tests/elements/stud_wall*.rs`. The docs site test report lists them all.

| Check | Result |
|---|---|
| Unsheathed and loaded walls vs the same frame built from `Truss`/`ElasticBeamColumn` elements | 1e-9 |
| Gravity takedowns (uniform, point, window, door, factored combinations, stacked storeys) vs hand calculations | 1e-9 (1e-6 for the continuous top plate) |
| One pier vs SDPWS/O86 bending + shear terms | 1e-6 (rigid stand-ins) |
| One pier, O86 law, several load levels vs O86 11.7.1.2 | 1e-6 |
| 2- and 5-storey walls vs O86 A.11.7.1 (code idealization) | 1e-5 |
| Pier M/jd vs the O86 hold-down force, every variant | 1e-6 |
| O86 law vs secant stiffness at the design load | 1e-8 (real plates) |
| Consistent tangent vs finite differences (O86, hysteretic) | 1e-6 |
| Cyclic `Steel01` sheathing vs the material driven directly | 1e-5 |
| Window wall with king studs and hold-downs at every pier end vs separate pier walls + drag strut (`Direct` and rim-stiff `Beam` top plates) | 1e-10 |
| FloorBand constant shear response, objectivity and rigid-body modes vs energy hand equations | 1e-12 |
| One-storey wall on FloorBand: base moment versus `V(H+d)` | 1e-9 |
| 2- and 5-storey walls with 254 mm floor bands vs an exact hand walk (A.11.7.1 terms, rigid band) | 1e-5 |
| Window wall under a floor band: pier moments vs `V_p(H+d)` (rim and top-plate drag couple) | 1% |
| Single-storey wall on `Anchorage` links vs O86 11.7.1.2 with `d_a = slack + T/k + crush(C)` | 1e-5 |
| Window wall with nonlinear `Anchorage` at all pier ends vs separate pier walls | 1e-10 |
| 3D walls (along X, Z, 30°: window wall, takedown, O86 law, L corner) vs the 2D wall | 1e-10 (1e-9 nonlinear) |

Realistic variants of the 2- and 5-storey walls, compared with O86 A.11.7.1:

| Variant | Deflection ÷ O86 | End-stud force ÷ M/jd |
|---|---|---|
| Interior studs, plate anchored along its length, no bearing | 1.000 | 1.00 |
| Interior studs bearing (foundation or stations) | 0.98–0.99 | 0.89–0.96 |
| Bottom plate anchored only at the corners | 1.01–1.05 | 1.00 |

The pier `chords` (M/jd) equal the O86 hold-down force in every variant.
