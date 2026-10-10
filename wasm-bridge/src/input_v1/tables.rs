//! Bulk per-(entity-kind) tables of the wire format (`docs/input-format.md`), one
//! set for both profiles (`header.ndm` 2 or 3). Dimension-agnostic entities
//! (nodes, trusses, zero-lengths, fibers, constraints, loads) have one table;
//! formulations that differ between 2D and 3D (beam-columns) have one table
//! each, suffixed `2d`/`3d`, and the table of the other profile must be empty
//! (`DecodeError::TableNotInProfile`).
//!
//! Node/element/pattern references inside these tables are already
//! resolved to dense 0-based row indices by the (not-yet-written)
//! `pysees` compiler — the "duplicate/missing tags and every reference
//! between entities" validation the compiler is responsible for
//! happens upstream of decode, so decode only ever sees plain array
//! indices, never sparse OpenSees-style tags.
//!
//! Each field is a plain `Vec`, not yet the transferable `Float64Array`/
//! `Uint8Array` values the wire format ultimately specifies for
//! `postMessage`: `boundary.rs`'s `serde-wasm-bindgen` decoding today
//! accepts an ordinary JS array for each of these (or a typed array,
//! copied element-by-element) rather than transferring one. Once real
//! transferable-array support lands, each `Vec<f64>`/`Vec<u32>` here is
//! exactly the shape a `Float64Array`/`Uint32Array` copies into. Every table
//! of [`super::CarapaceInputV1`] may be omitted, which is the same as empty.
//!
//! `#[serde(rename_all = "camelCase")]` throughout so the JS/TS shape
//! matches the PySees naming (`nodeI`, not `node_i`).

use serde::{Deserialize, Serialize};

/// Node table: `coords` stride `ndm` (x, y[, z]); `fixed` one bitmask byte per
/// node (bit `k` = DOF `k`: ux, uy, rz in 2D; ux, uy, uz, rx, ry, rz in 3D);
/// mass is sparse, addressed via a parallel node-index array since most nodes
/// carry none.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct NodeTable {
    pub coords: Vec<f64>,
    pub fixed: Vec<u8>,
    pub mass_node_index: Vec<u32>,
    /// Stride `ndf` (3 in 2D, 6 in 3D), parallel to `mass_node_index`.
    pub mass: Vec<f64>,
}

impl NodeTable {
    pub fn node_count(&self, ndm: usize) -> usize {
        self.coords.len() / ndm
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum TransformSpec {
    Linear,
    PDelta,
    Corotational,
}

/// `GeomTransf3`'s wire mirror. No `Corotational3` — spatial corotational
/// geometry doesn't exist in `core` yet (`GeomTransf3`'s own doc comment).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum TransformSpec3 {
    Linear3 { vec_xz: [f64; 3] },
    PDelta3 { vec_xz: [f64; 3] },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum IntegrationSpec {
    Legendre { points: u32 },
    Lobatto { points: u32 },
}

/// Same fields in 2D and 3D.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct TrussTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub area: Vec<f64>,
    /// Index into the material arena.
    pub material: Vec<u32>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ElasticBeamColumn2dTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub e: Vec<f64>,
    pub a: Vec<f64>,
    pub iz: Vec<f64>,
    pub transform: Vec<TransformSpec>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ElasticBeamColumn3dTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub e: Vec<f64>,
    pub g: Vec<f64>,
    pub a: Vec<f64>,
    /// Torsional constant.
    pub j: Vec<f64>,
    pub iy: Vec<f64>,
    pub iz: Vec<f64>,
    pub transform: Vec<TransformSpec3>,
    pub density: Vec<f64>,
}

/// Shared shape for 2D `DispBeamColumn` and `ForceBeamColumn` — both are one
/// prismatic fiber section (see [`FiberTable`]) replicated across
/// integration points by `core`'s own element constructors.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct FiberBeamColumn2dTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Index into `FiberTable::section_offsets`.
    pub fiber_section: Vec<u32>,
    pub integration: Vec<IntegrationSpec>,
    pub corotational: Vec<bool>,
    pub density: Vec<f64>,
}

/// Shared shape for `DispBeamColumn3` and `ForceBeamColumn3` — both are one
/// prismatic biaxial fiber section (see [`FiberTable`]) replicated across
/// integration points by `core`'s own element constructors. Unlike
/// [`FiberBeamColumn2dTable`], there is no `corotational` field: neither
/// spatial fiber element supports it yet, and no separate `transform` field —
/// `g`/`j`/`vec_xz` are passed directly to the constructor, since these
/// elements only ever use `Linear`-type geometry.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct FiberBeamColumn3dTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub g: Vec<f64>,
    /// Torsional constant — fiber sections don't carry torsion (see
    /// [`FiberTable`]'s doc comment), so it's supplied here as a decoupled
    /// elastic `G*J` term, same as `core::DispBeamColumn3`/`ForceBeamColumn3`.
    pub j: Vec<f64>,
    /// A vector not parallel to the member axis, fixing the local y/z
    /// orientation — see `core::GeomTransf3::Linear3`'s doc comment.
    pub vec_xz: Vec<[f64; 3]>,
    /// Index into `FiberTable::section_offsets`.
    pub fiber_section: Vec<u32>,
    pub integration: Vec<IntegrationSpec>,
    pub density: Vec<f64>,
}

/// One sparse `friction` row of [`ZeroLengthTable`]: the row's one `Friction`
/// (2D) / `Friction3` (3D) — see their doc comments for the field meanings.
/// `normal_dof` must already have a `materials` entry on the same row, and the
/// shear DOFs must not (checked by `core` via `debug_assert`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct FrictionRow {
    pub row: u32,
    pub normal_dof: u8,
    /// One shear DOF in 2D, two in 3D (coupled to the same normal force).
    pub shear_dofs: Vec<u8>,
    pub mu: f64,
    pub k0: f64,
    pub b: f64,
}

/// One sparse `orient` row — OpenSees' `-orient x1 x2 x3 [yp1 yp2 yp3]`: local
/// x is `x`; in 3D local z is `x × yp` and local y completes the frame; in 2D
/// `yp` is absent, `x[2]` must be 0 and local y is local x turned 90 degrees
/// counter-clockwise. Every DOF of the row is evaluated along those axes.
/// Rows without an entry use the global axes; at most one entry per row.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct OrientRow {
    pub row: u32,
    pub x: [f64; 3],
    /// Required in 3D, must be absent in 2D.
    #[serde(default)]
    #[tsify(optional)]
    pub yp: Option<[f64; 3]>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ZeroLengthTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Sparse `(zero_length row index, dof, material arena index)` — most
    /// DOFs on a given `ZeroLength` carry no material at all.
    pub materials: Vec<(u32, u8, u32)>,
    /// At most one entry per row (`core::ZeroLength` allows only one friction).
    #[serde(default)]
    #[tsify(optional)]
    pub friction: Vec<FrictionRow>,
    #[serde(default)]
    #[tsify(optional)]
    pub orient: Vec<OrientRow>,
}

/// A `ZeroLength` driven by a coupled `FiberSection` (axial + moment
/// response, `[ux, rz]` in 2D; axial + biaxial moment, `[ux, ry, rz]` in 3D)
/// instead of `ZeroLength`'s independent per-DOF materials — see
/// `core::ZeroLengthSection`'s doc comment.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ZeroLengthSectionTable {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Index into `FiberTable::section_offsets`.
    pub fiber_section: Vec<u32>,
    /// Sparse `(zero_length_section row index, dof, material arena index)` —
    /// an independent spring for a DOF the section has no resultant for: `uy`
    /// in 2D; `uy`/`uz` (shear) or `rx` (torsion) in 3D.
    pub materials: Vec<(u32, u8, u32)>,
    #[serde(default)]
    #[tsify(optional)]
    pub orient: Vec<OrientRow>,
}

/// Fibers for every `DispBeamColumn`/`ForceBeamColumn` section, flattened
/// and offset-indexed: section `k` occupies
/// `section_offsets[k]..section_offsets[k + 1]` in `y`/`z`/`area`/`material`.
/// `section_offsets` therefore has `num_sections + 1` entries. `z` is empty in
/// a 2D model and parallel to `y` in a 3D model (biaxial bending) — torsion is
/// deliberately excluded from the fiber loop in `core` (`FiberSection3`'s doc
/// comment) and supplied instead as each owning table's own `g`/`j` fields.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct FiberTable {
    pub section_offsets: Vec<u32>,
    pub y: Vec<f64>,
    #[serde(default)]
    #[tsify(optional)]
    pub z: Vec<f64>,
    pub area: Vec<f64>,
    /// Index into the material arena, parallel to `y`/`area`.
    pub material: Vec<u32>,
}

impl FiberTable {
    pub fn section_count(&self) -> usize {
        self.section_offsets.len().saturating_sub(1)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum TimeSeriesSpec {
    Constant,
    Linear { slope: f64 },
    Path { times: Vec<f64>, factors: Vec<f64> },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct LoadPatternTable {
    pub series: Vec<TimeSeriesSpec>,
    pub scale_factor: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct NodalLoadTable {
    /// Index into `LoadPatternTable`.
    pub pattern: Vec<u32>,
    /// Index into `NodeTable`.
    pub node: Vec<u32>,
    pub dof: Vec<u8>,
    pub value: Vec<f64>,
    /// Index into `SequenceSpec::stages`: this load is only registered on
    /// the `Domain` when that stage starts, not at decode time. This
    /// matters whenever an earlier stage's `LoadControl` shares one
    /// pseudo-time with every unfrozen pattern (core/src/analysis/
    /// integrator.rs) — a pattern's reference load must not exist yet if
    /// an earlier stage isn't meant to ramp it too (the same reason
    /// core/tests/elements/force_beam.rs's native two-phase test adds its
    /// lateral pattern's load only between phases, not upfront).
    pub stage: Vec<u32>,
}

/// A plane (2D continuum) material, referenced by `triangles`/`quads` rows
/// through the `planeMaterials` arena. Strain is `[eps_x, eps_y, gamma_xy]`
/// with engineering shear.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum PlaneMaterialSpec {
    /// Isotropic linear elastic.
    Isotropic {
        e: f64,
        nu: f64,
        state: PlaneStateSpec,
    },
    /// Orthotropic plane stress; `angle` is the counter-clockwise angle in
    /// radians from global x to material axis 1.
    Orthotropic {
        ex: f64,
        ey: f64,
        nu_xy: f64,
        g_xy: f64,
        angle: f64,
    },
    /// A symmetric positive-definite `D`, packed upper triangle row by row:
    /// `[d11, d12, d13, d22, d23, d33]`.
    ElasticMatrix { d: [f64; 6] },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum PlaneStateSpec {
    PlaneStress,
    PlaneStrain,
}

/// 3-node constant-strain triangles (2D models only). `nodeIds` has stride 3
/// (counter-clockwise); `material` indexes `planeMaterials`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct TriangleTable {
    pub node_ids: Vec<u32>,
    pub thickness: Vec<f64>,
    pub material: Vec<u32>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum Quad4FormulationSpec {
    Full,
    /// Wilson-Taylor incompatible modes; linear materials only.
    Enhanced,
}

/// 4-node bilinear quadrilaterals (2D models only). `nodeIds` has stride 4
/// (counter-clockwise); `material` indexes `planeMaterials`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct QuadTable {
    pub node_ids: Vec<u32>,
    pub thickness: Vec<f64>,
    pub material: Vec<u32>,
    pub density: Vec<f64>,
    pub formulation: Vec<Quad4FormulationSpec>,
}

/// A shell section, referenced by `shell4s` rows through the `shellSections` arena. Generalized
/// strain and resultants follow OpenSees' `ElasticMembranePlateSection` ordering.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ShellSectionSpec {
    /// Isotropic elastic membrane + plate bending + transverse shear: modulus, Poisson ratio,
    /// thickness and mass density (per unit volume).
    ElasticMembranePlate { e: f64, nu: f64, h: f64, rho: f64 },
}

/// 3-node DKT/Allman shells (3D models only). `nodeIds` has stride 3 (the node order sets the local normal by the
/// right-hand rule); `section` indexes `shellSections`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct Shell3Table {
    pub node_ids: Vec<u32>,
    pub section: Vec<u32>,
}

/// 4-node MITC4 shells (3D models only). `nodeIds` has stride 4 (the node order fixes the
/// local normal by the right-hand rule); `section` indexes `shellSections`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct Shell4Table {
    pub node_ids: Vec<u32>,
    pub section: Vec<u32>,
}

/// Every element formulation, 2D and 3D. A kind that does not belong to the
/// model's `ndm` is `DecodeError::ElementKindNotInProfile`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum ElementKind {
    Truss,
    ElasticBeamColumn2d,
    ElasticBeamColumn3d,
    DispBeamColumn2d,
    DispBeamColumn3d,
    ForceBeamColumn2d,
    ForceBeamColumn3d,
    ZeroLength,
    ZeroLengthSection,
    Tri3,
    Quad4,
    Shell3,
    Shell4,
}

impl ElementKind {
    pub const ALL: [ElementKind; 13] = [
        ElementKind::Truss,
        ElementKind::ElasticBeamColumn2d,
        ElementKind::ElasticBeamColumn3d,
        ElementKind::DispBeamColumn2d,
        ElementKind::DispBeamColumn3d,
        ElementKind::ForceBeamColumn2d,
        ElementKind::ForceBeamColumn3d,
        ElementKind::ZeroLength,
        ElementKind::ZeroLengthSection,
        ElementKind::Tri3,
        ElementKind::Quad4,
        ElementKind::Shell3,
        ElementKind::Shell4,
    ];

    /// The kind's table name as `DecodeError`s report it (snake_case).
    pub fn table(self) -> &'static str {
        match self {
            ElementKind::Truss => "trusses",
            ElementKind::ElasticBeamColumn2d => "elastic_beam_columns_2d",
            ElementKind::ElasticBeamColumn3d => "elastic_beam_columns_3d",
            ElementKind::DispBeamColumn2d => "disp_beam_columns_2d",
            ElementKind::DispBeamColumn3d => "disp_beam_columns_3d",
            ElementKind::ForceBeamColumn2d => "force_beam_columns_2d",
            ElementKind::ForceBeamColumn3d => "force_beam_columns_3d",
            ElementKind::ZeroLength => "zero_lengths",
            ElementKind::ZeroLengthSection => "zero_length_sections",
            ElementKind::Tri3 => "triangles",
            ElementKind::Quad4 => "quads",
            ElementKind::Shell3 => "shell3s",
            ElementKind::Shell4 => "shell4s",
        }
    }

    /// `Some(2)`/`Some(3)` for a formulation that exists only in that profile,
    /// `None` for one shared by both.
    pub fn only_ndm(self) -> Option<u8> {
        match self {
            ElementKind::Truss | ElementKind::ZeroLength | ElementKind::ZeroLengthSection => None,
            ElementKind::ElasticBeamColumn2d
            | ElementKind::DispBeamColumn2d
            | ElementKind::ForceBeamColumn2d
            | ElementKind::Tri3
            | ElementKind::Quad4 => Some(2),
            ElementKind::ElasticBeamColumn3d
            | ElementKind::DispBeamColumn3d
            | ElementKind::ForceBeamColumn3d
            | ElementKind::Shell3
            | ElementKind::Shell4 => Some(3),
        }
    }

    pub fn in_profile(self, ndm: u8) -> bool {
        self.only_ndm().is_none_or(|only| only == ndm)
    }

    /// The 2D continuum elements (they take body and edge loads).
    pub fn is_continuum(self) -> bool {
        matches!(self, ElementKind::Tri3 | ElementKind::Quad4)
    }

    /// The 3D shell elements (they take self-weight and pressure loads).
    pub fn is_shell(self) -> bool {
        matches!(self, ElementKind::Shell3 | ElementKind::Shell4)
    }

    /// Number of edges of a continuum element (zero for every other kind).
    pub fn edge_count(self) -> u8 {
        match self {
            ElementKind::Tri3 => 3,
            ElementKind::Quad4 => 4,
            _ => 0,
        }
    }

    /// Whether `core` applies a beam `Uniform` load to this kind (the beam-columns).
    pub fn accepts_uniform_load(self) -> bool {
        matches!(
            self,
            ElementKind::ElasticBeamColumn2d
                | ElementKind::ElasticBeamColumn3d
                | ElementKind::DispBeamColumn2d
                | ElementKind::DispBeamColumn3d
                | ElementKind::ForceBeamColumn2d
                | ElementKind::ForceBeamColumn3d
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ElementLoadSpec {
    /// Uniform load (force/length) in the element's *local* axes: `wx`
    /// along the member axis, `wy` transverse in local +y, and in 3D `wz`
    /// transverse in local +z (`vec_xz` defines `y`/`z`). `wz` is absent in a
    /// 2D model and defaults to 0 in a 3D one.
    ///
    /// Only the beam-column kinds currently honor this (`core`'s
    /// `Element::form_load_vector` falls through to zero for every other
    /// kind) — decode rejects any other `element_kind` explicitly rather
    /// than silently accepting a load that will never apply. Several loads
    /// on one element in one pattern sum.
    Uniform {
        wx: f64,
        wy: f64,
        #[serde(default)]
        #[tsify(optional)]
        wz: Option<f64>,
    },
    /// Body force per unit volume in global axes (continuum elements, 2D).
    Body { bx: f64, by: f64 },
    /// Traction per unit length of edge `edge` in global axes; edge `k` joins
    /// local node `k` to `k + 1` (continuum elements, 2D).
    EdgeTraction { edge: u8, tx: f64, ty: f64 },
    /// Pressure on edge `edge`, positive into the element along the inward
    /// normal (continuum elements, 2D).
    EdgePressure { edge: u8, pressure: f64 },
    /// Pressure per unit area on a shell, positive along its local normal `e3` (right-hand rule
    /// from the node order). Applied as consistent nodal forces.
    ShellPressure { pressure: f64 },
    /// Body acceleration on a shell in global axes (a gravity vector `g` points down); the force
    /// per unit area is `rho h b`. Note OpenSees' `ShellMITC4` `-selfWeight` has the opposite sign.
    ShellBody { bx: f64, by: f64, bz: f64 },
}

/// Identity multi-point constraints (`core::Domain::equal_dof`'s doc
/// comment): row `i` ties `constrained[i]`'s dofs listed in `dofs` exactly
/// to the same dofs of `retained[i]`. `dofs` is sparse per row — `(row,
/// dof)` pairs, mirroring `ZeroLengthTable::materials`'s own sparse
/// convention — since most ties only ever list one or two dofs, not every
/// dof a node has.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct EqualDofTable {
    pub retained: Vec<u32>,
    pub constrained: Vec<u32>,
    pub dofs: Vec<(u32, u8)>,
}

/// `core::Axis3`'s wire mirror — core's own type carries no `serde` derive
/// (`Axis3`'s doc comment: it's a plain `repr(usize)` enum for internal
/// indexing, not wire-facing), so this is decode's own translation, the
/// same reason `TransformSpec`/`TransformSpec3` mirror `GeomTransf`/
/// `GeomTransf3` rather than deriving on the `core` type directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub enum Axis3Spec {
    X,
    Y,
    Z,
}

/// Rigid diaphragm: row `i` ties every node listed against it in
/// `constrained` to `retained[i]` — in 2D the `ux` dof (`core::Domain::
/// rigid_diaphragm`'s doc comment); in 3D the two in-plane translational dofs
/// (perpendicular to `normal[i]`) plus the lever-arm rotation term
/// (`core::Domain3::rigid_diaphragm_about`). `constrained` is sparse per row —
/// `(row, node index)` pairs — since a diaphragm's node count varies.
/// `normal` is 3D only: omitted or empty means every row's normal is Y; a
/// 2D model must leave it empty.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct RigidDiaphragmTable {
    pub retained: Vec<u32>,
    #[serde(default)]
    #[tsify(optional)]
    pub normal: Vec<Axis3Spec>,
    pub constrained: Vec<(u32, u32)>,
}

/// Rigid links (`core::Domain::rigid_link`'s doc comment): `slave` moves with
/// `master` as a rigid body, small-rotation lever arm taken from the node
/// coordinates.
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct RigidLinkTable {
    pub master: Vec<u32>,
    pub slave: Vec<u32>,
}

/// General linear multi-point constraints (`core::Domain::add_constraint`):
/// row `i` defines `u[slaveNode[i], slaveDof[i]] = sum coeff * u[node, dof]`
/// over the terms `termOffsets[i]..termOffsets[i + 1]` of the flattened sparse
/// `termNode`/`termDof`/`termCoeff` arrays. `termOffsets` has `rows + 1`
/// entries (or none at all for an empty table).
#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct LinearConstraintTable {
    pub slave_node: Vec<u32>,
    pub slave_dof: Vec<u8>,
    pub term_offsets: Vec<u32>,
    pub term_node: Vec<u32>,
    pub term_dof: Vec<u8>,
    pub term_coeff: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct ElementLoadTable {
    pub pattern: Vec<u32>,
    pub element_kind: Vec<ElementKind>,
    /// Row index into the table named by the parallel `element_kind` entry.
    pub element_index: Vec<u32>,
    pub load: Vec<ElementLoadSpec>,
    /// See [`NodalLoadTable::stage`] — same per-stage registration timing.
    pub stage: Vec<u32>,
}
