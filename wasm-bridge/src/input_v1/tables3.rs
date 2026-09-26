//! Spatial ("space: 3") counterpart to `tables.rs` — same wire-format
//! conventions (dense 0-based row indices, every table always present,
//! possibly empty), just for `Domain3`/`Element3`'s six-DOF-per-node
//! profile instead of the planar three. `materials`/`LoadPatternTable`/
//! `NodalLoadTable`/`IntegrationSpec`/`TimeSeriesSpec` are entirely
//! dimension-agnostic (a material is a scalar strain->stress law regardless
//! of profile, a nodal load is just `(pattern, node, dof, value)`) and so
//! are reused as-is from `tables.rs` rather than duplicated here — only
//! element/geometry shapes that genuinely differ in 3D get their own type.

use serde::{Deserialize, Serialize};

/// Node table: `coords` stride 3 (x, y, z); `fixed` one bitmask byte per
/// node (bit 0..5 = ux, uy, uz, rx, ry, rz, matching `SpatialDof`'s order);
/// mass is sparse, addressed via a parallel node-index array since most
/// nodes carry none.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NodeTable3 {
    pub coords: Vec<f64>,
    pub fixed: Vec<u8>,
    pub mass_node_index: Vec<u32>,
    /// Stride 6 (mass_ux, mass_uy, mass_uz, mass_rx, mass_ry, mass_rz),
    /// parallel to `mass_node_index`.
    pub mass: Vec<f64>,
}

impl NodeTable3 {
    pub fn node_count(&self) -> usize {
        self.coords.len() / 3
    }
}

/// `GeomTransf3`'s wire mirror. No `Corotational3` — spatial corotational
/// geometry doesn't exist in `core` yet (`GeomTransf3`'s own doc comment).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum TransformSpec3 {
    Linear3 { vec_xz: [f64; 3] },
    PDelta3 { vec_xz: [f64; 3] },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TrussTable3 {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub area: Vec<f64>,
    /// Index into the material arena.
    pub material: Vec<u32>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ElasticBeamColumnTable3 {
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

/// Shared shape for `DispBeamColumn3` and `ForceBeamColumn3` — both are one
/// prismatic biaxial fiber section (see [`FiberTable3`]) replicated across
/// integration points by `core`'s own element constructors. Unlike
/// [`super::tables::FiberBeamColumnTable`], there is no `corotational`
/// field: neither spatial fiber element supports it yet (`core`'s
/// `DispBeamColumn3`/`ForceBeamColumn3` expose no `.with_corotational()`),
/// and no separate `transform` field — `g`/`j`/`vec_xz` are passed directly
/// to the constructor rather than wrapped in a `GeomTransf3`, since these
/// elements only ever use `Linear`-type geometry.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FiberBeamColumnTable3 {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    pub g: Vec<f64>,
    /// Torsional constant — fiber sections don't carry torsion (see
    /// [`FiberTable3`]'s doc comment), so it's supplied here as a decoupled
    /// elastic `G*J` term, same as `core::DispBeamColumn3`/`ForceBeamColumn3`.
    pub j: Vec<f64>,
    /// A vector not parallel to the member axis, fixing the local y/z
    /// orientation — see `core::GeomTransf3::Linear3`'s doc comment.
    pub vec_xz: Vec<[f64; 3]>,
    /// Index into `FiberTable3::section_offsets`.
    pub fiber_section: Vec<u32>,
    pub integration: Vec<super::tables::IntegrationSpec>,
    pub density: Vec<f64>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ZeroLengthTable3 {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Sparse `(zero_length row index, dof, material arena index)` — most
    /// DOFs on a given `ZeroLength3` carry no material at all.
    pub materials: Vec<(u32, u8, u32)>,
    /// Sparse `(zero_length row index, normal_dof, shear_dof_0, shear_dof_1,
    /// mu, k0, b)` — at most one entry per row (`core::ZeroLength3` allows
    /// only one `Friction3`, which independently couples *two* shear DOFs
    /// to the same normal force; see `Friction3`'s doc comment).
    #[allow(clippy::type_complexity)]
    pub friction: Vec<(u32, u8, u8, u8, f64, f64, f64)>,
}

/// A `ZeroLength3` driven by a coupled `FiberSection3` (axial + biaxial
/// moment response, `[ux, ry, rz]`) instead of `ZeroLength3`'s independent
/// per-DOF materials — see `core::ZeroLengthSection3`'s doc comment.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ZeroLengthSectionTable3 {
    pub node_i: Vec<u32>,
    pub node_j: Vec<u32>,
    /// Index into `FiberTable3::section_offsets`.
    pub fiber_section: Vec<u32>,
    /// Sparse `(zero_length_section row index, dof, material arena index)` —
    /// an independent spring for a DOF the section doesn't drive: `uy`/`uz`
    /// (shear, dofs 1/2) or `rx` (torsion, dof 3).
    pub materials: Vec<(u32, u8, u32)>,
}

/// Fibers for every `DispBeamColumn3`/`ForceBeamColumn3` section, flattened
/// and offset-indexed exactly like [`super::tables::FiberTable`], but with
/// both `y` and `z` coordinates (biaxial bending) — torsion is deliberately
/// excluded from the fiber loop in `core` (`FiberSection3`'s doc comment)
/// and supplied instead as each owning table's own `g`/`j` fields.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FiberTable3 {
    pub section_offsets: Vec<u32>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub area: Vec<f64>,
    /// Index into the material arena, parallel to `y`/`z`/`area`.
    pub material: Vec<u32>,
}

impl FiberTable3 {
    pub fn section_count(&self) -> usize {
        self.section_offsets.len().saturating_sub(1)
    }
}

/// `core::Axis3`'s wire mirror — core's own type carries no `serde` derive
/// (`Axis3`'s doc comment: it's a plain `repr(usize)` enum for internal
/// indexing, not wire-facing), so this is decode's own translation, the
/// same reason `TransformSpec`/`TransformSpec3` mirror `GeomTransf`/
/// `GeomTransf3` rather than deriving on the `core` type directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum Axis3Spec {
    X,
    Y,
    Z,
}

/// Spatial counterpart to [`super::tables::EqualDofTable`] — same shape,
/// dofs range `0..6` instead of `0..3`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EqualDofTable3 {
    pub retained: Vec<u32>,
    pub constrained: Vec<u32>,
    pub dofs: Vec<(u32, u8)>,
}

/// Spatial rigid diaphragm (`core::Domain3::rigid_diaphragm_about`'s doc
/// comment): row `i` ties every node listed against it in `constrained`'s
/// two in-plane translational dofs (perpendicular to `normal[i]`) to
/// `retained[i]`'s same in-plane translations plus the lever-arm rotation
/// term. `constrained` is sparse per row, same convention as
/// [`super::tables::RigidDiaphragmTable`].
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RigidDiaphragmTable3 {
    pub retained: Vec<u32>,
    pub normal: Vec<Axis3Spec>,
    pub constrained: Vec<(u32, u32)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum ElementKind3 {
    Truss,
    ElasticBeamColumn,
    DispBeamColumn,
    ForceBeamColumn,
    ZeroLength,
    ZeroLengthSection,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ElementLoadSpec3 {
    /// Only `ElasticBeamColumn3` currently honors this (mirroring planar
    /// `ElementLoadSpec::UniformTransverse`'s own restriction) — biaxial,
    /// since a spatial member bends in both local `y` and `z`.
    UniformTransverse { wy: f64, wz: f64 },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ElementLoadTable3 {
    pub pattern: Vec<u32>,
    pub element_kind: Vec<ElementKind3>,
    /// Row index into the table named by the parallel `element_kind` entry.
    pub element_index: Vec<u32>,
    pub load: Vec<ElementLoadSpec3>,
    /// See [`super::tables::NodalLoadTable::stage`] — same per-stage
    /// registration timing.
    pub stage: Vec<u32>,
}
