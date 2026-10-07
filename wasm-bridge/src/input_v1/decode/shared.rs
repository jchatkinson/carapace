//! Decoding shared by both profiles: everything that is bookkeeping around the
//! per-profile element constructors — nodes (generic over `<NDIM, NDOF>`),
//! fibers, equal-DOF/diaphragm/rigid-link/linear constraints, sparse per-row
//! grouping, load patterns, nodal and element loads, and recorder resolution
//! through an `(ElementKind, row) -> E::Id` map ([`ElementIds`]).

use std::collections::HashMap;

use carapace_core::model::{
    BeamIntegration, Domain, ElementLoadComponents, ElementOps, LoadPatternId, Material, Node,
    Orientation,
};
use slotmap::Key;

use super::stages::{compile_stages, load_series_of};
use crate::input_v1::error::DecodeError;
use crate::input_v1::sequence::RecorderSpec;
use crate::input_v1::session::{ModelSession, ResolvedRecorder};
use crate::input_v1::tables::{
    Axis3Spec, ElementKind, ElementLoadSpec, EqualDofTable, FiberTable, FrictionRow,
    IntegrationSpec, LinearConstraintTable, LoadPatternTable, NodeTable, OrientRow,
    RigidDiaphragmTable, RigidLinkTable,
};
use crate::input_v1::CarapaceInputV1;

pub(super) fn check_dof_within(
    table: &'static str,
    row: u32,
    dof: u8,
    ndof: u8,
) -> Result<(), DecodeError> {
    if dof < ndof {
        Ok(())
    } else {
        Err(DecodeError::InvalidDof { table, row, dof })
    }
}

fn invalid_row(table: &'static str, row: usize, reason: &'static str) -> DecodeError {
    DecodeError::InvalidRow {
        table,
        row: row as u32,
        reason,
    }
}

/// A node-table row to the profile's node handle.
pub(super) fn node_lookup<NId: Copy>(
    node_ids: &[NId],
) -> impl Fn(u32, &'static str) -> Result<NId, DecodeError> + '_ {
    move |row, table| {
        node_ids
            .get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownNodeIndex { table, row })
    }
}

pub(super) fn material_lookup(
    materials: &[Material],
) -> impl Fn(u32, &'static str) -> Result<Material, DecodeError> + '_ {
    move |row, table| {
        materials
            .get(row as usize)
            .cloned()
            .ok_or(DecodeError::UnknownMaterialIndex { table, row })
    }
}

/// Rejects a non-empty table that belongs to the other profile.
pub(super) fn reject_tables(tables: &[(&'static str, bool)]) -> Result<(), DecodeError> {
    for &(table, is_empty) in tables {
        if !is_empty {
            return Err(DecodeError::TableNotInProfile { table });
        }
    }
    Ok(())
}

/// `fibers.z`: empty in 2D, parallel to `y` in 3D.
pub(super) fn check_fibers(table: &FiberTable, ndm: u8) -> Result<(), DecodeError> {
    let expected = if ndm == 3 { table.y.len() } else { 0 };
    if table.z.len() != expected {
        return Err(invalid_row(
            "fibers",
            0,
            if ndm == 3 {
                "z must be parallel to y in a 3D model"
            } else {
                "z must be empty in a 2D model"
            },
        ));
    }
    Ok(())
}

pub(super) fn add_nodes<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, NId, E>,
    table: &NodeTable,
) -> Result<Vec<NId>, DecodeError>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    let node_count = table.node_count(NDIM);
    if table.coords.len() != node_count * NDIM {
        return Err(invalid_row(
            "nodes",
            0,
            "coords length is not a multiple of ndm",
        ));
    }
    if table.fixed.len() != node_count {
        return Err(invalid_row(
            "nodes",
            0,
            "fixed must have one entry per node",
        ));
    }
    if table.mass.len() != table.mass_node_index.len() * NDOF {
        return Err(invalid_row(
            "nodes",
            0,
            "mass must have ndf entries per massNodeIndex",
        ));
    }
    let mut mass_by_node: HashMap<u32, &[f64]> = HashMap::new();
    for (k, &node_index) in table.mass_node_index.iter().enumerate() {
        mass_by_node.insert(node_index, &table.mass[k * NDOF..(k + 1) * NDOF]);
    }

    let mut node_ids = Vec::with_capacity(node_count);
    for i in 0..node_count {
        let mut coords = [0.0; NDIM];
        coords.copy_from_slice(&table.coords[i * NDIM..(i + 1) * NDIM]);
        let mut node = Node::<NDIM, NDOF>::new(coords);
        let fixed_bits = table.fixed[i];
        for dof in 0..NDOF.min(8) {
            if fixed_bits & (1 << dof) != 0 {
                node = node.fix(dof);
            }
        }
        if let Some(mass) = mass_by_node.get(&(i as u32)) {
            for (dof, &value) in mass.iter().enumerate() {
                if value != 0.0 {
                    node = node.with_mass(dof, value);
                }
            }
        }
        node_ids.push(domain.add_node(node));
    }
    Ok(node_ids)
}

/// Decodes [`EqualDofTable`] against an already-populated `domain` (nodes
/// only — constraints don't reference elements or patterns), calling
/// `core::Domain::equal_dof` once per row.
pub(super) fn add_equal_dofs<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, NId, E>,
    table: &EqualDofTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
) -> Result<(), DecodeError>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    let mut dofs_by_row: Vec<Vec<usize>> = vec![Vec::new(); table.retained.len()];
    for &(row, dof) in &table.dofs {
        check_dof_within("equal_dofs", row, dof, NDOF as u8)?;
        dofs_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownConstraintRow {
                table: "equal_dofs",
                row,
            })?
            .push(dof as usize);
    }
    #[allow(clippy::needless_range_loop)] // parallel-indexes retained/constrained/dofs_by_row
    for i in 0..table.retained.len() {
        let retained = node_at(table.retained[i], "equal_dofs")?;
        let constrained = node_at(table.constrained[i], "equal_dofs")?;
        domain.equal_dof(retained, constrained, &dofs_by_row[i]);
    }
    Ok(())
}

/// Decodes [`RigidDiaphragmTable`], handing each row to `apply(retained,
/// constrained, normal)` — the 2D and 3D `Domain` constructors differ
/// (`rigid_diaphragm` vs `rigid_diaphragm_about`), so the domain call is the
/// caller's. `normal` defaults to Y; a 2D model must not give one.
pub(super) fn add_rigid_diaphragms<NId: Copy>(
    table: &RigidDiaphragmTable,
    ndm: u8,
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
    mut apply: impl FnMut(NId, &[NId], Axis3Spec),
) -> Result<(), DecodeError> {
    const TABLE: &str = "rigid_diaphragms";
    if ndm == 2 && !table.normal.is_empty() {
        return Err(invalid_row(TABLE, 0, "normal is only valid in a 3D model"));
    }
    if !table.normal.is_empty() && table.normal.len() != table.retained.len() {
        return Err(invalid_row(TABLE, 0, "normal must have one entry per row"));
    }
    let mut constrained_by_row: Vec<Vec<u32>> = vec![Vec::new(); table.retained.len()];
    for &(row, node) in &table.constrained {
        constrained_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownConstraintRow { table: TABLE, row })?
            .push(node);
    }
    #[allow(clippy::needless_range_loop)] // parallel-indexes retained/normal/constrained_by_row
    for i in 0..table.retained.len() {
        let retained = node_at(table.retained[i], TABLE)?;
        let constrained = constrained_by_row[i]
            .iter()
            .map(|&row| node_at(row, TABLE))
            .collect::<Result<Vec<_>, _>>()?;
        let normal = table.normal.get(i).copied().unwrap_or(Axis3Spec::Y);
        apply(retained, &constrained, normal);
    }
    Ok(())
}

/// Decodes [`RigidLinkTable`], handing each `(master, slave)` pair to `apply`
/// (`Domain::rigid_link` exists per profile).
pub(super) fn add_rigid_links<NId: Copy>(
    table: &RigidLinkTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
    mut apply: impl FnMut(NId, NId),
) -> Result<(), DecodeError> {
    const TABLE: &str = "rigid_links";
    if table.master.len() != table.slave.len() {
        return Err(invalid_row(TABLE, 0, "master and slave must be parallel"));
    }
    for i in 0..table.master.len() {
        apply(
            node_at(table.master[i], TABLE)?,
            node_at(table.slave[i], TABLE)?,
        );
    }
    Ok(())
}

/// Decodes [`LinearConstraintTable`] into `Domain::add_constraint` calls.
pub(super) fn add_linear_constraints<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, NId, E>,
    table: &LinearConstraintTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
) -> Result<(), DecodeError>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    const TABLE: &str = "linear_constraints";
    let rows = table.slave_node.len();
    let terms = table.term_node.len();
    if table.slave_dof.len() != rows {
        return Err(invalid_row(
            TABLE,
            0,
            "slaveNode and slaveDof must be parallel",
        ));
    }
    if table.term_dof.len() != terms || table.term_coeff.len() != terms {
        return Err(invalid_row(
            TABLE,
            0,
            "termNode, termDof and termCoeff must be parallel",
        ));
    }
    if rows == 0 && terms == 0 && table.term_offsets.is_empty() {
        return Ok(());
    }
    if table.term_offsets.len() != rows + 1
        || table.term_offsets.last().copied() != Some(terms as u32)
        || table.term_offsets.windows(2).any(|pair| pair[0] > pair[1])
    {
        return Err(invalid_row(
            TABLE,
            0,
            "termOffsets must be rows + 1 non-decreasing offsets ending at the term count",
        ));
    }
    for row in 0..rows {
        check_dof_within(TABLE, row as u32, table.slave_dof[row], NDOF as u8)?;
        let slave = node_at(table.slave_node[row], TABLE)?;
        let range = table.term_offsets[row] as usize..table.term_offsets[row + 1] as usize;
        let row_terms = range
            .map(|k| {
                check_dof_within(TABLE, row as u32, table.term_dof[k], NDOF as u8)?;
                Ok((
                    node_at(table.term_node[k], TABLE)?,
                    table.term_dof[k] as usize,
                    table.term_coeff[k],
                ))
            })
            .collect::<Result<Vec<_>, DecodeError>>()?;
        domain.add_constraint((slave, table.slave_dof[row] as usize), &row_terms);
    }
    Ok(())
}

/// One fiber section out of `table`, building each fiber with `make(i,
/// material)` (`Fiber::new(y, area, ..)` in 2D, `Fiber3::new(y, z, ..)` in 3D).
pub(super) fn build_fiber_section<F>(
    table: &FiberTable,
    section: u32,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
    make: impl Fn(usize, Material) -> F,
) -> Result<Vec<F>, DecodeError> {
    let idx = section as usize;
    let start = *table
        .section_offsets
        .get(idx)
        .ok_or(DecodeError::UnknownFiberSectionIndex { row: section })? as usize;
    let end = *table
        .section_offsets
        .get(idx + 1)
        .ok_or(DecodeError::UnknownFiberSectionIndex { row: section })? as usize;
    (start..end)
        .map(|i| Ok(make(i, material_at(table.material[i], "fibers")?)))
        .collect()
}

/// Sparse `(row, dof, material arena index)` entries grouped by row; a row past
/// the end of the owning table is `UnknownElementIndex`.
pub(super) fn materials_by_row(
    materials: &[(u32, u8, u32)],
    rows: usize,
    table: &'static str,
) -> Result<Vec<Vec<(u8, u32)>>, DecodeError> {
    let mut by_row: Vec<Vec<(u8, u32)>> = vec![Vec::new(); rows];
    for &(row, dof, material_index) in materials {
        by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex { table, row })?
            .push((dof, material_index));
    }
    Ok(by_row)
}

/// Per-row friction entry (at most one per row; a later entry replaces an
/// earlier one). `shear_dofs` must have `shear_count` entries.
pub(super) fn friction_by_row<'a>(
    friction: &'a [FrictionRow],
    rows: usize,
    shear_count: usize,
    table: &'static str,
) -> Result<Vec<Option<&'a FrictionRow>>, DecodeError> {
    let mut by_row = vec![None; rows];
    for entry in friction {
        let slot = by_row
            .get_mut(entry.row as usize)
            .ok_or(DecodeError::UnknownElementIndex {
                table,
                row: entry.row,
            })?;
        if entry.shear_dofs.len() != shear_count {
            return Err(invalid_row(
                table,
                entry.row as usize,
                if shear_count == 1 {
                    "friction needs 1 shear DOF in a 2D model"
                } else {
                    "friction needs 2 shear DOFs in a 3D model"
                },
            ));
        }
        *slot = Some(entry);
    }
    Ok(by_row)
}

/// Per-row `Orientation` from a sparse `orient` table (OpenSees' `-orient`).
/// `None` means "global axes". A row past the end of the owning table is
/// `UnknownElementIndex`; zero or parallel vectors (or, in 2D, `x[2] != 0`) are
/// `InvalidOrientation`; a `yp` in 2D, or none in 3D, is `InvalidRow`.
pub(super) fn orientations_by_row(
    orient: &[OrientRow],
    rows: usize,
    ndm: u8,
    table: &'static str,
) -> Result<Vec<Option<Orientation>>, DecodeError> {
    let mut by_row = vec![None; rows];
    for entry in orient {
        let row = entry.row;
        let slot = by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex { table, row })?;
        let invalid = DecodeError::InvalidOrientation { table, row };
        *slot = Some(if ndm == 2 {
            if entry.yp.is_some() {
                return Err(invalid_row(
                    table,
                    row as usize,
                    "yp is only valid in a 3D model",
                ));
            }
            if entry.x[2] != 0.0 {
                return Err(invalid);
            }
            Orientation::in_plane(entry.x[0], entry.x[1]).map_err(|_| invalid)?
        } else {
            let yp = entry
                .yp
                .ok_or_else(|| invalid_row(table, row as usize, "yp is required in a 3D model"))?;
            Orientation::new(entry.x, yp).map_err(|_| invalid)?
        });
    }
    Ok(by_row)
}

pub(super) fn integration_of(spec: IntegrationSpec) -> BeamIntegration {
    match spec {
        IntegrationSpec::Legendre { points } => BeamIntegration::Legendre {
            points: points as usize,
        },
        IntegrationSpec::Lobatto { points } => BeamIntegration::Lobatto {
            points: points as usize,
        },
    }
}

/// The decoded element handles of every kind, indexed by `(ElementKind, row)`.
pub(super) struct ElementIds<EId> {
    by_kind: [Vec<EId>; ElementKind::ALL.len()],
}

impl<EId: Copy> ElementIds<EId> {
    pub(super) fn new() -> Self {
        Self {
            by_kind: std::array::from_fn(|_| Vec::new()),
        }
    }

    pub(super) fn set(&mut self, kind: ElementKind, ids: Vec<EId>) {
        self.by_kind[kind as usize] = ids;
    }

    /// `ElementKindNotInProfile` for a kind of the other profile (`ndm`),
    /// `UnknownElementIndex` for a row past the end of the kind's table.
    pub(super) fn get(&self, kind: ElementKind, row: u32, ndm: u8) -> Result<EId, DecodeError> {
        let table = kind.table();
        if !kind.in_profile(ndm) {
            return Err(DecodeError::ElementKindNotInProfile { table });
        }
        self.by_kind[kind as usize]
            .get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownElementIndex { table, row })
    }
}

fn add_load_patterns<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &mut Domain<NDIM, NDOF, NId, E>,
    table: &LoadPatternTable,
) -> Vec<LoadPatternId>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    (0..table.series.len())
        .map(|i| {
            domain.add_load_pattern_scaled(load_series_of(&table.series[i]), table.scale_factor[i])
        })
        .collect()
}

/// Everything after the elements: load patterns, per-stage nodal and element
/// loads, stage compilation, recorder resolution, model validation, and the
/// session. `make_load(spec)` builds the profile's element load from a wire load, or
/// `None` when the profile has no such load (a 3D model's continuum loads).
pub(super) fn finish<const NDIM: usize, const NDOF: usize, NId, E>(
    mut domain: Domain<NDIM, NDOF, NId, E>,
    input: &CarapaceInputV1,
    node_ids: Vec<NId>,
    elements: ElementIds<E::Id>,
    make_load: impl Fn(&ElementLoadSpec) -> Option<E::Load>,
) -> Result<ModelSession<NDIM, NDOF, NId, E>, DecodeError>
where
    NId: Key + Copy,
    E: ElementOps<NDIM, NDOF, NId> + Clone,
    E::Load: Clone,
{
    let ndm = NDIM as u8;
    let ndof = NDOF as u8;
    let node_at = node_lookup(&node_ids);
    let pattern_ids = add_load_patterns(&mut domain, &input.load_patterns);
    let pattern_at = |row: u32| -> Result<LoadPatternId, DecodeError> {
        pattern_ids
            .get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownPatternIndex {
                table: "load_patterns",
                row,
            })
    };

    // Loads are *not* registered on the domain here: a stage's `LoadControl`
    // shares one pseudo-time across every unfrozen pattern
    // (core/src/analysis/integrator.rs), so a later stage's reference load
    // must not exist in the domain yet while an earlier stage ramps — see
    // `tables::NodalLoadTable::stage`. Each load is grouped by the stage
    // that registers it and only actually added in `ModelSession::
    // start_current_stage`, once that stage's `Domain` is current.
    let stage_count = input.sequence.stages.len();
    let mut nodal_loads_by_stage: Vec<Vec<(LoadPatternId, NId, usize, f64)>> =
        vec![Vec::new(); stage_count];
    let nodal = &input.nodal_loads;
    for i in 0..nodal.node.len() {
        let dof = nodal.dof[i];
        check_dof_within("nodal_loads", i as u32, dof, ndof)?;
        nodal_loads_by_stage
            .get_mut(nodal.stage[i] as usize)
            .ok_or(DecodeError::UnknownStageIndex {
                table: "nodal_loads",
                row: i as u32,
            })?
            .push((
                pattern_at(nodal.pattern[i])?,
                node_at(nodal.node[i], "nodal_loads")?,
                dof as usize,
                nodal.value[i],
            ));
    }

    let mut element_loads_by_stage: Vec<Vec<(LoadPatternId, E::Id, E::Load)>> =
        vec![Vec::new(); stage_count];
    let loads = &input.element_loads;
    for i in 0..loads.element_index.len() {
        let pattern = pattern_at(loads.pattern[i])?;
        let kind = loads.element_kind[i];
        if !kind.in_profile(ndm) {
            return Err(DecodeError::ElementKindNotInProfile {
                table: kind.table(),
            });
        }
        let spec = loads.load[i];
        let carried = match spec {
            ElementLoadSpec::Uniform { .. } => kind.accepts_uniform_load(),
            ElementLoadSpec::Body { .. }
            | ElementLoadSpec::EdgeTraction { .. }
            | ElementLoadSpec::EdgePressure { .. } => kind.is_continuum(),
        };
        if !carried {
            return Err(DecodeError::UnsupportedElementLoad {
                element_kind: kind.table(),
            });
        }
        let edge = match spec {
            ElementLoadSpec::EdgeTraction { edge, .. }
            | ElementLoadSpec::EdgePressure { edge, .. } => Some(edge),
            _ => None,
        };
        if edge.is_some_and(|edge| edge >= kind.edge_count()) {
            return Err(invalid_row(
                "element_loads",
                i,
                "edge is past the element's last edge",
            ));
        }
        let element_id = elements.get(kind, loads.element_index[i], ndm)?;
        if ndm == 2 {
            if let ElementLoadSpec::Uniform { wz: Some(wz), .. } = spec {
                if wz != 0.0 {
                    return Err(invalid_row(
                        "element_loads",
                        i,
                        "wz is only valid in a 3D model",
                    ));
                }
            }
        }
        let load = make_load(&spec).ok_or(DecodeError::UnsupportedElementLoad {
            element_kind: kind.table(),
        })?;
        element_loads_by_stage
            .get_mut(loads.stage[i] as usize)
            .ok_or(DecodeError::UnknownStageIndex {
                table: "element_loads",
                row: i as u32,
            })?
            .push((pattern, element_id, load));
    }

    let stages = compile_stages(
        &input.sequence.stages,
        &node_at,
        &pattern_at,
        ndof,
        ndm,
        nodal_loads_by_stage,
        element_loads_by_stage,
    )?;

    let recorders = input
        .sequence
        .recorders
        .iter()
        .enumerate()
        .map(|(index, recorder)| {
            let recorder_dof = |dof: u8| check_dof_within("recorders", index as u32, dof, ndof);
            let component_within = |component: u8, width: usize| {
                if (component as usize) < width {
                    Ok(())
                } else {
                    Err(DecodeError::InvalidRecorderComponent {
                        recorder: index as u32,
                        component,
                        width: width as u32,
                    })
                }
            };
            Ok(match *recorder {
                RecorderSpec::NodeDisp { node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::NodeDisp {
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::NodeVel { node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::NodeVel {
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::NodeAccel { node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::NodeAccel {
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::ElementForce {
                    element_kind,
                    element_index,
                    component,
                } => {
                    let element = elements.get(element_kind, element_index, ndm)?;
                    component_within(component, domain.element(element).local_force_width())?;
                    ResolvedRecorder::ElementForce { element, component }
                }
                RecorderSpec::ElementLoad {
                    element_kind,
                    element_index,
                    component,
                } => {
                    let element = elements.get(element_kind, element_index, ndm)?;
                    component_within(component, <E::Load as ElementLoadComponents>::COUNT)?;
                    ResolvedRecorder::ElementLoad { element, component }
                }
                RecorderSpec::ModeShape { mode, node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::ModeShape {
                        mode,
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::Reaction { node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::Reaction {
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::Fiber {
                    element_kind,
                    element_index,
                    point,
                    fiber,
                    quantity,
                } => ResolvedRecorder::Fiber {
                    element: elements.get(element_kind, element_index, ndm)?,
                    point,
                    fiber,
                    response: quantity,
                },
                RecorderSpec::GaussPoint {
                    element_kind,
                    element_index,
                    point,
                    quantity,
                    component,
                } => {
                    let element = elements.get(element_kind, element_index, ndm)?;
                    let count = domain.element_gauss_point_count(element);
                    if point as usize >= count {
                        return Err(DecodeError::InvalidGaussPoint {
                            recorder: index as u32,
                            point,
                            count: count as u32,
                        });
                    }
                    component_within(component, 3)?;
                    ResolvedRecorder::GaussPoint {
                        element,
                        point,
                        quantity,
                        component,
                    }
                }
            })
        })
        .collect::<Result<Vec<_>, DecodeError>>()?;

    domain
        .validate()
        .map_err(|error| DecodeError::InvalidModel {
            error: error.into(),
        })?;

    drop(node_at); // the opaque closure type keeps `node_ids` borrowed until dropped
    Ok(ModelSession::new(
        domain,
        node_ids,
        stages,
        recorders,
        input.header.record_initial,
    ))
}
