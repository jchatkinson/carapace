//! Decoding shared by both profiles: everything that is bookkeeping around the
//! per-profile element constructors — nodes (generic over `<NDIM, NDOF>`),
//! fibers, equal-DOF/diaphragm/rigid-link/linear constraints, sparse per-row
//! grouping, load patterns, nodal and element loads, and recorder resolution
//! through an `(ElementKind, row) -> E::Id` map ([`ElementIds`]).

use carapace_core::model::NodeKey;
use std::collections::HashMap;

use carapace_core::model::{
    transfer_floor, BeamIntegration, Domain, ElementLoadComponents, ElementOps, FloorBounds,
    FloorContinuity, FloorLoad, FloorSupport, LoadPatternId, Material, Node, Orientation,
    WallCaseLoad, WallLoad,
};

use super::stages::{compile_stages, load_series_of};
use crate::building::convert_wall_load;
use crate::input_v1::error::DecodeError;
use crate::input_v1::model_output::ModelOutput;
use crate::input_v1::sequence::{BulkQuantity, BulkTable, GaussQuantity, RecorderSpec};
use crate::input_v1::session::{BulkSource, ModelSession, ResolvedBulk, ResolvedRecorder};
use crate::input_v1::tables::{
    Axis3Spec, ElementKind, ElementLoadSpec, EqualDofTable, FiberTable, FrictionRow,
    IntegrationSpec, LinearConstraintTable, LoadPatternTable, NodeTable, OrientRow,
    RigidDiaphragmTable, RigidLinkTable, StudWallSpec,
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

pub(super) fn validate_stud_wall_geometry(
    nodes: &NodeTable,
    ndm: usize,
    wall: &StudWallSpec,
) -> Result<(), &'static str> {
    let point = |index: u32| -> [f64; 3] {
        let mut p = [0.; 3];
        for (axis, coord) in p.iter_mut().enumerate().take(ndm) {
            *coord = nodes.coords[index as usize * ndm + axis];
        }
        p
    };
    let p = wall.corners.map(point);
    let mut dx = [0.; 3];
    let mut dy = [0.; 3];
    for axis in 0..3 {
        dx[axis] = p[1][axis] - p[0][axis];
        dy[axis] = p[3][axis] - p[0][axis];
    }
    let norm = |v: [f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let (lx, ly) = (norm(dx), norm(dy));
    if !lx.is_finite() || !ly.is_finite() || lx <= f64::MIN_POSITIVE || ly <= f64::MIN_POSITIVE {
        return Err("wall axes are degenerate");
    }
    let dot = (0..3).map(|i| dx[i] * dy[i]).sum::<f64>();
    if !dot.is_finite() || dot.abs() > 1e-9 * lx * ly {
        return Err("wall axes are not orthogonal");
    }
    let ex = dx.map(|x| x / lx);
    let ey = dy.map(|x| x / ly);
    let cross = [
        ex[1] * ey[2] - ex[2] * ey[1],
        ex[2] * ey[0] - ex[0] * ey[2],
        ex[0] * ey[1] - ex[1] * ey[0],
    ];
    let length = wall.layout.length;
    // The exterior top row sits a floor-zone depth above the framing.
    let height = wall.layout.height + wall.layout.floor.map_or(0., |z| z.depth);
    let tol = 1e-9 * length.max(height).max(1.);
    let expected = [(0., 0.), (length, 0.), (length, height), (0., height)];
    let local = |q: [f64; 3]| {
        let r = [q[0] - p[0][0], q[1] - p[0][1], q[2] - p[0][2]];
        let x = (0..3).map(|i| r[i] * ex[i]).sum::<f64>();
        let y = (0..3).map(|i| r[i] * ey[i]).sum::<f64>();
        let off = (0..3).map(|i| r[i] * cross[i]).sum::<f64>().abs();
        (x, y, off)
    };
    for (i, q) in p.iter().enumerate() {
        let (x, y, off) = local(*q);
        if off
            .max((x - expected[i].0).abs())
            .max((y - expected[i].1).abs())
            > tol
        {
            return Err("corner geometry does not match the wall layout");
        }
    }
    for station in &wall.stations {
        let (x, y, off) = local(point(station.node));
        let target = match station.chord {
            crate::building::WallChordSpec::Top => height,
            crate::building::WallChordSpec::Bottom | crate::building::WallChordSpec::Base => 0.,
        };
        if (x - station.x).abs().max((y - target).abs()).max(off) > tol {
            return Err("station does not lie at its declared chord location");
        }
    }
    Ok(())
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

/// The wall loads of every pattern, gathered before the walls are built: each (pattern, wall) pair with loads is
/// one load case of that wall.
pub(super) struct WallLoadPlan {
    /// Per wall (row of `studWalls`): its load cases.
    pub cases: Vec<Vec<Vec<WallLoad>>>,
    /// One entry per case: `(pattern row, wall row, case index on that wall, stage)`.
    pub pairs: Vec<(u32, usize, u32, u32)>,
}

/// Groups `wallLoad` element loads and the reactions of `floorLoads` rows by (pattern, wall): all of a pattern's
/// loads on a wall become one case, so two patterns on a wall never merge and one pattern's loads never split. A
/// `wallLoad` row on another element kind, a wall row past the table or a floor row naming a missing pattern is an
/// error, as are `wallLoad` rows of one (pattern, wall) pair in different stages.
pub(super) fn wall_load_plan(
    input: &CarapaceInputV1,
    wall_count: usize,
) -> Result<WallLoadPlan, DecodeError> {
    use crate::input_v1::tables::{FloorTransferLoadSpec, WallTargetSpec};
    type Group = (Vec<WallLoad>, Option<u32>);
    let mut groups: std::collections::BTreeMap<(u32, usize), Group> = Default::default();
    let loads = &input.element_loads;
    for i in 0..loads.element_index.len() {
        let ElementLoadSpec::WallLoad { load } = loads.load[i] else {
            continue;
        };
        let kind = loads.element_kind[i];
        if !matches!(kind, ElementKind::StudWall) {
            return Err(DecodeError::UnsupportedElementLoad {
                element_kind: kind.table(),
            });
        }
        let wall = loads.element_index[i] as usize;
        if wall >= wall_count {
            return Err(DecodeError::UnknownElementIndex {
                table: kind.table(),
                row: wall as u32,
            });
        }
        let entry = groups.entry((loads.pattern[i], wall)).or_default();
        if entry.1.is_some_and(|stage| stage != loads.stage[i]) {
            return Err(invalid_row(
                "element_loads",
                i,
                "wallLoad rows of one pattern on one wall must share a stage",
            ));
        }
        entry.1 = Some(loads.stage[i]);
        entry.0.push(convert_wall_load(load));
    }
    for (fi, floor) in input.floor_loads.iter().enumerate() {
        if floor.pattern as usize >= input.load_patterns.series.len() {
            return Err(DecodeError::UnknownPatternIndex {
                table: "floor_loads",
                row: floor.pattern,
            });
        }
        let bounds = FloorBounds {
            s0: floor.bounds[0],
            s1: floor.bounds[1],
            t0: floor.bounds[2],
            t1: floor.bounds[3],
        };
        let supports = floor
            .supports
            .iter()
            .map(|s| FloorSupport {
                s: s.s,
                t0: s.t0,
                t1: s.t1,
            })
            .collect::<Vec<_>>();
        let floor_loads = floor
            .loads
            .iter()
            .map(|l| match *l {
                FloorTransferLoadSpec::Area { w, region } => FloorLoad::Area { slot: 0, w, region },
                FloorTransferLoadSpec::Line { w, s, t0, t1 } => FloorLoad::Line {
                    slot: 0,
                    w,
                    s,
                    t0,
                    t1,
                },
                FloorTransferLoadSpec::LineAlong { w, t, s0, s1 } => FloorLoad::LineAlong {
                    slot: 0,
                    w,
                    t,
                    s0,
                    s1,
                },
                FloorTransferLoadSpec::Point { p, s, t } => FloorLoad::Point { slot: 0, p, s, t },
            })
            .collect::<Vec<_>>();
        let continuity = match floor.continuity {
            crate::building::FloorContinuitySpec::Simple => FloorContinuity::Simple,
            crate::building::FloorContinuitySpec::Continuous => FloorContinuity::Continuous,
        };
        let reactions =
            transfer_floor(bounds, &supports, &floor_loads, continuity).map_err(|e| {
                DecodeError::InvalidFloorTransfer {
                    row: fi as u32,
                    reason: e.to_string(),
                }
            })?;
        for (si, support) in floor.supports.iter().enumerate() {
            let WallTargetSpec::StudWall { index } = support.wall;
            if index as usize >= wall_count {
                return Err(DecodeError::UnknownElementIndex {
                    table: "stud_walls",
                    row: index,
                });
            }
            let entry = groups.entry((floor.pattern, index as usize)).or_default();
            for loads in reactions.per_support[si]
                .wall_loads(support.x_at_t0, support.reversed)
                .into_values()
            {
                entry.0.extend(loads);
            }
        }
    }
    let mut plan = WallLoadPlan {
        cases: vec![Vec::new(); wall_count],
        pairs: Vec::new(),
    };
    for ((pattern, wall), (loads, stage)) in groups {
        plan.pairs.push((
            pattern,
            wall,
            plan.cases[wall].len() as u32,
            stage.unwrap_or(0),
        ));
        plan.cases[wall].push(loads);
    }
    Ok(plan)
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
    NId: NodeKey + Copy,
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
    NId: NodeKey + Copy,
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
    NId: NodeKey + Copy,
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

    /// Appends rows to `kind` after the ones its table registered.
    pub(super) fn extend(&mut self, kind: ElementKind, ids: impl IntoIterator<Item = EId>) {
        self.by_kind[kind as usize].extend(ids);
    }

    /// The number of rows registered for `kind`.
    pub(super) fn count(&self, kind: ElementKind) -> usize {
        self.by_kind[kind as usize].len()
    }

    /// Reverse map: each element id to its `(kind, row)` in the wire tables.
    pub(super) fn locations(&self) -> HashMap<EId, (ElementKind, u32)>
    where
        EId: std::hash::Hash + Eq,
    {
        let mut map = HashMap::new();
        for kind in ElementKind::ALL {
            for (row, &id) in self.by_kind[kind as usize].iter().enumerate() {
                map.insert(id, (kind, row as u32));
            }
        }
        map
    }

    /// Every element id of `kind`, in table order (`ElementKindNotInProfile` for a kind of the
    /// other profile).
    pub(super) fn all(&self, kind: ElementKind, ndm: u8) -> Result<&[EId], DecodeError> {
        if !kind.in_profile(ndm) {
            return Err(DecodeError::ElementKindNotInProfile {
                table: kind.table(),
            });
        }
        Ok(&self.by_kind[kind as usize])
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
    NId: NodeKey + Copy,
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
/// Per stage, the `(pattern, element, load)` triples that stage adds.
type StagedLoads<I, L> = Vec<Vec<(LoadPatternId, I, L)>>;

#[allow(clippy::too_many_arguments)]
/// The rows `rows` names of `all`, in that order (every row without it); `None` when one is past the end.
fn pick_rows<T: Copy>(all: &[T], rows: &Option<Vec<u32>>) -> Option<Vec<T>> {
    match rows {
        None => Some(all.to_vec()),
        Some(rows) => rows
            .iter()
            .map(|&row| all.get(row as usize).copied())
            .collect(),
    }
}

pub(super) fn finish<const NDIM: usize, const NDOF: usize, NId, E>(
    mut domain: Domain<NDIM, NDOF, NId, E>,
    input: &CarapaceInputV1,
    node_ids: Vec<NId>,
    elements: ElementIds<E::Id>,
    model_output: &ModelOutput,
    extra_nodal_loads: Vec<(u32, u32, NId, usize, f64)>,
    wall_pairs: &[(u32, usize, u32, u32)],
    make_load: impl Fn(&ElementLoadSpec) -> Option<E::Load>,
    make_wall_load: impl Fn(WallCaseLoad) -> E::Load,
) -> Result<ModelSession<NDIM, NDOF, NId, E>, DecodeError>
where
    NId: NodeKey + Copy,
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

    // Loads on nodes the decoder created itself (the rod end of a bar): (pattern row, stage, node, dof, value).
    for (pattern, stage, node, dof, value) in extra_nodal_loads {
        nodal_loads_by_stage
            .get_mut(stage as usize)
            .ok_or(DecodeError::UnknownStageIndex {
                table: "embedded_bars",
                row: 0,
            })?
            .push((pattern_at(pattern)?, node, dof, value));
    }

    let mut element_loads_by_stage: StagedLoads<E::Id, E::Load> = vec![Vec::new(); stage_count];
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
            ElementLoadSpec::Trapezoid { .. } => kind.accepts_trapezoid_load(),
            ElementLoadSpec::Body { .. }
            | ElementLoadSpec::EdgeTraction { .. }
            | ElementLoadSpec::EdgePressure { .. } => kind.is_continuum(),
            ElementLoadSpec::ShellPressure { .. } | ElementLoadSpec::ShellBody { .. } => {
                kind.is_shell()
            }
            ElementLoadSpec::WallLoad { .. } => matches!(kind, ElementKind::StudWall),
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
        if matches!(spec, ElementLoadSpec::WallLoad { .. }) {
            continue; // grouped into one load case per (pattern, wall) by `wall_load_plan`, added below
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
    for &(pattern, wall, case, stage) in wall_pairs {
        let element_id = elements.get(ElementKind::StudWall, wall as u32, ndm)?;
        let load = make_wall_load(WallCaseLoad { case, factor: 1.0 });
        element_loads_by_stage
            .get_mut(stage as usize)
            .ok_or(DecodeError::UnknownStageIndex {
                table: "element_loads",
                row: 0,
            })?
            .push((pattern_at(pattern)?, element_id, load));
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

    let mut bulks = Vec::new();
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
            Ok(match recorder.clone() {
                RecorderSpec::StudWallBase {
                    element_kind,
                    element_index,
                    component,
                } => {
                    if !matches!(element_kind, ElementKind::StudWall) {
                        return Err(DecodeError::RecorderNotApplicable {
                            recorder: index as u32,
                            reason: "studWallBase needs a stud-wall element",
                        });
                    }
                    ResolvedRecorder::StudWallBase {
                        element: elements.get(element_kind, element_index, ndm)?,
                        component,
                    }
                }
                RecorderSpec::DiaphragmBay {
                    element_index,
                    bay,
                    quantity,
                } => {
                    let element = elements.get(ElementKind::Diaphragm, element_index, ndm)?;
                    let count = domain
                        .element_diaphragm_response(element)
                        .map_or(0, |r| r.bays.len());
                    if bay as usize >= count {
                        return Err(DecodeError::InvalidRecorderBay {
                            recorder: index as u32,
                            bay,
                            count: count as u32,
                        });
                    }
                    ResolvedRecorder::DiaphragmBay {
                        element,
                        bay,
                        quantity,
                    }
                }
                RecorderSpec::StudWallPier {
                    element_kind,
                    element_index,
                    pier,
                    quantity,
                } => {
                    if !matches!(element_kind, ElementKind::StudWall) {
                        return Err(DecodeError::RecorderNotApplicable {
                            recorder: index as u32,
                            reason: "studWallPier needs a stud-wall element",
                        });
                    }
                    let element = elements.get(element_kind, element_index, ndm)?;
                    let count = domain
                        .element_stud_wall_response(element, 0.0)
                        .map_or(0, |r| r.piers.len());
                    if pier as usize >= count {
                        return Err(DecodeError::InvalidRecorderPier {
                            recorder: index as u32,
                            pier,
                            count: count as u32,
                        });
                    }
                    ResolvedRecorder::StudWallPier {
                        element,
                        pier,
                        quantity,
                    }
                }
                RecorderSpec::NodeDisp { node, dof } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::NodeDisp {
                        node: node_at(node, "recorders")?,
                        dof,
                    }
                }
                RecorderSpec::NodeDifference {
                    node,
                    other_node,
                    dof,
                } => {
                    recorder_dof(dof)?;
                    ResolvedRecorder::NodeDifference {
                        node: node_at(node, "recorders")?,
                        other_node: other_node.map(|n| node_at(n, "recorders")).transpose()?,
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
                    if quantity == GaussQuantity::Kc2 {
                        component_within(component, 1)?;
                        if domain.element_gauss_kc2(element, point as usize).is_none() {
                            return Err(DecodeError::RecorderNotApplicable {
                                recorder: index as u32,
                                reason: "kc2 needs a csfmConcrete element",
                            });
                        }
                    } else {
                        component_within(component, domain.element_gauss_component_count(element))?;
                    }
                    ResolvedRecorder::GaussPoint {
                        element,
                        point,
                        quantity,
                        component,
                    }
                }
                RecorderSpec::Bulk {
                    table,
                    quantity,
                    element_kind,
                    ref rows,
                } => {
                    let not_applicable = |reason| DecodeError::RecorderNotApplicable {
                        recorder: index as u32,
                        reason,
                    };
                    if (table == BulkTable::Elements) != element_kind.is_some() {
                        return Err(not_applicable(
                            "elementKind names the table of table: elements, and only of it",
                        ));
                    }
                    if quantity == BulkQuantity::MemberForce {
                        if table != BulkTable::Elements || element_kind != Some(ElementKind::StudWall) {
                            return Err(not_applicable(
                                "memberForce belongs to the elements table of stud walls, and only to it",
                            ));
                        }
                        let ids = pick_rows(elements.all(ElementKind::StudWall, ndm)?, rows)
                            .ok_or(not_applicable("a row is past the end of the table"))?;
                        let width = |id| {
                            domain.element_stud_wall_response(id, 0.).map(|r| {
                                3 + 3 * (r.studs.len() + r.jacks.len() + r.cripples.len())
                            })
                        };
                        let Some(&first) = ids.first() else {
                            return Err(not_applicable("the table is empty"));
                        };
                        let components = width(first)
                            .ok_or(not_applicable("memberForce needs a stud-wall element"))?;
                        if ids.iter().any(|&id| width(id) != Some(components)) {
                            return Err(not_applicable(
                                "memberForce rows must have as many members (record one wall per recorder)",
                            ));
                        }
                        bulks.push(ResolvedBulk {
                            recorder_index: index,
                            elements: ids,
                            nodes: Vec::new(),
                            points: 1,
                            components: components as u32,
                            quantity,
                            source: BulkSource::WallMembers,
                            emitted: 0,
                        });
                        return Ok(ResolvedRecorder::Bulk {
                            bulk: bulks.len() as u32 - 1,
                        });
                    }
                    let node_quantity = matches!(
                        quantity,
                        BulkQuantity::Displacement
                            | BulkQuantity::Reaction
                            | BulkQuantity::Velocity
                            | BulkQuantity::Acceleration
                    );
                    if (table == BulkTable::Nodes) != node_quantity {
                        return Err(not_applicable(
                            "displacement, reaction, velocity and acceleration belong to the nodes table, and only they",
                        ));
                    }
                    if table == BulkTable::Nodes {
                        let nodes = pick_rows(&node_ids, rows)
                            .ok_or(not_applicable("a row is past the end of the table"))?;
                        if nodes.is_empty() {
                            return Err(not_applicable("the table is empty"));
                        }
                        bulks.push(ResolvedBulk {
                            recorder_index: index,
                            elements: Vec::new(),
                            nodes,
                            points: 1,
                            components: ndof as u32,
                            quantity,
                            source: BulkSource::Nodes,
                            emitted: 0,
                        });
                        return Ok(ResolvedRecorder::Bulk {
                            bulk: bulks.len() as u32 - 1,
                        });
                    }
                    let zero_lengths = table == BulkTable::ZeroLengths;
                    if zero_lengths
                        != matches!(quantity, BulkQuantity::Deformation | BulkQuantity::Force)
                    {
                        return Err(not_applicable(
                            "deformation and force belong to the zeroLengths table, and only they",
                        ));
                    }
                    let fits_elements = matches!(
                        quantity,
                        BulkQuantity::EndForce | BulkQuantity::Load | BulkQuantity::Strain | BulkQuantity::Stress
                    );
                    if matches!(quantity, BulkQuantity::EndForce | BulkQuantity::Load) && table != BulkTable::Elements
                        || table == BulkTable::Elements && !fits_elements
                    {
                        return Err(not_applicable(
                            "the elements table records endForce, load, strain or stress; endForce and load only it",
                        ));
                    }
                    let kind = match table {
                        BulkTable::Elements => element_kind.expect("checked above"),
                        BulkTable::ZeroLengths => ElementKind::ZeroLength,
                        BulkTable::Triangles => ElementKind::Tri3,
                        BulkTable::Quads => ElementKind::Quad4,
                        _ => ElementKind::EmbeddedBar,
                    };
                    let bars = table == BulkTable::EmbeddedBars;
                    match quantity {
                        BulkQuantity::AxialForce | BulkQuantity::Ruptured if !bars => {
                            return Err(not_applicable(
                                "axialForce and ruptured belong to the embeddedBars table",
                            ))
                        }
                        BulkQuantity::PrincipalStrain
                        | BulkQuantity::PrincipalStress
                        | BulkQuantity::Kc2
                            if bars =>
                        {
                            return Err(not_applicable(
                                "principal values and kc2 belong to the triangles and quads tables",
                            ))
                        }
                        _ => {}
                    }
                    let ids = pick_rows(elements.all(kind, ndm)?, rows)
                        .ok_or(not_applicable("a row is past the end of the table"))?;
                    let Some(&first) = ids.first() else {
                        return Err(not_applicable("the table is empty"));
                    };
                    let shell = domain.element_shell_responses(first);
                    let gauss = domain.element_gauss_point_count(first) as u32;
                    let (source, points, components) = match quantity {
                        BulkQuantity::EndForce => {
                            let width = domain.element(first).local_force_width();
                            if width == 0 || ids.iter().any(|&id| domain.element(id).local_force_width() != width) {
                                return Err(not_applicable("the table's elements report no end forces"));
                            }
                            (BulkSource::EndForce, 1, width as u32)
                        }
                        BulkQuantity::Load => (BulkSource::Load, 1, ndm as u32),
                        _ if table == BulkTable::Elements => match &shell {
                            Some(points) => (BulkSource::Shell, points.len() as u32, domain.element_gauss_component_count(first) as u32),
                            None if gauss > 0 && !matches!(quantity, BulkQuantity::Kc2) => (BulkSource::Gauss, gauss, 3),
                            None => return Err(not_applicable("strain and stress need a shell or continuum table")),
                        },
                        _ if zero_lengths => (BulkSource::ZeroLengths, 1, 1),
                        _ if bars => (BulkSource::Bars, gauss, 1),
                        BulkQuantity::Kc2 => (BulkSource::Gauss, gauss, 1),
                        _ => (BulkSource::Gauss, gauss, 3),
                    };
                    if quantity == BulkQuantity::Kc2 {
                        let csfm = ids.iter().all(|&id| {
                            (0..points as usize).all(|p| domain.element_gauss_kc2(id, p).is_some())
                        });
                        if !csfm {
                            return Err(not_applicable("kc2 needs a csfmConcrete element"));
                        }
                    }
                    bulks.push(ResolvedBulk {
                        recorder_index: index,
                        elements: ids,
                        nodes: Vec::new(),
                        points,
                        components,
                        quantity,
                        source,
                        emitted: 0,
                    });
                    ResolvedRecorder::Bulk {
                        bulk: bulks.len() as u32 - 1,
                    }
                }
                RecorderSpec::EmbeddedBar {
                    bar,
                    segment,
                    quantity,
                } => {
                    if ndm == 3 {
                        return Err(DecodeError::ElementKindNotInProfile {
                            table: ElementKind::EmbeddedBar.table(),
                        });
                    }
                    let group = model_output.embedded_bars.get(bar as usize).ok_or(
                        DecodeError::UnknownElementIndex {
                            table: ElementKind::EmbeddedBar.table(),
                            row: bar,
                        },
                    )?;
                    if segment as usize >= group.segments.len() {
                        return Err(DecodeError::InvalidEmbeddedBarSegment {
                            recorder: index as u32,
                            bar,
                            segment,
                            count: group.segments.len() as u32,
                        });
                    }
                    ResolvedRecorder::EmbeddedBar {
                        element: elements.get(
                            ElementKind::EmbeddedBar,
                            group.first_segment + segment,
                            ndm,
                        )?,
                        quantity,
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
        bulks,
        elements.locations(),
        input.header.record_initial,
    ))
}

/// A stud-wall layout's bearing spec to the core type (the hold-down as an `Anchorage` material).
pub(super) fn decode_wall_bearing(
    l: &crate::building::WallLayoutSpec,
) -> Result<Option<carapace_core::model::WallBearing>, String> {
    use carapace_core::model::{AnchorageBearing, Material, WallBearing};
    let Some(b) = &l.bearing else {
        return Ok(None);
    };
    let hold_down = b
        .hold_down
        .map(|a| {
            let bearing = match a.bearing {
                crate::building::AnchorageBearingSpec::Linear { k } => {
                    AnchorageBearing::Linear { k }
                }
                crate::building::AnchorageBearingSpec::O86Crush { qr } => {
                    AnchorageBearing::O86Crush { qr }
                }
            };
            Material::anchorage(a.tension, a.slack, bearing)
                .map_err(|e| format!("bearing.holdDown: {e}"))
        })
        .transpose()?;
    Ok(Some(WallBearing {
        kv: b.kv,
        kh: b.kh,
        hold_down,
    }))
}

/// A station chord spec to the core type.
pub(super) fn decode_wall_chord(
    c: crate::building::WallChordSpec,
) -> carapace_core::model::WallChord {
    use crate::building::WallChordSpec as S;
    use carapace_core::model::WallChord;
    match c {
        S::Top => WallChord::Top,
        S::Bottom => WallChord::Bottom,
        S::Base => WallChord::Base,
    }
}
