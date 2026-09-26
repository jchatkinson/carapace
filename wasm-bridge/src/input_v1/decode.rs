//! Hand-written per-table decode: translates a [`CarapaceInputV1`] into a
//! [`Session`] by calling the same `Domain`/`Element`/`Material`
//! constructors native code calls — no generic deserialization into `core`
//! types, no reflection (pysees-handoff.md's "Decode and session model").
//! Panic-free: every lookup returns a [`DecodeError`] instead of indexing
//! past the end of a table, since this is defense in depth against
//! engine/schema version skew, not a re-check of what `pysees`'s compiler
//! is expected to have already validated.

use std::collections::HashMap;

use carapace_core::analysis::{Algorithm, ConvergenceTest, GroundMotion, Integrator, RayleighDamping};
use carapace_core::model::{
    BeamIntegration, DispBeamColumn, Domain, ElasticBeamColumn, Element, ElementId, ElementLoad,
    Fiber, FiberSection, ForceBeamColumn, Friction, GeomTransf, LoadPatternId, LoadSeries,
    Material, Node, NodeId, Truss, ZeroLength, ZeroLengthSection,
};

use super::error::DecodeError;
use super::materials::resolve_materials;
use super::sequence::{AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, StageSpec};
use super::session::{CompiledStage, CompiledStageKind, PlanarSession, ResolvedRecorder, Session};
use super::tables::{ElementKind, ElementLoadSpec, IntegrationSpec, TimeSeriesSpec, TransformSpec};
use super::CarapaceInputV1;

const PLANAR_SPACE: u8 = 2;
const SPATIAL_SPACE: u8 = 3;

/// Dispatches on `header.space`, chosen once here and never branched on
/// again downstream (`Session`'s own doc comment) — `decode_planar` and
/// `decode3::decode_spatial` are otherwise two independent per-table
/// decoders (only `compile_stages`/`check_dof_within` are shared, since
/// stage/integrator compilation doesn't touch element physics).
pub fn decode(input: CarapaceInputV1) -> Result<Session, DecodeError> {
    match input.header.space {
        PLANAR_SPACE => decode_planar(input).map(Session::Planar),
        SPATIAL_SPACE => super::decode3::decode_spatial(input).map(Session::Spatial),
        got => Err(DecodeError::UnsupportedSpace { got }),
    }
}

fn decode_planar(input: CarapaceInputV1) -> Result<PlanarSession, DecodeError> {
    let materials = resolve_materials(&input.materials)?;
    let material_at = |row: u32, table: &'static str| -> Result<Material, DecodeError> {
        materials
            .get(row as usize)
            .cloned()
            .ok_or(DecodeError::UnknownMaterialIndex { table, row })
    };

    let mut domain = Domain::new();
    let node_ids = add_nodes(&mut domain, &input.nodes);
    let node_at = |row: u32, table: &'static str| -> Result<NodeId, DecodeError> {
        node_ids
            .get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownNodeIndex { table, row })
    };
    add_equal_dofs(&mut domain, &input.equal_dofs, &node_at)?;
    add_rigid_diaphragms(&mut domain, &input.rigid_diaphragms, &node_at)?;
    let fiber_section_at = |section: u32| -> Result<Vec<Fiber>, DecodeError> {
        build_fiber_section(&input.fibers, section, &material_at)
    };

    let mut truss_ids = Vec::with_capacity(input.trusses.node_i.len());
    for i in 0..input.trusses.node_i.len() {
        let element = Truss::new(
            node_at(input.trusses.node_i[i], "trusses")?,
            node_at(input.trusses.node_j[i], "trusses")?,
            input.trusses.area[i],
            material_at(input.trusses.material[i], "trusses")?,
        )
        .with_density(input.trusses.density[i]);
        truss_ids.push(domain.add_element(Element::Truss(element)));
    }

    let mut elastic_beam_ids = Vec::with_capacity(input.elastic_beam_columns.node_i.len());
    for i in 0..input.elastic_beam_columns.node_i.len() {
        let table = &input.elastic_beam_columns;
        let element = ElasticBeamColumn::new(
            node_at(table.node_i[i], "elastic_beam_columns")?,
            node_at(table.node_j[i], "elastic_beam_columns")?,
            table.e[i],
            table.a[i],
            table.iz[i],
            transform_of(table.transform[i]),
        )
        .with_density(table.density[i]);
        elastic_beam_ids.push(domain.add_element(Element::ElasticBeamColumn(element)));
    }

    let mut disp_beam_ids = Vec::with_capacity(input.disp_beam_columns.node_i.len());
    for i in 0..input.disp_beam_columns.node_i.len() {
        let table = &input.disp_beam_columns;
        let mut element = DispBeamColumn::new(
            node_at(table.node_i[i], "disp_beam_columns")?,
            node_at(table.node_j[i], "disp_beam_columns")?,
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        if table.corotational[i] {
            element = element.with_corotational();
        }
        disp_beam_ids.push(domain.add_element(Element::DispBeamColumn(element)));
    }

    let mut force_beam_ids = Vec::with_capacity(input.force_beam_columns.node_i.len());
    for i in 0..input.force_beam_columns.node_i.len() {
        let table = &input.force_beam_columns;
        let mut element = ForceBeamColumn::new(
            node_at(table.node_i[i], "force_beam_columns")?,
            node_at(table.node_j[i], "force_beam_columns")?,
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        if table.corotational[i] {
            element = element.with_corotational();
        }
        force_beam_ids.push(domain.add_element(Element::ForceBeamColumn(element)));
    }

    let zero_length_ids = add_zero_lengths(&mut domain, &input.zero_lengths, &node_at, &material_at)?;
    let zero_length_section_ids = add_zero_length_sections(
        &mut domain,
        &input.zero_length_sections,
        &node_at,
        &material_at,
        &fiber_section_at,
    )?;

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
    // that registers it and only actually added in `PlanarSession::
    // start_current_stage`, once that stage's `Domain` is current.
    let stage_count = input.sequence.stages.len();
    let mut nodal_loads_by_stage: Vec<Vec<(LoadPatternId, NodeId, usize, f64)>> =
        vec![Vec::new(); stage_count];
    for i in 0..input.nodal_loads.node.len() {
        let dof = input.nodal_loads.dof[i];
        check_dof("nodal_loads", i as u32, dof)?;
        let stage = input.nodal_loads.stage[i];
        let entry =
            nodal_loads_by_stage
                .get_mut(stage as usize)
                .ok_or(DecodeError::UnknownStageIndex {
                    table: "nodal_loads",
                    row: i as u32,
                })?;
        entry.push((
            pattern_at(input.nodal_loads.pattern[i])?,
            node_at(input.nodal_loads.node[i], "nodal_loads")?,
            dof as usize,
            input.nodal_loads.value[i],
        ));
    }

    let mut element_loads_by_stage: Vec<Vec<(LoadPatternId, ElementId, ElementLoad)>> =
        vec![Vec::new(); stage_count];
    for i in 0..input.element_loads.element_index.len() {
        let pattern = pattern_at(input.element_loads.pattern[i])?;
        let kind = input.element_loads.element_kind[i];
        if kind != ElementKind::ElasticBeamColumn {
            return Err(DecodeError::UnsupportedElementLoad {
                element_kind: element_kind_name(kind),
            });
        }
        let row = input.element_loads.element_index[i];
        let element_id =
            *elastic_beam_ids
                .get(row as usize)
                .ok_or(DecodeError::UnknownElementIndex {
                    table: "elastic_beam_columns",
                    row,
                })?;
        let load = match input.element_loads.load[i] {
            ElementLoadSpec::UniformTransverse { w } => ElementLoad::UniformTransverse(w),
        };
        let stage = input.element_loads.stage[i];
        element_loads_by_stage
            .get_mut(stage as usize)
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
        3,
        2,
        nodal_loads_by_stage,
        element_loads_by_stage,
    )?;
    let element_at = |kind: ElementKind, row: u32| -> Result<ElementId, DecodeError> {
        let (ids, table) = match kind {
            ElementKind::Truss => (&truss_ids, "trusses"),
            ElementKind::ElasticBeamColumn => (&elastic_beam_ids, "elastic_beam_columns"),
            ElementKind::DispBeamColumn => (&disp_beam_ids, "disp_beam_columns"),
            ElementKind::ForceBeamColumn => (&force_beam_ids, "force_beam_columns"),
            ElementKind::ZeroLength => (&zero_length_ids, "zero_lengths"),
            ElementKind::ZeroLengthSection => (&zero_length_section_ids, "zero_length_sections"),
        };
        ids.get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownElementIndex { table, row })
    };
    let recorders = input
        .sequence
        .recorders
        .iter()
        .map(|recorder| {
            Ok(match *recorder {
                RecorderSpec::NodeDisp { node, dof } => ResolvedRecorder::NodeDisp {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec::NodeVel { node, dof } => ResolvedRecorder::NodeVel {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec::NodeAccel { node, dof } => ResolvedRecorder::NodeAccel {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec::ElementForce {
                    element_kind,
                    element_index,
                    component,
                } => ResolvedRecorder::ElementForce {
                    element: element_at(element_kind, element_index)?,
                    component,
                },
                RecorderSpec::ModeShape { mode, node, dof } => ResolvedRecorder::ModeShape {
                    mode,
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec::Reaction { node, dof } => ResolvedRecorder::Reaction {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec::Fiber {
                    element_kind,
                    element_index,
                    point,
                    fiber,
                    quantity,
                } => ResolvedRecorder::Fiber {
                    element: element_at(element_kind, element_index)?,
                    point,
                    fiber,
                    response: quantity,
                },
            })
        })
        .collect::<Result<Vec<_>, DecodeError>>()?;

    Ok(PlanarSession::new(domain, stages, recorders))
}

fn check_dof(table: &'static str, row: u32, dof: u8) -> Result<(), DecodeError> {
    check_dof_within(table, row, dof, 3)
}

/// `check_dof`, generalized to the caller's DOF-per-node count (`3` planar,
/// `6` spatial — see `decode3.rs`'s own `check_dof3`).
pub(super) fn check_dof_within(table: &'static str, row: u32, dof: u8, ndof: u8) -> Result<(), DecodeError> {
    if dof < ndof {
        Ok(())
    } else {
        Err(DecodeError::InvalidDof { table, row, dof })
    }
}

fn add_nodes(domain: &mut Domain, table: &super::tables::NodeTable) -> Vec<NodeId> {
    let mut mass_by_node: HashMap<u32, [f64; 3]> = HashMap::new();
    for (k, &node_index) in table.mass_node_index.iter().enumerate() {
        let mut mass = [0.0; 3];
        mass.copy_from_slice(&table.mass[k * 3..k * 3 + 3]);
        mass_by_node.insert(node_index, mass);
    }

    let node_count = table.node_count();
    let mut node_ids = Vec::with_capacity(node_count);
    for i in 0..node_count {
        let mut node = Node::new([table.coords[2 * i], table.coords[2 * i + 1]]);
        let fixed_bits = table.fixed[i];
        for dof in 0..3usize {
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
    node_ids
}

/// Decodes [`super::tables::EqualDofTable`] against an already-populated
/// `domain` (nodes only — constraints don't reference elements or
/// patterns), calling `core::Domain::equal_dof` once per row. Shared shape
/// and error handling with `decode3.rs`'s own `add_equal_dofs` — both take
/// a generic `node_at` closure, so this could be lifted to a shared
/// generic helper the way `compile_stages` is; kept separate for now since
/// each caller's `Domain`/`Domain3` type still differs.
fn add_equal_dofs(
    domain: &mut Domain,
    table: &super::tables::EqualDofTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
) -> Result<(), DecodeError> {
    let mut dofs_by_row: Vec<Vec<usize>> = vec![Vec::new(); table.retained.len()];
    for &(row, dof) in &table.dofs {
        check_dof("equal_dofs", row, dof)?;
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

/// Decodes [`super::tables::RigidDiaphragmTable`] against an
/// already-populated `domain` — see `add_equal_dofs`'s doc comment.
fn add_rigid_diaphragms(
    domain: &mut Domain,
    table: &super::tables::RigidDiaphragmTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
) -> Result<(), DecodeError> {
    let mut constrained_by_row: Vec<Vec<u32>> = vec![Vec::new(); table.retained.len()];
    for &(row, node) in &table.constrained {
        constrained_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownConstraintRow {
                table: "rigid_diaphragms",
                row,
            })?
            .push(node);
    }
    #[allow(clippy::needless_range_loop)] // parallel-indexes retained/constrained_by_row
    for i in 0..table.retained.len() {
        let retained = node_at(table.retained[i], "rigid_diaphragms")?;
        let constrained = constrained_by_row[i]
            .iter()
            .map(|&row| node_at(row, "rigid_diaphragms"))
            .collect::<Result<Vec<_>, _>>()?;
        domain.rigid_diaphragm(retained, &constrained);
    }
    Ok(())
}

fn build_fiber_section(
    table: &super::tables::FiberTable,
    section: u32,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<Fiber>, DecodeError> {
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
        .map(|i| {
            let material = material_at(table.material[i], "fibers")?;
            Ok(Fiber::new(table.y[i], table.area[i], material))
        })
        .collect()
}

/// `(normal_dof, shear_dof, mu, k0, b)` — see `ZeroLengthTable::friction`.
type FrictionRow = (u8, u8, f64, f64, f64);

fn add_zero_lengths(
    domain: &mut Domain,
    table: &super::tables::ZeroLengthTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<ElementId>, DecodeError> {
    let mut materials_by_row: Vec<Vec<(u8, u32)>> = vec![Vec::new(); table.node_i.len()];
    for &(row, dof, material_index) in &table.materials {
        materials_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex {
                table: "zero_lengths",
                row,
            })?
            .push((dof, material_index));
    }

    let mut friction_by_row: Vec<Option<FrictionRow>> = vec![None; table.node_i.len()];
    for &(row, normal_dof, shear_dof, mu, k0, b) in &table.friction {
        *friction_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex {
                table: "zero_lengths",
                row,
            })? = Some((normal_dof, shear_dof, mu, k0, b));
    }

    let mut ids = Vec::with_capacity(table.node_i.len());
    #[allow(clippy::needless_range_loop)] // parallel-indexes node_i/node_j/materials_by_row/friction_by_row
    for i in 0..table.node_i.len() {
        let mut element = ZeroLength::new(
            node_at(table.node_i[i], "zero_lengths")?,
            node_at(table.node_j[i], "zero_lengths")?,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element =
                element.with_material(dof as usize, material_at(material_index, "zero_lengths")?);
        }
        if let Some((normal_dof, shear_dof, mu, k0, b)) = friction_by_row[i] {
            element = element.with_friction(Friction::new(
                normal_dof as usize,
                shear_dof as usize,
                mu,
                k0,
                b,
            ));
        }
        ids.push(domain.add_element(Element::ZeroLength(element)));
    }
    Ok(ids)
}

fn add_zero_length_sections(
    domain: &mut Domain,
    table: &super::tables::ZeroLengthSectionTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
    fiber_section_at: &impl Fn(u32) -> Result<Vec<Fiber>, DecodeError>,
) -> Result<Vec<ElementId>, DecodeError> {
    let mut materials_by_row: Vec<Vec<(u8, u32)>> = vec![Vec::new(); table.node_i.len()];
    for &(row, dof, material_index) in &table.materials {
        materials_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex {
                table: "zero_length_sections",
                row,
            })?
            .push((dof, material_index));
    }

    let mut ids = Vec::with_capacity(table.node_i.len());
    #[allow(clippy::needless_range_loop)] // parallel-indexes node_i/node_j/fiber_section/materials_by_row
    for i in 0..table.node_i.len() {
        let section = FiberSection::new(fiber_section_at(table.fiber_section[i])?);
        let mut element = ZeroLengthSection::new(
            node_at(table.node_i[i], "zero_length_sections")?,
            node_at(table.node_j[i], "zero_length_sections")?,
            section,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element
                .with_material(dof as usize, material_at(material_index, "zero_length_sections")?);
        }
        ids.push(domain.add_element(Element::ZeroLengthSection(element)));
    }
    Ok(ids)
}

fn add_load_patterns(
    domain: &mut Domain,
    table: &super::tables::LoadPatternTable,
) -> Vec<LoadPatternId> {
    (0..table.series.len())
        .map(|i| domain.add_load_pattern_scaled(load_series_of(&table.series[i]), table.scale_factor[i]))
        .collect()
}

/// Shared by both profiles' decoders and by ground-motion decoding
/// (`compile_stages`): an accelerogram is exactly a `TimeSeriesSpec::Path`,
/// the same wire type a static `LoadPatternTable` entry uses (`core::
/// GroundMotion`'s own doc comment), so this one conversion serves both.
pub(super) fn load_series_of(spec: &TimeSeriesSpec) -> LoadSeries {
    match spec {
        TimeSeriesSpec::Constant => LoadSeries::Constant,
        TimeSeriesSpec::Linear { slope } => LoadSeries::Linear { slope: *slope },
        TimeSeriesSpec::Path { times, factors } => LoadSeries::Path {
            times: times.clone(),
            factors: factors.clone(),
        },
    }
}

fn transform_of(spec: TransformSpec) -> GeomTransf {
    match spec {
        TransformSpec::Linear => GeomTransf::Linear,
        TransformSpec::PDelta => GeomTransf::PDelta,
        TransformSpec::Corotational => GeomTransf::Corotational,
    }
}

fn integration_of(spec: IntegrationSpec) -> BeamIntegration {
    match spec {
        IntegrationSpec::Legendre { points } => BeamIntegration::Legendre {
            points: points as usize,
        },
        IntegrationSpec::Lobatto { points } => BeamIntegration::Lobatto {
            points: points as usize,
        },
    }
}

fn element_kind_name(kind: ElementKind) -> &'static str {
    match kind {
        ElementKind::Truss => "truss",
        ElementKind::ElasticBeamColumn => "elastic_beam_column",
        ElementKind::DispBeamColumn => "disp_beam_column",
        ElementKind::ForceBeamColumn => "force_beam_column",
        ElementKind::ZeroLength => "zero_length",
        ElementKind::ZeroLengthSection => "zero_length_section",
    }
}

/// Shared by both profiles' decoders (`decode3.rs` calls this too): the
/// stage-compilation logic (integrator/algorithm/convergence conversion,
/// held-pattern resolution) never touches element physics, only
/// `LoadPatternId`s and a caller-supplied `node_at`/`ndof`, so it's generic
/// over the node/element ID types and element-load type rather than
/// duplicated per profile.
/// Shared by both profiles' decoders (`decode3.rs` calls this too): stage
/// compilation for all three `StageSpec` kinds. `ndof` bounds a `Static`
/// stage's `DisplacementControl` dof (3 planar, 6 spatial); `ndim` bounds a
/// `Transient` stage's ground-motion `direction` (2 planar, 3 spatial —
/// translational axes only, never a rotation).
pub(super) fn compile_stages<NId: Copy, EId: Copy, Load: Clone>(
    stages: &[StageSpec],
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
    pattern_at: &impl Fn(u32) -> Result<LoadPatternId, DecodeError>,
    ndof: u8,
    ndim: u8,
    nodal_loads_by_stage: Vec<Vec<(LoadPatternId, NId, usize, f64)>>,
    element_loads_by_stage: Vec<Vec<(LoadPatternId, EId, Load)>>,
) -> Result<Vec<CompiledStage<NId, EId, Load>>, DecodeError> {
    stages
        .iter()
        .zip(nodal_loads_by_stage)
        .zip(element_loads_by_stage)
        .map(|((stage, pending_nodal_loads), pending_element_loads)| {
            let (steps, kind) = match stage {
                StageSpec::Static { steps, integrator, algorithm, convergence, hold_patterns_after, .. } => {
                    let integrator = match *integrator {
                        IntegratorSpec::LoadControl { increment } => Integrator::LoadControl { increment },
                        IntegratorSpec::DisplacementControl { node, dof, increment } => {
                            check_dof_within("sequence.integrator", node, dof, ndof)?;
                            Integrator::DisplacementControl { node: node_at(node, "sequence.integrator")?, dof: dof as usize, increment }
                        }
                    };
                    let algorithm = match algorithm {
                        AlgorithmSpec::Linear => Algorithm::Linear,
                        AlgorithmSpec::NewtonRaphson => Algorithm::NewtonRaphson,
                    };
                    let convergence = match convergence.unwrap_or(ConvergenceSpec::DEFAULT) {
                        ConvergenceSpec::NormUnbalance { tol, max_iter } => {
                            ConvergenceTest::NormUnbalance { tol, max_iter: max_iter as usize }
                        }
                        ConvergenceSpec::NormDispIncr { tol, max_iter } => {
                            ConvergenceTest::NormDispIncr { tol, max_iter: max_iter as usize }
                        }
                        ConvergenceSpec::EnergyIncr { tol, max_iter } => {
                            ConvergenceTest::EnergyIncr { tol, max_iter: max_iter as usize }
                        }
                    };
                    let hold_patterns_after = hold_patterns_after.iter().map(|&row| pattern_at(row)).collect::<Result<Vec<_>, _>>()?;
                    (*steps, CompiledStageKind::Static { integrator, algorithm, convergence, hold_patterns_after })
                }
                // A single eigensolve, not an iterative step loop — always
                // exactly one `advance()` step (`StageSpec::Modal`'s doc
                // comment).
                StageSpec::Modal { modes, .. } => (1, CompiledStageKind::Modal { num_modes: *modes as usize }),
                StageSpec::Transient { steps, dt, damping, ground_motions, .. } => {
                    let damping = RayleighDamping::new(damping.alpha_m, damping.beta_k);
                    let ground_motions = ground_motions
                        .iter()
                        .enumerate()
                        .map(|(row, gm)| {
                            check_dof_within("sequence.groundMotion", row as u32, gm.direction, ndim)?;
                            Ok(GroundMotion::new(gm.direction as usize, load_series_of(&gm.series))
                                .with_scale_factor(gm.scale_factor))
                        })
                        .collect::<Result<Vec<_>, DecodeError>>()?;
                    (*steps, CompiledStageKind::Transient { damping, dt: *dt, ground_motions })
                }
            };

            Ok(CompiledStage {
                id: stage.id().to_string(),
                steps,
                kind,
                pending_nodal_loads,
                pending_element_loads,
            })
        })
        .collect()
}
