//! Spatial ("space: 3") counterpart to `decode.rs` — same hand-written
//! per-table translation strategy (module doc comment on `decode.rs`), just
//! calling `Domain3`/`Element3`/`Friction3`/`FiberSection3` constructors
//! instead of their planar counterparts. Only `compile_stages`/
//! `check_dof_within` (stage/integrator compilation, which never touches
//! element physics) are shared with `decode.rs` — everything else here is
//! its own per-table loop, since the element constructors and table shapes
//! genuinely differ (see tables3.rs's module doc comment).

use std::collections::HashMap;

use carapace_core::model::{
    Axis3, BeamIntegration, DispBeamColumn3, Domain3, ElasticBeamColumn3, Element3, Element3Id,
    ElementLoad3, Fiber3, FiberSection3, ForceBeamColumn3, Friction3, GeomTransf3, LoadPatternId,
    Material, Node3, Node3Id, Truss3, ZeroLength3, ZeroLengthSection3,
};

use super::decode::{check_dof_within, compile_stages, load_series_of};
use super::error::DecodeError;
use super::materials::resolve_materials;
use super::sequence::RecorderSpec3;
use super::session::{ResolvedRecorder, SpatialSession};
use super::tables::IntegrationSpec;
use super::tables3::{Axis3Spec, ElementKind3, ElementLoadSpec3, TransformSpec3};
use super::CarapaceInputV1;

const NDOF3: u8 = 6;

pub(super) fn decode_spatial(input: CarapaceInputV1) -> Result<SpatialSession, DecodeError> {
    let materials = resolve_materials(&input.materials)?;
    let material_at = |row: u32, table: &'static str| -> Result<Material, DecodeError> {
        materials
            .get(row as usize)
            .cloned()
            .ok_or(DecodeError::UnknownMaterialIndex { table, row })
    };

    let mut domain = Domain3::new();
    let node_ids = add_nodes(&mut domain, &input.nodes3);
    let node_at = |row: u32, table: &'static str| -> Result<Node3Id, DecodeError> {
        node_ids
            .get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownNodeIndex { table, row })
    };
    add_equal_dofs(&mut domain, &input.equal_dofs3, &node_at)?;
    add_rigid_diaphragms(&mut domain, &input.rigid_diaphragms3, &node_at)?;
    let fiber_section_at = |section: u32| -> Result<Vec<Fiber3>, DecodeError> {
        build_fiber_section(&input.fibers3, section, &material_at)
    };

    let mut truss_ids = Vec::with_capacity(input.trusses3.node_i.len());
    for i in 0..input.trusses3.node_i.len() {
        let element = Truss3::new(
            node_at(input.trusses3.node_i[i], "trusses")?,
            node_at(input.trusses3.node_j[i], "trusses")?,
            input.trusses3.area[i],
            material_at(input.trusses3.material[i], "trusses")?,
        )
        .with_density(input.trusses3.density[i]);
        truss_ids.push(domain.add_element(Element3::Truss3(element)));
    }

    let mut elastic_beam_ids = Vec::with_capacity(input.elastic_beam_columns3.node_i.len());
    for i in 0..input.elastic_beam_columns3.node_i.len() {
        let table = &input.elastic_beam_columns3;
        let element = ElasticBeamColumn3::new(
            node_at(table.node_i[i], "elastic_beam_columns")?,
            node_at(table.node_j[i], "elastic_beam_columns")?,
            table.e[i],
            table.g[i],
            table.a[i],
            table.j[i],
            table.iy[i],
            table.iz[i],
            transform_of(table.transform[i]),
        )
        .with_density(table.density[i]);
        elastic_beam_ids.push(domain.add_element(Element3::ElasticBeamColumn3(element)));
    }

    let mut disp_beam_ids = Vec::with_capacity(input.disp_beam_columns3.node_i.len());
    for i in 0..input.disp_beam_columns3.node_i.len() {
        let table = &input.disp_beam_columns3;
        let element = DispBeamColumn3::new(
            node_at(table.node_i[i], "disp_beam_columns")?,
            node_at(table.node_j[i], "disp_beam_columns")?,
            table.g[i],
            table.j[i],
            table.vec_xz[i],
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        disp_beam_ids.push(domain.add_element(Element3::DispBeamColumn3(element)));
    }

    let mut force_beam_ids = Vec::with_capacity(input.force_beam_columns3.node_i.len());
    for i in 0..input.force_beam_columns3.node_i.len() {
        let table = &input.force_beam_columns3;
        let element = ForceBeamColumn3::new(
            node_at(table.node_i[i], "force_beam_columns")?,
            node_at(table.node_j[i], "force_beam_columns")?,
            table.g[i],
            table.j[i],
            table.vec_xz[i],
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        force_beam_ids.push(domain.add_element(Element3::ForceBeamColumn3(element)));
    }

    let zero_length_ids =
        add_zero_lengths(&mut domain, &input.zero_lengths3, &node_at, &material_at)?;
    let zero_length_section_ids = add_zero_length_sections(
        &mut domain,
        &input.zero_length_sections3,
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

    // See decode.rs's own comment on `nodal_loads_by_stage`/
    // `element_loads_by_stage` — identical per-stage registration timing,
    // just against `Node3Id`/`Element3Id`/`ElementLoad3`.
    let stage_count = input.sequence3.stages.len();
    let mut nodal_loads_by_stage: Vec<Vec<(LoadPatternId, Node3Id, usize, f64)>> =
        vec![Vec::new(); stage_count];
    for i in 0..input.nodal_loads3.node.len() {
        let dof = input.nodal_loads3.dof[i];
        check_dof_within("nodal_loads", i as u32, dof, NDOF3)?;
        let stage = input.nodal_loads3.stage[i];
        let entry =
            nodal_loads_by_stage
                .get_mut(stage as usize)
                .ok_or(DecodeError::UnknownStageIndex {
                    table: "nodal_loads",
                    row: i as u32,
                })?;
        entry.push((
            pattern_at(input.nodal_loads3.pattern[i])?,
            node_at(input.nodal_loads3.node[i], "nodal_loads")?,
            dof as usize,
            input.nodal_loads3.value[i],
        ));
    }

    let mut element_loads_by_stage: Vec<Vec<(LoadPatternId, Element3Id, ElementLoad3)>> =
        vec![Vec::new(); stage_count];
    for i in 0..input.element_loads3.element_index.len() {
        let pattern = pattern_at(input.element_loads3.pattern[i])?;
        let kind = input.element_loads3.element_kind[i];
        if kind != ElementKind3::ElasticBeamColumn {
            return Err(DecodeError::UnsupportedElementLoad {
                element_kind: element_kind_name(kind),
            });
        }
        let row = input.element_loads3.element_index[i];
        let element_id =
            *elastic_beam_ids
                .get(row as usize)
                .ok_or(DecodeError::UnknownElementIndex {
                    table: "elastic_beam_columns",
                    row,
                })?;
        let load = match input.element_loads3.load[i] {
            ElementLoadSpec3::UniformTransverse { wy, wz } => {
                ElementLoad3::UniformTransverse { wy, wz }
            }
        };
        let stage = input.element_loads3.stage[i];
        element_loads_by_stage
            .get_mut(stage as usize)
            .ok_or(DecodeError::UnknownStageIndex {
                table: "element_loads",
                row: i as u32,
            })?
            .push((pattern, element_id, load));
    }

    let stages = compile_stages(
        &input.sequence3.stages,
        &node_at,
        &pattern_at,
        NDOF3,
        3,
        nodal_loads_by_stage,
        element_loads_by_stage,
    )?;
    let element_at = |kind: ElementKind3, row: u32| -> Result<Element3Id, DecodeError> {
        let (ids, table) = match kind {
            ElementKind3::Truss => (&truss_ids, "trusses"),
            ElementKind3::ElasticBeamColumn => (&elastic_beam_ids, "elastic_beam_columns"),
            ElementKind3::DispBeamColumn => (&disp_beam_ids, "disp_beam_columns"),
            ElementKind3::ForceBeamColumn => (&force_beam_ids, "force_beam_columns"),
            ElementKind3::ZeroLength => (&zero_length_ids, "zero_lengths"),
            ElementKind3::ZeroLengthSection => (&zero_length_section_ids, "zero_length_sections"),
        };
        ids.get(row as usize)
            .copied()
            .ok_or(DecodeError::UnknownElementIndex { table, row })
    };
    let recorders = input
        .sequence3
        .recorders
        .iter()
        .map(|recorder| {
            Ok(match *recorder {
                RecorderSpec3::NodeDisp { node, dof } => ResolvedRecorder::NodeDisp {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec3::NodeVel { node, dof } => ResolvedRecorder::NodeVel {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec3::NodeAccel { node, dof } => ResolvedRecorder::NodeAccel {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec3::ElementForce {
                    element_kind,
                    element_index,
                    component,
                } => ResolvedRecorder::ElementForce {
                    element: element_at(element_kind, element_index)?,
                    component,
                },
                RecorderSpec3::ModeShape { mode, node, dof } => ResolvedRecorder::ModeShape {
                    mode,
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec3::Reaction { node, dof } => ResolvedRecorder::Reaction {
                    node: node_at(node, "recorders")?,
                    dof,
                },
                RecorderSpec3::Fiber {
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

    Ok(SpatialSession::new(domain, stages, recorders))
}

fn add_nodes(domain: &mut Domain3, table: &super::tables3::NodeTable3) -> Vec<Node3Id> {
    let mut mass_by_node: HashMap<u32, [f64; 6]> = HashMap::new();
    for (k, &node_index) in table.mass_node_index.iter().enumerate() {
        let mut mass = [0.0; 6];
        mass.copy_from_slice(&table.mass[k * 6..k * 6 + 6]);
        mass_by_node.insert(node_index, mass);
    }

    let node_count = table.node_count();
    let mut node_ids = Vec::with_capacity(node_count);
    for i in 0..node_count {
        let mut node = Node3::new([
            table.coords[3 * i],
            table.coords[3 * i + 1],
            table.coords[3 * i + 2],
        ]);
        let fixed_bits = table.fixed[i];
        for dof in 0..6usize {
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

/// See `decode.rs`'s own `add_equal_dofs`.
fn add_equal_dofs(
    domain: &mut Domain3,
    table: &super::tables3::EqualDofTable3,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
) -> Result<(), DecodeError> {
    let mut dofs_by_row: Vec<Vec<usize>> = vec![Vec::new(); table.retained.len()];
    for &(row, dof) in &table.dofs {
        check_dof_within("equal_dofs", row, dof, NDOF3)?;
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

/// See `decode.rs`'s own `add_rigid_diaphragms` — the spatial version also
/// converts each row's [`Axis3Spec`] to `core::Axis3` and calls
/// `Domain3::rigid_diaphragm_about` instead of the planar (single-dof,
/// normal-free) `Domain::rigid_diaphragm`.
fn add_rigid_diaphragms(
    domain: &mut Domain3,
    table: &super::tables3::RigidDiaphragmTable3,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
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
    #[allow(clippy::needless_range_loop)] // parallel-indexes retained/normal/constrained_by_row
    for i in 0..table.retained.len() {
        let retained = node_at(table.retained[i], "rigid_diaphragms")?;
        let constrained = constrained_by_row[i]
            .iter()
            .map(|&row| node_at(row, "rigid_diaphragms"))
            .collect::<Result<Vec<_>, _>>()?;
        let normal = match table.normal[i] {
            Axis3Spec::X => Axis3::X,
            Axis3Spec::Y => Axis3::Y,
            Axis3Spec::Z => Axis3::Z,
        };
        domain.rigid_diaphragm_about(retained, &constrained, normal);
    }
    Ok(())
}

fn build_fiber_section(
    table: &super::tables3::FiberTable3,
    section: u32,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<Fiber3>, DecodeError> {
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
            Ok(Fiber3::new(table.y[i], table.z[i], table.area[i], material))
        })
        .collect()
}

fn add_zero_lengths(
    domain: &mut Domain3,
    table: &super::tables3::ZeroLengthTable3,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<Element3Id>, DecodeError> {
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

    type FrictionRow3 = (u8, u8, u8, f64, f64, f64);
    let mut friction_by_row: Vec<Option<FrictionRow3>> = vec![None; table.node_i.len()];
    for &(row, normal_dof, shear_dof_0, shear_dof_1, mu, k0, b) in &table.friction {
        *friction_by_row
            .get_mut(row as usize)
            .ok_or(DecodeError::UnknownElementIndex {
                table: "zero_lengths",
                row,
            })? = Some((normal_dof, shear_dof_0, shear_dof_1, mu, k0, b));
    }

    let mut ids = Vec::with_capacity(table.node_i.len());
    #[allow(clippy::needless_range_loop)] // parallel-indexes node_i/node_j/materials_by_row/friction_by_row
    for i in 0..table.node_i.len() {
        let mut element = ZeroLength3::new(
            node_at(table.node_i[i], "zero_lengths")?,
            node_at(table.node_j[i], "zero_lengths")?,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element =
                element.with_material(dof as usize, material_at(material_index, "zero_lengths")?);
        }
        if let Some((normal_dof, shear_dof_0, shear_dof_1, mu, k0, b)) = friction_by_row[i] {
            element = element.with_friction(Friction3::new(
                normal_dof as usize,
                [shear_dof_0 as usize, shear_dof_1 as usize],
                mu,
                k0,
                b,
            ));
        }
        ids.push(domain.add_element(Element3::ZeroLength3(element)));
    }
    Ok(ids)
}

fn add_zero_length_sections(
    domain: &mut Domain3,
    table: &super::tables3::ZeroLengthSectionTable3,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
    fiber_section_at: &impl Fn(u32) -> Result<Vec<Fiber3>, DecodeError>,
) -> Result<Vec<Element3Id>, DecodeError> {
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
        let section = FiberSection3::new(fiber_section_at(table.fiber_section[i])?);
        let mut element = ZeroLengthSection3::new(
            node_at(table.node_i[i], "zero_length_sections")?,
            node_at(table.node_j[i], "zero_length_sections")?,
            section,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element
                .with_material(dof as usize, material_at(material_index, "zero_length_sections")?);
        }
        ids.push(domain.add_element(Element3::ZeroLengthSection3(element)));
    }
    Ok(ids)
}

fn add_load_patterns(
    domain: &mut Domain3,
    table: &super::tables::LoadPatternTable,
) -> Vec<LoadPatternId> {
    (0..table.series.len())
        .map(|i| domain.add_load_pattern_scaled(load_series_of(&table.series[i]), table.scale_factor[i]))
        .collect()
}

fn transform_of(spec: TransformSpec3) -> GeomTransf3 {
    match spec {
        TransformSpec3::Linear3 { vec_xz } => GeomTransf3::linear(vec_xz),
        TransformSpec3::PDelta3 { vec_xz } => GeomTransf3::p_delta(vec_xz),
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

fn element_kind_name(kind: ElementKind3) -> &'static str {
    match kind {
        ElementKind3::Truss => "truss",
        ElementKind3::ElasticBeamColumn => "elastic_beam_column",
        ElementKind3::DispBeamColumn => "disp_beam_column",
        ElementKind3::ForceBeamColumn => "force_beam_column",
        ElementKind3::ZeroLength => "zero_length",
        ElementKind3::ZeroLengthSection => "zero_length_section",
    }
}
