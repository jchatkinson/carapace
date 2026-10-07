//! The 2D (`ndm: 2`) element decoder: `Domain`/`Element`/`Friction`/
//! `FiberSection` constructors from the per-formulation tables, then the
//! shared tail (`shared::finish`). Same hand-written per-table translation
//! strategy as `elements_3d.rs`; only the element constructors and table
//! shapes differ (see `tables.rs`).

use carapace_core::model::{
    DispBeamColumn, Domain, ElasticBeamColumn, Element, ElementId, ElementLoad, Fiber,
    FiberSection, ForceBeamColumn, Friction, GeomTransf, Material, NodeId, PlaneMaterial, Quad4,
    Quad4Formulation, Tri3, Truss, ZeroLength, ZeroLengthSection,
};

use super::shared::{
    add_equal_dofs, add_linear_constraints, add_nodes, add_rigid_diaphragms, add_rigid_links,
    build_fiber_section, check_fibers, finish, friction_by_row, integration_of, material_lookup,
    materials_by_row, node_lookup, orientations_by_row, reject_tables, ElementIds,
};
use crate::input_v1::error::DecodeError;
use crate::input_v1::materials::resolve_materials;
use crate::input_v1::session::Session2;
use crate::input_v1::tables::{
    ElementKind, ElementLoadSpec, PlaneMaterialSpec, PlaneStateSpec, Quad4FormulationSpec,
    QuadTable, TransformSpec, TriangleTable, ZeroLengthSectionTable, ZeroLengthTable,
};
use crate::input_v1::CarapaceInputV1;

const NDM: u8 = 2;

pub(super) fn decode(input: CarapaceInputV1) -> Result<Session2, DecodeError> {
    reject_tables(&[
        (
            "elastic_beam_columns_3d",
            input.elastic_beam_columns_3d.node_i.is_empty(),
        ),
        (
            "disp_beam_columns_3d",
            input.disp_beam_columns_3d.node_i.is_empty(),
        ),
        (
            "force_beam_columns_3d",
            input.force_beam_columns_3d.node_i.is_empty(),
        ),
    ])?;
    check_fibers(&input.fibers, NDM)?;
    let materials = resolve_materials(&input.materials)?;
    let material_at = material_lookup(&materials);

    let mut domain: Domain = Domain::new();
    let node_ids = add_nodes(&mut domain, &input.nodes)?;
    let node_at = node_lookup(&node_ids);
    add_equal_dofs(&mut domain, &input.equal_dofs, &node_at)?;
    add_rigid_diaphragms(
        &input.rigid_diaphragms,
        NDM,
        &node_at,
        |retained, constrained, _| domain.rigid_diaphragm(retained, constrained),
    )?;
    add_rigid_links(&input.rigid_links, &node_at, |master, slave| {
        domain.rigid_link(master, slave)
    })?;
    add_linear_constraints(&mut domain, &input.linear_constraints, &node_at)?;
    let fiber_section_at = |section: u32| -> Result<Vec<Fiber>, DecodeError> {
        build_fiber_section(&input.fibers, section, &material_at, |i, material| {
            Fiber::new(input.fibers.y[i], input.fibers.area[i], material)
        })
    };

    let mut elements = ElementIds::new();

    let table = &input.trusses;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = Truss::new(
            node_at(table.node_i[i], "trusses")?,
            node_at(table.node_j[i], "trusses")?,
            table.area[i],
            material_at(table.material[i], "trusses")?,
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element::Truss(element)));
    }
    elements.set(ElementKind::Truss, ids);

    let table = &input.elastic_beam_columns_2d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = ElasticBeamColumn::new(
            node_at(table.node_i[i], "elastic_beam_columns_2d")?,
            node_at(table.node_j[i], "elastic_beam_columns_2d")?,
            table.e[i],
            table.a[i],
            table.iz[i],
            transform_of(table.transform[i]),
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element::ElasticBeamColumn(element)));
    }
    elements.set(ElementKind::ElasticBeamColumn2d, ids);

    let table = &input.disp_beam_columns_2d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let mut element = DispBeamColumn::new(
            node_at(table.node_i[i], "disp_beam_columns_2d")?,
            node_at(table.node_j[i], "disp_beam_columns_2d")?,
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        if table.corotational[i] {
            element = element.with_corotational();
        }
        ids.push(domain.add_element(Element::DispBeamColumn(element)));
    }
    elements.set(ElementKind::DispBeamColumn2d, ids);

    let table = &input.force_beam_columns_2d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let mut element = ForceBeamColumn::new(
            node_at(table.node_i[i], "force_beam_columns_2d")?,
            node_at(table.node_j[i], "force_beam_columns_2d")?,
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        if table.corotational[i] {
            element = element.with_corotational();
        }
        ids.push(domain.add_element(Element::ForceBeamColumn(element)));
    }
    elements.set(ElementKind::ForceBeamColumn2d, ids);

    elements.set(
        ElementKind::ZeroLength,
        add_zero_lengths(&mut domain, &input.zero_lengths, &node_at, &material_at)?,
    );
    elements.set(
        ElementKind::ZeroLengthSection,
        add_zero_length_sections(
            &mut domain,
            &input.zero_length_sections,
            &node_at,
            &material_at,
            &fiber_section_at,
        )?,
    );

    let plane_materials = resolve_plane_materials(&input.plane_materials)?;
    elements.set(
        ElementKind::Tri3,
        add_triangles(&mut domain, &input.triangles, &node_at, &plane_materials)?,
    );
    elements.set(
        ElementKind::Quad4,
        add_quads(&mut domain, &input.quads, &node_at, &plane_materials)?,
    );

    drop(node_at); // the opaque closure type keeps `node_ids` borrowed until dropped
    finish(domain, &input, node_ids, elements, |spec| {
        Some(match *spec {
            ElementLoadSpec::Uniform { wx, wy, .. } => ElementLoad::uniform(wx, wy),
            ElementLoadSpec::Body { bx, by } => ElementLoad::body(bx, by),
            ElementLoadSpec::EdgeTraction { edge, tx, ty } => {
                ElementLoad::edge_traction(edge as usize, tx, ty)
            }
            ElementLoadSpec::EdgePressure { edge, pressure } => {
                ElementLoad::edge_pressure(edge as usize, pressure)
            }
        })
    })
}

fn resolve_plane_materials(specs: &[PlaneMaterialSpec]) -> Result<Vec<PlaneMaterial>, DecodeError> {
    specs
        .iter()
        .enumerate()
        .map(|(index, spec)| {
            match *spec {
                PlaneMaterialSpec::Isotropic {
                    e,
                    nu,
                    state: PlaneStateSpec::PlaneStress,
                } => PlaneMaterial::plane_stress(e, nu),
                PlaneMaterialSpec::Isotropic {
                    e,
                    nu,
                    state: PlaneStateSpec::PlaneStrain,
                } => PlaneMaterial::plane_strain(e, nu),
                PlaneMaterialSpec::Orthotropic {
                    ex,
                    ey,
                    nu_xy,
                    g_xy,
                    angle,
                } => PlaneMaterial::orthotropic(ex, ey, nu_xy, g_xy, angle),
                PlaneMaterialSpec::ElasticMatrix { d } => PlaneMaterial::elastic_matrix(d),
            }
            .map_err(|error| DecodeError::InvalidPlaneMaterial {
                index: index as u32,
                reason: error.message(),
            })
        })
        .collect()
}

/// Checks that every per-row array of a continuum table has `rows` entries
/// (`node_ids` has `rows * nodes_per_row`).
fn check_continuum_lengths(
    table: &'static str,
    rows: usize,
    lengths: &[usize],
) -> Result<(), DecodeError> {
    if lengths.iter().any(|&len| len != rows) {
        return Err(DecodeError::InvalidRow {
            table,
            row: 0,
            reason: "every per-row array must have one entry per element",
        });
    }
    Ok(())
}

fn plane_material_at<'a>(
    materials: &'a [PlaneMaterial],
    table: &'static str,
    row: u32,
) -> Result<&'a PlaneMaterial, DecodeError> {
    materials
        .get(row as usize)
        .ok_or(DecodeError::UnknownMaterialIndex { table, row })
}

fn add_triangles(
    domain: &mut Domain,
    table: &TriangleTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    materials: &[PlaneMaterial],
) -> Result<Vec<ElementId>, DecodeError> {
    const TABLE: &str = "triangles";
    let rows = table.thickness.len();
    if table.node_ids.len() != 3 * rows {
        return Err(DecodeError::InvalidRow {
            table: TABLE,
            row: 0,
            reason: "nodeIds must have stride 3",
        });
    }
    check_continuum_lengths(TABLE, rows, &[table.material.len(), table.density.len()])?;
    (0..rows)
        .map(|i| {
            let n = |k: usize| node_at(table.node_ids[3 * i + k], TABLE);
            let element = Tri3::new(
                n(0)?,
                n(1)?,
                n(2)?,
                table.thickness[i],
                plane_material_at(materials, TABLE, table.material[i])?.clone(),
            )
            .with_density(table.density[i]);
            Ok(domain.add_element(Element::Tri3(element)))
        })
        .collect()
}

fn add_quads(
    domain: &mut Domain,
    table: &QuadTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    materials: &[PlaneMaterial],
) -> Result<Vec<ElementId>, DecodeError> {
    const TABLE: &str = "quads";
    let rows = table.thickness.len();
    if table.node_ids.len() != 4 * rows {
        return Err(DecodeError::InvalidRow {
            table: TABLE,
            row: 0,
            reason: "nodeIds must have stride 4",
        });
    }
    check_continuum_lengths(
        TABLE,
        rows,
        &[
            table.material.len(),
            table.density.len(),
            table.formulation.len(),
        ],
    )?;
    (0..rows)
        .map(|i| {
            let n = |k: usize| node_at(table.node_ids[4 * i + k], TABLE);
            let formulation = match table.formulation[i] {
                Quad4FormulationSpec::Full => Quad4Formulation::Full,
                Quad4FormulationSpec::Enhanced => Quad4Formulation::Enhanced,
            };
            let element = Quad4::new(
                [n(0)?, n(1)?, n(2)?, n(3)?],
                table.thickness[i],
                plane_material_at(materials, TABLE, table.material[i])?.clone(),
            )
            .with_density(table.density[i])
            .with_formulation(formulation);
            Ok(domain.add_element(Element::Quad4(element)))
        })
        .collect()
}

fn add_zero_lengths(
    domain: &mut Domain,
    table: &ZeroLengthTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<ElementId>, DecodeError> {
    const TABLE: &str = "zero_lengths";
    let rows = table.node_i.len();
    let materials_by_row = materials_by_row(&table.materials, rows, TABLE)?;
    let friction_by_row = friction_by_row(&table.friction, rows, 1, TABLE)?;
    let orientation_by_row = orientations_by_row(&table.orient, rows, NDM, TABLE)?;

    let mut ids = Vec::with_capacity(rows);
    #[allow(clippy::needless_range_loop)]
    // parallel-indexes node_i/node_j/materials_by_row/friction_by_row/orientation_by_row
    for i in 0..rows {
        let mut element = ZeroLength::new(
            node_at(table.node_i[i], TABLE)?,
            node_at(table.node_j[i], TABLE)?,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element.with_material(dof as usize, material_at(material_index, TABLE)?);
        }
        if let Some(friction) = friction_by_row[i] {
            element = element.with_friction(Friction::new(
                friction.normal_dof as usize,
                friction.shear_dofs[0] as usize,
                friction.mu,
                friction.k0,
                friction.b,
            ));
        }
        if let Some(orientation) = orientation_by_row[i] {
            element = element.with_orientation(orientation).map_err(|_| {
                DecodeError::InvalidOrientation {
                    table: TABLE,
                    row: i as u32,
                }
            })?;
        }
        ids.push(domain.add_element(Element::ZeroLength(element)));
    }
    Ok(ids)
}

fn add_zero_length_sections(
    domain: &mut Domain,
    table: &ZeroLengthSectionTable,
    node_at: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
    fiber_section_at: &impl Fn(u32) -> Result<Vec<Fiber>, DecodeError>,
) -> Result<Vec<ElementId>, DecodeError> {
    const TABLE: &str = "zero_length_sections";
    let rows = table.node_i.len();
    let materials_by_row = materials_by_row(&table.materials, rows, TABLE)?;
    let orientation_by_row = orientations_by_row(&table.orient, rows, NDM, TABLE)?;

    let mut ids = Vec::with_capacity(rows);
    #[allow(clippy::needless_range_loop)]
    // parallel-indexes node_i/node_j/fiber_section/materials_by_row/orientation_by_row
    for i in 0..rows {
        let section = FiberSection::new(fiber_section_at(table.fiber_section[i])?);
        let mut element = ZeroLengthSection::new(
            node_at(table.node_i[i], TABLE)?,
            node_at(table.node_j[i], TABLE)?,
            section,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element.with_material(dof as usize, material_at(material_index, TABLE)?);
        }
        if let Some(orientation) = orientation_by_row[i] {
            element = element.with_orientation(orientation).map_err(|_| {
                DecodeError::InvalidOrientation {
                    table: TABLE,
                    row: i as u32,
                }
            })?;
        }
        ids.push(domain.add_element(Element::ZeroLengthSection(element)));
    }
    Ok(ids)
}

fn transform_of(spec: TransformSpec) -> GeomTransf {
    match spec {
        TransformSpec::Linear => GeomTransf::Linear,
        TransformSpec::PDelta => GeomTransf::PDelta,
        TransformSpec::Corotational => GeomTransf::Corotational,
    }
}
