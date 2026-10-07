//! The 2D (`ndm: 2`) element decoder: `Domain`/`Element`/`Friction`/
//! `FiberSection` constructors from the per-formulation tables, then the
//! shared tail (`shared::finish`). Same hand-written per-table translation
//! strategy as `elements_3d.rs`; only the element constructors and table
//! shapes differ (see `tables.rs`).

use carapace_core::model::{
    DispBeamColumn, Domain, ElasticBeamColumn, Element, ElementId, ElementLoad, Fiber,
    FiberSection, ForceBeamColumn, Friction, GeomTransf, Material, NodeId, Truss, ZeroLength,
    ZeroLengthSection,
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
    ElementKind, TransformSpec, ZeroLengthSectionTable, ZeroLengthTable,
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

    drop(node_at); // the opaque closure type keeps `node_ids` borrowed until dropped
    finish(domain, &input, node_ids, elements, |wx, wy, _| {
        ElementLoad::Uniform { wx, wy }
    })
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
