//! The 3D (`ndm: 3`) element decoder — see `elements_2d.rs`; this one calls
//! the `Domain3`/`Element3`/`Friction3`/`FiberSection3` constructors.

use carapace_core::model::{
    Axis3, ChordDistribution, Diaphragm, DiaphragmLayout, DispBeamColumn3, Domain3,
    ElasticBeamColumn3, Element3, Element3Id, ElementLoad3, Fiber3, FiberSection3, FloorZone,
    ForceBeamColumn3, Friction3, GeomTransf3, Material, Node3Id, Shell3, Shell4, ShellModifiers,
    ShellSection, StudWall, StudWallLayout, Truss3, WallChord, WallLoad, WallOpening, WallShearLaw,
    WallSheathing, WallStation, WallStud, ZeroLength3, ZeroLengthSection3,
};

use super::shared::{
    add_equal_dofs, add_linear_constraints, add_nodes, add_rigid_diaphragms, add_rigid_links,
    build_fiber_section, check_fibers, finish, friction_by_row, integration_of, material_lookup,
    materials_by_row, node_lookup, orientations_by_row, reject_tables, validate_stud_wall_geometry,
    wall_load_plan, ElementIds,
};
use crate::building::{ChordDistributionSpec, WallChordSpec, WallShearLawSpec};
use crate::input_v1::error::DecodeError;
use crate::input_v1::materials::resolve_materials;
use crate::input_v1::model_output::ModelOutput;
use crate::input_v1::session::Session3;
use crate::input_v1::tables::{
    Axis3Spec, ElementKind, ElementLoadSpec, Shell3Table, Shell4Table, ShellSectionSpec,
    TransformSpec3, ZeroLengthSectionTable, ZeroLengthTable,
};
use crate::input_v1::CarapaceInputV1;

const NDM: u8 = 3;

pub(super) fn decode(input: CarapaceInputV1) -> Result<Session3, DecodeError> {
    reject_tables(&[
        (
            "elastic_beam_columns_2d",
            input.elastic_beam_columns_2d.node_i.is_empty(),
        ),
        (
            "disp_beam_columns_2d",
            input.disp_beam_columns_2d.node_i.is_empty(),
        ),
        (
            "force_beam_columns_2d",
            input.force_beam_columns_2d.node_i.is_empty(),
        ),
        ("plane_materials", input.plane_materials.is_empty()),
        ("embedded_bars", input.embedded_bars.is_empty()),
        (
            "triangles",
            input.triangles.thickness.is_empty() && input.triangles.node_ids.is_empty(),
        ),
        (
            "quads",
            input.quads.thickness.is_empty() && input.quads.node_ids.is_empty(),
        ),
    ])?;
    check_fibers(&input.fibers, NDM)?;
    let materials = resolve_materials(&input.materials)?;
    let material_at = material_lookup(&materials);

    let mut domain = Domain3::new();
    let node_ids = add_nodes(&mut domain, &input.nodes)?;
    let node_at = node_lookup(&node_ids);
    add_equal_dofs(&mut domain, &input.equal_dofs, &node_at)?;
    add_rigid_diaphragms(
        &input.rigid_diaphragms,
        NDM,
        &node_at,
        |retained, constrained, normal| {
            let normal = match normal {
                Axis3Spec::X => Axis3::X,
                Axis3Spec::Y => Axis3::Y,
                Axis3Spec::Z => Axis3::Z,
            };
            domain.rigid_diaphragm_about(retained, constrained, normal)
        },
    )?;
    add_rigid_links(&input.rigid_links, &node_at, |master, slave| {
        domain.rigid_link(master, slave)
    })?;
    add_linear_constraints(&mut domain, &input.linear_constraints, &node_at)?;
    let fiber_section_at = |section: u32| -> Result<Vec<Fiber3>, DecodeError> {
        build_fiber_section(&input.fibers, section, &material_at, |i, material| {
            Fiber3::new(
                input.fibers.y[i],
                input.fibers.z[i],
                input.fibers.area[i],
                material,
            )
        })
    };

    let mut elements = ElementIds::new();

    let table = &input.trusses;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = Truss3::new(
            node_at(table.node_i[i], "trusses")?,
            node_at(table.node_j[i], "trusses")?,
            table.area[i],
            material_at(table.material[i], "trusses")?,
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element3::Truss3(element)));
    }
    elements.set(ElementKind::Truss, ids);

    let table = &input.elastic_beam_columns_3d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = ElasticBeamColumn3::new(
            node_at(table.node_i[i], "elastic_beam_columns_3d")?,
            node_at(table.node_j[i], "elastic_beam_columns_3d")?,
            table.e[i],
            table.g[i],
            table.a[i],
            table.j[i],
            table.iy[i],
            table.iz[i],
            transform_of(table.transform[i]),
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element3::ElasticBeamColumn3(element)));
    }
    elements.set(ElementKind::ElasticBeamColumn3d, ids);

    let table = &input.disp_beam_columns_3d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = DispBeamColumn3::new(
            node_at(table.node_i[i], "disp_beam_columns_3d")?,
            node_at(table.node_j[i], "disp_beam_columns_3d")?,
            table.g[i],
            table.j[i],
            table.vec_xz[i],
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element3::DispBeamColumn3(element)));
    }
    elements.set(ElementKind::DispBeamColumn3d, ids);

    let table = &input.force_beam_columns_3d;
    let mut ids = Vec::with_capacity(table.node_i.len());
    for i in 0..table.node_i.len() {
        let element = ForceBeamColumn3::new(
            node_at(table.node_i[i], "force_beam_columns_3d")?,
            node_at(table.node_j[i], "force_beam_columns_3d")?,
            table.g[i],
            table.j[i],
            table.vec_xz[i],
            fiber_section_at(table.fiber_section[i])?,
            integration_of(table.integration[i]),
        )
        .with_density(table.density[i]);
        ids.push(domain.add_element(Element3::ForceBeamColumn3(element)));
    }
    elements.set(ElementKind::ForceBeamColumn3d, ids);

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

    let shell_sections = resolve_shell_sections(&input.shell_sections)?;
    elements.set(
        ElementKind::Shell3,
        add_shell3s(&mut domain, &input.shell3s, &node_at, &shell_sections)?,
    );
    elements.set(
        ElementKind::Shell4,
        add_shell4s(&mut domain, &input.shell4s, &node_at, &shell_sections)?,
    );

    let wall_plan = wall_load_plan(&input, input.stud_walls.len())?;
    elements.set(
        ElementKind::StudWall,
        decode_stud_walls3(
            &mut domain,
            &input.stud_walls,
            &input.nodes,
            &node_at,
            &wall_plan.cases,
        )?,
    );
    elements.set(
        ElementKind::Diaphragm,
        decode_diaphragms(&mut domain, &input.diaphragms, &node_at)?,
    );

    drop(node_at); // the opaque closure type keeps `node_ids` borrowed until dropped
    finish(
        domain,
        &input,
        node_ids,
        elements,
        &ModelOutput::default(),
        Vec::new(),
        &wall_plan.pairs,
        |spec| match *spec {
            ElementLoadSpec::Uniform { wx, wy, wz } => {
                Some(ElementLoad3::uniform(wx, wy, wz.unwrap_or(0.0)))
            }
            ElementLoadSpec::ShellPressure { pressure } => Some(ElementLoad3::pressure(pressure)),
            ElementLoadSpec::ShellBody { bx, by, bz } => Some(ElementLoad3::body(bx, by, bz)),
            _ => None,
        },
        |case| ElementLoad3 {
            wall: Some(case),
            ..Default::default()
        },
    )
}

fn layout3(l: &crate::building::WallLayoutSpec) -> StudWallLayout {
    StudWallLayout {
        length: l.length,
        height: l.height,
        studs: l
            .studs
            .iter()
            .map(|s| WallStud { x: s.x, ea: s.ea })
            .collect(),
        top_ea: l.top_ea,
        top_ei: l.top_ei,
        top_distribution: match l.top_distribution {
            ChordDistributionSpec::Beam => ChordDistribution::Beam,
            ChordDistributionSpec::Direct => ChordDistribution::Direct,
        },
        bottom_ea: l.bottom_ea,
        bottom_ei: l.bottom_ei,
        openings: l
            .openings
            .iter()
            .map(|o| WallOpening {
                x0: o.x0,
                x1: o.x1,
                head_height: o.head_height,
                header_ei: o.header_ei,
                jack_ea: o.jack_ea,
                door: o.door,
            })
            .collect(),
        sheathing: l.sheathing.as_ref().map(|s| WallSheathing {
            law: match s.law {
                WallShearLawSpec::Linear { ga } => WallShearLaw::Linear { ga },
                WallShearLawSpec::O86 {
                    bv,
                    nail_spacing,
                    nail_diameter,
                } => WallShearLaw::O86 {
                    bv,
                    nail_spacing,
                    nail_diameter,
                },
            },
        }),
        foundation: l.foundation,
        floor: l.floor.map(|z| FloorZone {
            depth: z.depth,
            ga: z.ga,
            rim_ea: z.rim_ea,
        }),
    }
}
fn decode_stud_walls3(
    d: &mut Domain3,
    t: &[crate::input_v1::tables::StudWallSpec],
    node_table: &crate::input_v1::tables::NodeTable,
    node: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    cases: &[Vec<Vec<WallLoad>>],
) -> Result<Vec<Element3Id>, DecodeError> {
    let mut out = Vec::with_capacity(t.len());
    for (i, r) in t.iter().enumerate() {
        let c = r
            .corners
            .map(|x| node(x, "stud_walls"))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .unwrap();
        let st = r
            .stations
            .iter()
            .map(|s| {
                Ok(WallStation {
                    chord: match s.chord {
                        WallChordSpec::Top => WallChord::Top,
                        WallChordSpec::Bottom => WallChord::Bottom,
                    },
                    x: s.x,
                    node: node(s.node, "stud_walls")?,
                })
            })
            .collect::<Result<Vec<_>, DecodeError>>()?;
        validate_stud_wall_geometry(node_table, 3, r).map_err(|reason| {
            DecodeError::InvalidMacroElement {
                table: "stud_walls",
                row: i as u32,
                reason: reason.into(),
            }
        })?;
        let w = StudWall::new3(c, st, layout3(&r.layout), cases[i].clone()).map_err(|e| {
            DecodeError::InvalidMacroElement {
                table: "stud_walls",
                row: i as u32,
                reason: e.to_string(),
            }
        })?;
        out.push(d.add_element(Element3::StudWall(w)));
    }
    Ok(out)
}
fn decode_diaphragms(
    d: &mut Domain3,
    t: &[crate::input_v1::tables::DiaphragmSpec],
    node: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
) -> Result<Vec<Element3Id>, DecodeError> {
    let mut out = Vec::with_capacity(t.len());
    for (i, r) in t.iter().enumerate() {
        let l = &r.layout;
        let laws = l
            .bays
            .iter()
            .map(|law| {
                law.map(|x| match x {
                    WallShearLawSpec::Linear { ga } => {
                        carapace_core::model::WallShearLaw::Linear { ga }
                    }
                    WallShearLawSpec::O86 {
                        bv,
                        nail_spacing,
                        nail_diameter,
                    } => carapace_core::model::WallShearLaw::O86 {
                        bv,
                        nail_spacing,
                        nail_diameter,
                    },
                })
            })
            .collect();
        let nodes = r
            .nodes
            .iter()
            .map(|&x| node(x, "diaphragms"))
            .collect::<Result<Vec<_>, _>>()?;
        let e = Diaphragm::new(
            nodes,
            DiaphragmLayout {
                u: l.u.clone(),
                v: l.v.clone(),
                level: l.level,
                bays: laws,
                ea_u: l.ea_u.clone(),
                ea_v: l.ea_v.clone(),
            },
        )
        .map_err(|reason| DecodeError::InvalidMacroElement {
            table: "diaphragms",
            row: i as u32,
            reason,
        })?;
        out.push(d.add_element(Element3::Diaphragm(e)));
    }
    Ok(out)
}

fn resolve_shell_sections(specs: &[ShellSectionSpec]) -> Result<Vec<ShellSection>, DecodeError> {
    specs
        .iter()
        .enumerate()
        .map(|(index, spec)| {
            let ShellSectionSpec::ElasticMembranePlate {
                e,
                nu,
                h,
                rho,
                modifiers,
            } = *spec;
            let spec = modifiers.unwrap_or_default();
            let modifiers = ShellModifiers {
                f11: spec.f11.unwrap_or(1.0),
                f22: spec.f22.unwrap_or(1.0),
                f12: spec.f12.unwrap_or(1.0),
                m11: spec.m11.unwrap_or(1.0),
                m22: spec.m22.unwrap_or(1.0),
                m12: spec.m12.unwrap_or(1.0),
                v13: spec.v13.unwrap_or(1.0),
                v23: spec.v23.unwrap_or(1.0),
            };
            ShellSection::elastic_membrane_plate_modified(e, nu, h, rho, modifiers).map_err(
                |error| DecodeError::InvalidShellSection {
                    index: index as u32,
                    reason: error.message(),
                },
            )
        })
        .collect()
}

fn add_shell3s(
    domain: &mut Domain3,
    table: &Shell3Table,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    sections: &[ShellSection],
) -> Result<Vec<Element3Id>, DecodeError> {
    const TABLE: &str = "shell3s";
    let rows = table.section.len();
    if table.node_ids.len() != 3 * rows {
        return Err(DecodeError::InvalidRow {
            table: TABLE,
            row: 0,
            reason: "nodeIds must have stride 3",
        });
    }
    (0..rows)
        .map(|i| {
            let n = |k: usize| node_at(table.node_ids[3 * i + k], TABLE);
            let section = sections.get(table.section[i] as usize).ok_or(
                DecodeError::UnknownMaterialIndex {
                    table: TABLE,
                    row: table.section[i],
                },
            )?;
            let element = Shell3::new([n(0)?, n(1)?, n(2)?], section.clone());
            Ok(domain.add_element(Element3::Shell3(element)))
        })
        .collect()
}

fn add_shell4s(
    domain: &mut Domain3,
    table: &Shell4Table,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    sections: &[ShellSection],
) -> Result<Vec<Element3Id>, DecodeError> {
    const TABLE: &str = "shell4s";
    let rows = table.section.len();
    if table.node_ids.len() != 4 * rows {
        return Err(DecodeError::InvalidRow {
            table: TABLE,
            row: 0,
            reason: "nodeIds must have stride 4",
        });
    }
    (0..rows)
        .map(|i| {
            let n = |k: usize| node_at(table.node_ids[4 * i + k], TABLE);
            let section = sections.get(table.section[i] as usize).ok_or(
                DecodeError::UnknownMaterialIndex {
                    table: TABLE,
                    row: table.section[i],
                },
            )?;
            let element = Shell4::new([n(0)?, n(1)?, n(2)?, n(3)?], section.clone());
            Ok(domain.add_element(Element3::Shell4(element)))
        })
        .collect()
}

fn add_zero_lengths(
    domain: &mut Domain3,
    table: &ZeroLengthTable,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
) -> Result<Vec<Element3Id>, DecodeError> {
    const TABLE: &str = "zero_lengths";
    let rows = table.node_i.len();
    let materials_by_row = materials_by_row(&table.materials, rows, TABLE)?;
    let friction_by_row = friction_by_row(&table.friction, rows, 2, TABLE)?;
    let orientation_by_row = orientations_by_row(&table.orient, rows, NDM, TABLE)?;

    let mut ids = Vec::with_capacity(rows);
    #[allow(clippy::needless_range_loop)]
    // parallel-indexes node_i/node_j/materials_by_row/friction_by_row/orientation_by_row
    for i in 0..rows {
        let mut element = ZeroLength3::new(
            node_at(table.node_i[i], TABLE)?,
            node_at(table.node_j[i], TABLE)?,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element.with_material(dof as usize, material_at(material_index, TABLE)?);
        }
        if let Some(friction) = friction_by_row[i] {
            element = element.with_friction(Friction3::new(
                friction.normal_dof as usize,
                [
                    friction.shear_dofs[0] as usize,
                    friction.shear_dofs[1] as usize,
                ],
                friction.mu,
                friction.k0,
                friction.b,
            ));
        }
        if let Some(orientation) = orientation_by_row[i] {
            element = element.with_orientation(orientation);
        }
        ids.push(domain.add_element(Element3::ZeroLength3(element)));
    }
    Ok(ids)
}

fn add_zero_length_sections(
    domain: &mut Domain3,
    table: &ZeroLengthSectionTable,
    node_at: &impl Fn(u32, &'static str) -> Result<Node3Id, DecodeError>,
    material_at: &impl Fn(u32, &'static str) -> Result<Material, DecodeError>,
    fiber_section_at: &impl Fn(u32) -> Result<Vec<Fiber3>, DecodeError>,
) -> Result<Vec<Element3Id>, DecodeError> {
    const TABLE: &str = "zero_length_sections";
    let rows = table.node_i.len();
    let materials_by_row = materials_by_row(&table.materials, rows, TABLE)?;
    let orientation_by_row = orientations_by_row(&table.orient, rows, NDM, TABLE)?;

    let mut ids = Vec::with_capacity(rows);
    #[allow(clippy::needless_range_loop)]
    // parallel-indexes node_i/node_j/fiber_section/materials_by_row/orientation_by_row
    for i in 0..rows {
        let section = FiberSection3::new(fiber_section_at(table.fiber_section[i])?);
        let mut element = ZeroLengthSection3::new(
            node_at(table.node_i[i], TABLE)?,
            node_at(table.node_j[i], TABLE)?,
            section,
        );
        for &(dof, material_index) in &materials_by_row[i] {
            element = element.with_material(dof as usize, material_at(material_index, TABLE)?);
        }
        if let Some(orientation) = orientation_by_row[i] {
            element = element.with_orientation(orientation);
        }
        ids.push(domain.add_element(Element3::ZeroLengthSection3(element)));
    }
    Ok(ids)
}

fn transform_of(spec: TransformSpec3) -> GeomTransf3 {
    match spec {
        TransformSpec3::Linear3 { vec_xz } => GeomTransf3::linear(vec_xz),
        TransformSpec3::PDelta3 { vec_xz } => GeomTransf3::p_delta(vec_xz),
    }
}
