//! The 2D (`ndm: 2`) element decoder: `Domain`/`Element`/`Friction`/
//! `FiberSection` constructors from the per-formulation tables, then the
//! shared tail (`shared::finish`). Same hand-written per-table translation
//! strategy as `elements_3d.rs`; only the element constructors and table
//! shapes differ (see `tables.rs`).

use carapace_core::model::continuum::{bond_polyline, embed_polyline_in, BondSpec, EndSpring};
use carapace_core::model::{
    ChordDistribution, CsfmConcreteLaw, CsfmRebarRouting, DispBeamColumn, Domain,
    ElasticBeamColumn, Element, ElementId, ElementLoad, EmbeddedBar, EmbeddedBarHost, Fiber,
    FiberSection, FloorBand, FloorBandLayout, ForceBeamColumn, Friction, GeomTransf, Material,
    NodeId, PlaneMaterial, Quad4, Quad4Formulation, StudWall, StudWallLayout, Tri3, Truss,
    WallChord, WallLoad, WallOpening, WallShearLaw, WallSheathing, WallStation, WallStud,
    ZeroLength, ZeroLengthSection,
};

use super::shared::{
    add_equal_dofs, add_linear_constraints, add_nodes, add_rigid_diaphragms, add_rigid_links,
    build_fiber_section, check_fibers, finish, floor_wall_loads, friction_by_row, integration_of,
    material_lookup, materials_by_row, node_lookup, orientations_by_row, reject_tables,
    validate_stud_wall_geometry, ElementIds,
};
use crate::building::{ChordDistributionSpec, WallChordSpec, WallLoadSpec, WallShearLawSpec};
use crate::input_v1::error::DecodeError;
use crate::input_v1::materials::resolve_materials;
use crate::input_v1::model_output::{
    BondGroup, BondPointInfo, EmbeddedBarGroup, EmbeddedBarSegment, ModelOutput, RebarRouting,
};
use crate::input_v1::session::Session2;
use crate::input_v1::tables::{
    CsfmConcreteLawSpec, ElementKind, ElementLoadSpec, EmbeddedBarHostsSpec, EmbeddedBarSpec,
    PlaneMaterialSpec, PlaneStateSpec, Quad4FormulationSpec, QuadTable, TransformSpec,
    TriangleTable, ZeroLengthSectionTable, ZeroLengthTable,
};
use crate::input_v1::CarapaceInputV1;

const NDM: u8 = 2;

pub(super) fn decode(input: CarapaceInputV1) -> Result<Session2, DecodeError> {
    reject_tables(&[
        ("diaphragms", input.diaphragms.is_empty()),
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
        ("shell_sections", input.shell_sections.is_empty()),
        (
            "shell3s",
            input.shell3s.node_ids.is_empty() && input.shell3s.section.is_empty(),
        ),
        (
            "shell4s",
            input.shell4s.node_ids.is_empty() && input.shell4s.section.is_empty(),
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

    let floor_wall_loads = floor_wall_loads(&input.floor_loads, input.stud_walls.len())?;
    elements.set(
        ElementKind::StudWall,
        decode_stud_walls(
            &mut domain,
            &input.stud_walls,
            &input.nodes,
            &node_at,
            &floor_wall_loads,
        )?,
    );
    elements.set(
        ElementKind::FloorBand,
        decode_floor_bands(&mut domain, &input.floor_bands, &node_at)?,
    );

    let mut bar_end_loads = Vec::new();
    let model_output = add_embedded_bars(
        &mut domain,
        &input.embedded_bars,
        &materials,
        &mut elements,
        &node_ids,
        &mut bar_end_loads,
    )?;

    drop(node_at); // the opaque closure type keeps `node_ids` borrowed until dropped
    let mut session = finish(
        domain,
        &input,
        node_ids,
        elements,
        &model_output,
        bar_end_loads,
        |spec| {
            Some(match *spec {
                ElementLoadSpec::Uniform { wx, wy, .. } => ElementLoad::uniform(wx, wy),
                ElementLoadSpec::Trapezoid {
                    wx_i,
                    wy_i,
                    wx_j,
                    wy_j,
                } => ElementLoad::trapezoid(wx_i, wy_i, wx_j, wy_j),
                ElementLoadSpec::Body { bx, by } => ElementLoad::body(bx, by),
                ElementLoadSpec::EdgeTraction { edge, tx, ty } => {
                    ElementLoad::edge_traction(edge as usize, tx, ty)
                }
                ElementLoadSpec::EdgePressure { edge, pressure } => {
                    ElementLoad::edge_pressure(edge as usize, pressure)
                }
                ElementLoadSpec::ShellPressure { .. } | ElementLoadSpec::ShellBody { .. } => {
                    return None
                }
                ElementLoadSpec::WallSlots { factors } => ElementLoad {
                    wall: factors,
                    ..Default::default()
                },
            })
        },
    )?;
    session.set_model_output(model_output);
    Ok(session)
}

fn decode_layout(l: &crate::building::WallLayoutSpec) -> StudWallLayout {
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
    }
}
fn decode_wall_load(l: WallLoadSpec) -> WallLoad {
    match l {
        WallLoadSpec::TopPoint { x, fx, fy } => WallLoad::TopPoint { x, fx, fy },
        WallLoadSpec::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        } => WallLoad::TopLine {
            x0,
            x1,
            wx0,
            wx1,
            wy0,
            wy1,
        },
    }
}
fn decode_stud_walls(
    domain: &mut Domain,
    t: &[crate::input_v1::tables::StudWallSpec],
    node_table: &crate::input_v1::tables::NodeTable,
    node: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
    extra: &[[Vec<WallLoad>; 8]],
) -> Result<Vec<ElementId>, DecodeError> {
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
        validate_stud_wall_geometry(node_table, 2, r).map_err(|reason| {
            DecodeError::InvalidMacroElement {
                table: "stud_walls",
                row: i as u32,
                reason: reason.into(),
            }
        })?;
        let mut loads: [Vec<WallLoad>; 8] = std::array::from_fn(|_| Vec::new());
        for x in &r.loads {
            if x.slot >= 8 {
                return Err(DecodeError::InvalidWallSlot {
                    table: "stud_walls",
                    row: i as u32,
                    slot: x.slot as u32,
                });
            }
            loads[x.slot].push(decode_wall_load(x.load));
        }
        for k in 0..8 {
            loads[k].extend(extra[i][k].iter().copied());
        }
        let wall = StudWall::new(c, st, decode_layout(&r.layout), loads).map_err(|e| {
            DecodeError::InvalidMacroElement {
                table: "stud_walls",
                row: i as u32,
                reason: e.to_string(),
            }
        })?;
        out.push(domain.add_element(Element::StudWall(wall)));
    }
    Ok(out)
}
fn decode_floor_bands(
    domain: &mut Domain,
    t: &[crate::input_v1::tables::FloorBandSpec],
    node: &impl Fn(u32, &'static str) -> Result<NodeId, DecodeError>,
) -> Result<Vec<ElementId>, DecodeError> {
    let mut out = Vec::with_capacity(t.len());
    for (i, r) in t.iter().enumerate() {
        let b = r
            .bottom
            .iter()
            .map(|&x| node(x, "floor_bands"))
            .collect::<Result<Vec<_>, _>>()?;
        let top = r
            .top
            .iter()
            .map(|&x| node(x, "floor_bands"))
            .collect::<Result<Vec<_>, _>>()?;
        let l = &r.layout;
        let band = FloorBand::new(
            b,
            top,
            FloorBandLayout {
                length: l.length,
                depth: l.depth,
                xs: l.xs.clone(),
                bays: l.bays.clone(),
                rim_ea: l.rim_ea,
            },
        )
        .map_err(|e| DecodeError::InvalidMacroElement {
            table: "floor_bands",
            row: i as u32,
            reason: e.to_string(),
        })?;
        out.push(domain.add_element(Element::FloorBand(band)));
    }
    Ok(out)
}

/// Discretizes each `embeddedBars` row through `embed_polyline_in` and registers its segments
/// as `EmbeddedBar` elements, in row order (segment indices are flat across bars).
fn add_embedded_bars(
    domain: &mut Domain,
    table: &[EmbeddedBarSpec],
    materials: &[Material],
    elements: &mut ElementIds<ElementId>,
    node_ids: &[NodeId],
    bar_end_loads: &mut Vec<(u32, u32, NodeId, usize, f64)>,
) -> Result<ModelOutput, DecodeError> {
    const TABLE: &str = "embedded_bars";
    let invalid = |row: usize, reason: &'static str| DecodeError::InvalidRow {
        table: TABLE,
        row: row as u32,
        reason,
    };
    let node_row = |id: NodeId| node_ids.iter().position(|&n| n == id).unwrap_or(0) as u32;
    let mut all_segments = Vec::new();
    let mut groups = Vec::with_capacity(table.len());
    for (row, bar) in table.iter().enumerate() {
        if bar.polyline.len() < 4 || bar.polyline.len() % 2 != 0 {
            return Err(invalid(
                row,
                "polyline must hold at least two (x, y) vertices",
            ));
        }
        if !(bar.area.is_finite() && bar.area > 0.0) || bar.polyline.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                row,
                "area must be finite and positive, and the polyline finite",
            ));
        }
        let material =
            materials
                .get(bar.material as usize)
                .ok_or(DecodeError::UnknownMaterialIndex {
                    table: TABLE,
                    row: bar.material,
                })?;
        let allowed = match &bar.hosts {
            None | Some(EmbeddedBarHostsSpec::All(_)) => None,
            Some(EmbeddedBarHostsSpec::List(list)) => Some(
                list.iter()
                    .map(|host| {
                        if !matches!(host.kind, ElementKind::Tri3 | ElementKind::Quad4) {
                            return Err(invalid(row, "a host must be a tri3 or quad4 element"));
                        }
                        elements.get(host.kind, host.index, NDM)
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            ),
        };
        let polyline: Vec<[f64; 2]> = bar.polyline.chunks(2).map(|v| [v[0], v[1]]).collect();
        let segments = embed_polyline_in(domain, &polyline, bar.area, material, allowed.as_deref())
            .map_err(|error| DecodeError::InvalidModel {
                error: error.into(),
            })?;

        let first_segment = all_segments.len() as u32;
        let mut reported = Vec::with_capacity(segments.len());
        let mut bond_group = None;
        let mut rods: Vec<EmbeddedBar> = Vec::new();
        let restrained = [
            bar.anchored_start.unwrap_or(false),
            bar.anchored_end.unwrap_or(false),
        ];
        // A bar with a bond law, or perfectly bonded with an end anchored beyond the mesh, is built from rod nodes.
        let rebuilt = bar.bond.is_some() || restrained.iter().any(|&r| r);
        if rebuilt {
            let spec = match &bar.bond {
                Some(wire) => BondSpec {
                    perimeter: wire.perimeter,
                    fbd: wire.fbd,
                    gb: wire.gb,
                    rf: wire.rf.unwrap_or(1e-5),
                    slip_limit: wire.slip_limit,
                    lock_angle: wire.lock_angle.unwrap_or(15.0).to_radians(),
                    ends: [wire.start, wire.end].map(|end| {
                        end.map(|e| EndSpring {
                            fau: e.fau,
                            ka: e.ka,
                            ksn: e.ksn.unwrap_or(1e-2),
                        })
                    }),
                    perfect: false,
                    restrained,
                },
                None => BondSpec {
                    perimeter: 1.0,
                    fbd: 1.0,
                    gb: 1.0,
                    rf: 1e-5,
                    slip_limit: 1.0,
                    lock_angle: 15f64.to_radians(),
                    ends: [None, None],
                    perfect: true,
                    restrained,
                },
            };
            spec.validate().map_err(|reason| invalid(row, reason))?;
            let built =
                bond_polyline(domain, &segments, bar.area, material, &spec).map_err(|error| {
                    DecodeError::InvalidModel {
                        error: error.into(),
                    }
                })?;
            for (end, load) in [(0, bar.start_load), (built.points.len() - 1, bar.end_load)] {
                let Some(load) = load else { continue };
                let point = &built.points[end];
                let (Some(node), false) = (point.rod, restrained[usize::from(end != 0)]) else {
                    return Err(invalid(
                        row,
                        "a bar end load needs a bond spec or an anchored other end, and an end that is not anchored",
                    ));
                };
                // Positive `force` is tension: outward along the bar, against the tangent at the first end.
                let sign = if end == 0 { -1.0 } else { 1.0 };
                let stage = load.stage.unwrap_or(0);
                for (dof, component) in point.tangent.iter().enumerate() {
                    if *component != 0.0 {
                        bar_end_loads.push((
                            load.pattern,
                            stage,
                            node,
                            dof,
                            sign * load.force * component,
                        ));
                    }
                }
            }
            let mut points: Vec<BondPointInfo> = built
                .points
                .iter()
                .map(|p| BondPointInfo {
                    position: p.coords,
                    tangent: p.tangent,
                    locked: p.rod.is_none(),
                    tributary: 0.0,
                    spring: None,
                    end_spring: None,
                })
                .collect();
            let first_row = elements.count(ElementKind::ZeroLength) as u32;
            let mut ids = Vec::with_capacity(built.springs.len());
            for (offset, spring) in built.springs.into_iter().enumerate() {
                let point = &mut points[spring.point];
                if spring.end {
                    point.end_spring = Some(first_row + offset as u32);
                } else {
                    point.spring = Some(first_row + offset as u32);
                    point.tributary = spring.tributary;
                }
                ids.push(domain.add_element(Element::ZeroLength(spring.element)));
            }
            elements.extend(ElementKind::ZeroLength, ids);
            if bar.bond.is_some() {
                bond_group = Some(BondGroup { points });
            }
            rods = built.rods;
        }
        for (index, segment) in segments.into_iter().enumerate() {
            let [start, end] = segment
                .end_points()
                .expect("embed_polyline sets end points");
            let (host_kind, host_nodes) = match segment.host() {
                EmbeddedBarHost::Tri3 { nodes, .. } => (
                    ElementKind::Tri3,
                    nodes.iter().map(|&n| node_row(n)).collect(),
                ),
                EmbeddedBarHost::Quad4 { nodes, .. } => (
                    ElementKind::Quad4,
                    nodes.iter().map(|&n| node_row(n)).collect(),
                ),
                EmbeddedBarHost::Rod { .. } => {
                    unreachable!("embed_polyline binds every segment to a host")
                }
            };
            reported.push(EmbeddedBarSegment {
                start,
                end,
                length: segment.length(),
                host_kind,
                host_nodes,
            });
            let element = if rebuilt {
                rods[index].clone()
            } else {
                segment
            };
            all_segments.push(domain.add_element(Element::EmbeddedBar(element)));
        }
        let (routing, sr) = match material {
            Material::CsfmRebar(fields) => (
                Some(match fields.routing {
                    CsfmRebarRouting::Tcm => RebarRouting::Tcm,
                    CsfmRebarRouting::Pom => RebarRouting::Pom,
                }),
                Some(fields.sr),
            ),
            _ => (None, None),
        };
        groups.push(EmbeddedBarGroup {
            first_segment,
            segments: reported,
            routing,
            sr,
            bond: bond_group,
        });
    }
    elements.set(ElementKind::EmbeddedBar, all_segments);
    Ok(ModelOutput {
        embedded_bars: groups,
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
                PlaneMaterialSpec::CsfmConcrete {
                    fck,
                    ecm,
                    law,
                    softening,
                    tensile_regularization,
                } => PlaneMaterial::csfm_concrete(
                    fck,
                    ecm,
                    match law {
                        CsfmConcreteLawSpec::ParabolaRectangle => {
                            CsfmConcreteLaw::ParabolaRectangle
                        }
                        CsfmConcreteLawSpec::ElasticIdeallyPlastic => {
                            CsfmConcreteLaw::ElasticIdeallyPlastic
                        }
                    },
                    softening,
                    tensile_regularization,
                ),
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
