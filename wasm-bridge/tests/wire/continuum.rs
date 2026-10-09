//! the `planeMaterials`, `triangles` and `quads` tables, body and
//! edge element loads, and Gauss-point recorders, decoded and run end to end
//! against the same models built natively in `carapace-core`.

use crate::common::empty_input;
use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, Element, ElementLoad, LoadSeries, Node, PlaneMaterial, Quad4, Quad4Formulation, Tri3,
};
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, GaussQuantity, IntegratorSpec, RecorderSpec, SequenceSpec,
    StageSpec,
};
use carapace_wasm::input_v1::tables::*;
use carapace_wasm::input_v1::{decode, CarapaceInputV1, DecodeError, Header, StepOutcome};

fn one_step(recorders: Vec<RecorderSpec>) -> SequenceSpec {
    SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 5,
            }),
            hold_patterns_after: vec![],
        }],
        recorders,
    }
}

fn sample(outcome: &StepOutcome, recorder: usize) -> f64 {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder)
        .and_then(|batch| batch.samples.last())
        .map(|&(_, value)| value)
        .unwrap_or_else(|| panic!("recorder {recorder} recorded nothing"))
}

/// Two quads side by side, the left one enhanced, plus a triangle capping the right edge:
///
/// ```text
/// 3 --- 4 --- 5
/// |  Q0 |  Q1 |      nodes 0..6 on y = 0 / 1 (row-major), triangle (1,2,4)? no:
/// 0 --- 1 --- 2      the triangle is (2, 5, 6) hanging off the right edge at node 6.
/// ```
fn panel() -> CarapaceInputV1 {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![
            0.0, 0.0, 2.0, 0.0, 4.0, 0.2, //
            0.0, 1.5, 2.1, 1.6, 4.0, 1.4, //
            5.0, 0.7,
        ],
        fixed: vec![0b011, 0, 0, 0b011, 0, 0, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.plane_materials = vec![
        PlaneMaterialSpec::Isotropic {
            e: 3000.0,
            nu: 0.25,
            state: PlaneStateSpec::PlaneStress,
        },
        PlaneMaterialSpec::Orthotropic {
            ex: 3000.0,
            ey: 1500.0,
            nu_xy: 0.2,
            g_xy: 700.0,
            angle: 0.4,
        },
    ];
    input.quads = QuadTable {
        node_ids: vec![0, 1, 4, 3, 1, 2, 5, 4],
        thickness: vec![0.5, 0.5],
        material: vec![0, 1],
        density: vec![0.0, 0.0],
        formulation: vec![Quad4FormulationSpec::Enhanced, Quad4FormulationSpec::Full],
    };
    input.triangles = TriangleTable {
        node_ids: vec![2, 6, 5],
        thickness: vec![0.5],
        material: vec![0],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Linear { slope: 1.0 },
            TimeSeriesSpec::Constant,
        ],
        scale_factor: vec![1.0, 0.5],
    };
    input
}

/// Loads on the panel: a body force and an edge pressure through pattern 0, a second body
/// force through pattern 1 (scale 0.5), and an edge traction on the triangle.
fn panel_loads(input: &mut CarapaceInputV1) {
    input.element_loads = ElementLoadTable {
        pattern: vec![0, 0, 1, 0],
        element_kind: vec![
            ElementKind::Quad4,
            ElementKind::Quad4,
            ElementKind::Quad4,
            ElementKind::Tri3,
        ],
        element_index: vec![0, 1, 0, 0],
        load: vec![
            ElementLoadSpec::Body { bx: 0.5, by: -2.0 },
            ElementLoadSpec::EdgePressure {
                edge: 1,
                pressure: 3.0,
            },
            ElementLoadSpec::Body { bx: 1.0, by: 4.0 },
            ElementLoadSpec::EdgeTraction {
                edge: 0,
                tx: 1.0,
                ty: -1.5,
            },
        ],
        stage: vec![0, 0, 0, 0],
    };
}

/// The same panel built natively: returns node displacements (7 nodes x 2) and the Gauss
/// strain/stress of every element.
fn native_panel() -> (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut domain = Domain::new();
    let coords = [
        [0.0, 0.0],
        [2.0, 0.0],
        [4.0, 0.2],
        [0.0, 1.5],
        [2.1, 1.6],
        [4.0, 1.4],
        [5.0, 0.7],
    ];
    let ids: Vec<_> = coords
        .iter()
        .enumerate()
        .map(|(k, &c)| {
            let mut node = Node::new(c);
            if k == 0 || k == 3 {
                node = node.fix(0).fix(1);
            }
            domain.add_node(node)
        })
        .collect();
    let m0 = PlaneMaterial::plane_stress(3000.0, 0.25).unwrap();
    let m1 = PlaneMaterial::orthotropic(3000.0, 1500.0, 0.2, 700.0, 0.4).unwrap();
    let q0 = domain.add_element(Element::Quad4(
        Quad4::new([ids[0], ids[1], ids[4], ids[3]], 0.5, m0.clone())
            .with_formulation(Quad4Formulation::Enhanced),
    ));
    let q1 = domain.add_element(Element::Quad4(Quad4::new(
        [ids[1], ids[2], ids[5], ids[4]],
        0.5,
        m1,
    )));
    let t0 = domain.add_element(Element::Tri3(Tri3::new(ids[2], ids[6], ids[5], 0.5, m0)));
    let p0 = domain.add_load_pattern_scaled(LoadSeries::Linear { slope: 1.0 }, 1.0);
    let p1 = domain.add_load_pattern_scaled(LoadSeries::Constant, 0.5);
    domain.add_element_load(p0, q0, ElementLoad::body(0.5, -2.0));
    domain.add_element_load(p0, q1, ElementLoad::edge_pressure(1, 3.0));
    domain.add_element_load(p1, q0, ElementLoad::body(1.0, 4.0));
    domain.add_element_load(p0, t0, ElementLoad::edge_traction(0, 1.0, -1.5));
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis.step().unwrap();
    let d = analysis.domain();
    let displacements = ids
        .iter()
        .flat_map(|&id| [d.node(id).displacement[0], d.node(id).displacement[1]])
        .collect();
    let mut strains = Vec::new();
    let mut stresses = Vec::new();
    for id in [q0, q1, t0] {
        for point in d.element_gauss_responses(id).unwrap() {
            strains.push(point.strain);
            stresses.push(point.stress);
        }
    }
    (displacements, strains, stresses)
}

/// A small panel (two Quad4 and a Tri3) decoded from the wire format and run through a session
/// gives the same nodal displacements and Gauss-point strains and stresses as the same model built
/// directly with `carapace-core`.
#[test]
fn a_panel_runs_end_to_end_identically_to_the_same_model_built_natively() {
    let mut input = panel();
    panel_loads(&mut input);
    let mut recorders = Vec::new();
    for node in 0..7 {
        for dof in 0..2 {
            recorders.push(RecorderSpec::NodeDisp { node, dof });
        }
    }
    // Gauss points: quad 0 (4), quad 1 (4), triangle (1); strain then stress, 3 components each.
    let points = [
        (ElementKind::Quad4, 0, 4),
        (ElementKind::Quad4, 1, 4),
        (ElementKind::Tri3, 0, 1),
    ];
    for quantity in [GaussQuantity::Strain, GaussQuantity::Stress] {
        for (element_kind, element_index, count) in points {
            for point in 0..count {
                for component in 0..3 {
                    recorders.push(RecorderSpec::GaussPoint {
                        element_kind,
                        element_index,
                        point,
                        quantity,
                        component,
                    });
                }
            }
        }
    }
    input.sequence = one_step(recorders);
    let mut session = decode(input).expect("panel decodes");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");

    let (displacements, strains, stresses) = native_panel();
    for (k, expected) in displacements.iter().enumerate() {
        assert_eq!(sample(&outcome, k), *expected, "displacement {k}");
    }
    assert!(displacements.iter().any(|d| d.abs() > 1e-6));
    let base = 14;
    let per_quantity = strains.len() * 3;
    for (g, (strain, stress)) in strains.iter().zip(&stresses).enumerate() {
        for c in 0..3 {
            assert_eq!(
                sample(&outcome, base + 3 * g + c),
                strain[c],
                "strain {g}/{c}"
            );
            assert_eq!(
                sample(&outcome, base + per_quantity + 3 * g + c),
                stress[c],
                "stress {g}/{c}"
            );
        }
    }
}

/// Wire-format version of the body-force summation: several body loads on one Quad4 across patterns
/// with different series and scale factors sum to the same displacement as one pattern with the
/// total load.
#[test]
fn several_body_force_patterns_through_the_wire_sum_with_their_factors() {
    let run = |bodies: &[(u32, f64, f64)], scale_factor: Vec<f64>| {
        let mut input = panel();
        input.load_patterns.scale_factor = scale_factor;
        input.element_loads = ElementLoadTable {
            pattern: bodies.iter().map(|b| b.0).collect(),
            element_kind: vec![ElementKind::Quad4; bodies.len()],
            element_index: vec![0; bodies.len()],
            load: bodies
                .iter()
                .map(|b| ElementLoadSpec::Body { bx: b.1, by: b.2 })
                .collect(),
            stage: vec![0; bodies.len()],
        };
        input.sequence = one_step(vec![RecorderSpec::NodeDisp { node: 4, dof: 1 }]);
        let mut session = decode(input).unwrap();
        sample(&session.advance(10), 0)
    };
    // Pattern 0 is Linear (factor 1 at the end), pattern 1 Constant: totals (1 + 0.5 * 1 * ...)
    // two rows on pattern 0 add, and pattern 1 enters with its scale factor 2.
    let split = run(
        &[(0, 0.0, -1.0), (0, 0.0, -2.0), (1, 0.0, -0.5)],
        vec![1.0, 2.0],
    );
    let single = run(&[(0, 0.0, -4.0)], vec![1.0, 2.0]);
    assert!(
        (split - single).abs() < 1e-14 * single.abs(),
        "{split} vs {single}"
    );
    assert!(single != 0.0);
}

fn decode_error(input: CarapaceInputV1) -> DecodeError {
    decode(input).err().expect("expected a decode error")
}

/// The continuum tables, loads and Gauss-point recorders survive a serialize/deserialize round trip
/// unchanged, with the expected camelCase JSON shape (`planeStress`, `edgePressure`, `gaussPoint`),
/// and still decode.
#[test]
fn the_tables_round_trip_through_serde() {
    let mut input = panel();
    panel_loads(&mut input);
    input.sequence = one_step(vec![RecorderSpec::GaussPoint {
        element_kind: ElementKind::Quad4,
        element_index: 1,
        point: 3,
        quantity: GaussQuantity::Stress,
        component: 2,
    }]);
    let value = serde_json::to_value(&input).unwrap();
    let back: CarapaceInputV1 = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(serde_json::to_value(&back).unwrap(), value);
    assert_eq!(value["quads"]["formulation"][0], "enhanced");
    assert_eq!(value["planeMaterials"][0]["state"], "planeStress");
    assert_eq!(value["planeMaterials"][1]["kind"], "orthotropic");
    assert_eq!(value["elementLoads"]["load"][1]["kind"], "edgePressure");
    assert_eq!(value["sequence"]["recorders"][0]["response"], "gaussPoint");
    assert!(decode(back).is_ok());
}

/// Malformed continuum tables are structured decode errors: unknown material or node indices, rows
/// with the wrong number of nodes, and other shape problems, each reported with the table and row.
#[test]
fn bad_references_and_shapes_are_decode_errors() {
    let mut input = panel();
    input.quads.material[1] = 7;
    assert_eq!(
        decode_error(input),
        DecodeError::UnknownMaterialIndex {
            table: "quads",
            row: 7
        }
    );
    let mut input = panel();
    input.triangles.node_ids[1] = 99;
    assert_eq!(
        decode_error(input),
        DecodeError::UnknownNodeIndex {
            table: "triangles",
            row: 99
        }
    );
    let mut input = panel();
    input.quads.node_ids.pop();
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidRow { table: "quads", .. }
    ));
    let mut input = panel();
    input.plane_materials[0] = PlaneMaterialSpec::Isotropic {
        e: -1.0,
        nu: 0.3,
        state: PlaneStateSpec::PlaneStress,
    };
    assert_eq!(
        decode_error(input),
        DecodeError::InvalidPlaneMaterial {
            index: 0,
            reason: "a modulus is not finite and positive"
        }
    );
    let mut input = panel();
    input.plane_materials[0] = PlaneMaterialSpec::Isotropic {
        e: 1.0,
        nu: 0.5,
        state: PlaneStateSpec::PlaneStrain,
    };
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidPlaneMaterial { index: 0, .. }
    ));
    // A clockwise quad is a model error from core's validation.
    let mut input = panel();
    input.quads.node_ids = vec![0, 3, 4, 1, 1, 2, 5, 4];
    assert!(matches!(
        decode_error(input),
        DecodeError::InvalidModel { .. }
    ));
}

/// Continuum tables (`quads`, `triangles`, `planeMaterials`) and continuum recorders are rejected
/// in a 3D input with `TableNotInProfile`, since those elements are 2D only.
#[test]
fn continuum_tables_and_kinds_are_not_in_a_3d_model() {
    let mut input = empty_input(3);
    input.quads = panel().quads;
    assert_eq!(
        decode_error(input),
        DecodeError::TableNotInProfile { table: "quads" }
    );
    let mut input = empty_input(3);
    input.triangles = panel().triangles;
    assert_eq!(
        decode_error(input),
        DecodeError::TableNotInProfile { table: "triangles" }
    );
    let mut input = empty_input(3);
    input.plane_materials = panel().plane_materials;
    assert_eq!(
        decode_error(input),
        DecodeError::TableNotInProfile {
            table: "plane_materials"
        }
    );
    let mut input = empty_input(3);
    input.sequence = one_step(vec![RecorderSpec::GaussPoint {
        element_kind: ElementKind::Quad4,
        element_index: 0,
        point: 0,
        quantity: GaussQuantity::Stress,
        component: 0,
    }]);
    assert_eq!(
        decode_error(input),
        DecodeError::ElementKindNotInProfile { table: "quads" }
    );
}

/// Element loads of the wrong kind for their element (a beam load on a continuum element, continuum
/// loads on a beam) and out-of-range edge indices are decode errors, and nothing is accumulated.
#[test]
fn wrong_kind_loads_and_out_of_range_edges_are_rejected_without_accumulating() {
    let load_on = |kind, load| {
        let mut input = panel();
        input.elastic_beam_columns_2d = ElasticBeamColumn2dTable {
            node_i: vec![0],
            node_j: vec![1],
            e: vec![1.0],
            a: vec![1.0],
            iz: vec![1.0],
            transform: vec![TransformSpec::Linear],
            density: vec![0.0],
        };
        input.element_loads = ElementLoadTable {
            pattern: vec![0],
            element_kind: vec![kind],
            element_index: vec![0],
            load: vec![load],
            stage: vec![0],
        };
        input.sequence = one_step(vec![]);
        decode_error(input)
    };
    let unsupported = |kind: ElementKind| DecodeError::UnsupportedElementLoad {
        element_kind: kind.table(),
    };
    // A beam load on a continuum element, continuum loads on a beam.
    let uniform = ElementLoadSpec::Uniform {
        wx: 1.0,
        wy: 0.0,
        wz: None,
    };
    assert_eq!(
        load_on(ElementKind::Quad4, uniform),
        unsupported(ElementKind::Quad4)
    );
    assert_eq!(
        load_on(ElementKind::Tri3, uniform),
        unsupported(ElementKind::Tri3)
    );
    let body = ElementLoadSpec::Body { bx: 0.0, by: 1.0 };
    assert_eq!(
        load_on(ElementKind::ElasticBeamColumn2d, body),
        unsupported(ElementKind::ElasticBeamColumn2d)
    );
    let pressure = |edge| ElementLoadSpec::EdgePressure {
        edge,
        pressure: 1.0,
    };
    assert_eq!(
        load_on(ElementKind::ElasticBeamColumn2d, pressure(0)),
        unsupported(ElementKind::ElasticBeamColumn2d)
    );
    // Edge indices: a triangle has three, a quad four.
    assert!(matches!(
        load_on(ElementKind::Tri3, pressure(3)),
        DecodeError::InvalidRow {
            table: "element_loads",
            ..
        }
    ));
    assert!(matches!(
        load_on(ElementKind::Quad4, pressure(4)),
        DecodeError::InvalidRow {
            table: "element_loads",
            ..
        }
    ));
}

/// Gauss-point recorders are checked against the element: Quad4 has 4 points and Tri3 has 1, and an
/// out-of-range point or component is `InvalidGaussPoint` or an invalid-component error rather than
/// a silent empty recording.
#[test]
fn gauss_point_recorders_are_bounds_checked_against_the_element() {
    let with = |element_kind, element_index, point, component| {
        let mut input = panel();
        input.sequence = one_step(vec![RecorderSpec::GaussPoint {
            element_kind,
            element_index,
            point,
            quantity: GaussQuantity::Strain,
            component,
        }]);
        decode(input).err()
    };
    assert_eq!(with(ElementKind::Quad4, 0, 3, 2), None);
    assert_eq!(with(ElementKind::Tri3, 0, 0, 2), None);
    assert_eq!(
        with(ElementKind::Quad4, 0, 4, 0),
        Some(DecodeError::InvalidGaussPoint {
            recorder: 0,
            point: 4,
            count: 4
        })
    );
    assert_eq!(
        with(ElementKind::Tri3, 0, 1, 0),
        Some(DecodeError::InvalidGaussPoint {
            recorder: 0,
            point: 1,
            count: 1
        })
    );
    assert_eq!(
        with(ElementKind::Quad4, 1, 0, 3),
        Some(DecodeError::InvalidRecorderComponent {
            recorder: 0,
            component: 3,
            width: 3
        })
    );
    assert_eq!(
        with(ElementKind::Quad4, 5, 0, 0),
        Some(DecodeError::UnknownElementIndex {
            table: "quads",
            row: 5
        })
    );
    // A kind that is not a continuum element has no Gauss points.
    let mut input = panel();
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.0],
        material: vec![0],
        density: vec![0.0],
    };
    input.sequence = one_step(vec![RecorderSpec::GaussPoint {
        element_kind: ElementKind::Truss,
        element_index: 0,
        point: 0,
        quantity: GaussQuantity::Strain,
        component: 0,
    }]);
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::InvalidGaussPoint {
            recorder: 0,
            point: 0,
            count: 0
        })
    );
}

#[allow(dead_code)]
fn _header_is_public(header: Header) -> u8 {
    header.ndm
}
