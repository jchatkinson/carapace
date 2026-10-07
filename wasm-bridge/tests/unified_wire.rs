//! Phase 2 (docs/domain-refactor-plan.md §2.5): the unified `CarapaceInputV1`.
//! Profile rejection, omitted tables, serde round trip, `ndm` validation,
//! the `rigidLinks`/`linearConstraints` tables against the same models built
//! natively in `core`, and recorder component bounds.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    Domain, Domain3, Element, Element3, Material, Node, Node3, ZeroLength, ZeroLength3,
};
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::*;
use carapace_wasm::input_v1::{decode, CarapaceInputV1, DecodeError, Header, StepOutcome};

fn empty_input(ndm: u8) -> CarapaceInputV1 {
    CarapaceInputV1 {
        header: Header {
            schema_version: 1,
            ndm,
            engine_version: "test".to_string(),
            record_initial: false,
        },
        nodes: Default::default(),
        materials: Vec::new(),
        fibers: Default::default(),
        trusses: Default::default(),
        elastic_beam_columns_2d: Default::default(),
        elastic_beam_columns_3d: Default::default(),
        disp_beam_columns_2d: Default::default(),
        disp_beam_columns_3d: Default::default(),
        force_beam_columns_2d: Default::default(),
        force_beam_columns_3d: Default::default(),
        zero_lengths: Default::default(),
        zero_length_sections: Default::default(),
        equal_dofs: Default::default(),
        rigid_diaphragms: Default::default(),
        rigid_links: Default::default(),
        linear_constraints: Default::default(),
        load_patterns: Default::default(),
        nodal_loads: Default::default(),
        element_loads: Default::default(),
        sequence: Default::default(),
    }
}

fn last_sample(outcome: &StepOutcome, recorder_index: usize) -> Option<f64> {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder_index)
        .and_then(|batch| batch.samples.last())
        .map(|&(_, value)| value)
}

/// One linear `LoadControl` step of a unit-slope pattern.
fn one_linear_step(recorders: Vec<RecorderSpec>) -> SequenceSpec {
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

fn unit_pattern(input: &mut CarapaceInputV1) {
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
}

fn nodal_loads(input: &mut CarapaceInputV1, node: u32, values: &[f64]) {
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0; values.len()],
        node: vec![node; values.len()],
        dof: (0..values.len() as u8).collect(),
        value: values.to_vec(),
        stage: vec![0; values.len()],
    };
}

fn run(input: CarapaceInputV1) -> StepOutcome {
    let mut session = decode(input).expect("model decodes");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    outcome
}

// ---------------------------------------------------------------------------
// Profile rejection, omitted tables, serde, ndm.
// ---------------------------------------------------------------------------

fn beam_tables_3d() -> (ElasticBeamColumn3dTable, FiberBeamColumn3dTable) {
    (
        ElasticBeamColumn3dTable {
            node_i: vec![0],
            ..Default::default()
        },
        FiberBeamColumn3dTable {
            node_i: vec![0],
            ..Default::default()
        },
    )
}

#[test]
fn a_table_that_belongs_to_the_other_profile_is_rejected() {
    // 3D tables in a 2D model.
    let (elastic, fiber) = beam_tables_3d();
    let mut input = empty_input(2);
    input.elastic_beam_columns_3d = elastic.clone();
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "elastic_beam_columns_3d"
        })
    );
    let mut input = empty_input(2);
    input.disp_beam_columns_3d = fiber.clone();
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "disp_beam_columns_3d"
        })
    );
    let mut input = empty_input(2);
    input.force_beam_columns_3d = fiber;
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "force_beam_columns_3d"
        })
    );

    // 2D tables in a 3D model.
    let mut input = empty_input(3);
    input.elastic_beam_columns_2d.node_i = vec![0];
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "elastic_beam_columns_2d"
        })
    );
    let mut input = empty_input(3);
    input.disp_beam_columns_2d.node_i = vec![0];
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "disp_beam_columns_2d"
        })
    );
    let mut input = empty_input(3);
    input.force_beam_columns_2d.node_i = vec![0];
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::TableNotInProfile {
            table: "force_beam_columns_2d"
        })
    );
}

#[test]
fn an_element_kind_of_the_other_profile_is_rejected_in_loads_and_recorders() {
    let mut input = empty_input(2);
    unit_pattern(&mut input);
    input.sequence = one_linear_step(vec![]);
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ElasticBeamColumn3d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy: 1.0,
            wz: None,
        }],
        stage: vec![0],
    };
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::ElementKindNotInProfile {
            table: "elastic_beam_columns_3d"
        })
    );

    let mut input = empty_input(3);
    input.sequence = one_linear_step(vec![RecorderSpec::ElementForce {
        element_kind: ElementKind::ForceBeamColumn2d,
        element_index: 0,
        component: 0,
    }]);
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::ElementKindNotInProfile {
            table: "force_beam_columns_2d"
        })
    );
}

#[test]
fn unknown_ndm_is_rejected() {
    for got in [0, 1, 4, 255] {
        assert_eq!(
            decode(empty_input(got)).err(),
            Some(DecodeError::UnsupportedNdm { got })
        );
    }
}

#[test]
fn omitted_tables_decode_as_empty() {
    // Only the header: every table (and the sequence) is omitted.
    for ndm in [2, 3] {
        let json = serde_json::json!({
            "header": { "schemaVersion": 1, "ndm": ndm, "engineVersion": "t" }
        });
        let input: CarapaceInputV1 = serde_json::from_value(json).expect("header-only input");
        let mut session = decode(input).expect("empty model decodes");
        assert!(session.advance(1).done);
    }

    // A truss carried by the shared tables alone, with every other table omitted.
    let json = serde_json::json!({
        "header": { "schemaVersion": 1, "ndm": 2, "engineVersion": "t" },
        "nodes": { "coords": [0.0, 0.0, 100.0, 0.0], "fixed": [7, 6], "massNodeIndex": [], "mass": [] },
        "materials": [{ "kind": "elastic", "e": 1000.0 }],
        "trusses": { "nodeI": [0], "nodeJ": [1], "area": [2.0], "material": [0], "density": [0.0] },
        "loadPatterns": { "series": [{ "kind": "linear", "slope": 1.0 }], "scaleFactor": [1.0] },
        "nodalLoads": { "pattern": [0], "node": [1], "dof": [0], "value": [50.0], "stage": [0] },
        "sequence": {
            "stages": [{
                "kind": "static", "id": "only", "steps": 1,
                "integrator": { "kind": "loadControl", "increment": 1.0 },
                "algorithm": "linear", "holdPatternsAfter": []
            }],
            "recorders": [{ "response": "nodeDisp", "node": 1, "dof": 0 }]
        }
    });
    let input: CarapaceInputV1 = serde_json::from_value(json).expect("sparse input");
    let outcome = run(input);
    let u = last_sample(&outcome, 0).expect("sample");
    assert!((u - 50.0 * 100.0 / (2.0 * 1000.0)).abs() < 1e-12, "u = {u}");
}

/// A 3D model touching the tables with the new unified shapes.
fn populated_3d() -> CarapaceInputV1 {
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0],
        fixed: vec![0b111111, 0, 0],
        mass_node_index: vec![1],
        mass: vec![1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 10.0 }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 2],
        y: vec![1.0, -1.0],
        z: vec![0.5, -0.5],
        area: vec![1.0, 1.0],
        material: vec![0, 0],
    };
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![FrictionRow {
            row: 0,
            normal_dof: 0,
            shear_dofs: vec![1, 2],
            mu: 0.3,
            k0: 5.0,
            b: 0.01,
        }],
        orient: vec![OrientRow {
            row: 0,
            x: [1.0, 0.0, 0.0],
            yp: Some([0.0, 1.0, 0.0]),
        }],
    };
    input.rigid_diaphragms = RigidDiaphragmTable {
        retained: vec![1],
        normal: vec![Axis3Spec::Z],
        constrained: vec![(0, 2)],
    };
    input.rigid_links = RigidLinkTable {
        master: vec![0],
        slave: vec![1],
    };
    input.linear_constraints = LinearConstraintTable {
        slave_node: vec![2],
        slave_dof: vec![2],
        term_offsets: vec![0, 2],
        term_node: vec![0, 1],
        term_dof: vec![2, 2],
        term_coeff: vec![0.5, 0.5],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ElasticBeamColumn3d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy: 1.0,
            wz: Some(2.0),
        }],
        stage: vec![0],
    };
    input.sequence = one_linear_step(vec![RecorderSpec::ElementForce {
        element_kind: ElementKind::ElasticBeamColumn3d,
        element_index: 0,
        component: 11,
    }]);
    input
}

#[test]
fn the_unified_input_round_trips_through_serde() {
    let mut two_d = empty_input(2);
    two_d.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![7, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    two_d.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![FrictionRow {
            row: 0,
            normal_dof: 0,
            shear_dofs: vec![1],
            mu: 0.3,
            k0: 5.0,
            b: 0.01,
        }],
        orient: vec![OrientRow {
            row: 0,
            x: [0.0, 1.0, 0.0],
            yp: None,
        }],
    };
    two_d.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::DispBeamColumn2d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 1.0,
            wy: 2.0,
            wz: None,
        }],
        stage: vec![0],
    };
    for input in [two_d, populated_3d()] {
        let value = serde_json::to_value(&input).expect("serializes");
        let back: CarapaceInputV1 = serde_json::from_value(value.clone()).expect("deserializes");
        assert_eq!(
            serde_json::to_value(&back).expect("serializes again"),
            value
        );
    }
    // The wire names the plan fixes.
    let value = serde_json::to_value(populated_3d()).unwrap();
    assert_eq!(value["header"]["ndm"], 3);
    for key in [
        "elasticBeamColumns3d",
        "elasticBeamColumns2d",
        "dispBeamColumns2d",
        "forceBeamColumns3d",
        "rigidLinks",
        "linearConstraints",
    ] {
        assert!(value.get(key).is_some(), "missing table {key}");
    }
    assert_eq!(value["zeroLengths"]["friction"][0]["shearDofs"][1], 2);
    assert_eq!(value["zeroLengths"]["orient"][0]["yp"][1], 1.0);
    assert_eq!(value["linearConstraints"]["termOffsets"][1], 2);
}

// ---------------------------------------------------------------------------
// rigidLinks / linearConstraints against the same model built natively.
// ---------------------------------------------------------------------------

const K2: [f64; 3] = [400.0, 250.0, 1000.0];
const LOAD2: [f64; 3] = [10.0, -6.0, 20.0];
const D2: (f64, f64) = (3.0, 4.0);

/// 2D: master node 0 (loaded), slave node 1 at offset `D2`, fixed ground node 2
/// with `K2` springs to the slave.
fn rigid_link_input_2d() -> CarapaceInputV1 {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, D2.0, D2.1, 0.0, 0.0],
        fixed: vec![0, 0, 0b111],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = K2.iter().map(|&e| MaterialSpec::Elastic { e }).collect();
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![2],
        node_j: vec![1],
        materials: vec![(0, 0, 0), (0, 1, 1), (0, 2, 2)],
        friction: vec![],
        orient: vec![],
    };
    unit_pattern(&mut input);
    nodal_loads(&mut input, 0, &LOAD2);
    input.sequence = one_linear_step(vec![
        RecorderSpec::NodeDisp { node: 0, dof: 0 },
        RecorderSpec::NodeDisp { node: 0, dof: 1 },
        RecorderSpec::NodeDisp { node: 0, dof: 2 },
        RecorderSpec::NodeDisp { node: 1, dof: 0 },
        RecorderSpec::NodeDisp { node: 1, dof: 1 },
        RecorderSpec::NodeDisp { node: 1, dof: 2 },
    ]);
    input
}

fn native_2d(
    link: impl FnOnce(&mut Domain, carapace_core::model::NodeId, carapace_core::model::NodeId),
) -> [f64; 6] {
    let mut domain: Domain = Domain::new();
    let master = domain.add_node(Node::new([0.0, 0.0]));
    let slave = domain.add_node(Node::new([D2.0, D2.1]));
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, slave)
            .with_material(0, Material::Elastic { e: K2[0] })
            .with_material(1, Material::Elastic { e: K2[1] })
            .with_material(2, Material::Elastic { e: K2[2] }),
    ));
    link(&mut domain, master, slave);
    for (dof, value) in LOAD2.into_iter().enumerate() {
        domain.load_node(master, dof, value);
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis.step().expect("native model solves");
    let (m, s) = (
        analysis.domain().node(master).displacement,
        analysis.domain().node(slave).displacement,
    );
    [m[0], m[1], m[2], s[0], s[1], s[2]]
}

fn assert_identical(outcome: &StepOutcome, native: [f64; 6]) {
    for (i, expected) in native.into_iter().enumerate() {
        let got = last_sample(outcome, i).expect("recorded");
        assert_eq!(
            got, expected,
            "recorder {i}: wire {got} vs native {expected}"
        );
    }
}

#[test]
fn a_rigid_links_table_matches_the_native_2d_rigid_link() {
    let mut input = rigid_link_input_2d();
    input.rigid_links = RigidLinkTable {
        master: vec![0],
        slave: vec![1],
    };
    let outcome = run(input);
    assert_identical(&outcome, native_2d(|d, m, s| d.rigid_link(m, s)));
    // The slave really moved with the master (not a vacuous comparison).
    let (m2, s0) = (
        last_sample(&outcome, 2).unwrap(),
        last_sample(&outcome, 3).unwrap(),
    );
    assert!(m2 != 0.0 && (s0 - (last_sample(&outcome, 0).unwrap() - m2 * D2.1)).abs() < 1e-12);
}

#[test]
fn a_linear_constraints_table_matches_the_native_add_constraint() {
    // The 2D rigid link written out as three general constraints.
    let (dx, dy) = D2;
    let mut input = rigid_link_input_2d();
    input.linear_constraints = LinearConstraintTable {
        slave_node: vec![1, 1, 1],
        slave_dof: vec![0, 1, 2],
        term_offsets: vec![0, 2, 4, 5],
        term_node: vec![0, 0, 0, 0, 0],
        term_dof: vec![0, 2, 1, 2, 2],
        term_coeff: vec![1.0, -dy, 1.0, dx, 1.0],
    };
    let outcome = run(input);
    let native = native_2d(|d, m, s| {
        d.add_constraint((s, 0), &[(m, 0, 1.0), (m, 2, -dy)]);
        d.add_constraint((s, 1), &[(m, 1, 1.0), (m, 2, dx)]);
        d.add_constraint((s, 2), &[(m, 2, 1.0)]);
    });
    assert_identical(&outcome, native);
    // And it is the rigid link, to rounding.
    let link = native_2d(|d, m, s| d.rigid_link(m, s));
    for (a, b) in native.iter().zip(link) {
        assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0));
    }
}

#[test]
fn a_rigid_links_table_matches_the_native_3d_rigid_link() {
    let k = [300.0, 200.0, 100.0, 500.0, 600.0, 700.0];
    let d = [2.0, -1.0, 3.0];
    let load = [5.0, -3.0, 8.0, 1.0, 2.0, -4.0];

    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, d[0], d[1], d[2], 0.0, 0.0, 0.0],
        fixed: vec![0, 0, 0b111111],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = k.iter().map(|&e| MaterialSpec::Elastic { e }).collect();
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![2],
        node_j: vec![1],
        materials: (0..6).map(|dof| (0, dof as u8, dof)).collect(),
        friction: vec![],
        orient: vec![],
    };
    input.rigid_links = RigidLinkTable {
        master: vec![0],
        slave: vec![1],
    };
    unit_pattern(&mut input);
    nodal_loads(&mut input, 0, &load);
    input.sequence = one_linear_step(
        (0..6)
            .map(|dof| RecorderSpec::NodeDisp { node: 0, dof })
            .chain((0..6).map(|dof| RecorderSpec::NodeDisp { node: 1, dof }))
            .collect(),
    );
    let outcome = run(input);

    let mut domain = Domain3::new();
    let master = domain.add_node(Node3::new([0.0; 3]));
    let slave = domain.add_node(Node3::new(d));
    let ground = domain.add_node(
        Node3::new([0.0; 3])
            .fix(0)
            .fix(1)
            .fix(2)
            .fix(3)
            .fix(4)
            .fix(5),
    );
    let mut spring = ZeroLength3::new(ground, slave);
    for (dof, &e) in k.iter().enumerate() {
        spring = spring.with_material(dof, Material::Elastic { e });
    }
    domain.add_element(Element3::ZeroLength3(spring));
    domain.rigid_link(master, slave);
    for (dof, value) in load.into_iter().enumerate() {
        domain.load_node(master, dof, value);
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis.step().expect("native model solves");
    for dof in 0..6 {
        let m = analysis.domain().node(master).displacement[dof];
        let s = analysis.domain().node(slave).displacement[dof];
        assert_eq!(last_sample(&outcome, dof).unwrap(), m, "master dof {dof}");
        assert_eq!(
            last_sample(&outcome, 6 + dof).unwrap(),
            s,
            "slave dof {dof}"
        );
    }
}

// ---------------------------------------------------------------------------
// Recorder component bounds.
// ---------------------------------------------------------------------------

fn truss_input(ndm: u8) -> CarapaceInputV1 {
    let mut input = empty_input(ndm);
    let (coords, fixed) = if ndm == 2 {
        (vec![0.0, 0.0, 1.0, 0.0], vec![0b111, 0b110])
    } else {
        (vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0], vec![0b111111, 0b111110])
    };
    input.nodes = NodeTable {
        coords,
        fixed,
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 10.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.0],
        material: vec![0],
        density: vec![0.0],
    };
    input
}

#[test]
fn an_out_of_range_recorder_component_is_a_decode_error() {
    for (ndm, force_width, load_width) in [(2u8, 6u8, 16u8), (3, 12, 3)] {
        let element_kind = ElementKind::Truss;
        // The last valid component decodes; one past it is `InvalidRecorderComponent`.
        for (kind_is_load, width) in [(false, force_width), (true, load_width)] {
            let recorder = |component| {
                if kind_is_load {
                    RecorderSpec::ElementLoad {
                        element_kind,
                        element_index: 0,
                        component,
                    }
                } else {
                    RecorderSpec::ElementForce {
                        element_kind,
                        element_index: 0,
                        component,
                    }
                }
            };
            let mut input = truss_input(ndm);
            input.sequence = one_linear_step(vec![recorder(width - 1)]);
            assert!(decode(input).is_ok(), "ndm {ndm} width {width}");

            let mut input = truss_input(ndm);
            input.sequence = one_linear_step(vec![
                RecorderSpec::NodeDisp { node: 0, dof: 0 },
                recorder(width),
            ]);
            assert_eq!(
                decode(input).err(),
                Some(DecodeError::InvalidRecorderComponent {
                    recorder: 1,
                    component: width,
                    width: width as u32,
                }),
                "ndm {ndm} load {kind_is_load}"
            );
        }
    }
}

#[test]
fn recorder_dofs_are_bounded_by_the_profile() {
    let mut input = truss_input(2);
    input.sequence = one_linear_step(vec![RecorderSpec::NodeDisp { node: 1, dof: 3 }]);
    assert_eq!(
        decode(input).err(),
        Some(DecodeError::InvalidDof {
            table: "recorders",
            row: 0,
            dof: 3
        })
    );
    let mut input = truss_input(3);
    input.sequence = one_linear_step(vec![RecorderSpec::Reaction { node: 1, dof: 5 }]);
    assert!(decode(input).is_ok());
}

// ---------------------------------------------------------------------------
// Row shapes that depend on the profile.
// ---------------------------------------------------------------------------

#[test]
fn profile_dependent_row_shapes_are_checked() {
    let invalid = |table, reason| DecodeError::InvalidRow {
        table,
        row: 0,
        reason,
    };
    let friction = |shear_dofs: Vec<u8>| ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![FrictionRow {
            row: 0,
            normal_dof: 0,
            shear_dofs,
            mu: 0.3,
            k0: 5.0,
            b: 0.01,
        }],
        orient: vec![],
    };

    // Friction: one shear DOF in 2D, two in 3D.
    let mut input = truss_input(2);
    input.zero_lengths = friction(vec![1, 2]);
    assert_eq!(
        decode(input).err(),
        Some(invalid(
            "zero_lengths",
            "friction needs 1 shear DOF in a 2D model"
        ))
    );
    let mut input = truss_input(3);
    input.zero_lengths = friction(vec![1]);
    assert_eq!(
        decode(input).err(),
        Some(invalid(
            "zero_lengths",
            "friction needs 2 shear DOFs in a 3D model"
        ))
    );

    // Orientation: `yp` only in 3D.
    let orient = |yp| ZeroLengthSectionTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
        orient: vec![OrientRow {
            row: 0,
            x: [1.0, 0.0, 0.0],
            yp,
        }],
    };
    for (ndm, yp, reason) in [
        (2u8, Some([0.0, 1.0, 0.0]), "yp is only valid in a 3D model"),
        (3, None, "yp is required in a 3D model"),
    ] {
        let mut input = truss_input(ndm);
        input.fibers = FiberTable {
            section_offsets: vec![0, 1],
            y: vec![0.0],
            z: if ndm == 3 { vec![0.0] } else { vec![] },
            area: vec![1.0],
            material: vec![0],
        };
        input.zero_length_sections = orient(yp);
        assert_eq!(
            decode(input).err(),
            Some(invalid("zero_length_sections", reason))
        );
    }

    // Fibers: `z` empty in 2D, parallel to `y` in 3D.
    for (ndm, z) in [(2u8, vec![0.0]), (3, vec![])] {
        let mut input = truss_input(ndm);
        input.fibers = FiberTable {
            section_offsets: vec![0, 1],
            y: vec![0.0],
            z,
            area: vec![1.0],
            material: vec![0],
        };
        assert!(matches!(
            decode(input).err(),
            Some(DecodeError::InvalidRow {
                table: "fibers",
                ..
            })
        ));
    }

    // A diaphragm normal is 3D only; a 2D element load takes no `wz`.
    let mut input = truss_input(2);
    input.rigid_diaphragms = RigidDiaphragmTable {
        retained: vec![0],
        normal: vec![Axis3Spec::Y],
        constrained: vec![(0, 1)],
    };
    assert_eq!(
        decode(input).err(),
        Some(invalid(
            "rigid_diaphragms",
            "normal is only valid in a 3D model"
        ))
    );
    let mut input = truss_input(2);
    input.elastic_beam_columns_2d = ElasticBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![1.0],
        a: vec![1.0],
        iz: vec![1.0],
        transform: vec![TransformSpec::Linear],
        density: vec![0.0],
    };
    unit_pattern(&mut input);
    input.sequence = one_linear_step(vec![]);
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ElasticBeamColumn2d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy: 1.0,
            wz: Some(1.0),
        }],
        stage: vec![0],
    };
    assert_eq!(
        decode(input).err(),
        Some(invalid("element_loads", "wz is only valid in a 3D model"))
    );

    // Malformed linear-constraint offsets.
    let mut input = truss_input(2);
    input.linear_constraints = LinearConstraintTable {
        slave_node: vec![1],
        slave_dof: vec![0],
        term_offsets: vec![0, 2],
        term_node: vec![0],
        term_dof: vec![0],
        term_coeff: vec![1.0],
    };
    assert!(matches!(
        decode(input).err(),
        Some(DecodeError::InvalidRow {
            table: "linear_constraints",
            ..
        })
    ));
}
