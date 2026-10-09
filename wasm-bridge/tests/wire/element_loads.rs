//! Uniform element loads through the wire format, and rejection on unsupported elements.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::decode;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::*;

/// A cantilever `ElasticBeamColumn` under a local-axis uniform load (`wx` axial + `wy`
/// transverse) driven through the wire tables: tip displacements match the closed forms, and the
/// element-force recorders report the member end forces (statics: the fixed end carries the whole
/// load, the free end none), which only holds if recorders account for the element load.
#[test]
fn decodes_a_uniform_element_load_and_records_loaded_member_end_forces() {
    let (length, e, area, iz) = (100.0, 30_000.0, 10.0, 1000.0);
    let (wx, wy) = (2.0, -1.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.elastic_beam_columns_2d = ElasticBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![e],
        a: vec![area],
        iz: vec![iz],
        transform: vec![TransformSpec::Linear],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    // Two rows on the same element: they must sum to (wx, wy).
    input.element_loads = ElementLoadTable {
        pattern: vec![0, 0],
        element_kind: vec![ElementKind::ElasticBeamColumn2d; 2],
        element_index: vec![0, 0],
        load: vec![
            ElementLoadSpec::Uniform {
                wx,
                wy: 0.0,
                wz: None,
            },
            ElementLoadSpec::Uniform {
                wx: 0.0,
                wy,
                wz: None,
            },
        ],
        stage: vec![0, 0],
    };
    let element_force = |component| RecorderSpec::ElementForce {
        element_kind: ElementKind::ElasticBeamColumn2d,
        element_index: 0,
        component,
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
            element_force(0), // axial at i
            element_force(1), // shear at i
            element_force(2), // moment at i
            element_force(3), // axial at j
            element_force(4), // shear at j
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = [
        wx * length * length / (2.0 * e * area),
        wy * length.powi(4) / (8.0 * e * iz),
        -wx * length,
        -wy * length,
        -wy * length * length / 2.0,
        0.0,
        0.0,
    ];
    for (index, expected) in expected.iter().enumerate() {
        let (_, value) = last_sample(&outcome, index).expect("recorder sample");
        assert!(
            (value - expected).abs() < 1e-6,
            "recorder {index}: expected {expected}, got {value}"
        );
    }
}

/// Element loads on a kind `core` doesn't apply them to are rejected up front rather than
/// silently dropped.
#[test]
fn rejects_an_element_load_on_an_unsupported_element_kind() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.0],
        material: vec![0],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Constant],
        scale_factor: vec![1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::Truss],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy: 1.0,
            wz: None,
        }],
        stage: vec![0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
            hold_patterns_after: vec![],
        }],
        recorders: vec![],
    };
    assert!(decode(input).is_err());
}

/// `elementKind: "dispBeamColumn"` accepts a uniform element load: a two-fiber elastic section
/// reproduces `E·A`/`E·I`, so the tip deflection is the closed form and the recorded member end
/// forces are plain statics (fixed end carries the load, free end nothing).
#[test]
fn decodes_a_uniform_load_on_a_disp_beam_column() {
    let (length, e, area, iz, wy): (f64, f64, f64, f64, f64) = (100.0, 30_000.0, 2.0, 1000.0, -1.0);
    let h = (iz / area).sqrt();
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 2],
        z: vec![],
        y: vec![h, -h],
        area: vec![area / 2.0; 2],
        material: vec![0, 0],
    };
    input.disp_beam_columns_2d = FiberBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        integration: vec![IntegrationSpec::Lobatto { points: 5 }],
        corotational: vec![false],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::DispBeamColumn2d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy,
            wz: None,
        }],
        stage: vec![0],
    };
    let element_force = |component| RecorderSpec::ElementForce {
        element_kind: ElementKind::DispBeamColumn2d,
        element_index: 0,
        component,
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
            element_force(1), // shear at i
            element_force(2), // moment at i
            element_force(4), // shear at j
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    let expected = [
        wy * length.powi(4) / (8.0 * e * iz),
        -wy * length,
        -wy * length * length / 2.0,
        0.0,
    ];
    for (index, expected) in expected.iter().enumerate() {
        let (_, value) = last_sample(&outcome, index).expect("recorder sample");
        assert!(
            (value - expected).abs() < 1e-6,
            "recorder {index}: expected {expected}, got {value}"
        );
    }
}

/// `elementKind: "forceBeamColumn"` accepts a uniform element load: a two-fiber elastic section
/// reproduces `E·A`/`E·I`, so the tip deflection is the closed form and the recorded member end
/// forces are plain statics (fixed end carries the load, free end nothing).
#[test]
fn decodes_a_uniform_load_on_a_force_beam_column() {
    let (length, e, area, iz, wy): (f64, f64, f64, f64, f64) = (100.0, 30_000.0, 2.0, 1000.0, -1.0);
    let h = (iz / area).sqrt();
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 2],
        z: vec![],
        y: vec![h, -h],
        area: vec![area / 2.0; 2],
        material: vec![0, 0],
    };
    input.force_beam_columns_2d = FiberBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        integration: vec![IntegrationSpec::Lobatto { points: 5 }],
        corotational: vec![false],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ForceBeamColumn2d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform {
            wx: 0.0,
            wy,
            wz: None,
        }],
        stage: vec![0],
    };
    let element_force = |component| RecorderSpec::ElementForce {
        element_kind: ElementKind::ForceBeamColumn2d,
        element_index: 0,
        component,
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
            element_force(1), // shear at i
            element_force(2), // moment at i
            element_force(4), // shear at j
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    let expected = [
        wy * length.powi(4) / (8.0 * e * iz),
        -wy * length,
        -wy * length * length / 2.0,
        0.0,
    ];
    for (index, expected) in expected.iter().enumerate() {
        let (_, value) = last_sample(&outcome, index).expect("recorder sample");
        assert!(
            (value - expected).abs() < 1e-6,
            "recorder {index}: expected {expected}, got {value}"
        );
    }
}
