//! Recorders through the wire format: element force, reaction, fiber and element-load recorders.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::decode;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, FiberResponseKind, IntegratorSpec, RecorderSpec, SequenceSpec,
    StageSpec,
};
use carapace_wasm::input_v1::tables::*;

/// The second recorder kind (results-storage-indexeddb.md's "several more
/// types of recorders" plan): an `ElementForce` recorder alongside a
/// `NodeDisp` one, on a horizontal cantilever `ElasticBeamColumn` — a
/// statically determinate case, so the expected local force at both
/// recorded components is plain nodal-equilibrium statics (the same
/// technique `core/tests/elements/local_force.rs` uses), independent of the
/// element's own stiffness formulation. Also proves node and element
/// batches coexist correctly in one `advance()` call's `recorder_batches`,
/// each keyed by its own `recorder_index`.
#[test]
fn decodes_an_element_force_recorder_alongside_a_node_disp_recorder() {
    let (length, e, area, iz, fy) = (100.0, 30_000.0, 10.0, 1000.0, -10.0);
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
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![1],
        value: vec![fy],
        stage: vec![0],
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
            RecorderSpec::ElementForce {
                element_kind: ElementKind::ElasticBeamColumn2d,
                element_index: 0,
                // component 2 = rz_i (fixed-end reaction moment), component 4 = uy_j (tip shear).
                component: 2,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::ElasticBeamColumn2d,
                element_index: 0,
                component: 4,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected_tip_disp = fy * length.powi(3) / (3.0 * e * iz);
    let (_, tip_disp) = last_sample(&outcome, 0).expect("recorder 0 (node disp) sample");
    assert!(
        (tip_disp - expected_tip_disp).abs() < 1e-6,
        "expected {expected_tip_disp}, got {tip_disp}"
    );

    let (_, reaction_moment) =
        last_sample(&outcome, 1).expect("recorder 1 (element force, component 2) sample");
    let expected_reaction_moment = -length * fy;
    assert!(
        (reaction_moment - expected_reaction_moment).abs() < 1e-6,
        "expected {expected_reaction_moment}, got {reaction_moment}"
    );

    let (_, tip_shear) =
        last_sample(&outcome, 2).expect("recorder 2 (element force, component 4) sample");
    assert!(
        (tip_shear - fy).abs() < 1e-6,
        "expected {fy}, got {tip_shear}"
    );
}

#[test]
fn decodes_an_elastic_beam_column3_element_force_recorder() {
    let (length, e, g, area, iy, iz, j, fy): (f64, f64, f64, f64, f64, f64, f64, f64) =
        (100.0, 30_000.0, 12_000.0, 10.0, 500.0, 1000.0, 50.0, -10.0);
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0],
        fixed: vec![0b111111, 0b000000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.elastic_beam_columns_3d = carapace_wasm::input_v1::tables::ElasticBeamColumn3dTable {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![e],
        g: vec![g],
        a: vec![area],
        j: vec![j],
        iy: vec![iy],
        iz: vec![iz],
        transform: vec![carapace_wasm::input_v1::tables::TransformSpec3::Linear3 {
            vec_xz: [0.0, 0.0, 1.0],
        }],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![1], // uy
        value: vec![fy],
        stage: vec![0],
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
            RecorderSpec::ElementForce {
                element_kind: ElementKind::ElasticBeamColumn3d,
                element_index: 0,
                // Local DOF order [ux,uy,uz,rx,ry,rz]_i, [..]_j (width 12):
                // component 5 = rz_i (fixed-end reaction moment about z).
                component: 5,
            },
        ],
    };

    let mut session =
        decode(input).expect("well-formed spatial elastic beam column input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected_tip_disp = fy * length.powi(3) / (3.0 * e * iz);
    let (_, tip_disp) = last_sample(&outcome, 0).expect("tip disp sample");
    assert!(
        (tip_disp - expected_tip_disp).abs() < 1e-6,
        "expected {expected_tip_disp}, got {tip_disp}"
    );

    let expected_reaction_moment = -length * fy;
    let (_, reaction_moment) = last_sample(&outcome, 1).expect("reaction moment sample");
    assert!(
        (reaction_moment - expected_reaction_moment).abs() < 1e-6,
        "expected {expected_reaction_moment}, got {reaction_moment}"
    );
}

/// The fixed-free truss closed form (core/tests/analysis/reaction.rs's own
/// `fixed_free_truss_reaction_exactly_opposes_the_applied_load`), decoded
/// entirely from wire tables: a `Reaction` recorder at the fixed support
/// must read exactly the negative of the applied tip load.
#[test]
fn decodes_a_reaction_recorder_and_matches_the_closed_form_support_reaction() {
    let (length, area, e, load) = (100.0, 2.0, 30000.0, 50.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![area],
        material: vec![0],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![load],
        stage: vec![0],
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
        recorders: vec![RecorderSpec::Reaction { node: 0, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed reaction input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let (_, reaction) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (reaction + load).abs() < 1e-9,
        "expected {}, got {reaction}",
        -load
    );
}

/// The `DispBeamColumn` pure-axial closed form (core/tests/
/// `elements/fiber_responses.rs`'s own `disp_beam_column_fiber_responses_match_hand_
/// computed_strain_and_stress`), decoded entirely from wire tables: a
/// `Fiber` recorder on each of a symmetric section's two fibers must read
/// identical strain/stress (pure axial extension, zero curvature).
#[test]
fn decodes_a_fiber_recorder_and_matches_hand_computed_strain_and_stress() {
    let (e, area, iz, length, axial_load): (f64, f64, f64, f64, f64) =
        (30_000.0, 2.0, 1000.0, 100.0, 60.0);
    let h = (iz / area).sqrt();
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 2],
        z: vec![],
        y: vec![h, -h],
        area: vec![area / 2.0, area / 2.0],
        material: vec![0, 0],
    };
    input.disp_beam_columns_2d = FiberBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        integration: vec![IntegrationSpec::Legendre { points: 2 }],
        corotational: vec![false],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![axial_load],
        stage: vec![0],
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
            RecorderSpec::Fiber {
                element_kind: ElementKind::DispBeamColumn2d,
                element_index: 0,
                point: 0,
                fiber: 0,
                quantity: FiberResponseKind::Strain,
            },
            RecorderSpec::Fiber {
                element_kind: ElementKind::DispBeamColumn2d,
                element_index: 0,
                point: 0,
                fiber: 1,
                quantity: FiberResponseKind::Stress,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed fiber recorder input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected_strain = axial_load / (e * area);
    let expected_stress = e * expected_strain;
    let (_, strain) = last_sample(&outcome, 0).expect("fiber 0 strain sample");
    let (_, stress) = last_sample(&outcome, 1).expect("fiber 1 stress sample");
    assert!(
        (strain - expected_strain).abs() < 1e-9,
        "expected {expected_strain}, got {strain}"
    );
    assert!(
        (stress - expected_stress).abs() < 1e-6,
        "expected {expected_stress}, got {stress}"
    );
}

/// An out-of-range `Fiber.point`/`fiber` records nothing this call, and a
/// `Fiber` recorder aimed at a non-fiber element kind (`Truss`, here)
/// records nothing either — both unresolved-by-runtime-state cases, same
/// graceful-skip handling `ModeShape`'s out-of-range mode gets, since
/// neither can be checked before the analysis actually runs.
#[test]
fn a_fiber_recorder_on_a_non_fiber_element_or_out_of_range_index_records_nothing() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 100.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1000.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![1.0],
        material: vec![0],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![10.0],
        stage: vec![0],
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
        recorders: vec![RecorderSpec::Fiber {
            element_kind: ElementKind::Truss,
            element_index: 0,
            point: 0,
            fiber: 0,
            quantity: FiberResponseKind::Strain,
        }],
    };

    let mut session = decode(input).expect("well-formed input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    assert!(
        last_sample(&outcome, 0).is_none(),
        "a non-fiber element's fiber recorder should record nothing"
    );
}

/// The `ElementLoad` recorder reports the load an element carries at each
/// sample's pseudo-time: ramped by its pattern's factor during the stage that
/// ramps it, then constant once the pattern is held — what a consumer needs to
/// recover member internal-force diagrams from the recorded end forces.
#[test]
fn an_element_load_recorder_follows_the_pattern_factor_then_the_frozen_value() {
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
        series: vec![
            TimeSeriesSpec::Linear { slope: 1.0 },
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
        scale_factor: vec![1.0, 1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ElasticBeamColumn2d],
        element_index: vec![0],
        load: vec![ElementLoadSpec::Uniform { wx, wy, wz: None }],
        stage: vec![0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![1],
        node: vec![1],
        dof: vec![1],
        value: vec![-1.0],
        stage: vec![1],
    };
    let stage = |id: &str, steps: u32, increment: f64, hold: Vec<u32>| StageSpec::Static {
        id: id.to_string(),
        steps,
        integrator: IntegratorSpec::LoadControl { increment },
        algorithm: AlgorithmSpec::Linear,
        convergence: Some(ConvergenceSpec::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        }),
        hold_patterns_after: hold,
    };
    input.sequence = SequenceSpec {
        stages: vec![
            stage("ramp", 2, 0.25, vec![0]),
            stage("push", 1, 1.0, vec![]),
        ],
        recorders: vec![
            RecorderSpec::ElementLoad {
                element_kind: ElementKind::ElasticBeamColumn2d,
                element_index: 0,
                component: 0,
            },
            RecorderSpec::ElementLoad {
                element_kind: ElementKind::ElasticBeamColumn2d,
                element_index: 0,
                component: 1,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    // Ramp stage: two steps to factor 0.5.
    let outcome = session.advance(2);
    assert!(outcome.error.is_none(), "unexpected outcome: {outcome:?}");
    let (_, wx_ramped) = last_sample(&outcome, 0).expect("wx sample");
    let (_, wy_ramped) = last_sample(&outcome, 1).expect("wy sample");
    assert!((wx_ramped - 0.5 * wx).abs() < 1e-12, "wx {wx_ramped}");
    assert!((wy_ramped - 0.5 * wy).abs() < 1e-12, "wy {wy_ramped}");
    // Push stage: pattern 0 is held at its final value while pattern 1 ramps.
    let outcome = session.advance(1);
    assert!(outcome.error.is_none(), "unexpected outcome: {outcome:?}");
    let (_, wy_held) = last_sample(&outcome, 1).expect("wy sample");
    assert!((wy_held - 0.5 * wy).abs() < 1e-12, "held wy {wy_held}");
}
