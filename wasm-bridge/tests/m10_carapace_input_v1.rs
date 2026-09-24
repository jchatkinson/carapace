//! M10 (implementation-plan.md / docs/pysees-handoff.md): exercises
//! `carapace_wasm::input_v1::decode` end to end by hand-building a
//! `CarapaceInputV1` value the way `pysees`'s (not-yet-written) compiler
//! eventually will, then driving the resulting `Session` with `advance`.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::{
    ElasticBeamColumnTable, ElementLoadTable, FiberBeamColumnTable, FiberTable, IntegrationSpec,
    LoadPatternTable, NodalLoadTable, NodeTable, TimeSeriesSpec, TrussTable, ZeroLengthTable,
};
use carapace_wasm::input_v1::{decode, CarapaceInputV1, Header};

fn empty_input(header_space: u8) -> CarapaceInputV1 {
    CarapaceInputV1 {
        header: Header {
            schema_version: 1,
            space: header_space,
            engine_version: "test".to_string(),
        },
        nodes: NodeTable::default(),
        materials: Vec::new(),
        fibers: FiberTable::default(),
        trusses: TrussTable::default(),
        elastic_beam_columns: ElasticBeamColumnTable::default(),
        disp_beam_columns: FiberBeamColumnTable::default(),
        force_beam_columns: FiberBeamColumnTable::default(),
        zero_lengths: ZeroLengthTable::default(),
        load_patterns: LoadPatternTable::default(),
        nodal_loads: NodalLoadTable::default(),
        element_loads: ElementLoadTable::default(),
        sequence: SequenceSpec::default(),
    }
}

/// Smoke test for the whole node/material/element/pattern/sequence
/// pipeline: the M1 truss case (`core/tests/m1_truss.rs`'s closed form),
/// built entirely from wire-format tables instead of direct `core` calls.
#[test]
fn decodes_a_single_truss_and_matches_the_closed_form_displacement() {
    let (length, area, e, load) = (100.0, 2.0, 30000.0, 50.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b110], // node 0 fully fixed, node 1 fixed except ux
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
        recorders: vec![RecorderSpec { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = load * length / (area * e);
    let (_, got) = *session
        .recorder_samples(0)
        .last()
        .expect("one recorded sample");
    assert!(
        (got - expected).abs() < 1e-9,
        "expected {expected}, got {got}"
    );
}

#[test]
fn rejects_the_spatial_profile_before_constructing_a_session() {
    let err = decode(empty_input(3)).unwrap_err();
    assert_eq!(
        err,
        carapace_wasm::input_v1::DecodeError::UnsupportedSpace { got: 3 }
    );
}

#[test]
fn rejects_a_modal_stage_with_a_named_diagnostic_not_a_panic() {
    let mut input = empty_input(2);
    input.sequence.stages = vec![StageSpec::Modal {
        id: "modes".to_string(),
        modes: 2,
    }];
    let err = decode(input).unwrap_err();
    assert_eq!(
        err,
        carapace_wasm::input_v1::DecodeError::UnsupportedStage {
            stage_id: "modes".to_string(),
            kind: "modal"
        }
    );
}

/// The handoff's required first acceptance case: a 2D fiber-section
/// cantilever runs gravity (`LoadControl`), freezes it, then a
/// displacement-controlled lateral pushover — built the same way
/// `core/tests/m8_force_beam_column.rs` builds and checks a
/// `ForceBeamColumn` cantilever, but assembled entirely from decoded wire
/// tables. Elastic-range and post-yield behavior must agree with that
/// native acceptance test, and gravity's frozen axial displacement must
/// survive the pushover stage untouched (`core/tests/m_load_patterns.rs`'s
/// same check, here across a decoded multi-stage `AnalysisSequence`).
#[test]
fn gravity_then_pushover_matches_native_force_beam_column_behavior() {
    let (e, length) = (30000.0_f64, 100.0);
    let area = 0.25;
    let ys = [20.0_f64, 15.0, 10.0, 5.0];
    let iz = area * ys.iter().map(|y| y * y * 2.0).sum::<f64>();
    let eyp = 0.001;
    let fy = e * eyp;
    let my = fy * iz / ys[0];
    let p_yield = my / length;

    // Fiber section: 8 fibers (±20, ±15, ±10, ±5), ElasticPP.
    let mut fiber_y = Vec::new();
    let mut fiber_area = Vec::new();
    let mut fiber_material = Vec::new();
    for y in ys.iter().flat_map(|&y| [y, -y]) {
        fiber_y.push(y);
        fiber_area.push(area);
        fiber_material.push(0u32); // single shared ElasticPP material entry
    }

    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b111, 0b000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::ElasticPp { e, eyp }];
    input.fibers = FiberTable {
        section_offsets: vec![0, fiber_y.len() as u32],
        y: fiber_y,
        area: fiber_area,
        material: fiber_material,
    };
    input.force_beam_columns = FiberBeamColumnTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        integration: vec![IntegrationSpec::Lobatto { points: 3 }],
        corotational: vec![false],
        density: vec![0.0],
    };
    // Pattern 0: gravity (axial, dof 0), registered on stage 0. Pattern 1:
    // lateral reference load (transverse, dof 1), registered only on stage
    // 1 — *not* upfront, since stage 0's `LoadControl` shares one
    // pseudo-time across every unfrozen pattern, so an already-registered
    // pattern 1 would ramp during gravity too (the same reason
    // core/tests/m8_force_beam_column.rs's native two-phase test adds its
    // lateral load only between phases). `DisplacementControl` needs this
    // reference load's nonzero sensitivity on the controlled dof, same as
    // core/tests/m4_analysis.rs.
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Linear { slope: 1.0 },
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
        scale_factor: vec![1.0, 1.0],
    };
    // Small relative to this section's axial yield capacity (area_total *
    // fy = 2.0 * 30.0 = 60.0) and, more importantly, small enough that the
    // uniform axial strain it adds doesn't push the elastic-range pushover
    // check's outer fiber (already at half its yield curvature by design,
    // mirroring core/tests/m8_force_beam_column.rs) past `eyp`.
    let gravity_axial = -10.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 1],
        node: vec![1, 1],
        dof: vec![0, 1],
        value: vec![gravity_axial, -1.0],
        stage: vec![0, 1],
    };
    input.sequence = SequenceSpec {
        stages: vec![
            StageSpec::Static {
                id: "gravity".to_string(),
                steps: 2,
                integrator: IntegratorSpec::LoadControl { increment: 0.5 },
                algorithm: AlgorithmSpec::NewtonRaphson,
                convergence: Some(ConvergenceSpec::NormUnbalance {
                    tol: 1e-10,
                    max_iter: 30,
                }),
                hold_patterns_after: vec![0],
            },
            StageSpec::Static {
                id: "pushover".to_string(),
                steps: 3,
                integrator: IntegratorSpec::DisplacementControl {
                    node: 1,
                    dof: 1,
                    increment: -0.5 * p_yield * length.powi(3) / (3.0 * e * iz),
                },
                algorithm: AlgorithmSpec::NewtonRaphson,
                convergence: Some(ConvergenceSpec::NormUnbalance {
                    tol: 1e-10,
                    max_iter: 30,
                }),
                hold_patterns_after: vec![],
            },
        ],
        recorders: vec![
            RecorderSpec { node: 1, dof: 0 },
            RecorderSpec { node: 1, dof: 1 },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");

    // Drain the gravity stage first so we can read its frozen axial value.
    let mut outcome = session.advance(2);
    assert!(
        outcome.stage_complete && outcome.error.is_none(),
        "gravity stage should complete cleanly: {outcome:?}"
    );
    let axial_after_gravity = session.recorder_samples(0).last().unwrap().1;
    assert!(
        axial_after_gravity.abs() > 1e-9,
        "gravity should produce nonzero axial displacement"
    );

    // First pushover step stays within the elastic range: must match the
    // closed-form cantilever deflection exactly, same as m8's elastic
    // check. The section's fiber layout is symmetric and still uniformly
    // elastic at this point, so — same reasoning as
    // core/tests/m_load_patterns.rs's frozen-pattern check, here for a
    // fiber section rather than `ElasticBeamColumn` — gravity's frozen
    // axial displacement must also still be completely untouched.
    outcome = session.advance(1);
    assert!(
        outcome.error.is_none(),
        "elastic-range pushover step should converge: {outcome:?}"
    );
    let (_, transverse_1) = *session.recorder_samples(1).last().unwrap();
    let (_, axial_1) = *session.recorder_samples(0).last().unwrap();
    let expected_elastic = -0.5 * p_yield * length.powi(3) / (3.0 * e * iz);
    assert!(
        (transverse_1 - expected_elastic).abs() < 1e-9,
        "expected {expected_elastic}, got {transverse_1}"
    );
    assert!(
        (axial_1 - axial_after_gravity).abs() < 1e-9,
        "frozen gravity pattern's axial contribution must be unchanged while the section is still elastic: {axial_after_gravity} vs {axial_1}"
    );

    // Remaining two steps push well past yield (cumulative target 1.5x the
    // elastic-range increment, past the outer fibers' `eyp`). Checking the
    // controlled dof's own displacement here would be circular —
    // `DisplacementControl` drives it to exactly the requested cumulative
    // increment regardless of material state, converged or not. Instead
    // check the *force* (`load_factor`, since the reference lateral load is
    // a unit load) needed to reach that displacement: a softened section
    // needs measurably less force than pure elastic stiffness would predict
    // for the same displacement, the load-controlled equivalent of m8's
    // "softer than elastic" check.
    //
    // Axial/bending decoupling itself does *not* survive into the
    // post-yield region here (unlike m_load_patterns.rs's `ElasticBeamColumn`
    // case): once the fiber layout's yielding is asymmetric, the section's
    // *tangent* couples axial and bending even though its geometry is
    // symmetric — a real property of fiber sections, not a decode bug — so
    // this test doesn't assert axial invariance past this point.
    outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "pushover stage should finish cleanly: {outcome:?}"
    );

    let (_, transverse_final) = *session.recorder_samples(1).last().unwrap();
    let elastic_force_for_final_displacement =
        transverse_final.abs() * 3.0 * e * iz / length.powi(3);
    assert!(
        outcome.load_factor.abs() < elastic_force_for_final_displacement * 0.99,
        "softened section should need measurably less force than elastic stiffness predicts: needed {}, elastic would need {elastic_force_for_final_displacement}",
        outcome.load_factor.abs()
    );
}
