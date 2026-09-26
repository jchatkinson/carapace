//! M10 (implementation-plan.md / docs/pysees-handoff.md): exercises
//! `carapace_wasm::input_v1::decode` end to end by hand-building a
//! `CarapaceInputV1` value the way `pysees`'s (not-yet-written) compiler
//! eventually will, then driving the resulting `Session` with `advance`.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, FiberResponseKind, IntegratorSpec, RecorderSpec, SequenceSpec,
    SequenceSpec3, StageSpec,
};
use carapace_wasm::input_v1::tables::{
    ElasticBeamColumnTable, ElementKind, ElementLoadTable, EqualDofTable, FiberBeamColumnTable,
    FiberTable, IntegrationSpec, LoadPatternTable, NodalLoadTable, NodeTable, RigidDiaphragmTable,
    TimeSeriesSpec, TransformSpec, TrussTable, ZeroLengthSectionTable, ZeroLengthTable,
};
use carapace_wasm::input_v1::tables3::{
    ElasticBeamColumnTable3, ElementLoadTable3, EqualDofTable3, FiberBeamColumnTable3,
    FiberTable3, NodeTable3, RigidDiaphragmTable3, TrussTable3, ZeroLengthSectionTable3,
    ZeroLengthTable3,
};
use carapace_wasm::input_v1::{decode, CarapaceInputV1, Header, StepOutcome};

/// The last sample a given recorder produced in one `advance()` call's outcome, or `None` if
/// that recorder didn't record this call (e.g. the call took no steps).
fn last_sample(outcome: &StepOutcome, recorder_index: usize) -> Option<(f64, f64)> {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder_index)
        .and_then(|batch| batch.samples.last())
        .copied()
}

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
        zero_length_sections: ZeroLengthSectionTable::default(),
        equal_dofs: EqualDofTable::default(),
        rigid_diaphragms: RigidDiaphragmTable::default(),
        load_patterns: LoadPatternTable::default(),
        nodal_loads: NodalLoadTable::default(),
        element_loads: ElementLoadTable::default(),
        sequence: SequenceSpec::default(),

        nodes3: NodeTable3::default(),
        fibers3: FiberTable3::default(),
        trusses3: TrussTable3::default(),
        elastic_beam_columns3: ElasticBeamColumnTable3::default(),
        disp_beam_columns3: FiberBeamColumnTable3::default(),
        force_beam_columns3: FiberBeamColumnTable3::default(),
        zero_lengths3: ZeroLengthTable3::default(),
        zero_length_sections3: ZeroLengthSectionTable3::default(),
        equal_dofs3: EqualDofTable3::default(),
        rigid_diaphragms3: RigidDiaphragmTable3::default(),
        nodal_loads3: NodalLoadTable::default(),
        element_loads3: ElementLoadTable3::default(),
        sequence3: SequenceSpec3::default(),
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
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = load * length / (area * e);
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (got - expected).abs() < 1e-9,
        "expected {expected}, got {got}"
    );

    let batch = outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == 0)
        .expect("recorder 0 batch");
    assert_eq!(batch.first_sample, 0);
    assert_eq!(batch.stage_index, 0);
    assert_eq!(batch.samples.len(), 1);
}

/// The second recorder kind (results-storage-indexeddb.md's "several more
/// types of recorders" plan): an `ElementForce` recorder alongside a
/// `NodeDisp` one, on a horizontal cantilever `ElasticBeamColumn` — a
/// statically determinate case, so the expected local force at both
/// recorded components is plain nodal-equilibrium statics (the same
/// technique `core/tests/element_local_force.rs` uses), independent of the
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
    input.elastic_beam_columns = ElasticBeamColumnTable {
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::ElasticBeamColumn,
                element_index: 0,
                // component 2 = rz_i (fixed-end reaction moment), component 4 = uy_j (tip shear).
                component: 2,
            },
            RecorderSpec::ElementForce {
                element_kind: ElementKind::ElasticBeamColumn,
                element_index: 0,
                component: 4,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected_tip_disp = fy * length.powi(3) / (3.0 * e * iz);
    let (_, tip_disp) = last_sample(&outcome, 0).expect("recorder 0 (node disp) sample");
    assert!((tip_disp - expected_tip_disp).abs() < 1e-6, "expected {expected_tip_disp}, got {tip_disp}");

    let (_, reaction_moment) = last_sample(&outcome, 1).expect("recorder 1 (element force, component 2) sample");
    let expected_reaction_moment = -length * fy;
    assert!(
        (reaction_moment - expected_reaction_moment).abs() < 1e-6,
        "expected {expected_reaction_moment}, got {reaction_moment}"
    );

    let (_, tip_shear) = last_sample(&outcome, 2).expect("recorder 2 (element force, component 4) sample");
    assert!((tip_shear - fy).abs() < 1e-6, "expected {fy}, got {tip_shear}");
}

/// `ZeroLength`'s friction coupling (`ZeroLengthTable::friction`), driven
/// through the decoder rather than directly against `core::Friction` (see
/// `core/src/model/elements/zero_length.rs`'s own
/// `friction_sticks_below_yield_surface`/`friction_slides_past_yield_surface`
/// unit tests, which this exercises end to end from wire tables): a normal
/// spring in compression plus a friction-coupled shear DOF, checked both
/// below and past the Coulomb yield surface.
#[test]
fn decodes_a_zero_length_with_friction_coupling() {
    let (k_normal, mu, k0, b) = (1000.0, 0.3, 500.0, 0.01);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111, 0b100], // node 1 free in ux/uy (normal/shear), rz fixed
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: k_normal }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)], // dof 0 (normal) uses material arena index 0
        friction: vec![(0, 0, 1, mu, k0, b)], // normal_dof=0, shear_dof=1
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Constant, TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0, 1.0],
    };
    // Pattern 0: a sustained compressive normal displacement, held constant
    // via a `Constant` series so the LoadControl integrator's ramp doesn't
    // also scale it. Pattern 1: ramps the shear direction.
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 1],
        node: vec![1, 1],
        dof: vec![0, 1],
        value: vec![-1.0, 5.0],
        stage: vec![0, 0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 1 }],
    };

    let mut session = decode(input).expect("well-formed friction input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    // A `LoadControl` integrator with `increment: 1.0` reaches the full
    // nodal-load values directly, so this is a plain equilibrium check, not
    // a linear analysis proportionality check: normal_force = -1000 ->
    // yield_force = 300; trial = k0*5.0 = 2500, well past yield, so the
    // shear DOF displaces past its own reference load's linear response.
    let (_, shear_disp) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(shear_disp.is_finite() && shear_disp != 0.0, "shear DOF should have displaced: {shear_disp}");
}

/// `ZeroLengthSection` (`core::ZeroLengthSection`, wired via
/// `ZeroLengthSectionTable`): reproduces the same closed-form axial/moment
/// stiffness `core`'s own
/// `zero_length_section_reproduces_symmetric_ea_ei` unit test checks
/// directly against `FiberSection`, here decoded entirely from wire tables
/// with a two-fiber symmetric section and no independent `uy` spring.
#[test]
fn decodes_a_zero_length_section_and_matches_closed_form_axial_stiffness() {
    let (e, area, iz, load): (f64, f64, f64, f64) = (30_000.0, 2.0, 1000.0, 50.0);
    let h = (iz / area).sqrt();
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111, 0b110], // node 1 free only in ux
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 2],
        y: vec![h, -h],
        area: vec![area / 2.0, area / 2.0],
        material: vec![0, 0],
    };
    input.zero_length_sections = ZeroLengthSectionTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed zero-length-section input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected = load / (e * area);
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((got - expected).abs() < 1e-9, "expected {expected}, got {got}");
}

/// `Material::Hysteretic`/`Material::Pinching4` (`materials.rs`'s newly
/// added arena entries): both decode without error and, for the
/// symmetric-envelope well-within-elastic-range case checked here, behave
/// like a plain elastic spring — a decode-plumbing smoke test, not a
/// re-check of either material's own hysteresis logic (that lives in
/// `core/src/model/materials/{hysteretic,pinching4}.rs`'s unit tests).
#[test]
fn decodes_hysteretic_and_pinching4_materials_in_the_elastic_range() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111, 0b110, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![
        MaterialSpec::Hysteretic {
            mom1p: 10.0,
            rot1p: 0.01,
            mom2p: 15.0,
            rot2p: 0.02,
            mom3p: 15.0,
            rot3p: 0.03,
            mom1n: -10.0,
            rot1n: -0.01,
            mom2n: -15.0,
            rot2n: -0.02,
            mom3n: -15.0,
            rot3n: -0.03,
            pinch_x: 0.5,
            pinch_y: 0.5,
            damfc1: 0.0,
            damfc2: 0.0,
            beta: 0.0,
        },
        MaterialSpec::Pinching4 {
            stress1p: 10.0,
            strain1p: 0.01,
            stress2p: 15.0,
            strain2p: 0.02,
            stress3p: 17.0,
            strain3p: 0.03,
            stress4p: 10.0,
            strain4p: 0.04,
            stress1n: -10.0,
            strain1n: -0.01,
            stress2n: -15.0,
            strain2n: -0.02,
            stress3n: -17.0,
            strain3n: -0.03,
            stress4n: -10.0,
            strain4n: -0.04,
            r_disp_p: 0.5,
            r_force_p: 0.25,
            u_force_p: 0.05,
            r_disp_n: 0.5,
            r_force_n: 0.25,
            u_force_n: 0.05,
            gamma_k_params: [0.0; 4],
            gamma_k_limit: 0.0,
            gamma_d_params: [0.0; 4],
            gamma_d_limit: 0.0,
            gamma_f_params: [0.0; 4],
            gamma_f_limit: 0.0,
            gamma_e: 10.0,
            dmg_cyc: carapace_wasm::input_v1::materials::Pinching4DmgCycSpec::EnergyBased,
        },
    ];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 0],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (1, 0, 1)],
        friction: vec![],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 0],
        node: vec![1, 2],
        dof: vec![0, 0],
        value: vec![1.0, 1.0],
        stage: vec![0, 0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::NewtonRaphson,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 20 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed Hysteretic/Pinching4 input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let (_, hysteretic_disp) = last_sample(&outcome, 0).expect("hysteretic sample");
    let (_, pinching4_disp) = last_sample(&outcome, 1).expect("pinching4 sample");
    assert!(hysteretic_disp.is_finite() && hysteretic_disp > 0.0, "got {hysteretic_disp}");
    assert!(pinching4_disp.is_finite() && pinching4_disp > 0.0, "got {pinching4_disp}");
}

#[test]
fn rejects_an_unrecognized_space_before_constructing_a_session() {
    let err = decode(empty_input(4)).unwrap_err();
    assert_eq!(
        err,
        carapace_wasm::input_v1::DecodeError::UnsupportedSpace { got: 4 }
    );
}

/// The M5 golden-ratio closed form (core/tests/m5_modal.rs), decoded
/// entirely from wire tables: a `Modal` stage's mode frequencies come back
/// via `ModeShape` recorders, whose "time" component (`StepOutcome::
/// load_factor`'s doc comment) is the mode's own frequency rather than a
/// pseudo-time or load factor.
#[test]
fn decodes_a_modal_stage_and_matches_the_golden_ratio_closed_form() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0, 2.0, 0.0],
        fixed: vec![0b111, 0b110, 0b110],
        mass_node_index: vec![1, 2],
        mass: vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0, 1],
        node_j: vec![1, 2],
        materials: vec![(0, 0, 0), (1, 0, 0)],
        friction: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 2,
        }],
        recorders: vec![
            RecorderSpec::ModeShape { mode: 0, node: 1, dof: 0 },
            RecorderSpec::ModeShape { mode: 1, node: 1, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mode0_freq, _) = last_sample(&outcome, 0).expect("mode 0 sample");
    let (mode1_freq, _) = last_sample(&outcome, 1).expect("mode 1 sample");
    assert!((mode0_freq - 1.0 / phi).abs() < 1e-9, "expected {}, got {mode0_freq}", 1.0 / phi);
    assert!((mode1_freq - phi).abs() < 1e-9, "expected {phi}, got {mode1_freq}");
}

/// An out-of-range `ModeShape.mode` records nothing this call, rather than
/// erroring — `decode` can't validate it upfront (a mode index isn't a
/// table row; it's only meaningful once `modal_analysis` has actually run
/// and produced however many modes it was asked for).
#[test]
fn a_mode_shape_recorder_past_the_computed_mode_count_records_nothing() {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![1],
        mass: vec![1.0, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal { id: "modes".to_string(), modes: 1 }],
        recorders: vec![RecorderSpec::ModeShape { mode: 5, node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");
    assert!(last_sample(&outcome, 0).is_none(), "out-of-range mode should record nothing");
}

/// Rayleigh-damped free vibration (core/tests/m6_dynamics.rs's SDOF case),
/// decoded from a `Transient` stage: matches the closed-form damped
/// free-vibration displacement after a known number of Newmark steps.
#[test]
fn decodes_a_transient_stage_and_matches_damped_free_vibration_closed_form() {
    let (mass, k, dt, steps) = (1.0_f64, 100.0, 0.001, 500);
    let alpha_m = 0.2; // Rayleigh mass-proportional damping coefficient
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 1.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![1],
        mass: vec![mass, 0.0, 0.0],
    };
    input.materials = vec![MaterialSpec::Elastic { e: k }];
    input.zero_lengths = ZeroLengthTable {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)],
        friction: vec![],
    };
    // The wire format has no initial-displacement field yet, so this
    // exercises the Newmark integration a different (still closed-form)
    // way: a `Constant` reference load, already active when the transient
    // stage starts, is a suddenly-applied constant force from rest — the
    // classic damped step response, converging to static equilibrium
    // `u_static = force/k` with decaying oscillation, rather than free
    // vibration decaying to zero from a nonzero initial offset.
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Constant],
        scale_factor: vec![1.0],
    };
    let applied_force = 10.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![applied_force],
        stage: vec![0],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Transient {
            id: "dynamic".to_string(),
            steps,
            dt,
            damping: carapace_wasm::input_v1::sequence::DampingSpec { alpha_m, beta_k: 0.0 },
            ground_motions: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed transient input should decode");
    let outcome = session.advance(steps);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    // Closed-form damped step response from rest (Chopra, "Dynamics of
    // Structures"): u(t) = u_static * (1 - e^{-zeta*omega*t} * (cos(omega_d
    // t) + (zeta*omega/omega_d) * sin(omega_d t))), damping ratio
    // `zeta = alpha_m/(2*omega)` (mass-proportional-only Rayleigh damping).
    let omega = (k / mass).sqrt();
    let zeta = alpha_m / (2.0 * omega);
    let omega_d = omega * (1.0 - zeta * zeta).sqrt();
    let u_static = applied_force / k;
    let t = dt * steps as f64;
    let expected = u_static
        * (1.0
            - (-zeta * omega * t).exp()
                * ((omega_d * t).cos() + (zeta * omega / omega_d) * (omega_d * t).sin()));

    // Newmark discretization error, not exact — same tolerance
    // core/tests/m6_dynamics.rs's own closed-form comparison uses.
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((got - expected).abs() < 1e-4, "expected {expected}, got {got}");
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
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
        ],
    };

    let mut session = decode(input).expect("well-formed planar input should decode");

    // Drain the gravity stage first so we can read its frozen axial value.
    let mut outcome = session.advance(2);
    assert!(
        outcome.stage_complete && outcome.error.is_none(),
        "gravity stage should complete cleanly: {outcome:?}"
    );
    let axial_after_gravity = last_sample(&outcome, 0).unwrap().1;
    assert!(
        axial_after_gravity.abs() > 1e-9,
        "gravity should produce nonzero axial displacement"
    );
    // Gravity took both budgeted steps in this one advance() call, so recorder 0's batch should
    // cover samples [0, 2) of stage 0 — the "exact sample/time/value ordering across multiple
    // advance() calls and stage transitions" this bounded-batch slice needs to prove.
    let gravity_batch = outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == 0)
        .expect("recorder 0 batch");
    assert_eq!(gravity_batch.first_sample, 0);
    assert_eq!(gravity_batch.stage_index, 0);
    assert_eq!(gravity_batch.samples.len(), 2);

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
    let (_, transverse_1) = last_sample(&outcome, 1).unwrap();
    let (_, axial_1) = last_sample(&outcome, 0).unwrap();
    let expected_elastic = -0.5 * p_yield * length.powi(3) / (3.0 * e * iz);
    assert!(
        (transverse_1 - expected_elastic).abs() < 1e-9,
        "expected {expected_elastic}, got {transverse_1}"
    );
    assert!(
        (axial_1 - axial_after_gravity).abs() < 1e-9,
        "frozen gravity pattern's axial contribution must be unchanged while the section is still elastic: {axial_after_gravity} vs {axial_1}"
    );
    // This call crossed into the pushover stage: its one sample continues the running count
    // from gravity's two (first_sample == 2), tagged with the new stage index.
    let pushover_batch = outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == 1)
        .expect("recorder 1 batch");
    assert_eq!(pushover_batch.first_sample, 2);
    assert_eq!(pushover_batch.stage_index, 1);
    assert_eq!(pushover_batch.samples.len(), 1);

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

    let (_, transverse_final) = last_sample(&outcome, 1).unwrap();
    let elastic_force_for_final_displacement =
        transverse_final.abs() * 3.0 * e * iz / length.powi(3);
    assert!(
        outcome.load_factor.abs() < elastic_force_for_final_displacement * 0.99,
        "softened section should need measurably less force than elastic stiffness predicts: needed {}, elastic would need {elastic_force_for_final_displacement}",
        outcome.load_factor.abs()
    );
}

/// The M1 fixed-free truss closed form (core/tests/reaction.rs's own
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::Reaction { node: 0, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed reaction input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let (_, reaction) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((reaction + load).abs() < 1e-9, "expected {}, got {reaction}", -load);
}

/// The `DispBeamColumn` pure-axial closed form (core/tests/
/// fiber_recorder.rs's own `disp_beam_column_fiber_responses_match_hand_
/// computed_strain_and_stress`), decoded entirely from wire tables: a
/// `Fiber` recorder on each of a symmetric section's two fibers must read
/// identical strain/stress (pure axial extension, zero curvature).
#[test]
fn decodes_a_fiber_recorder_and_matches_hand_computed_strain_and_stress() {
    let (e, area, iz, length, axial_load): (f64, f64, f64, f64, f64) = (30_000.0, 2.0, 1000.0, 100.0, 60.0);
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
        y: vec![h, -h],
        area: vec![area / 2.0, area / 2.0],
        material: vec![0, 0],
    };
    input.disp_beam_columns = FiberBeamColumnTable {
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::Fiber {
                element_kind: ElementKind::DispBeamColumn,
                element_index: 0,
                point: 0,
                fiber: 0,
                quantity: FiberResponseKind::Strain,
            },
            RecorderSpec::Fiber {
                element_kind: ElementKind::DispBeamColumn,
                element_index: 0,
                point: 0,
                fiber: 1,
                quantity: FiberResponseKind::Stress,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed fiber recorder input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected_strain = axial_load / (e * area);
    let expected_stress = e * expected_strain;
    let (_, strain) = last_sample(&outcome, 0).expect("fiber 0 strain sample");
    let (_, stress) = last_sample(&outcome, 1).expect("fiber 1 stress sample");
    assert!((strain - expected_strain).abs() < 1e-9, "expected {expected_strain}, got {strain}");
    assert!((stress - expected_stress).abs() < 1e-6, "expected {expected_stress}, got {stress}");
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
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
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");
    assert!(last_sample(&outcome, 0).is_none(), "a non-fiber element's fiber recorder should record nothing");
}

/// `core::Domain::equal_dof` (`EqualDofTable`, `ConstraintHandler::
/// Transformation` picked automatically by `session::start_current_stage`
/// once `Domain::has_mp_constraints()` is true): node 2 carries no element
/// of its own and is free only in `ux`, tied via `equal_dof` to node 1's
/// `ux` — so its recorded displacement must exactly equal node 1's, which
/// itself matches the plain M1 closed-form truss elongation.
#[test]
fn decodes_an_equal_dof_constraint_tying_one_nodes_ux_to_another() {
    let (length, area, e, load) = (100.0, 2.0, 30_000.0, 50.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0, 2.0 * length, 0.0],
        // node 0 fully fixed; node 1 free only in ux; node 2 (no element of
        // its own) also free only in ux, tied to node 1's via equal_dof.
        fixed: vec![0b111, 0b110, 0b110],
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
    input.equal_dofs = EqualDofTable {
        retained: vec![1],
        constrained: vec![2],
        dofs: vec![(0, 0)],
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed equal_dof input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected = load * length / (area * e);
    let (_, retained_disp) = last_sample(&outcome, 0).expect("retained node sample");
    let (_, constrained_disp) = last_sample(&outcome, 1).expect("constrained node sample");
    assert!(
        (retained_disp - expected).abs() < 1e-9,
        "expected {expected}, got {retained_disp}"
    );
    assert!(
        (constrained_disp - retained_disp).abs() < 1e-9,
        "equal_dof-tied node should exactly match the retained node's ux: {constrained_disp} vs {retained_disp}"
    );
}

/// `core::Domain::rigid_diaphragm` (`RigidDiaphragmTable`) — same setup as
/// the `equal_dof` test above, but tying *two* otherwise-elementless nodes'
/// `ux` to the retained node in one call, proving the table's sparse
/// per-row `constrained` list decodes correctly.
#[test]
fn decodes_a_rigid_diaphragm_tying_two_nodes_ux_to_the_retained_node() {
    let (length, area, e, load) = (100.0, 2.0, 30_000.0, 50.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0, 2.0 * length, 0.0, 3.0 * length, 0.0],
        fixed: vec![0b111, 0b110, 0b110, 0b110],
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
    input.rigid_diaphragms = RigidDiaphragmTable {
        retained: vec![1],
        constrained: vec![(0, 2), (0, 3)],
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
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
            RecorderSpec::NodeDisp { node: 3, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed rigid_diaphragm input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected = load * length / (area * e);
    let (_, retained_disp) = last_sample(&outcome, 0).expect("retained node sample");
    assert!((retained_disp - expected).abs() < 1e-9, "expected {expected}, got {retained_disp}");
    for recorder_index in [1, 2] {
        let (_, disp) = last_sample(&outcome, recorder_index).expect("diaphragm-tied node sample");
        assert!(
            (disp - retained_disp).abs() < 1e-9,
            "diaphragm-tied node should exactly match the retained node's ux: {disp} vs {retained_disp}"
        );
    }
}
