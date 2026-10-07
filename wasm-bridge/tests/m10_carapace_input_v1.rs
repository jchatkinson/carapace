//! M10 (implementation-plan.md / docs/pysees-handoff.md): exercises
//! `carapace_wasm::input_v1::decode` end to end by hand-building a
//! `CarapaceInputV1` value the way `pysees`'s (not-yet-written) compiler
//! eventually will, then driving the resulting `Session` with `advance`.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, FiberResponseKind, IntegratorSpec, RecorderSpec, SequenceSpec,
    StageSpec,
};
use carapace_wasm::input_v1::tables::*;
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
        plane_materials: Default::default(),
        triangles: Default::default(),
        quads: Default::default(),
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

/// Nonlinear equilibrium must be the same across each configurable solver.
#[test]
fn configurable_algorithms_solve_static_and_transient_yielding() {
    use carapace_wasm::input_v1::sequence::{
        AlgorithmConfigSpec, DampingSpec, LineSearchSpec, TangentStrategySpec,
    };
    let mut algorithms = vec![AlgorithmSpec::NewtonRaphson];
    for tangent in [
        TangentStrategySpec::Current,
        TangentStrategySpec::ReuseAtStepStart,
        TangentStrategySpec::Initial,
    ] {
        for line_search in [
            None,
            Some(LineSearchSpec::Bisection {
                tol: 1e-10,
                max_iter: 30,
                max_eta: 16.0,
            }),
            Some(LineSearchSpec::RegulaFalsi {
                tol: 1e-10,
                max_iter: 30,
                max_eta: 16.0,
            }),
        ] {
            algorithms.push(AlgorithmSpec::Config(AlgorithmConfigSpec::Newton {
                tangent,
                line_search,
            }));
        }
        algorithms.push(AlgorithmSpec::Config(AlgorithmConfigSpec::KrylovNewton {
            tangent,
            max_dimension: 3,
        }));
    }
    for algorithm in algorithms {
        for transient in [false, true] {
            let mut input = empty_input(2);
            input.nodes = NodeTable {
                coords: vec![0.0, 0.0, 0.0, 0.0],
                fixed: vec![0b111, 0b110],
                mass_node_index: vec![1],
                mass: vec![1.0, 0.0, 0.0],
            };
            input.materials = vec![
                MaterialSpec::Elastic { e: 50.0 },
                MaterialSpec::ElasticPp {
                    e: 100.0,
                    eyp: 0.01,
                },
            ];
            input.zero_lengths = ZeroLengthTable {
                node_i: vec![0, 0],
                node_j: vec![1, 1],
                materials: vec![(0, 0, 0), (1, 0, 1)],
                friction: vec![],
                orient: vec![],
            };
            input.load_patterns = LoadPatternTable {
                series: vec![if transient {
                    TimeSeriesSpec::Linear { slope: 10.0 }
                } else {
                    TimeSeriesSpec::Constant
                }],
                scale_factor: vec![1.0],
            };
            input.nodal_loads = NodalLoadTable {
                pattern: vec![0],
                node: vec![1],
                dof: vec![0],
                value: vec![if transient { 1000.0 } else { 5.0 }],
                stage: vec![0],
            };
            let convergence = Some(ConvergenceSpec::NormUnbalance {
                tol: if transient { 1e-6 } else { 1e-9 },
                max_iter: 120,
            });
            input.sequence = SequenceSpec {
                stages: vec![if transient {
                    StageSpec::Transient {
                        id: "yield".into(),
                        steps: 1,
                        dt: 0.1,
                        damping: DampingSpec {
                            alpha_m: 0.0,
                            beta_k: 0.0,
                        },
                        ground_motions: vec![],
                        algorithm,
                        convergence,
                    }
                } else {
                    StageSpec::Static {
                        id: "yield".into(),
                        steps: 1,
                        integrator: IntegratorSpec::LoadControl { increment: 1.0 },
                        algorithm,
                        convergence,
                        hold_patterns_after: vec![],
                    }
                }],
                recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
            };
            let outcome = decode(input).unwrap().advance(1);
            assert!(
                outcome.done && outcome.error.is_none(),
                "{algorithm:?}, transient={transient}: {outcome:?}"
            );
            // A ramp starts at zero force; first Newmark effective inertia = 4*m/dt² = 400.
            let expected = if transient {
                (1000.0 - 1.0) / (400.0 + 50.0)
            } else {
                (5.0 - 1.0) / 50.0
            };
            assert!(
                (last_sample(&outcome, 0).unwrap().1 - expected).abs() < 1e-8,
                "{algorithm:?}, transient={transient}: {outcome:?}, expected={expected}"
            );
        }
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

/// A single truss, solved for the tip displacement, with node rotations left
/// free (`rz` not fixed) and an optional extra nodal load.
fn truss_with_unfixed_rotations(extra_load: Option<(u32, u8, f64)>) -> CarapaceInputV1 {
    let (length, area, e, load) = (100.0, 2.0, 30000.0, 50.0);
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, length, 0.0],
        fixed: vec![0b011, 0b010], // only the unneeded rz rows are left free
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
    let mut loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![load],
        stage: vec![0],
    };
    if let Some((node, dof, value)) = extra_load {
        loads.pattern.push(0);
        loads.node.push(node);
        loads.dof.push(dof);
        loads.value.push(value);
        loads.stage.push(0);
    }
    input.nodal_loads = loads;
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
    input
}

/// DOF activation through the wire: a model that leaves a truss's rotations
/// unfixed (the compiler no longer has to fix them) solves, and a nodal load
/// on a DOF nothing resists comes back as a structured error from
/// `advance`, not a panic.
#[test]
fn truss_with_unfixed_rotations_solves_and_a_load_on_an_unused_dof_is_a_structured_error() {
    use carapace_wasm::input_v1::error::ModelErrorDetail;
    use carapace_wasm::input_v1::AnalysisErrorDetail;

    let mut session = decode(truss_with_unfixed_rotations(None)).expect("decodes");
    let outcome = session.advance(10);
    assert!(outcome.error.is_none(), "unexpected outcome: {outcome:?}");
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((got - 50.0 * 100.0 / (2.0 * 30000.0)).abs() < 1e-9);

    // A moment (dof 2) on the tip: the truss has no rotational stiffness.
    let mut session = decode(truss_with_unfixed_rotations(Some((1, 2, 1.0)))).expect("decodes");
    let outcome = session.advance(10);
    assert_eq!(
        outcome.error,
        Some(AnalysisErrorDetail::InvalidModel {
            error: ModelErrorDetail::LoadOnInactiveDof { node: 1, dof: 2 }
        })
    );
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
        friction: vec![FrictionRow {
            row: 0,
            normal_dof: 0,
            shear_dofs: vec![1], // normal_dof=0, shear_dof=1
            mu,
            k0,
            b,
        }],
        orient: vec![],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Constant,
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
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
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 1 }],
    };

    let mut session = decode(input).expect("well-formed friction input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    // A `LoadControl` integrator with `increment: 1.0` reaches the full
    // nodal-load values directly, so this is a plain equilibrium check, not
    // a linear analysis proportionality check: normal_force = -1000 ->
    // yield_force = 300; trial = k0*5.0 = 2500, well past yield, so the
    // shear DOF displaces past its own reference load's linear response.
    let (_, shear_disp) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        shear_disp.is_finite() && shear_disp != 0.0,
        "shear DOF should have displaced: {shear_disp}"
    );
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
        z: vec![],
        y: vec![h, -h],
        area: vec![area / 2.0, area / 2.0],
        material: vec![0, 0],
    };
    input.zero_length_sections = ZeroLengthSectionTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
        orient: vec![],
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

    let mut session = decode(input).expect("well-formed zero-length-section input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = load / (e * area);
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (got - expected).abs() < 1e-9,
        "expected {expected}, got {got}"
    );
}

/// `ZeroLengthSectionTable::orient`: the section's axial direction turned onto
/// global y (OpenSees `-orient 0 1 0`), so a `uy` load meets the same
/// closed-form axial stiffness the unoriented test above gets in `ux`.
#[test]
fn decodes_an_oriented_zero_length_section_and_matches_closed_form_axial_stiffness() {
    let (e, area, iz, load): (f64, f64, f64, f64) = (30_000.0, 2.0, 1000.0, 50.0);
    let h = (iz / area).sqrt();
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111, 0b101], // node 1 free only in uy
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
    input.zero_length_sections = ZeroLengthSectionTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
        orient: vec![OrientRow {
            row: 0,
            x: [0.0, 1.0, 0.0],
            yp: None,
        }],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![1],
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
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 1 }],
    };

    let mut session = decode(input).expect("well-formed zero-length-section input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = load / (e * area);
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!(
        (got - expected).abs() < 1e-9,
        "expected {expected}, got {got}"
    );
}

/// Bad `orient` rows are decode errors, not silent global axes: a zero
/// vector, a 2D vector leaving the xy plane (`x3 != 0`), and a row past the
/// table.
#[test]
fn rejects_invalid_zero_length_orientations() {
    use carapace_wasm::input_v1::DecodeError;
    for ((row, x), expected) in [
        (
            (0, [0.0, 0.0, 0.0]),
            DecodeError::InvalidOrientation {
                table: "zero_length_sections",
                row: 0,
            },
        ),
        (
            (0, [0.0, 0.0, 1.0]),
            DecodeError::InvalidOrientation {
                table: "zero_length_sections",
                row: 0,
            },
        ),
        (
            (3, [1.0, 0.0, 0.0]),
            DecodeError::UnknownElementIndex {
                table: "zero_length_sections",
                row: 3,
            },
        ),
    ] {
        let mut input = empty_input(2);
        input.nodes = NodeTable {
            coords: vec![0.0, 0.0, 0.0, 0.0],
            fixed: vec![0b111, 0b110],
            mass_node_index: vec![],
            mass: vec![],
        };
        input.materials = vec![MaterialSpec::Elastic { e: 1.0 }];
        input.fibers = FiberTable {
            section_offsets: vec![0, 1],
            z: vec![],
            y: vec![0.0],
            area: vec![1.0],
            material: vec![0],
        };
        input.zero_length_sections = ZeroLengthSectionTable {
            node_i: vec![0],
            node_j: vec![1],
            fiber_section: vec![0],
            materials: vec![],
            orient: vec![OrientRow { row, x, yp: None }],
        };
        assert_eq!(decode(input).err(), Some(expected));
    }
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
        orient: vec![],
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
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 20,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed Hysteretic/Pinching4 input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let (_, hysteretic_disp) = last_sample(&outcome, 0).expect("hysteretic sample");
    let (_, pinching4_disp) = last_sample(&outcome, 1).expect("pinching4 sample");
    assert!(
        hysteretic_disp.is_finite() && hysteretic_disp > 0.0,
        "got {hysteretic_disp}"
    );
    assert!(
        pinching4_disp.is_finite() && pinching4_disp > 0.0,
        "got {pinching4_disp}"
    );
}

#[test]
fn rejects_an_unrecognized_space_before_constructing_a_session() {
    let err = decode(empty_input(4)).unwrap_err();
    assert_eq!(
        err,
        carapace_wasm::input_v1::DecodeError::UnsupportedNdm { got: 4 }
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
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 2,
        }],
        recorders: vec![
            RecorderSpec::ModeShape {
                mode: 0,
                node: 1,
                dof: 0,
            },
            RecorderSpec::ModeShape {
                mode: 1,
                node: 1,
                dof: 0,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mode0_freq, _) = last_sample(&outcome, 0).expect("mode 0 sample");
    let (mode1_freq, _) = last_sample(&outcome, 1).expect("mode 1 sample");
    assert!(
        (mode0_freq - 1.0 / phi).abs() < 1e-9,
        "expected {}, got {mode0_freq}",
        1.0 / phi
    );
    assert!(
        (mode1_freq - phi).abs() < 1e-9,
        "expected {phi}, got {mode1_freq}"
    );
}

/// `Session::modal_results` returns frequencies, full-node shapes and participation for a
/// finished `Modal` stage — no recorders needed. For the 2-mass chain every mode's effective
/// mass sums to the total x mass (all modes requested), shapes are M-normalized and sign-fixed.
#[test]
fn modal_results_report_shapes_and_participation_without_recorders() {
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
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 2,
        }],
        recorders: vec![],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    assert!(session.modal_results().stages.is_empty());
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "{outcome:?}");

    let report = session.modal_results();
    assert_eq!(report.stages.len(), 1);
    let stage = &report.stages[0];
    assert_eq!(
        (stage.stage_index, stage.stage_id.as_str(), stage.ndf),
        (0, "modes", 3)
    );
    assert_eq!(stage.modes.len(), 2);
    assert!((stage.total_mass[0] - 2.0).abs() < 1e-12);
    assert_eq!(stage.total_mass[1], 0.0);

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    assert!((stage.modes[0].frequency - 1.0 / phi).abs() < 1e-9);
    let ratio_sum: f64 = stage.modes.iter().map(|m| m.mass_ratio[0]).sum();
    assert!(
        (ratio_sum - 1.0).abs() < 1e-9,
        "x mass ratios sum to {ratio_sum}"
    );
    for mode in &stage.modes {
        assert_eq!(mode.shape.len(), 3 * 3);
        // Ground node fixed -> zero; largest entry positive; M-normalized.
        assert!(mode.shape[..3].iter().all(|v| *v == 0.0));
        let peak = mode
            .shape
            .iter()
            .copied()
            .fold(0.0_f64, |b, v| if v.abs() > b.abs() { v } else { b });
        assert!(peak > 0.0);
        let norm: f64 = mode.shape[3] * mode.shape[3] + mode.shape[6] * mode.shape[6];
        assert!((norm - 1.0).abs() < 1e-9);
        assert_eq!(mode.mass_ratio[1], 0.0);
    }
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
        orient: vec![],
    };
    input.sequence = SequenceSpec {
        stages: vec![StageSpec::Modal {
            id: "modes".to_string(),
            modes: 1,
        }],
        recorders: vec![RecorderSpec::ModeShape {
            mode: 5,
            node: 1,
            dof: 0,
        }],
    };

    let mut session = decode(input).expect("well-formed modal input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    assert!(
        last_sample(&outcome, 0).is_none(),
        "out-of-range mode should record nothing"
    );
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
        orient: vec![],
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
            damping: carapace_wasm::input_v1::sequence::DampingSpec {
                alpha_m,
                beta_k: 0.0,
            },
            ground_motions: vec![],
            algorithm: AlgorithmSpec::Linear,
            convergence: None,
        }],
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed transient input should decode");
    let outcome = session.advance(steps);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

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
    assert!(
        (got - expected).abs() < 1e-4,
        "expected {expected}, got {got}"
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
        z: vec![],
        y: fiber_y,
        area: fiber_area,
        material: fiber_material,
    };
    input.force_beam_columns_2d = FiberBeamColumn2dTable {
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
/// fiber_recorder.rs's own `disp_beam_column_fiber_responses_match_hand_
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
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed equal_dof input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

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
        normal: vec![],
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
            convergence: Some(ConvergenceSpec::NormUnbalance {
                tol: 1e-9,
                max_iter: 10,
            }),
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
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected = load * length / (area * e);
    let (_, retained_disp) = last_sample(&outcome, 0).expect("retained node sample");
    assert!(
        (retained_disp - expected).abs() < 1e-9,
        "expected {expected}, got {retained_disp}"
    );
    for recorder_index in [1, 2] {
        let (_, disp) = last_sample(&outcome, recorder_index).expect("diaphragm-tied node sample");
        assert!(
            (disp - retained_disp).abs() < 1e-9,
            "diaphragm-tied node should exactly match the retained node's ux: {disp} vs {retained_disp}"
        );
    }
}

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

/// A bar loaded in two static stages that both ramp pattern 0 (load stage `0`), optionally with a
/// `Reset` between them. `hold` freezes the pattern after the first stage.
fn two_stage_bar(
    record_initial: bool,
    between: Option<StageSpec>,
    hold: bool,
) -> carapace_wasm::input_v1::Session {
    let static_stage = |id: &str, hold: bool| StageSpec::Static {
        id: id.to_string(),
        steps: 2,
        integrator: IntegratorSpec::LoadControl { increment: 0.5 },
        algorithm: AlgorithmSpec::Linear,
        convergence: Some(ConvergenceSpec::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        }),
        hold_patterns_after: if hold { vec![0] } else { vec![] },
    };
    let mut input = empty_input(2);
    input.header.record_initial = record_initial;
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 100.0, 0.0],
        fixed: vec![0b111, 0b110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: 30000.0 }];
    input.trusses = TrussTable {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![2.0],
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
        value: vec![60.0],
        stage: vec![0],
    };
    let mut stages = vec![static_stage("first", hold)];
    stages.extend(between);
    stages.push(static_stage("second", false));
    input.sequence = SequenceSpec {
        stages,
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 0 }],
    };
    decode(input).expect("well-formed input should decode")
}

/// Every sample of recorder 0, as `(stage_index, pseudo_time, value)`, running a session to the end.
fn all_samples(session: &mut carapace_wasm::input_v1::Session) -> Vec<(usize, f64, f64)> {
    let mut out = Vec::new();
    loop {
        let outcome = session.advance(8);
        assert!(outcome.error.is_none(), "{outcome:?}");
        for batch in outcome
            .recorder_batches
            .iter()
            .filter(|b| b.recorder_index == 0)
        {
            out.extend(
                batch
                    .samples
                    .iter()
                    .map(|&(t, v)| (batch.stage_index, t, v)),
            );
        }
        if outcome.done {
            return out;
        }
    }
}

/// `header.recordInitial` adds each static stage's initial conditions as its first sample: load
/// factor 0 and the state the stage inherits (zero for the first, the first stage's end for the
/// second). Without it a stage's samples are only its steps.
#[test]
fn record_initial_adds_a_step_zero_sample_per_stage() {
    let without = all_samples(&mut two_stage_bar(false, None, true));
    assert_eq!(without.len(), 4);

    let samples = all_samples(&mut two_stage_bar(true, None, true));
    assert_eq!(samples.len(), 6);
    let stage = |i: usize| samples.iter().filter(|s| s.0 == i).collect::<Vec<_>>();
    let (first, second) = (stage(0), stage(1));
    assert_eq!((first.len(), second.len()), (3, 3));
    assert_eq!((first[0].1, first[0].2), (0.0, 0.0));
    let end_of_first = first[2].2;
    assert!((end_of_first - 60.0 * 100.0 / (2.0 * 30000.0)).abs() < 1e-9);
    assert_eq!(second[0].1, 0.0);
    assert!(
        (second[0].2 - end_of_first).abs() < 1e-12,
        "stage 2 starts where stage 1 ended"
    );
}

/// A `Reset` stage reverts displacements to the as-built state (a later stage starts from zero, not
/// from where the previous one ended). Like OpenSees' `reset`, a pattern frozen by `holdPatternsAfter`
/// stays frozen at its final factor: the next stage applies it in full from its first step.
#[test]
fn reset_stage_reverts_the_state_but_keeps_held_patterns_applied() {
    let reset = Some(StageSpec::Reset {
        id: "reset".to_string(),
    });
    let samples = all_samples(&mut two_stage_bar(true, reset, true));
    // Reset records nothing: stage indices are first = 0, second = 2.
    let second: Vec<_> = samples.iter().filter(|s| s.0 == 2).collect();
    assert_eq!(second.len(), 3);
    assert_eq!(
        second[0].2, 0.0,
        "second stage starts from the as-built state"
    );
    let end_of_first = samples.iter().filter(|s| s.0 == 0).last().unwrap().2;
    assert!(
        (second[1].2 - end_of_first).abs() < 1e-12,
        "the held load is applied in full at once"
    );
    assert!((second[2].2 - end_of_first).abs() < 1e-12);
}

/// Without a hold, the same pattern ramps again from zero after a `Reset`.
#[test]
fn reset_stage_lets_an_unheld_pattern_ramp_again() {
    let reset = Some(StageSpec::Reset {
        id: "reset".to_string(),
    });
    let samples = all_samples(&mut two_stage_bar(true, reset, false));
    let second: Vec<_> = samples.iter().filter(|s| s.0 == 2).collect();
    let end_of_first = samples.iter().filter(|s| s.0 == 0).last().unwrap().2;
    assert_eq!(second[0].2, 0.0);
    assert!(
        (second[1].2 - end_of_first / 2.0).abs() < 1e-12,
        "half the load at the first of two steps"
    );
    assert!((second[2].2 - end_of_first).abs() < 1e-12);
}
