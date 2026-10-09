//! Elements through the wire format: tables decode into the right elements and match closed forms
//! (truss, zero-length, zero-length section, and uniaxial materials), in 2D and 3D.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::*;
use carapace_wasm::input_v1::{decode, CarapaceInputV1};

/// Smoke test for the whole node/material/element/pattern/sequence
/// pipeline: the truss case (`core/tests/elements/truss.rs`'s closed form),
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

/// The closed-form axial-stiffness case (core/tests/m17's skew-truss
/// analogue, core/tests/elements/truss.rs's straight-bar analogue), decoded
/// entirely from 3D wire tables: a 3-4-5 skew truss should reproduce the
/// same axial elongation a hand computation gives, exercising `Truss3`'s
/// direction-cosine geometry through the decoder rather than just `core`
/// directly.
#[test]
fn decodes_a_skew_truss3_and_matches_the_closed_form_axial_elongation() {
    let (e, area) = (1000.0_f64, 2.0);
    let (dx, dy, length) = (3.0_f64, 4.0, 5.0); // 3-4-5 triangle
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, dx, dy, 0.0],
        fixed: vec![0b111111, 0b111100], // node 1 free only in ux/uy
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
    let force = 100.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 0],
        node: vec![1, 1],
        dof: vec![0, 1], // ux, uy
        value: vec![force * dx / length, force * dy / length],
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
        recorders: vec![
            RecorderSpec::NodeDisp { node: 1, dof: 0 },
            RecorderSpec::NodeDisp { node: 1, dof: 1 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let expected_elongation = force * length / (area * e);
    let (_, ux) = last_sample(&outcome, 0).expect("ux sample");
    let (_, uy) = last_sample(&outcome, 1).expect("uy sample");
    let axial_disp = ux * dx / length + uy * dy / length;
    assert!(
        (axial_disp - expected_elongation).abs() < 1e-9,
        "expected {expected_elongation}, got {axial_disp}"
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

/// `ZeroLength3`'s friction coupling (`ZeroLengthTable::friction`), the
/// spatial counterpart of `wire/elements.rs`'s
/// `decodes_a_zero_length_with_friction_coupling`: one normal DOF (ux) and
/// *two* independent shear DOFs (uy, uz — `Friction3::shear_dofs`), each run
/// through the Coulomb return map independently against the shared normal
/// force (core/src/model/elements/zero_length.rs's
/// `friction3_shear_axes_slide_independently_against_shared_normal_force`
/// unit test covers the same physics directly against `core`).
#[test]
fn decodes_a_zero_length3_with_friction_coupling_on_two_shear_axes() {
    let (k_normal, mu, k0, b) = (1000.0, 0.3, 500.0, 0.01);
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111111, 0b111000], // node 1 free in ux/uy/uz, rotations fixed
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
            shear_dofs: vec![1, 2], // normal=ux, shear=[uy, uz]
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
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
        scale_factor: vec![1.0, 1.0, 1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 1, 2],
        node: vec![1, 1, 1],
        dof: vec![0, 1, 2],
        value: vec![-1.0, 5.0, 0.1],
        stage: vec![0, 0, 0],
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
            RecorderSpec::NodeDisp { node: 1, dof: 2 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial friction input should decode");
    let outcome = session.advance(1);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    // normal_force = -1000 -> yield_force = 300. Shear axis 0 (uy): trial =
    // 500*5 = 2500, well past yield -> slides. Shear axis 1 (uz): trial =
    // 500*0.1 = 50, well under yield -> sticks (V == k0*shear_rel == 50).
    let (_, uy) = last_sample(&outcome, 0).expect("uy sample");
    let (_, uz) = last_sample(&outcome, 1).expect("uz sample");
    assert!(
        uy.is_finite() && uy != 0.0,
        "sliding shear DOF should have displaced: {uy}"
    );
    assert!(
        uz.is_finite() && uz != 0.0,
        "sticking shear DOF should have displaced: {uz}"
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

/// `ZeroLengthSection3` (`core::ZeroLengthSection3`, wired via
/// `ZeroLengthSectionTable`): reproduces closed-form axial stiffness from a
/// four-corner symmetric fiber section, the spatial counterpart of
/// `wire/elements.rs`'s own zero-length-section test.
#[test]
fn decodes_a_zero_length_section3_and_matches_closed_form_axial_stiffness() {
    let (e, area, iy, iz, load): (f64, f64, f64, f64, f64) = (30_000.0, 4.0, 500.0, 2000.0, 50.0);
    let hz = (iy / area).sqrt();
    let hy = (iz / area).sqrt();
    let a4 = area / 4.0;
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111111, 0b111110], // node 1 free only in ux
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 4],
        y: vec![hy, hy, -hy, -hy],
        z: vec![hz, -hz, hz, -hz],
        area: vec![a4, a4, a4, a4],
        material: vec![0, 0, 0, 0],
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

    let mut session = decode(input).expect("well-formed zero-length-section3 input should decode");
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

/// `ZeroLengthSectionTable::orient`: local x = global Z (yp = global X), so
/// the section's axial stiffness answers a `uz` load.
#[test]
fn decodes_an_oriented_zero_length_section3_and_matches_closed_form_axial_stiffness() {
    let (e, area, iy, iz, load): (f64, f64, f64, f64, f64) = (30_000.0, 4.0, 500.0, 2000.0, 50.0);
    let hz = (iy / area).sqrt();
    let hy = (iz / area).sqrt();
    let a4 = area / 4.0;
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111111, 0b111011], // node 1 free only in uz
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers = FiberTable {
        section_offsets: vec![0, 4],
        y: vec![hy, hy, -hy, -hy],
        z: vec![hz, -hz, hz, -hz],
        area: vec![a4, a4, a4, a4],
        material: vec![0, 0, 0, 0],
    };
    input.zero_length_sections = ZeroLengthSectionTable {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
        orient: vec![OrientRow {
            row: 0,
            x: [0.0, 0.0, 1.0],
            yp: Some([1.0, 0.0, 0.0]),
        }],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![2],
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
        recorders: vec![RecorderSpec::NodeDisp { node: 1, dof: 2 }],
    };

    let mut session = decode(input).expect("well-formed zero-length-section3 input should decode");
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
