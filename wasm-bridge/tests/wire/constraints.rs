//! `equalDofs` and `rigidDiaphragms` tables through the wire format, in 2D and 3D.

use crate::common::{empty_input, last_sample};
use carapace_wasm::input_v1::decode;
use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_wasm::input_v1::tables::*;

/// `core::Domain::equal_dof` (`EqualDofTable`, `ConstraintHandler::
/// Transformation` picked automatically by `session::start_current_stage`
/// once `Domain::has_mp_constraints()` is true): node 2 carries no element
/// of its own and is free only in `ux`, tied via `equal_dof` to node 1's
/// `ux` — so its recorded displacement must exactly equal node 1's, which
/// itself matches the plain closed-form truss elongation.
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

/// `core::Domain::equal_dof` through the spatial path (`EqualDofTable`) —
/// same shape as `wire/elements.rs`'s own `equal_dof` test, just
/// against `Node3Id`s: node 2 carries no element of its own and is free
/// only in `ux`, tied to node 1's `ux`.
#[test]
fn decodes_an_equal_dof3_constraint_tying_one_nodes_ux_to_another() {
    let (length, area, e, load) = (100.0, 2.0, 30_000.0, 50.0);
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0, 2.0 * length, 0.0, 0.0],
        // node 0 fully fixed; node 1 free only in ux; node 2 (no element of
        // its own) also free only in ux, tied to node 1's via equal_dof.
        fixed: vec![0b111111, 0b111110, 0b111110],
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

    let mut session = decode(input).expect("well-formed spatial equal_dof input should decode");
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

/// `core::Domain3::rigid_diaphragm_about` through the spatial path
/// (`RigidDiaphragmTable`): reuses the skew-truss closed form from
/// `decodes_a_skew_truss3_and_matches_the_closed_form_axial_elongation`
/// to give the retained node a known, nonzero translation in *both*
/// in-plane axes (`x`/`z`, perpendicular to `normal = Y`) while its
/// rotation about `normal` stays fixed at zero — so the lever-arm term in
/// `Domain3::rigid_diaphragm_about`'s affine tie vanishes and the
/// constrained node (which carries no element of its own) must reproduce
/// the retained node's `ux`/`uz` exactly.
#[test]
fn decodes_a_rigid_diaphragm3_and_ties_translation_with_no_lever_arm_when_untwisted() {
    let (e, area) = (1000.0_f64, 2.0);
    let (dx, dz, length) = (3.0_f64, 4.0, 5.0); // 3-4-5 triangle, in the x-z plane
    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, dx, 0.0, dz, 2.0 * dx, 0.0, 2.0 * dz],
        // node 0 fully fixed; node 1 free only in ux/uz (rotation about
        // normal=Y stays fixed at zero, so the diaphragm's lever-arm term
        // is exactly zero); node 2 (no element of its own) is the
        // rigid-diaphragm-constrained node, also free only in ux/uz.
        fixed: vec![0b111111, 0b111010, 0b111010],
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
        normal: vec![Axis3Spec::Y],
        constrained: vec![(0, 2)],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    let force = 100.0;
    input.nodal_loads = NodalLoadTable {
        pattern: vec![0, 0],
        node: vec![1, 1],
        dof: vec![0, 2], // ux, uz
        value: vec![force * dx / length, force * dz / length],
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
            RecorderSpec::NodeDisp { node: 1, dof: 2 },
            RecorderSpec::NodeDisp { node: 2, dof: 0 },
            RecorderSpec::NodeDisp { node: 2, dof: 2 },
        ],
    };

    let mut session =
        decode(input).expect("well-formed spatial rigid_diaphragm input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );

    let (_, retained_ux) = last_sample(&outcome, 0).expect("retained ux sample");
    let (_, retained_uz) = last_sample(&outcome, 1).expect("retained uz sample");
    let (_, constrained_ux) = last_sample(&outcome, 2).expect("constrained ux sample");
    let (_, constrained_uz) = last_sample(&outcome, 3).expect("constrained uz sample");

    let axial_disp = retained_ux * dx / length + retained_uz * dz / length;
    let expected_elongation = force * length / (area * e);
    assert!(
        (axial_disp - expected_elongation).abs() < 1e-9,
        "expected {expected_elongation}, got {axial_disp}"
    );
    assert!(
        (constrained_ux - retained_ux).abs() < 1e-9,
        "diaphragm-tied node should exactly match the retained node's ux with no rotation: {constrained_ux} vs {retained_ux}"
    );
    assert!(
        (constrained_uz - retained_uz).abs() < 1e-9,
        "diaphragm-tied node should exactly match the retained node's uz with no rotation: {constrained_uz} vs {retained_uz}"
    );
}
