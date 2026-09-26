//! Spatial ("space: 3") counterpart to `m10_carapace_input_v1.rs`: exercises
//! `carapace_wasm::input_v1::decode`'s spatial path (`decode3.rs`) end to
//! end, the same way that file exercises the planar path — hand-building a
//! `CarapaceInputV1` with its `*3` tables populated and driving the
//! resulting `Session::Spatial` with `advance`.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec3, SequenceSpec3, StageSpec,
};
use carapace_wasm::input_v1::tables::{
    ElasticBeamColumnTable, ElementLoadTable, EqualDofTable, FiberBeamColumnTable, FiberTable,
    LoadPatternTable, NodalLoadTable, NodeTable, RigidDiaphragmTable, TimeSeriesSpec, TrussTable,
    ZeroLengthSectionTable, ZeroLengthTable,
};
use carapace_wasm::input_v1::tables3::{
    Axis3Spec, ElementKind3, EqualDofTable3, FiberTable3, NodeTable3, RigidDiaphragmTable3,
    TrussTable3, ZeroLengthSectionTable3, ZeroLengthTable3,
};
use carapace_wasm::input_v1::{decode, CarapaceInputV1, Header, StepOutcome};

fn last_sample(outcome: &StepOutcome, recorder_index: usize) -> Option<(f64, f64)> {
    outcome
        .recorder_batches
        .iter()
        .find(|batch| batch.recorder_index == recorder_index)
        .and_then(|batch| batch.samples.last())
        .copied()
}

fn empty_input() -> CarapaceInputV1 {
    CarapaceInputV1 {
        header: Header {
            schema_version: 1,
            space: 3,
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
        sequence: Default::default(),

        nodes3: NodeTable3::default(),
        fibers3: FiberTable3::default(),
        trusses3: TrussTable3::default(),
        elastic_beam_columns3: Default::default(),
        disp_beam_columns3: Default::default(),
        force_beam_columns3: Default::default(),
        zero_lengths3: ZeroLengthTable3::default(),
        zero_length_sections3: ZeroLengthSectionTable3::default(),
        equal_dofs3: EqualDofTable3::default(),
        rigid_diaphragms3: RigidDiaphragmTable3::default(),
        nodal_loads3: NodalLoadTable::default(),
        element_loads3: Default::default(),
        sequence3: SequenceSpec3::default(),
    }
}

/// The M17/M1 closed-form axial-stiffness case (core/tests/m17's skew-truss
/// analogue, core/tests/m1_truss.rs's straight-bar analogue), decoded
/// entirely from `*3` wire tables: a 3-4-5 skew truss should reproduce the
/// same axial elongation a hand computation gives, exercising `Truss3`'s
/// direction-cosine geometry through the decoder rather than just `core`
/// directly.
#[test]
fn decodes_a_skew_truss3_and_matches_the_closed_form_axial_elongation() {
    let (e, area) = (1000.0_f64, 2.0);
    let (dx, dy, length) = (3.0_f64, 4.0, 5.0); // 3-4-5 triangle
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, dx, dy, 0.0],
        fixed: vec![0b111111, 0b111100], // node 1 free only in ux/uy
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.trusses3 = TrussTable3 {
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
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0, 0],
        node: vec![1, 1],
        dof: vec![0, 1], // ux, uy
        value: vec![force * dx / length, force * dy / length],
        stage: vec![0, 0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec3::NodeDisp { node: 1, dof: 0 },
            RecorderSpec3::NodeDisp { node: 1, dof: 1 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected_elongation = force * length / (area * e);
    let (_, ux) = last_sample(&outcome, 0).expect("ux sample");
    let (_, uy) = last_sample(&outcome, 1).expect("uy sample");
    let axial_disp = ux * dx / length + uy * dy / length;
    assert!(
        (axial_disp - expected_elongation).abs() < 1e-9,
        "expected {expected_elongation}, got {axial_disp}"
    );
}

/// `ZeroLength3`'s friction coupling (`ZeroLengthTable3::friction`), the
/// spatial counterpart of `m10_carapace_input_v1.rs`'s
/// `decodes_a_zero_length_with_friction_coupling`: one normal DOF (ux) and
/// *two* independent shear DOFs (uy, uz — `Friction3::shear_dofs`), each run
/// through the Coulomb return map independently against the shared normal
/// force (core/src/model/elements/zero_length.rs's
/// `friction3_shear_axes_slide_independently_against_shared_normal_force`
/// unit test covers the same physics directly against `core`).
#[test]
fn decodes_a_zero_length3_with_friction_coupling_on_two_shear_axes() {
    let (k_normal, mu, k0, b) = (1000.0, 0.3, 500.0, 0.01);
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111111, 0b111000], // node 1 free in ux/uy/uz, rotations fixed
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e: k_normal }];
    input.zero_lengths3 = ZeroLengthTable3 {
        node_i: vec![0],
        node_j: vec![1],
        materials: vec![(0, 0, 0)], // dof 0 (normal) uses material arena index 0
        friction: vec![(0, 0, 1, 2, mu, k0, b)], // normal=ux, shear=[uy, uz]
    };
    input.load_patterns = LoadPatternTable {
        series: vec![
            TimeSeriesSpec::Constant,
            TimeSeriesSpec::Linear { slope: 1.0 },
            TimeSeriesSpec::Linear { slope: 1.0 },
        ],
        scale_factor: vec![1.0, 1.0, 1.0],
    };
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0, 1, 2],
        node: vec![1, 1, 1],
        dof: vec![0, 1, 2],
        value: vec![-1.0, 5.0, 0.1],
        stage: vec![0, 0, 0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec3::NodeDisp { node: 1, dof: 1 },
            RecorderSpec3::NodeDisp { node: 1, dof: 2 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial friction input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    // normal_force = -1000 -> yield_force = 300. Shear axis 0 (uy): trial =
    // 500*5 = 2500, well past yield -> slides. Shear axis 1 (uz): trial =
    // 500*0.1 = 50, well under yield -> sticks (V == k0*shear_rel == 50).
    let (_, uy) = last_sample(&outcome, 0).expect("uy sample");
    let (_, uz) = last_sample(&outcome, 1).expect("uz sample");
    assert!(uy.is_finite() && uy != 0.0, "sliding shear DOF should have displaced: {uy}");
    assert!(uz.is_finite() && uz != 0.0, "sticking shear DOF should have displaced: {uz}");
}

/// `ZeroLengthSection3` (`core::ZeroLengthSection3`, wired via
/// `ZeroLengthSectionTable3`): reproduces closed-form axial stiffness from a
/// four-corner symmetric fiber section, the spatial counterpart of
/// `m10_carapace_input_v1.rs`'s own zero-length-section test.
#[test]
fn decodes_a_zero_length_section3_and_matches_closed_form_axial_stiffness() {
    let (e, area, iy, iz, load): (f64, f64, f64, f64, f64) = (30_000.0, 4.0, 500.0, 2000.0, 50.0);
    let hz = (iy / area).sqrt();
    let hy = (iz / area).sqrt();
    let a4 = area / 4.0;
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        fixed: vec![0b111111, 0b111110], // node 1 free only in ux
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.fibers3 = FiberTable3 {
        section_offsets: vec![0, 4],
        y: vec![hy, hy, -hy, -hy],
        z: vec![hz, -hz, hz, -hz],
        area: vec![a4, a4, a4, a4],
        material: vec![0, 0, 0, 0],
    };
    input.zero_length_sections3 = ZeroLengthSectionTable3 {
        node_i: vec![0],
        node_j: vec![1],
        fiber_section: vec![0],
        materials: vec![],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![load],
        stage: vec![0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![RecorderSpec3::NodeDisp { node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed zero-length-section3 input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected = load / (e * area);
    let (_, got) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((got - expected).abs() < 1e-9, "expected {expected}, got {got}");
}

#[test]
fn decodes_an_elastic_beam_column3_element_force_recorder() {
    let (length, e, g, area, iy, iz, j, fy): (f64, f64, f64, f64, f64, f64, f64, f64) =
        (100.0, 30_000.0, 12_000.0, 10.0, 500.0, 1000.0, 50.0, -10.0);
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0],
        fixed: vec![0b111111, 0b000000],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.elastic_beam_columns3 = carapace_wasm::input_v1::tables3::ElasticBeamColumnTable3 {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![e],
        g: vec![g],
        a: vec![area],
        j: vec![j],
        iy: vec![iy],
        iz: vec![iz],
        transform: vec![carapace_wasm::input_v1::tables3::TransformSpec3::Linear3 {
            vec_xz: [0.0, 0.0, 1.0],
        }],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![1], // uy
        value: vec![fy],
        stage: vec![0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec3::NodeDisp { node: 1, dof: 1 },
            RecorderSpec3::ElementForce {
                element_kind: ElementKind3::ElasticBeamColumn,
                element_index: 0,
                // Local DOF order [ux,uy,uz,rx,ry,rz]_i, [..]_j (width 12):
                // component 5 = rz_i (fixed-end reaction moment about z).
                component: 5,
            },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial elastic beam column input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected_tip_disp = fy * length.powi(3) / (3.0 * e * iz);
    let (_, tip_disp) = last_sample(&outcome, 0).expect("tip disp sample");
    assert!((tip_disp - expected_tip_disp).abs() < 1e-6, "expected {expected_tip_disp}, got {tip_disp}");

    let expected_reaction_moment = -length * fy;
    let (_, reaction_moment) = last_sample(&outcome, 1).expect("reaction moment sample");
    assert!(
        (reaction_moment - expected_reaction_moment).abs() < 1e-6,
        "expected {expected_reaction_moment}, got {reaction_moment}"
    );
}

/// The spatial SDOF-truss closed form (core/tests/m19_spatial_dynamics.rs's
/// `spatial_truss_mass_matches_sdof_closed_form_frequency_via_modal_analysis3`),
/// decoded entirely from `*3` wire tables — proves `Modal`/`ModeShape` work
/// through `decode3.rs`'s spatial path, not just the planar one already
/// checked in `m10_carapace_input_v1.rs`.
#[test]
fn decodes_a_modal_stage3_and_matches_the_sdof_truss_closed_form() {
    let (e, area, length, density) = (30_000.0_f64, 2.0, 100.0, 0.5);
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0],
        fixed: vec![0b111111, 0b111110], // node 1 free only in ux
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.trusses3 = TrussTable3 {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![area],
        material: vec![0],
        density: vec![density],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Modal { id: "modes".to_string(), modes: 1 }],
        recorders: vec![RecorderSpec3::ModeShape { mode: 0, node: 1, dof: 0 }],
    };

    let mut session = decode(input).expect("well-formed spatial modal input should decode");
    let outcome = session.advance(1);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let k = e * area / length;
    let m = density * area * length / 2.0;
    let expected = (k / m).sqrt();
    let (got_freq, _) = last_sample(&outcome, 0).expect("one recorded sample");
    assert!((got_freq - expected).abs() < 1e-9, "expected omega={expected}, got {got_freq}");
}

/// `core::Domain::equal_dof` through the spatial path (`EqualDofTable3`) —
/// same shape as `m10_carapace_input_v1.rs`'s own `equal_dof` test, just
/// against `Node3Id`s: node 2 carries no element of its own and is free
/// only in `ux`, tied to node 1's `ux`.
#[test]
fn decodes_an_equal_dof3_constraint_tying_one_nodes_ux_to_another() {
    let (length, area, e, load) = (100.0, 2.0, 30_000.0, 50.0);
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
        coords: vec![0.0, 0.0, 0.0, length, 0.0, 0.0, 2.0 * length, 0.0, 0.0],
        // node 0 fully fixed; node 1 free only in ux; node 2 (no element of
        // its own) also free only in ux, tied to node 1's via equal_dof.
        fixed: vec![0b111111, 0b111110, 0b111110],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.materials = vec![MaterialSpec::Elastic { e }];
    input.trusses3 = TrussTable3 {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![area],
        material: vec![0],
        density: vec![0.0],
    };
    input.equal_dofs3 = EqualDofTable3 {
        retained: vec![1],
        constrained: vec![2],
        dofs: vec![(0, 0)],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0],
        node: vec![1],
        dof: vec![0],
        value: vec![load],
        stage: vec![0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec3::NodeDisp { node: 1, dof: 0 },
            RecorderSpec3::NodeDisp { node: 2, dof: 0 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial equal_dof input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

    let expected = load * length / (area * e);
    let (_, retained_disp) = last_sample(&outcome, 0).expect("retained node sample");
    let (_, constrained_disp) = last_sample(&outcome, 1).expect("constrained node sample");
    assert!((retained_disp - expected).abs() < 1e-9, "expected {expected}, got {retained_disp}");
    assert!(
        (constrained_disp - retained_disp).abs() < 1e-9,
        "equal_dof-tied node should exactly match the retained node's ux: {constrained_disp} vs {retained_disp}"
    );
}

/// `core::Domain3::rigid_diaphragm_about` through the spatial path
/// (`RigidDiaphragmTable3`): reuses the skew-truss closed form from
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
    let mut input = empty_input();
    input.nodes3 = NodeTable3 {
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
    input.trusses3 = TrussTable3 {
        node_i: vec![0],
        node_j: vec![1],
        area: vec![area],
        material: vec![0],
        density: vec![0.0],
    };
    input.rigid_diaphragms3 = RigidDiaphragmTable3 {
        retained: vec![1],
        normal: vec![Axis3Spec::Y],
        constrained: vec![(0, 2)],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    let force = 100.0;
    input.nodal_loads3 = NodalLoadTable {
        pattern: vec![0, 0],
        node: vec![1, 1],
        dof: vec![0, 2], // ux, uz
        value: vec![force * dx / length, force * dz / length],
        stage: vec![0, 0],
    };
    input.sequence3 = SequenceSpec3 {
        stages: vec![StageSpec::Static {
            id: "only".to_string(),
            steps: 1,
            integrator: IntegratorSpec::LoadControl { increment: 1.0 },
            algorithm: AlgorithmSpec::Linear,
            convergence: Some(ConvergenceSpec::NormUnbalance { tol: 1e-9, max_iter: 10 }),
            hold_patterns_after: vec![],
        }],
        recorders: vec![
            RecorderSpec3::NodeDisp { node: 1, dof: 0 },
            RecorderSpec3::NodeDisp { node: 1, dof: 2 },
            RecorderSpec3::NodeDisp { node: 2, dof: 0 },
            RecorderSpec3::NodeDisp { node: 2, dof: 2 },
        ],
    };

    let mut session = decode(input).expect("well-formed spatial rigid_diaphragm input should decode");
    let outcome = session.advance(10);
    assert!(outcome.done && outcome.error.is_none(), "unexpected outcome: {outcome:?}");

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
