//! Trapezoidal (linearly varying) element loads on `elasticBeamColumn2d` through the wire format.
//!
//! References are closed forms for a single prismatic member: tip displacements of a cantilever
//! (Euler-Bernoulli, `EI` and `EA` constant) and the fixed-end forces of a fixed-fixed beam under a
//! triangular load, plus consistency checks between equivalent load descriptions. All use a linear
//! algorithm, so tolerances are round-off level (1e-9 relative).
//! Rejection on element kinds that cannot carry a varying load is checked at the end.

use crate::common::{empty_input, last_sample};
use carapace_input::input_v1::decode;
use carapace_input::input_v1::materials::MaterialSpec;
use carapace_input::input_v1::sequence::{
    AlgorithmSpec, ConvergenceSpec, IntegratorSpec, RecorderSpec, SequenceSpec, StageSpec,
};
use carapace_input::input_v1::tables::*;

const L: f64 = 120.0;
const E: f64 = 29_000.0;
const A: f64 = 12.0;
const IZ: f64 = 850.0;

fn trapezoid(wx_i: f64, wy_i: f64, wx_j: f64, wy_j: f64) -> ElementLoadSpec {
    ElementLoadSpec::Trapezoid {
        wx_i,
        wy_i,
        wx_j,
        wy_j,
    }
}

fn uniform(wx: f64, wy: f64) -> ElementLoadSpec {
    ElementLoadSpec::Uniform { wx, wy, wz: None }
}

/// One elastic beam from `(0,0)` to `end`, node i fixed and node j fixed per `fixed_j`, carrying `loads`
/// (all in one pattern). Returns the node-j `ux`, `uy`, `rz` followed by the six element end forces.
fn solve(end: [f64; 2], fixed_j: u8, loads: Vec<ElementLoadSpec>) -> Vec<f64> {
    let mut input = empty_input(2);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, end[0], end[1]],
        fixed: vec![0b111, fixed_j],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.elastic_beam_columns_2d = ElasticBeamColumn2dTable {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![E],
        a: vec![A],
        iz: vec![IZ],
        transform: vec![TransformSpec::Linear],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Linear { slope: 1.0 }],
        scale_factor: vec![1.0],
    };
    let count = loads.len();
    input.element_loads = ElementLoadTable {
        pattern: vec![0; count],
        element_kind: vec![ElementKind::ElasticBeamColumn2d; count],
        element_index: vec![0; count],
        load: loads,
        stage: vec![0; count],
    };
    let mut recorders: Vec<RecorderSpec> = (0..3)
        .map(|dof| RecorderSpec::NodeDisp { node: 1, dof })
        .collect();
    recorders.extend((0..6).map(|component| RecorderSpec::ElementForce {
        element_kind: ElementKind::ElasticBeamColumn2d,
        element_index: 0,
        component,
    }));
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
        recorders,
    };
    let mut session = decode(input).expect("well-formed planar input should decode");
    let outcome = session.advance(10);
    assert!(
        outcome.done && outcome.error.is_none(),
        "unexpected outcome: {outcome:?}"
    );
    (0..9)
        .map(|index| last_sample(&outcome, index).expect("recorder sample").1)
        .collect()
}

fn assert_close(actual: &[f64], expected: &[f64], what: &str) {
    for (index, (a, e)) in actual.iter().zip(expected).enumerate() {
        let scale = e.abs().max(
            actual
                .iter()
                .chain(expected)
                .fold(0.0_f64, |m, v| m.max(v.abs()))
                * 1e-6,
        );
        assert!(
            (a - e).abs() <= 1e-9 * scale.max(1e-12),
            "{what}, value {index}: expected {e}, got {a}"
        );
    }
}

/// A cantilever under a transverse trapezoid: the tip deflection and rotation are the uniform part
/// (`wI`) plus the triangle that grows to `wJ - wI` at the free end, `11 d L^4 / 120 EI` and
/// `d L^3 / 8 EI`; the root carries the resultant and its moment.
#[test]
fn cantilever_transverse_trapezoid_matches_closed_form() {
    let (wi, wj) = (-1.5, -4.0);
    let d = wj - wi;
    let got = solve([L, 0.0], 0b000, vec![trapezoid(0.0, wi, 0.0, wj)]);
    let uy = wi * L.powi(4) / (8.0 * E * IZ) + 11.0 * d * L.powi(4) / (120.0 * E * IZ);
    let rz = wi * L.powi(3) / (6.0 * E * IZ) + d * L.powi(3) / (8.0 * E * IZ);
    let shear = -(wi * L + d * L / 2.0);
    let moment = -(wi * L * L / 2.0 + d * L * L / 3.0);
    assert_close(&got[..3], &[0.0, uy, rz], "tip displacement");
    assert_close(
        &got[3..],
        &[0.0, shear, moment, 0.0, 0.0, 0.0],
        "root shear and moment, free end unloaded",
    );
}

/// A trapezoid that decreases to zero at the free end (the mirrored triangle) uses `d = -w`:
/// tip deflection `w L^4 / 30 EI` for a triangle peaking at the root.
#[test]
fn cantilever_triangle_peaking_at_the_root() {
    let w = -2.0;
    let got = solve([L, 0.0], 0b000, vec![trapezoid(0.0, w, 0.0, 0.0)]);
    assert_close(
        &got[..3],
        &[
            0.0,
            w * L.powi(4) / (30.0 * E * IZ),
            w * L.powi(3) / (24.0 * E * IZ),
        ],
        "tip displacement",
    );
}

/// Axial trapezoid on a cantilever: tip elongation `(wI L^2 / 2 + d L^2 / 3) / EA`, root force the resultant.
#[test]
fn cantilever_axial_trapezoid_matches_closed_form() {
    let (wi, wj) = (3.0, -1.0);
    let d = wj - wi;
    let got = solve([L, 0.0], 0b000, vec![trapezoid(wi, 0.0, wj, 0.0)]);
    let ux = (wi * L * L / 2.0 + d * L * L / 3.0) / (E * A);
    assert_close(&got[..3], &[ux, 0.0, 0.0], "tip displacement");
    assert_close(&got[3..4], &[-(wi * L + d * L / 2.0)], "root axial force");
}

/// A fixed-fixed beam under a triangle that is zero at node i: fixed-end shears `3wL/20`, `7wL/20` and
/// moments `wL^2/30`, `wL^2/20` (the larger at the loaded end), and no displacement.
#[test]
fn fixed_fixed_triangular_load_gives_the_fixed_end_forces() {
    let w = -6.0;
    let got = solve([L, 0.0], 0b111, vec![trapezoid(0.0, 0.0, 0.0, w)]);
    assert_close(&got[..3], &[0.0; 3], "displacement");
    let expected = [
        0.0,
        -3.0 * w * L / 20.0,
        -w * L * L / 30.0,
        0.0,
        -7.0 * w * L / 20.0,
        w * L * L / 20.0,
    ];
    assert_close(&got[3..], &expected, "end forces");
}

/// Reversing the taper mirrors the end forces: peak at node i swaps the shears and the moments.
#[test]
fn reversed_triangle_mirrors_the_fixed_end_forces() {
    let w = -6.0;
    let got = solve([L, 0.0], 0b111, vec![trapezoid(0.0, w, 0.0, 0.0)]);
    let expected = [
        0.0,
        -7.0 * w * L / 20.0,
        -w * L * L / 20.0,
        0.0,
        -3.0 * w * L / 20.0,
        w * L * L / 30.0,
    ];
    assert_close(&got[3..], &expected, "end forces");
}

/// Equal end intensities are a uniform load, and a trapezoid is its uniform part plus a triangle.
#[test]
fn trapezoid_agrees_with_uniform_and_triangle_decompositions() {
    for fixed_j in [0b000, 0b111, 0b011] {
        let a = solve([L, 0.0], fixed_j, vec![trapezoid(1.0, -2.0, 1.0, -2.0)]);
        let b = solve([L, 0.0], fixed_j, vec![uniform(1.0, -2.0)]);
        assert_close(&a, &b, "equal ends vs uniform");

        let c = solve([L, 0.0], fixed_j, vec![trapezoid(0.5, -2.0, -1.0, -5.0)]);
        let d = solve(
            [L, 0.0],
            fixed_j,
            vec![uniform(0.5, -2.0), trapezoid(0.0, 0.0, -1.5, -3.0)],
        );
        assert_close(&c, &d, "trapezoid vs uniform + triangle");

        let zero = solve([L, 0.0], fixed_j, vec![trapezoid(0.0, 0.0, 0.0, 0.0)]);
        assert_close(
            &zero,
            &solve([L, 0.0], fixed_j, vec![]),
            "zero trapezoid vs no load",
        );
    }
}

/// The load is in the member's local axes: a vertical cantilever (local x = global +y, local y = global -x)
/// deflects in global x by the negative of the local closed-form tip deflection.
#[test]
fn trapezoid_follows_the_member_axes() {
    let (wi, wj) = (1.0, 3.0);
    let d = wj - wi;
    let got = solve([0.0, L], 0b000, vec![trapezoid(0.0, wi, 0.0, wj)]);
    let local = wi * L.powi(4) / (8.0 * E * IZ) + 11.0 * d * L.powi(4) / (120.0 * E * IZ);
    assert_close(&got[..2], &[-local, 0.0], "global tip displacement");
}

/// A varying load is only integrated for `elasticBeamColumn2d`: other kinds and the 3D profile reject it
/// instead of silently applying a uniform load or dropping it.
#[test]
fn rejects_a_trapezoid_on_other_element_kinds() {
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
        load: vec![trapezoid(0.0, 1.0, 0.0, 2.0)],
        stage: vec![0],
    };
    assert!(
        decode(input).is_err(),
        "a truss must reject a trapezoid load"
    );

    let mut input = empty_input(3);
    input.nodes = NodeTable {
        coords: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        fixed: vec![0b111111, 0],
        mass_node_index: vec![],
        mass: vec![],
    };
    input.elastic_beam_columns_3d = ElasticBeamColumn3dTable {
        node_i: vec![0],
        node_j: vec![1],
        e: vec![E],
        g: vec![E / 2.6],
        a: vec![A],
        j: vec![1.0],
        iy: vec![IZ],
        iz: vec![IZ],
        transform: vec![TransformSpec3::Linear3 {
            vec_xz: [0.0, 0.0, 1.0],
        }],
        density: vec![0.0],
    };
    input.load_patterns = LoadPatternTable {
        series: vec![TimeSeriesSpec::Constant],
        scale_factor: vec![1.0],
    };
    input.element_loads = ElementLoadTable {
        pattern: vec![0],
        element_kind: vec![ElementKind::ElasticBeamColumn3d],
        element_index: vec![0],
        load: vec![trapezoid(0.0, 1.0, 0.0, 2.0)],
        stage: vec![0],
    };
    assert!(
        decode(input).is_err(),
        "a 3D beam must reject a trapezoid load"
    );
}
