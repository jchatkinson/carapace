//! Bulk recorders against the scalar recorders they replace: the same values, sample for sample, over
//! the modal, static, transient and reset stages of 3D frame and shell models compiled by PySees.

use carapace_input::input_v1::sequence::{BulkQuantity, BulkTable, GaussQuantity, RecorderSpec};
use carapace_input::input_v1::tables::ElementKind;
use carapace_input::input_v1::{decode, CarapaceInputV1};

const NDF: usize = 6;

fn fixture(name: &str) -> CarapaceInputV1 {
    let path = format!(
        "{}/tests/fixtures/bulk/{name}.json",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// One bulk recorder and, per value of its frame, the scalar recorder that must agree with it.
struct Check {
    bulk: usize,
    scalars: Vec<usize>,
    /// Velocity and acceleration: scalars record only in transient stages, the bulk reads `0` elsewhere.
    dynamic: bool,
}

fn bulk(
    table: BulkTable,
    quantity: BulkQuantity,
    element_kind: Option<ElementKind>,
    rows: Option<Vec<u32>>,
) -> RecorderSpec {
    RecorderSpec::Bulk {
        table,
        quantity,
        element_kind,
        rows,
    }
}

fn check(name: &str) {
    let mut input = fixture(name);
    let nodes = (input.nodes.coords.len() / 3) as u32;
    let beams = input.elastic_beam_columns_3d.node_i.len() as u32;
    let shells = (input.shell4s.node_ids.len() / 4) as u32;
    let mut recorders = Vec::new();
    let mut checks = Vec::new();
    let mut add = |recorders: &mut Vec<RecorderSpec>,
                   spec: RecorderSpec,
                   scalars: Vec<RecorderSpec>,
                   dynamic: bool| {
        let first = recorders.len();
        recorders.extend(scalars);
        let scalars = (first..recorders.len()).collect();
        recorders.push(spec);
        checks.push(Check {
            bulk: recorders.len() - 1,
            scalars,
            dynamic,
        });
    };
    // Every other node, in reverse: `rows` picks and orders.
    let rows: Vec<u32> = (0..nodes).rev().step_by(2).collect();
    for (quantity, dynamic) in [
        (BulkQuantity::Displacement, false),
        (BulkQuantity::Reaction, false),
        (BulkQuantity::Velocity, true),
        (BulkQuantity::Acceleration, true),
    ] {
        let scalars = rows
            .iter()
            .flat_map(|&node| {
                (0..NDF as u8).map(move |dof| match quantity {
                    BulkQuantity::Displacement => RecorderSpec::NodeDisp { node, dof },
                    BulkQuantity::Reaction => RecorderSpec::Reaction { node, dof },
                    BulkQuantity::Velocity => RecorderSpec::NodeVel { node, dof },
                    _ => RecorderSpec::NodeAccel { node, dof },
                })
            })
            .collect();
        add(
            &mut recorders,
            bulk(BulkTable::Nodes, quantity, None, Some(rows.clone())),
            scalars,
            dynamic,
        );
    }
    for (kind, count, width) in [
        (ElementKind::ElasticBeamColumn3d, beams, 12),
        (ElementKind::Shell4, shells, 24),
    ] {
        if count == 0 {
            continue;
        }
        let rows: Vec<u32> = (0..count).rev().collect();
        let force = rows
            .iter()
            .flat_map(|&element_index| {
                (0..width).map(move |component| RecorderSpec::ElementForce {
                    element_kind: kind,
                    element_index,
                    component,
                })
            })
            .collect();
        add(
            &mut recorders,
            bulk(
                BulkTable::Elements,
                BulkQuantity::EndForce,
                Some(kind),
                Some(rows.clone()),
            ),
            force,
            false,
        );
        let load = rows
            .iter()
            .flat_map(|&element_index| {
                (0..3).map(move |component| RecorderSpec::ElementLoad {
                    element_kind: kind,
                    element_index,
                    component,
                })
            })
            .collect();
        add(
            &mut recorders,
            bulk(
                BulkTable::Elements,
                BulkQuantity::Load,
                Some(kind),
                Some(rows.clone()),
            ),
            load,
            false,
        );
        if kind == ElementKind::Shell4 {
            for (quantity, gauss) in [
                (BulkQuantity::Stress, GaussQuantity::Stress),
                (BulkQuantity::Strain, GaussQuantity::Strain),
            ] {
                let scalars = rows
                    .iter()
                    .flat_map(|&element_index| {
                        (0..4).flat_map(move |point| {
                            (0..8).map(move |component| RecorderSpec::GaussPoint {
                                element_kind: kind,
                                element_index,
                                point,
                                quantity: gauss,
                                component,
                            })
                        })
                    })
                    .collect();
                add(
                    &mut recorders,
                    bulk(
                        BulkTable::Elements,
                        quantity,
                        Some(kind),
                        Some(rows.clone()),
                    ),
                    scalars,
                    false,
                );
            }
        }
    }
    input.sequence.recorders = recorders;
    let kinds: Vec<bool> = input
        .sequence
        .stages
        .iter()
        .map(|s| {
            matches!(
                s,
                carapace_input::input_v1::sequence::StageSpec::Transient { .. }
            )
        })
        .collect();

    let count = input.sequence.recorders.len();
    let mut session = decode(input).expect("fixture decodes");
    let mut scalar: Vec<Vec<(f64, f64)>> = vec![Vec::new(); count];
    let mut frames: Vec<Vec<(f64, bool)>> = vec![Vec::new(); count];
    let mut values: Vec<Vec<f64>> = vec![Vec::new(); count];
    loop {
        let outcome = session.advance(1000);
        assert!(outcome.error.is_none(), "{name}: {:?}", outcome.error);
        for b in &outcome.recorder_batches {
            scalar[b.recorder_index].extend(&b.samples);
        }
        for b in &outcome.bulk_batches {
            frames[b.recorder_index].extend(b.progress.iter().map(|&p| (p, kinds[b.stage_index])));
            values[b.recorder_index].extend(&b.values);
        }
        if outcome.done {
            break;
        }
    }
    for c in &checks {
        let width = c.scalars.len();
        let frames = &frames[c.bulk];
        assert!(
            !frames.is_empty(),
            "{name}: bulk {} recorded nothing",
            c.bulk
        );
        assert_eq!(
            values[c.bulk].len(),
            frames.len() * width,
            "{name}: bulk {} frame width",
            c.bulk
        );
        // Relative to the largest value of the whole recorder: a free DOF's reaction is rounding noise.
        let scale = c
            .scalars
            .iter()
            .flat_map(|&s| &scalar[s])
            .fold(1e-300_f64, |m, &(_, v)| m.max(v.abs()));
        for (k, &s) in c.scalars.iter().enumerate() {
            let mut samples = scalar[s].iter();
            for (f, &(progress, transient)) in frames.iter().enumerate() {
                let value = values[c.bulk][f * width + k];
                if c.dynamic && !transient {
                    assert_eq!(
                        value, 0.0,
                        "{name}: bulk {} reads a static velocity",
                        c.bulk
                    );
                    continue;
                }
                let &(t, expected) = samples
                    .next()
                    .unwrap_or_else(|| panic!("{name}: scalar {s} has fewer samples"));
                assert_eq!(t, progress, "{name}: scalar {s} frame {f}");
                assert!(
                    (value - expected).abs() <= 1e-9 * scale,
                    "{name}: bulk {} value {k} frame {f}: {value} vs {expected}",
                    c.bulk
                );
            }
            assert!(
                samples.next().is_none(),
                "{name}: scalar {s} has more samples"
            );
        }
    }
}

#[test]
fn bulk_recorders_match_scalar_recorders_on_a_3d_frame() {
    check("frame3d");
}

#[test]
fn bulk_recorders_match_scalar_recorders_on_a_shell_tower() {
    check("tower");
}
