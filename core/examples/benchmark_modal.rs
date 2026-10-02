//! Modal-analysis benchmark on an elastic planar frame.
//! Run with `cargo run --release -p carapace-core --example benchmark_modal`.
//! Optional positional arguments: stories, bays, modes, repetitions.

use std::time::Instant;

use carapace_core::analysis::modal_analysis;
use carapace_core::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Node};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let arg = |index: usize, default: usize| {
        args.get(index)
            .map(|value| value.parse::<usize>().expect("arguments must be integers"))
            .unwrap_or(default)
    };
    let (stories, bays, num_modes, repetitions) = (arg(0, 10), arg(1, 3), arg(2, 6), arg(3, 21));
    let free_dofs = stories * (bays + 1) * 3;
    assert!(stories > 0 && bays > 0 && repetitions > 0);
    assert!(num_modes > 0 && num_modes <= free_dofs);

    let mut domain = Domain::new();
    let mut nodes = Vec::new();
    for story in 0..=stories {
        let mut row = Vec::new();
        for bay in 0..=bays {
            let mut node = Node::new([bay as f64 * 6000.0, story as f64 * 3000.0]);
            if story == 0 {
                node = node.fix(0).fix(1).fix(2);
            } else {
                node = node.with_mass(0, 1.0).with_mass(1, 1.0).with_mass(2, 1.0e6);
            }
            row.push(domain.add_node(node));
        }
        nodes.push(row);
    }
    for story in 1..=stories {
        for bay in 0..=bays {
            domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                nodes[story - 1][bay],
                nodes[story][bay],
                200000.0,
                8000.0,
                2.0e8,
                GeomTransf::Linear,
            )));
        }
        for bay in 0..bays {
            domain.add_element(Element::ElasticBeamColumn(ElasticBeamColumn::new(
                nodes[story][bay],
                nodes[story][bay + 1],
                200000.0,
                6000.0,
                3.0e8,
                GeomTransf::Linear,
            )));
        }
    }

    // Warm up without including model construction in the measured time.
    for _ in 0..3 {
        std::hint::black_box(modal_analysis(&mut domain, num_modes).expect("modal solve failed"));
    }
    let mut times_ms = Vec::with_capacity(repetitions);
    let mut frequencies = Vec::new();
    for _ in 0..repetitions {
        let start = Instant::now();
        let modes = modal_analysis(&mut domain, num_modes).expect("modal solve failed");
        times_ms.push(start.elapsed().as_secs_f64() * 1000.0);
        frequencies = modes.iter().map(|mode| mode.frequency).collect();
        std::hint::black_box(modes);
    }
    times_ms.sort_by(f64::total_cmp);
    let median_ms = if repetitions % 2 == 0 {
        (times_ms[repetitions / 2 - 1] + times_ms[repetitions / 2]) / 2.0
    } else {
        times_ms[repetitions / 2]
    };
    println!(
        "{{\"stories\":{stories},\"bays\":{bays},\"free_dofs\":{free_dofs},\"modes\":{num_modes},\"repetitions\":{repetitions},\"median_ms\":{median_ms},\"min_ms\":{},\"max_ms\":{},\"frequencies\":{frequencies:?}}}",
        times_ms[0], times_ms[repetitions - 1],
    );
}
