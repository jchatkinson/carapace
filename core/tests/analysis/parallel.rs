//! The solver thread count (`set_parallelism`, `parallel` feature) must not change results: a frame large
//! enough to take the supernodal factorization path is solved serially and on 4 threads, and every nodal
//! displacement must be bit-identical.

use carapace_core::analysis::{
    parallelism, set_parallelism, Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest,
    Integrator, Parallelism,
};
use carapace_core::model::{Domain, ElasticBeamColumn, Element, GeomTransf, Node};

/// Displacements of a `bays x stories` rectangular frame under a lateral roof-line load.
fn solve_frame(bays: usize, stories: usize) -> Vec<f64> {
    let mut domain = Domain::new();
    let mut ids = Vec::new();
    for s in 0..=stories {
        for b in 0..=bays {
            let mut node = Node::new([b as f64 * 6000.0, s as f64 * 3000.0]);
            if s == 0 {
                node = node.fix(0).fix(1).fix(2);
            }
            ids.push(domain.add_node(node));
        }
    }
    let at = |s: usize, b: usize| ids[s * (bays + 1) + b];
    let beam = |i, j| {
        Element::ElasticBeamColumn(ElasticBeamColumn::new(
            i,
            j,
            200000.0,
            8000.0,
            2.0e8,
            GeomTransf::Linear,
        ))
    };
    for s in 0..=stories {
        for b in 0..=bays {
            if s < stories {
                domain.add_element(beam(at(s, b), at(s + 1, b)));
            }
            if s > 0 && b < bays {
                domain.add_element(beam(at(s, b), at(s, b + 1)));
            }
        }
    }
    for s in 1..=stories {
        domain.load_node(at(s, 0), 0, 10.0 * s as f64);
    }
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-6,
            max_iter: 10,
        })
        .build(domain);
    analysis.step().expect("frame should solve");
    ids.iter()
        .flat_map(|&id| analysis.domain().node(id).displacement)
        .collect()
}

#[test]
fn thread_count_does_not_change_the_solution() {
    assert!(set_parallelism(Parallelism::Serial));
    let serial = solve_frame(20, 40);
    assert!(set_parallelism(Parallelism::Threads(4)));
    assert_eq!(parallelism(), Parallelism::Threads(4));
    let threaded = solve_frame(20, 40);
    set_parallelism(Parallelism::Serial);
    assert!(serial.iter().any(|&v| v != 0.0));
    assert_eq!(serial, threaded);
}

#[allow(dead_code)]
fn element_and_domain_are_thread_safe() {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<Domain>();
    send_sync::<carapace_core::model::Domain3>();
}
