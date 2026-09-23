//! Shared helpers for driving a `DisplacementControl` cyclic protocol on a
//! single controlled DOF — used by `tests/m_cyclic_protocol.rs`'s pattern
//! and `examples/cyclic_material.rs` (the openseespy comparison harness in
//! `comparison/`).

use crate::analysis::{Analysis, Integrator};
use crate::model::NodeId;

/// One `DisplacementControl` leg via `Analysis::set_integrator` (no
/// rebuild) targeting an *absolute* displacement — the increment passed to
/// `DisplacementControl` is computed fresh from wherever the controlled DOF
/// actually is right now, not assumed from the previous leg's target.
pub fn leg(analysis: &mut Analysis, node: NodeId, dof: usize, target: f64) -> f64 {
    let current = analysis.domain().node(node).displacement[dof];
    analysis.set_integrator(Integrator::DisplacementControl {
        node,
        dof,
        increment: target - current,
    });
    analysis.step().expect("cyclic protocol leg should converge").load_factor
}

/// Walks the controlled DOF from wherever it is to `target` in steps no
/// larger than `step`, via repeated `leg` calls. Calls `on_point` after
/// every leg with `(displacement, load_factor)` so callers can record a
/// full trace, not just the endpoint.
pub fn ramp_to(
    analysis: &mut Analysis,
    node: NodeId,
    dof: usize,
    target: f64,
    step: f64,
    mut on_point: impl FnMut(f64, f64),
) {
    loop {
        let current = analysis.domain().node(node).displacement[dof];
        let remaining = target - current;
        if remaining.abs() < 1e-12 {
            return;
        }
        let next = current + remaining.clamp(-step, step);
        let load_factor = leg(analysis, node, dof, next);
        on_point(next, load_factor);
    }
}

/// Runs a reverse-cyclic protocol: for each `(amplitude, cycles)` tier, in
/// order, ramps 0 -> +amplitude -> -amplitude, `cycles` times, then back to
/// 0 before the next tier. Calls `on_point` after every leg — this is the
/// full recorded trace a stress-strain (or force-displacement) plot is
/// built from.
pub fn run_cyclic_protocol(
    analysis: &mut Analysis,
    node: NodeId,
    dof: usize,
    tiers: &[(f64, usize)],
    step: f64,
    mut on_point: impl FnMut(f64, f64),
) {
    for &(amplitude, cycles) in tiers {
        for _ in 0..cycles {
            ramp_to(analysis, node, dof, amplitude, step, &mut on_point);
            ramp_to(analysis, node, dof, -amplitude, step, &mut on_point);
        }
        ramp_to(analysis, node, dof, 0.0, step, &mut on_point);
    }
}
