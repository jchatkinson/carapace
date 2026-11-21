//! Runs independent models concurrently (parameter sweeps, load cases, many drawings): one analysis per
//! task, so results are exactly what one-at-a-time runs give and do not depend on the thread count.

use rayon::prelude::*;

use carapace_core::analysis::{parallelism, set_parallelism, Parallelism};

use super::{decode, CarapaceInputV1, DecodeError, RunOptions, RunReport};

/// Decodes and runs every input on up to `threads` threads (`0` = all cores), returning one entry per
/// input in input order.
///
/// The solver's own threading is switched off for the duration (two levels of parallelism would
/// oversubscribe the machine) and restored afterwards; it is a process-wide setting, so do not change it
/// from another thread meanwhile.
pub fn run_many(
    inputs: Vec<CarapaceInputV1>,
    threads: usize,
    options: &RunOptions<'_>,
) -> Vec<Result<RunReport, DecodeError>> {
    let previous = parallelism();
    set_parallelism(Parallelism::Serial);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("thread pool");
    let results = pool.install(|| {
        inputs
            .into_par_iter()
            .map(|input| decode(input).map(|mut session| session.run(options, |_| {})))
            .collect()
    });
    set_parallelism(previous);
    results
}
