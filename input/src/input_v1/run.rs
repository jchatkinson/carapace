//! Run-to-completion for native callers: drives [`Session::advance`] in budgeted slices, reports
//! progress and honors a cancellation flag between slices, and merges every slice's batches into one
//! [`RunReport`]. The wasm boundary keeps stepping itself (a worker interleaves its own messages).

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};

use serde::Serialize;

use super::session::{
    AnalysisErrorDetail, BulkBatch, LimitSample, ModalResultsReport, Session, StepDiagnostic,
    StepOutcome,
};

/// How [`Session::run`] slices the analysis.
#[derive(Debug, Clone, Copy)]
pub struct RunOptions<'a> {
    /// Steps per `advance` call; a progress callback and the cancellation check run between calls.
    pub step_budget: u32,
    /// Checked between slices; once set the run stops with `cancelled: true`.
    pub cancel: Option<&'a AtomicBool>,
}

impl Default for RunOptions<'_> {
    fn default() -> Self {
        Self {
            step_budget: 64,
            cancel: None,
        }
    }
}

/// One scalar recorder's samples over a stage: `(pseudo_time, value)` pairs.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RecorderSeries {
    pub recorder_index: usize,
    pub stage_index: usize,
    pub samples: Vec<(f64, f64)>,
}

/// Everything a finished (or failed, or cancelled) run produced.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunReport {
    /// Every stage ran to its end with no error and no cancellation.
    pub completed: bool,
    pub cancelled: bool,
    /// The first analysis error, which stops the sequence; results before it are still reported.
    pub error: Option<AnalysisErrorDetail>,
    /// The stage `error` belongs to (or the last stage reached).
    pub stage_id: Option<String>,
    pub total_steps: u64,
    pub recorders: Vec<RecorderSeries>,
    pub bulk: Vec<BulkBatch>,
    pub limit_history: Vec<LimitSample>,
    pub step_diagnostics: Vec<StepDiagnostic>,
    pub modal: ModalResultsReport,
}

#[derive(Default)]
struct Accumulator {
    steps: u64,
    series: BTreeMap<(usize, usize), Vec<(f64, f64)>>,
    bulk: BTreeMap<(usize, usize), BulkBatch>,
    limits: Vec<LimitSample>,
    diagnostics: Vec<StepDiagnostic>,
}

impl Accumulator {
    fn take(&mut self, outcome: StepOutcome) {
        self.steps += u64::from(outcome.steps_taken);
        for batch in outcome.recorder_batches {
            self.series
                .entry((batch.stage_index, batch.recorder_index))
                .or_default()
                .extend(batch.samples);
        }
        for batch in outcome.bulk_batches {
            match self.bulk.entry((batch.stage_index, batch.recorder_index)) {
                std::collections::btree_map::Entry::Occupied(mut slot) => {
                    let merged = slot.get_mut();
                    merged.progress.extend(batch.progress);
                    merged.values.extend(batch.values);
                }
                std::collections::btree_map::Entry::Vacant(slot) => {
                    slot.insert(batch);
                }
            }
        }
        self.limits.extend(outcome.limit_history);
        self.diagnostics.extend(outcome.step_diagnostics);
    }
}

impl Session {
    /// Runs every remaining stage, calling `progress` after each slice with that slice's outcome.
    pub fn run(
        &mut self,
        options: &RunOptions<'_>,
        mut progress: impl FnMut(&StepOutcome),
    ) -> RunReport {
        let mut acc = Accumulator::default();
        let mut error = None;
        let mut cancelled = false;
        let mut done = false;
        let budget = options.step_budget.max(1);
        while !done {
            if options.cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                cancelled = true;
                break;
            }
            let outcome = self.advance(budget);
            progress(&outcome);
            done = outcome.done;
            if outcome.error.is_some() {
                error = outcome.error.clone();
            }
            acc.take(outcome);
        }
        RunReport {
            completed: done && error.is_none(),
            cancelled,
            stage_id: self.current_stage_id().map(str::to_owned),
            error,
            total_steps: acc.steps,
            recorders: acc
                .series
                .into_iter()
                .map(|((stage_index, recorder_index), samples)| RecorderSeries {
                    recorder_index,
                    stage_index,
                    samples,
                })
                .collect(),
            bulk: acc.bulk.into_values().collect(),
            limit_history: acc.limits,
            step_diagnostics: acc.diagnostics,
            modal: self.modal_results(),
        }
    }
}
