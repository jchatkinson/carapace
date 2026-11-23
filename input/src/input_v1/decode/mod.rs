//! Hand-written per-table decode: translates a [`CarapaceInputV1`] into a
//! [`Session`] by calling the same `Domain`/`Element`/`Material`
//! constructors native code calls — no generic deserialization into `core`
//! types, no reflection.
//! Panic-free: every lookup returns a [`DecodeError`] instead of indexing
//! past the end of a table, since this is defense in depth against
//! engine/schema version skew, not a re-check of what `pysees`'s compiler
//! is expected to have already validated.
//!
//! The profile is chosen once from `header.ndm` and never branched on again
//! downstream. `shared` holds the bookkeeping common to both (nodes,
//! constraints, loads, recorders), `stages` the stage compilation, and
//! `elements_2d`/`elements_3d` the two small per-profile element decoders.

mod elements_2d;
mod elements_3d;
mod shared;
mod stages;

use super::error::DecodeError;
use super::session::Session;
use super::CarapaceInputV1;

pub fn decode(input: CarapaceInputV1) -> Result<Session, DecodeError> {
    // Validate all path data before it reaches interpolation, including inertial excitations.
    let check = |series: &super::tables::TimeSeriesSpec| -> Result<(), DecodeError> {
        use super::tables::TimeSeriesSpec;
        if let TimeSeriesSpec::Path { times, factors }
        | TimeSeriesSpec::BoundedPath { times, factors, .. } = series
        {
            if times.is_empty()
                || times.len() != factors.len()
                || times.iter().any(|t| !t.is_finite())
                || factors.iter().any(|v| !v.is_finite())
                || times.windows(2).any(|w| w[1] <= w[0])
            {
                return Err(DecodeError::InvalidAnalysisOption {
                    stage: "timeSeries".into(),
                    field: "path",
                });
            }
        }
        Ok(())
    };
    for series in &input.load_patterns.series {
        check(series)?;
    }
    for stage in &input.sequence.stages {
        if let super::sequence::StageSpec::Transient {
            ground_motions,
            damping,
            steps,
            ..
        } = stage
        {
            if *steps == 0
                || !damping.alpha_m.is_finite()
                || !damping.beta_k.is_finite()
                || damping.alpha_m < 0.0
                || damping.beta_k < 0.0
                || damping.modal_anchors.is_some_and(|a| {
                    a.mode1 == 0
                        || a.mode2 == 0
                        || !a.ratio.is_finite()
                        || !(0.0..=1.0).contains(&a.ratio)
                })
            {
                return Err(DecodeError::InvalidAnalysisOption {
                    stage: stage.id().into(),
                    field: "damping/steps",
                });
            }
            for motion in ground_motions {
                check(&motion.series)?;
                if !motion.scale_factor.is_finite() {
                    return Err(DecodeError::InvalidAnalysisOption {
                        stage: stage.id().into(),
                        field: "groundMotion.scaleFactor",
                    });
                }
            }
        }
    }
    match input.header.ndm {
        2 => elements_2d::decode(input).map(Session::D2),
        3 => elements_3d::decode(input).map(Session::D3),
        got => Err(DecodeError::UnsupportedNdm { got }),
    }
}
