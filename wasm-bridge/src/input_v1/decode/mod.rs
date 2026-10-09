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
    match input.header.ndm {
        2 => elements_2d::decode(input).map(Session::D2),
        3 => elements_3d::decode(input).map(Session::D3),
        got => Err(DecodeError::UnsupportedNdm { got }),
    }
}
