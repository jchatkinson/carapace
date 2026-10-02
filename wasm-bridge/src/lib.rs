//! JavaScript interface to the Carapace engine.
//! Binding glue and TypeScript definitions are generated during the wasm build.

pub mod input_v1;

mod boundary;
pub use boundary::{decode_input, WasmSession};

// Preserve existing milestone exports for verification callers.
pub mod verification;
pub use verification::*;
