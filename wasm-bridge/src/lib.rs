//! JavaScript interface to the Carapace engine.
//! Binding glue and TypeScript definitions are generated during the wasm build.

pub mod input_v1;
pub mod material_probe;

mod boundary;
pub use boundary::{create_material_probe, decode_input, WasmMaterialProbe, WasmSession};

// Preserve existing milestone exports for verification callers.
pub mod verification;
pub use verification::*;
