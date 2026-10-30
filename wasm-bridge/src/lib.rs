//! JavaScript interface to the Carapace engine.
//! Binding glue and TypeScript definitions are generated during the wasm build.

pub mod input_v1;
pub mod material_probe;
pub mod section_probe;

mod boundary;
pub use boundary::{
    create_material_probe, create_section_probe, decode_input, WasmMaterialProbe, WasmSectionProbe,
    WasmSession,
};
