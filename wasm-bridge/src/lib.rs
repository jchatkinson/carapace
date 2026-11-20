//! JavaScript interface to the Carapace engine.
//! Binding glue and TypeScript definitions are generated during the wasm build.
//! The model layer itself (wire format, decode, sessions, probes) lives in `carapace-input`.

mod boundary;
pub use boundary::{
    create_material_probe, create_section_probe, decode_input, WasmMaterialProbe, WasmSectionProbe,
    WasmSession,
};
