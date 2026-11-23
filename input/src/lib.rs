//! Model-from-data layer of the Carapace engine: the `CarapaceInputV1` wire format, its decoder, the
//! steppable `Session`, recorders and result streaming, and the section and material probes.
//! Native-first and free of wasm dependencies; the `tsify` feature only adds TypeScript definition
//! derives for `carapace-wasm`.

pub mod building;
pub mod input_v1;
pub mod material_probe;
pub mod section_probe;
