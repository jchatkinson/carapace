#[cfg(all(feature = "parallel", target_arch = "wasm32"))]
compile_error!("the `parallel` feature is native-only and must not be enabled for wasm32");

pub mod analysis;
pub mod model;
pub mod testkit;
