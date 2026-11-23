//! Native-first builder for stacked 3D stud-wall buildings. See `docs/building.md`.

mod build;
mod spec;

pub use build::{
    build, BuildError, Building, BuildingWall, HoldDownForce, LevelHandle, StoreyDrift,
};
pub use spec::*;
