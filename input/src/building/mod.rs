//! Native-first builder for stacked 3D stud-wall buildings. See `docs/building.md`.

mod build;
mod expand;
mod spec;

pub use build::{
    build, BuildError, Building, BuildingWall, HoldDownForce, LevelHandle, StoreyDrift,
};
pub use expand::{
    expand, ExpandOptions, Expansion, ExpansionMap, PierLink, StationNode, StoreyMap, UpAxis,
    WallMap,
};
pub use spec::*;
