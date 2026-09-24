//! `DecodeError` — decode's own diagnostic type, in `AnalysisError`'s style
//! (implementation-plan.md §2.8): tagged variants carrying context, not
//! sentinel codes. This is defense in depth against what `pysees`'s
//! compiler can't see (engine/schema version skew, a stale snapshot run
//! against a newer/older wasm build) — the compiler is expected to catch
//! entity-level problems (missing tags, unsupported materials, ...) first.

#[derive(Debug, Clone, PartialEq)]
pub enum DecodeError {
    /// `header.space` is neither the planar profile this decoder
    /// implements nor (yet) any other recognized value.
    UnsupportedSpace {
        got: u8,
    },
    /// A well-formed stage kind this decoder doesn't implement yet
    /// (`modal`/`transient` today).
    UnsupportedStage {
        stage_id: String,
        kind: &'static str,
    },
    UnknownNodeIndex {
        table: &'static str,
        row: u32,
    },
    UnknownMaterialIndex {
        table: &'static str,
        row: u32,
    },
    CyclicMaterialReference {
        index: u32,
    },
    UnknownPatternIndex {
        table: &'static str,
        row: u32,
    },
    UnknownFiberSectionIndex {
        row: u32,
    },
    UnknownElementIndex {
        kind: &'static str,
        row: u32,
    },
    UnknownStageIndex {
        table: &'static str,
        row: u32,
    },
    /// An element load referenced an element kind `core` doesn't apply
    /// that load to (today, only `ElasticBeamColumn` honors
    /// `UniformTransverse`).
    UnsupportedElementLoad {
        element_kind: &'static str,
    },
    InvalidDof {
        table: &'static str,
        row: u32,
        dof: u8,
    },
}
