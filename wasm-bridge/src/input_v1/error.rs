//! `DecodeError` — decode's own diagnostic type, in `AnalysisError`'s style
//! (implementation-plan.md §2.8): tagged variants carrying context, not
//! sentinel codes. This is defense in depth against what `pysees`'s
//! compiler can't see (engine/schema version skew, a stale snapshot run
//! against a newer/older wasm build) — the compiler is expected to catch
//! entity-level problems (missing tags, unsupported materials, ...) first.

use carapace_core::model::ModelError;
use serde::Serialize;

/// `Serialize`, not `Deserialize` — a `DecodeError` only ever flows *out*
/// to JS (`boundary.rs`), as a `{ kind: "...", ... }`-shaped object.
#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum DecodeError {
    /// An analysis parameter is non-finite or outside its valid range.
    InvalidAnalysisOption {
        stage: String,
        field: &'static str,
    },
    /// `header.space` is neither the planar profile this decoder
    /// implements nor (yet) any other recognized value.
    UnsupportedSpace {
        got: u8,
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
        table: &'static str,
        row: u32,
    },
    UnknownStageIndex {
        table: &'static str,
        row: u32,
    },
    /// An element load referenced an element kind `core` doesn't apply
    /// that load to (today, only the beam-column kinds honor
    /// `Uniform`).
    UnsupportedElementLoad {
        element_kind: &'static str,
    },
    InvalidDof {
        table: &'static str,
        row: u32,
        dof: u8,
    },
    /// A sparse per-row entry (`EqualDofTable::dofs`/`RigidDiaphragmTable::
    /// constrained` and their spatial counterparts) named a row past the
    /// end of that table's own dense `retained` list.
    UnknownConstraintRow {
        table: &'static str,
        row: u32,
    },
    /// A zero-length element's `orient` vectors are zero or parallel, or (for
    /// a 2D element) leave the xy plane, so no local frame exists.
    InvalidOrientation {
        table: &'static str,
        row: u32,
    },
    /// The assembled model failed `core`'s `Domain::validate` (a nodal load
    /// on a DOF nothing uses, an element rejecting its own geometry, ...).
    InvalidModel {
        error: ModelErrorDetail,
    },
}

/// `core::ModelError`, restated for the wire. Node and element indices are
/// the rows of the node table and of the element tables in insertion order
/// (all element kinds in decode order), so a caller can map them back.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ModelErrorDetail {
    LoadOnInactiveDof {
        node: usize,
        dof: usize,
    },
    DuplicateSlave {
        node: usize,
        dof: usize,
    },
    SlaveIsFixed {
        node: usize,
        dof: usize,
    },
    ConstraintCycle {
        node: usize,
        dof: usize,
    },
    ConstraintOnPrescribedDof {
        node: usize,
        dof: usize,
    },
    InconsistentInitialState {
        node: usize,
        dof: usize,
    },
    MassOnConstrainedDof {
        node: usize,
        dof: usize,
    },
    InvalidElement {
        element: usize,
        reason: &'static str,
    },
}

impl From<ModelError> for ModelErrorDetail {
    fn from(error: ModelError) -> Self {
        match error {
            ModelError::LoadOnInactiveDof { node, dof } => {
                ModelErrorDetail::LoadOnInactiveDof { node, dof }
            }
            ModelError::DuplicateSlave { node, dof } => {
                ModelErrorDetail::DuplicateSlave { node, dof }
            }
            ModelError::SlaveIsFixed { node, dof } => ModelErrorDetail::SlaveIsFixed { node, dof },
            ModelError::ConstraintCycle { node, dof } => {
                ModelErrorDetail::ConstraintCycle { node, dof }
            }
            ModelError::ConstraintOnPrescribedDof { node, dof } => {
                ModelErrorDetail::ConstraintOnPrescribedDof { node, dof }
            }
            ModelError::InconsistentInitialState { node, dof } => {
                ModelErrorDetail::InconsistentInitialState { node, dof }
            }
            ModelError::MassOnConstrainedDof { node, dof } => {
                ModelErrorDetail::MassOnConstrainedDof { node, dof }
            }
            ModelError::InvalidElement { element, reason } => {
                ModelErrorDetail::InvalidElement { element, reason }
            }
        }
    }
}
