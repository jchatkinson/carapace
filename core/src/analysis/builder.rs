use slotmap::Key;

use crate::model::{Domain, ElementOps, NodeId};

use super::{Algorithm, Analysis, ConstraintHandler, ConvergenceTest, Integrator, SparseSolver};

/// Typestate analysis composition (§2.5): each stage exposes only the next
/// piece that must be wired, and only `AnalysisBuilder<Ready>` exposes
/// `.build()`. Illegal sequencing (e.g. calling `.build()` before an
/// algorithm is set) is a compile error, not a runtime `setLinks` ordering
/// bug like Xara's `BasicAnalysisBuilder`.
///
/// `AnalysisBuilder<S>` itself isn't generic over the kinematic profile —
/// only its state markers are (`WithIntegrator<NId>` etc., defaulted to the
/// planar profile), since an `Integrator<NId>` only commits to which node ID
/// type it controls. The rest of the profile (`NDIM`/`NDOF`/`ELEMENT_DOF`/
/// `E::Id`/`E`) is inferred at `.build(domain)` time from `domain` itself.
pub struct AnalysisBuilder<S = Unwired> {
    state: S,
}

pub struct Unwired;

pub struct WithConstraintHandler {
    constraint_handler: ConstraintHandler,
}

pub struct WithIntegrator<NId = NodeId> {
    constraint_handler: ConstraintHandler,
    integrator: Integrator<NId>,
}

pub struct WithAlgorithm<NId = NodeId> {
    constraint_handler: ConstraintHandler,
    integrator: Integrator<NId>,
    algorithm: Algorithm,
}

pub struct Ready<NId = NodeId> {
    constraint_handler: ConstraintHandler,
    integrator: Integrator<NId>,
    algorithm: Algorithm,
    test: ConvergenceTest,
}

impl Default for AnalysisBuilder<Unwired> {
    fn default() -> Self {
        Self::new()
    }
}

impl AnalysisBuilder<Unwired> {
    pub fn new() -> Self {
        AnalysisBuilder { state: Unwired }
    }

    pub fn constraint_handler(self, constraint_handler: ConstraintHandler) -> AnalysisBuilder<WithConstraintHandler> {
        AnalysisBuilder {
            state: WithConstraintHandler { constraint_handler },
        }
    }
}

impl AnalysisBuilder<WithConstraintHandler> {
    pub fn integrator<NId>(self, integrator: Integrator<NId>) -> AnalysisBuilder<WithIntegrator<NId>> {
        AnalysisBuilder {
            state: WithIntegrator {
                constraint_handler: self.state.constraint_handler,
                integrator,
            },
        }
    }
}

impl<NId> AnalysisBuilder<WithIntegrator<NId>> {
    pub fn algorithm(self, algorithm: Algorithm) -> AnalysisBuilder<WithAlgorithm<NId>> {
        AnalysisBuilder {
            state: WithAlgorithm {
                constraint_handler: self.state.constraint_handler,
                integrator: self.state.integrator,
                algorithm,
            },
        }
    }
}

impl<NId> AnalysisBuilder<WithAlgorithm<NId>> {
    pub fn test(self, test: ConvergenceTest) -> AnalysisBuilder<Ready<NId>> {
        AnalysisBuilder {
            state: Ready {
                constraint_handler: self.state.constraint_handler,
                integrator: self.state.integrator,
                algorithm: self.state.algorithm,
                test,
            },
        }
    }
}

impl<NId: Copy> AnalysisBuilder<Ready<NId>> {
    /// Numbers DOFs and builds the sparsity/solver setup once, here — not
    /// re-checked every step (§4.4; no live re-solve loop per §1). Generic
    /// over the rest of the profile (`NDIM`/`NDOF`/`ELEMENT_DOF`/`E::Id`/`E`),
    /// inferred from `domain`'s type — this is what lets one `.build()`
    /// serve both `Domain`/`Analysis` and `Domain3`/`Analysis3`.
    ///
    /// # Panics
    /// If `domain` has any `equal_dof`/`rigid_diaphragm` constraint but
    /// `ConstraintHandler::Plain` was selected — `Plain` can't resolve
    /// multi-point constraints (see its doc comment); use `Transformation`.
    /// This is a model-construction error, not a runtime condition, so it's
    /// caught here rather than threaded through `Result` (§2.8 is about
    /// real runtime failure, not misuse of the builder).
    pub fn build<const NDIM: usize, const NDOF: usize, const ELEMENT_DOF: usize, E>(
        self,
        mut domain: Domain<NDIM, NDOF, ELEMENT_DOF, NId, E>,
    ) -> Analysis<NDIM, NDOF, ELEMENT_DOF, NId, E>
    where
        NId: Key,
        E: ElementOps<NDIM, NDOF, ELEMENT_DOF, NId>,
    {
        assert!(
            !(matches!(self.state.constraint_handler, ConstraintHandler::Plain) && domain.has_mp_constraints()),
            "ConstraintHandler::Plain can't resolve multi-point constraints (equal_dof/rigid_diaphragm) — use ConstraintHandler::Transformation"
        );
        domain.number_dofs();
        Analysis {
            domain,
            constraint_handler: self.state.constraint_handler,
            integrator: self.state.integrator,
            algorithm: self.state.algorithm,
            test: self.state.test,
            solver: SparseSolver::new(),
            step_count: 0,
            load_factor: 0.0,
        }
    }
}
