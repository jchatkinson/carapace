use crate::model::Domain;

use super::{Algorithm, Analysis, ConstraintHandler, ConvergenceTest, Integrator, SparseSolver};

/// Typestate analysis composition (§2.5): each stage exposes only the next
/// piece that must be wired, and only `AnalysisBuilder<Ready>` exposes
/// `.build()`. Illegal sequencing (e.g. calling `.build()` before an
/// algorithm is set) is a compile error, not a runtime `setLinks` ordering
/// bug like Xara's `BasicAnalysisBuilder`.
pub struct AnalysisBuilder<S> {
    state: S,
}

pub struct Unwired;

pub struct WithConstraintHandler {
    constraint_handler: ConstraintHandler,
}

pub struct WithIntegrator {
    constraint_handler: ConstraintHandler,
    integrator: Integrator,
}

pub struct WithAlgorithm {
    constraint_handler: ConstraintHandler,
    integrator: Integrator,
    algorithm: Algorithm,
}

pub struct Ready {
    constraint_handler: ConstraintHandler,
    integrator: Integrator,
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
    pub fn integrator(self, integrator: Integrator) -> AnalysisBuilder<WithIntegrator> {
        AnalysisBuilder {
            state: WithIntegrator {
                constraint_handler: self.state.constraint_handler,
                integrator,
            },
        }
    }
}

impl AnalysisBuilder<WithIntegrator> {
    pub fn algorithm(self, algorithm: Algorithm) -> AnalysisBuilder<WithAlgorithm> {
        AnalysisBuilder {
            state: WithAlgorithm {
                constraint_handler: self.state.constraint_handler,
                integrator: self.state.integrator,
                algorithm,
            },
        }
    }
}

impl AnalysisBuilder<WithAlgorithm> {
    pub fn test(self, test: ConvergenceTest) -> AnalysisBuilder<Ready> {
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

impl AnalysisBuilder<Ready> {
    /// Numbers DOFs and builds the sparsity/solver setup once, here — not
    /// re-checked every step (§4.4; no live re-solve loop per §1).
    ///
    /// # Panics
    /// If `domain` has any `equal_dof`/`rigid_diaphragm` constraint but
    /// `ConstraintHandler::Plain` was selected — `Plain` can't resolve
    /// multi-point constraints (see its doc comment); use `Transformation`.
    /// This is a model-construction error, not a runtime condition, so it's
    /// caught here rather than threaded through `Result` (§2.8 is about
    /// real runtime failure, not misuse of the builder).
    pub fn build(self, mut domain: Domain) -> Analysis {
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
