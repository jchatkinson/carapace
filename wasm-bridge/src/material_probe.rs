//! Uniaxial material probe: one persistent unit zero-length spring driven by
//! prescribed displacement (`Analysis::step_prescribed`), the equivalent of
//! OpenSees' `sp`-constrained zeroLength model. Strain is the relative
//! displacement `node2.x - node1.x` and stress is the material's tensile
//! resisting force, `Domain::reaction` at node 2 (OpenSees' `-nodeReaction(1, 1)`).
//! Unit area and gauge length are implicit.

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, AnalysisError, ConstraintHandler, ConvergenceTest,
    Integrator, TangentStrategy,
};
use carapace_core::model::{Domain, Element, Material, Node, NodeId, ZeroLength};
use serde::{Deserialize, Serialize};

use crate::input_v1::materials::{resolve_materials, MaterialSpec};

#[derive(Debug, Clone, Deserialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct MaterialProbeConfig {
    pub material: MaterialSpec,
    #[serde(default)]
    pub initial_strain: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, tsify::Tsify)]
#[serde(rename_all = "camelCase")]
pub struct MaterialProbeResponse {
    pub strain: f64,
    pub stress: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, tsify::Tsify)]
#[serde(tag = "kind", rename_all = "camelCase", rename_all_fields = "camelCase")]
pub enum MaterialProbeError {
    /// The material kind isn't in the probe's supported set.
    UnsupportedMaterial { material: &'static str },
    /// The target strain isn't finite.
    InvalidTarget { target: f64 },
    InvalidProbe { reason: String },
    /// The nonlinear solve failed; the last committed state is intact.
    SolverFailure { target: f64, detail: String },
}

pub struct MaterialProbe {
    material: Material,
    initial_strain: f64,
    analysis: Analysis,
    node: NodeId,
    strain: f64,
}

fn kind_name(spec: &MaterialSpec) -> Option<&'static str> {
    Some(match spec {
        MaterialSpec::Elastic { .. } | MaterialSpec::ElasticPp { .. } | MaterialSpec::Ent { .. } => {
            return None
        }
        MaterialSpec::Steel01 { .. } | MaterialSpec::Concrete01 { .. } => return None,
        MaterialSpec::Gap { .. } => "gap",
        MaterialSpec::Steel02 { .. } => "steel02",
        MaterialSpec::Concrete02 { .. } => "concrete02",
        MaterialSpec::Parallel { .. } => "parallel",
        MaterialSpec::Series { .. } => "series",
        MaterialSpec::MinMax { .. } => "minMax",
        MaterialSpec::Hysteretic { .. } => "hysteretic",
        MaterialSpec::Pinching4 { .. } => "pinching4",
    })
}

fn build(material: &Material) -> (Analysis, NodeId) {
    let mut domain = Domain::new();
    let ground = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let node = domain.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    domain.add_element(Element::ZeroLength(
        ZeroLength::new(ground, node).with_material(0, material.clone()),
    ));
    // No free DOFs: the solve is a trivial equilibrium check that still commits
    // material state, and a zero-tangent branch has nothing to divide by.
    let analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 0.0 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 10,
        })
        .build(domain);
    (analysis, node)
}

impl MaterialProbe {
    pub fn new(config: &MaterialProbeConfig) -> Result<Self, MaterialProbeError> {
        if let Some(material) = kind_name(&config.material) {
            return Err(MaterialProbeError::UnsupportedMaterial { material });
        }
        let initial_strain = config.initial_strain.unwrap_or(0.0);
        if !initial_strain.is_finite() {
            return Err(MaterialProbeError::InvalidTarget {
                target: initial_strain,
            });
        }
        let material = resolve_materials(std::slice::from_ref(&config.material))
            .map_err(|error| MaterialProbeError::InvalidProbe {
                reason: format!("{error:?}"),
            })?
            .remove(0);
        let mut probe = Self::fresh(material, initial_strain);
        probe.apply_strain(initial_strain)?;
        Ok(probe)
    }

    fn fresh(material: Material, initial_strain: f64) -> Self {
        let (analysis, node) = build(&material);
        MaterialProbe {
            material,
            initial_strain,
            analysis,
            node,
            strain: 0.0,
        }
    }

    /// Imposes `target` and returns the response. A repeat of the current
    /// strain reports the committed response without stepping, so it never
    /// advances material state. A failed target leaves the last committed
    /// state intact.
    pub fn apply_strain(&mut self, target: f64) -> Result<MaterialProbeResponse, MaterialProbeError> {
        if !target.is_finite() {
            return Err(MaterialProbeError::InvalidTarget { target });
        }
        if target != self.strain {
            self.analysis
                .step_prescribed(&[(self.node, 0, target)])
                .map_err(|error: AnalysisError| MaterialProbeError::SolverFailure {
                    target,
                    detail: format!("{error:?}"),
                })?;
            self.strain = target;
        }
        Ok(MaterialProbeResponse {
            strain: self.strain,
            stress: self.analysis.domain().reaction(self.node, 0, 0.0),
        })
    }

    /// Back to the freshly constructed state: committed history discarded.
    pub fn reset(&mut self) -> Result<(), MaterialProbeError> {
        *self = Self::fresh(self.material.clone(), self.initial_strain);
        self.apply_strain(self.initial_strain).map(|_| ())
    }
}
