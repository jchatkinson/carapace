//! Stage compilation: `StageSpec` -> `CompiledStage` for all three stage kinds
//! (integrator/algorithm/convergence conversion, held-pattern resolution).
//! Never touches element physics, only `LoadPatternId`s and a caller-supplied
//! `node_at`/`ndof`, so it is generic over the node/element ID types and the
//! element-load type.

use carapace_core::analysis::{
    Algorithm, AnalysisError, ArcDirection, ArcLength, ArcPredictor, ArcScales, ArcSeed,
    Backtracking, ConvergenceTest, DisplacementTarget, GroundMotion, Integrator, LineSearch,
    LoadFactorTarget, RayleighDamping, TangentStrategy,
};
use carapace_core::model::{LoadPatternId, LoadSeries};

use super::shared::check_dof_within;
use crate::input_v1::error::DecodeError;
use crate::input_v1::sequence::{
    AlgorithmConfigSpec, AlgorithmSpec, ArcDirectionSpec, ArcLengthSpec, ArcPredictorSpec,
    ArcScalesSpec, ConvergenceSpec, IntegratorSpec, LegacyAlgorithmSpec, LineSearchSpec, StageSpec,
    TangentStrategySpec,
};
use crate::input_v1::session::{CompiledStage, CompiledStageKind};
use crate::input_v1::tables::TimeSeriesSpec;

/// Used by load patterns and by ground-motion decoding
/// (`compile_stages`): an accelerogram is exactly a `TimeSeriesSpec::Path`,
/// the same wire type a static `LoadPatternTable` entry uses (`core::
/// GroundMotion`'s own doc comment), so this one conversion serves both.
pub(super) fn load_series_of(spec: &TimeSeriesSpec) -> LoadSeries {
    match spec {
        TimeSeriesSpec::Constant => LoadSeries::Constant,
        TimeSeriesSpec::Linear { slope } => LoadSeries::Linear { slope: *slope },
        TimeSeriesSpec::Path { times, factors } => LoadSeries::Path {
            times: times.clone(),
            factors: factors.clone(),
        },
    }
}

/// Stage compilation for all three `StageSpec` kinds. `ndof` bounds a `Static`
/// stage's `DisplacementControl` dof (3 planar, 6 spatial); `ndim` bounds a
/// `Transient` stage's ground-motion `direction` (2 planar, 3 spatial —
/// translational axes only, never a rotation).
pub(super) fn compile_stages<NId: Copy, EId: Copy, Load: Clone>(
    stages: &[StageSpec],
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
    pattern_at: &impl Fn(u32) -> Result<LoadPatternId, DecodeError>,
    ndof: u8,
    ndim: u8,
    nodal_loads_by_stage: Vec<Vec<(LoadPatternId, NId, usize, f64)>>,
    element_loads_by_stage: Vec<Vec<(LoadPatternId, EId, Load)>>,
) -> Result<Vec<CompiledStage<NId, EId, Load>>, DecodeError> {
    stages
        .iter()
        .zip(nodal_loads_by_stage)
        .zip(element_loads_by_stage)
        .map(|((stage, pending_nodal_loads), pending_element_loads)| {
            let (steps, kind) = match stage {
                StageSpec::Static {
                    steps,
                    integrator,
                    algorithm,
                    convergence,
                    hold_patterns_after,
                    ..
                } => {
                    let is_arc_length = matches!(integrator, IntegratorSpec::ArcLength(_));
                    let algorithm = algorithm_of(*algorithm, stage.id())?;
                    let convergence = convergence_of(*convergence, stage.id(), is_arc_length)?;
                    let integrator = match integrator {
                        IntegratorSpec::LoadControl { increment } => Integrator::LoadControl {
                            increment: *increment,
                        },
                        IntegratorSpec::DisplacementControl {
                            node,
                            dof,
                            increment,
                        } => {
                            check_dof_within("sequence.integrator", *node, *dof, ndof)?;
                            Integrator::DisplacementControl {
                                node: node_at(*node, "sequence.integrator")?,
                                dof: *dof as usize,
                                increment: *increment,
                            }
                        }
                        IntegratorSpec::ArcLength(spec) => {
                            let config = arc_length_of(spec, node_at, ndof)?;
                            config
                                .validate_with(&algorithm, &convergence)
                                .map_err(|error| analysis_option(stage.id(), error))?;
                            Integrator::ArcLength(config)
                        }
                    };
                    let hold_patterns_after = hold_patterns_after
                        .iter()
                        .map(|&row| pattern_at(row))
                        .collect::<Result<Vec<_>, _>>()?;
                    (
                        *steps,
                        CompiledStageKind::Static {
                            integrator,
                            algorithm,
                            convergence,
                            hold_patterns_after,
                        },
                    )
                }
                // A single eigensolve, not an iterative step loop — always
                // exactly one `advance()` step (`StageSpec::Modal`'s doc
                // comment).
                StageSpec::Modal { modes, .. } => (
                    1,
                    CompiledStageKind::Modal {
                        num_modes: *modes as usize,
                    },
                ),
                // Takes no steps: the domain is swapped for its pristine copy when the stage starts.
                StageSpec::Reset { .. } => (0, CompiledStageKind::Reset),
                StageSpec::Transient {
                    steps,
                    dt,
                    damping,
                    ground_motions,
                    algorithm,
                    convergence,
                    ..
                } => {
                    if !dt.is_finite() || *dt <= 0.0 {
                        return Err(invalid_option(stage.id(), "dt"));
                    }
                    let damping = RayleighDamping::new(damping.alpha_m, damping.beta_k);
                    let ground_motions =
                        ground_motions
                            .iter()
                            .enumerate()
                            .map(|(row, gm)| {
                                check_dof_within(
                                    "sequence.groundMotion",
                                    row as u32,
                                    gm.direction,
                                    ndim,
                                )?;
                                Ok(GroundMotion::new(
                                    gm.direction as usize,
                                    load_series_of(&gm.series),
                                )
                                .with_scale_factor(gm.scale_factor))
                            })
                            .collect::<Result<Vec<_>, DecodeError>>()?;
                    (
                        *steps,
                        CompiledStageKind::Transient {
                            damping,
                            algorithm: algorithm_of(*algorithm, stage.id())?,
                            convergence: convergence_of(*convergence, stage.id(), false)?,
                            dt: *dt,
                            ground_motions,
                        },
                    )
                }
            };

            Ok(CompiledStage {
                id: stage.id().to_string(),
                steps,
                kind,
                pending_nodal_loads,
                pending_element_loads,
            })
        })
        .collect()
}

// Shared by both profiles and both iterative stage kinds.
fn invalid_option(stage: &str, field: &'static str) -> DecodeError {
    DecodeError::InvalidAnalysisOption {
        stage: stage.to_owned(),
        field,
    }
}

fn tangent_of(spec: TangentStrategySpec) -> TangentStrategy {
    match spec {
        TangentStrategySpec::Current => TangentStrategy::Current,
        TangentStrategySpec::ReuseAtStepStart => TangentStrategy::ReuseAtStepStart,
        TangentStrategySpec::Initial => TangentStrategy::Initial,
    }
}

fn algorithm_of(spec: AlgorithmSpec, stage: &str) -> Result<Algorithm, DecodeError> {
    Ok(match spec {
        AlgorithmSpec::Legacy(LegacyAlgorithmSpec::Linear)
        | AlgorithmSpec::Config(AlgorithmConfigSpec::Linear) => Algorithm::Linear,
        AlgorithmSpec::Legacy(LegacyAlgorithmSpec::NewtonRaphson) => Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        },
        AlgorithmSpec::Config(AlgorithmConfigSpec::Newton {
            tangent,
            line_search,
        }) => {
            let line_search = line_search
                .map(|search| {
                    let (tol, max_iter, max_eta) = match search {
                        LineSearchSpec::Bisection {
                            tol,
                            max_iter,
                            max_eta,
                        }
                        | LineSearchSpec::RegulaFalsi {
                            tol,
                            max_iter,
                            max_eta,
                        } => (tol, max_iter, max_eta),
                    };
                    if !tol.is_finite() || tol <= 0.0 {
                        return Err(invalid_option(stage, "algorithm.lineSearch.tol"));
                    }
                    if max_iter == 0 {
                        return Err(invalid_option(stage, "algorithm.lineSearch.maxIter"));
                    }
                    if !max_eta.is_finite() || max_eta < 1.0 {
                        return Err(invalid_option(stage, "algorithm.lineSearch.maxEta"));
                    }
                    Ok(match search {
                        LineSearchSpec::Bisection { .. } => LineSearch::Bisection {
                            tol,
                            max_iter: max_iter as usize,
                            max_eta,
                        },
                        LineSearchSpec::RegulaFalsi { .. } => LineSearch::RegulaFalsi {
                            tol,
                            max_iter: max_iter as usize,
                            max_eta,
                        },
                    })
                })
                .transpose()?;
            Algorithm::Newton {
                tangent: tangent_of(tangent),
                line_search,
            }
        }
        AlgorithmSpec::Config(AlgorithmConfigSpec::KrylovNewton {
            tangent,
            max_dimension,
        }) => {
            if max_dimension == 0 {
                return Err(invalid_option(stage, "algorithm.maxDimension"));
            }
            Algorithm::KrylovNewton {
                tangent: tangent_of(tangent),
                max_dimension: max_dimension as usize,
            }
        }
    })
}

fn convergence_of(
    spec: Option<ConvergenceSpec>,
    stage: &str,
    arc_length: bool,
) -> Result<ConvergenceTest, DecodeError> {
    let spec = spec.unwrap_or(if arc_length {
        ConvergenceSpec::ARC_LENGTH_DEFAULT
    } else {
        ConvergenceSpec::DEFAULT
    });
    let test = match spec {
        ConvergenceSpec::NormUnbalance { tol, max_iter } => ConvergenceTest::NormUnbalance {
            tol,
            max_iter: max_iter as usize,
        },
        ConvergenceSpec::NormDispIncr { tol, max_iter } => ConvergenceTest::NormDispIncr {
            tol,
            max_iter: max_iter as usize,
        },
        ConvergenceSpec::EnergyIncr { tol, max_iter } => ConvergenceTest::EnergyIncr {
            tol,
            max_iter: max_iter as usize,
        },
        ConvergenceSpec::Combined {
            force_tol,
            moment_tol,
            relative_tol,
            displacement_tol,
            max_iter,
        } => ConvergenceTest::Combined {
            force_tol,
            moment_tol: moment_tol.unwrap_or(force_tol),
            relative_tol: relative_tol.unwrap_or(1e-6),
            displacement_tol,
            max_iter: max_iter as usize,
        },
    };
    test.validate()
        .map_err(|error| analysis_option(stage, error))?;
    Ok(test)
}

/// `core`'s own validation reports `AnalysisError::InvalidOption` with the
/// wire field name; restate it as this decoder's error for `stage`.
fn analysis_option(stage: &str, error: AnalysisError) -> DecodeError {
    match error {
        AnalysisError::InvalidOption { field } => invalid_option(stage, field),
        _ => invalid_option(stage, "integrator"),
    }
}

/// `ArcLengthSpec` -> `core::ArcLength`, resolving node indices. Numeric
/// ranges are checked afterwards by `ArcLength::validate_with`.
fn arc_length_of<NId: Copy>(
    spec: &ArcLengthSpec,
    node_at: &impl Fn(u32, &'static str) -> Result<NId, DecodeError>,
    ndof: u8,
) -> Result<ArcLength<NId>, DecodeError> {
    let scales = match spec.scales {
        ArcScalesSpec::Explicit {
            displacement,
            rotation,
            load,
        } => ArcScales::Explicit {
            displacement,
            rotation,
            load,
        },
        ArcScalesSpec::Auto { load } => ArcScales::Auto { load },
    };
    let mut config = ArcLength::adaptive(
        spec.initial_radius,
        spec.min_radius.unwrap_or(spec.initial_radius),
        spec.max_radius.unwrap_or(spec.initial_radius),
        scales,
    );
    if let Some(target) = spec.target_iterations {
        config.target_iterations = target as usize;
    }
    if let Some(retries) = spec.max_retries {
        config.max_retries = retries as usize;
    }
    config.direction = match spec.direction {
        ArcDirectionSpec::Increasing => ArcDirection::Increasing,
        ArcDirectionSpec::Decreasing => ArcDirection::Decreasing,
    };
    config.predictor = match spec.predictor {
        ArcPredictorSpec::Secant => ArcPredictor::Secant,
        ArcPredictorSpec::Tangent => ArcPredictor::Tangent,
    };
    if let Some(seed) = &spec.seed {
        let table = "sequence.integrator.seed";
        let components = seed
            .components
            .iter()
            .map(|c| {
                check_dof_within(table, c.node, c.dof, ndof)?;
                Ok((node_at(c.node, table)?, c.dof as usize, c.value))
            })
            .collect::<Result<Vec<_>, DecodeError>>()?;
        config.seed = Some(ArcSeed {
            components,
            load: seed.load,
        });
    }
    if let Some(tolerance) = spec.arc_tolerance {
        config.arc_tolerance = tolerance;
    }
    config.correction_tolerance = spec.correction_tolerance;
    if let Some(backtracking) = spec.backtracking {
        config.backtracking = Backtracking {
            armijo: backtracking.armijo,
            min_step: backtracking.min_step,
        };
    }
    if let Some(stop) = spec.stop {
        if let Some(target) = stop.displacement {
            let table = "sequence.integrator.stop";
            check_dof_within(table, target.node, target.dof, ndof)?;
            config.stop.displacement = Some(DisplacementTarget {
                node: node_at(target.node, table)?,
                dof: target.dof as usize,
                value: target.value,
                exact: target.exact,
            });
        }
        config.stop.load_factor = stop.load_factor.map(|target| LoadFactorTarget {
            value: target.value,
            exact: target.exact,
        });
        config.stop.load_factor_zero_crossing = stop.load_factor_zero_crossing;
        config.stop.max_chord_length = stop.max_chord_length;
    }
    Ok(config)
}
