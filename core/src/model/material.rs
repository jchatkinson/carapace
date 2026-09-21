/// Uniaxial material catalog. Closed enum, not a trait object (§3.1 of the
/// implementation plan) — the set is fixed at compile time. `Steel01`,
/// `Concrete01` and the `Parallel`/`Series`/`MinMax` composites land at M7
/// stage 2.
///
/// # Trial vs. commit
///
/// A `Material` value holds only *committed* state (e.g. `ElasticPP`'s
/// plastic strain `ep`) — there is no OpenSees-style parallel `T`/`C`
/// field pair. Instead there are two pure (`&self`, non-mutating) views
/// over that committed state:
///
/// - [`trial_stress_tangent`](Material::trial_stress_tangent) — called
///   every Newton iteration, as many times as needed, always relative to
///   the same fixed committed state. Since it can't mutate anything,
///   there's no way for one iteration's (possibly discarded) trial to leak
///   into the next — a real correctness requirement for path-dependent
///   materials, not just a convenience: `setTrialStrain` in OpenSees is
///   `&mut self`, but its `T*` fields are always recomputed from the `C*`
///   fields, i.e. it's *this* pure-function property, not the mutation
///   itself, that OpenSees's parallel field pair exists to guarantee.
/// - [`commit`](Material::commit) — called once, only after a step's
///   Newton iteration actually converges (or, for `Algorithm::Linear`,
///   after its one unconditional solve), producing the new committed
///   `Material`. `Domain::commit` is what calls this and writes the
///   result back into the model; nothing else ever mutates a `Material`.
///
/// This gets rollback "for free": since nothing is mutated until
/// `commit`, a step that fails partway through a Newton loop has left no
/// material in an inconsistent state to revert — only nodal displacement
/// (mutated eagerly every iteration, which is correct Newton behavior)
/// needs `Analysis`/`TransientAnalysis`'s snapshot-and-restore-on-failure.
///
/// `Material` is `Clone`, not `Copy`: `Parallel`/`Series`/`MinMax` hold a
/// `Vec`/`Box` of child materials, so every call site that used to rely on
/// an implicit copy (array-repeat initializers, `match *self`, ...) needs
/// an explicit `.clone()` instead.
#[derive(Debug, Clone)]
pub enum Material {
    Elastic {
        e: f64,
    },
    /// Elastic-perfectly-plastic, with real permanent-set tracking via the
    /// committed plastic strain `ep` — standard 1D return-mapping
    /// plasticity. Construct via [`Material::elastic_pp`], not the struct
    /// literal: `ep` is internal history, not a material parameter, and
    /// must start at zero.
    ElasticPP {
        e: f64,
        eyp: f64,
        ep: f64,
    },
    /// One-sided gap: zero stiffness until the strain closes the gap
    /// (`strain + gap < 0`, i.e. compressive strain past the opening), then
    /// linear elastic. Stateless — a genuine function of current strain
    /// only (no permanent set to track), unlike `ElasticPP`.
    Gap {
        e: f64,
        gap: f64,
    },
    /// "No tension": zero stiffness in tension (`strain >= 0`), linear
    /// elastic in compression. Equivalent to `Gap { e, gap: 0.0 }`, kept as
    /// its own variant since OpenSees exposes it separately and it's the
    /// simplest possible nonlinear material (useful as the base case for M2
    /// wiring). Stateless, same reasoning as `Gap`.
    Ent {
        e: f64,
    },
    /// Bilinear kinematic hardening with isotropic hardening on load
    /// reversal (Filippou 1983), ported from OpenSees' `Steel01`. Construct
    /// via [`Material::steel01`].
    ///
    /// `min_strain`/`max_strain` track the largest excursion ever committed
    /// in each direction (used only to update `shift_p`/`shift_n` on a
    /// reversal, so they lag one reversal behind, exactly like OpenSees'
    /// `CminStrain`/`CmaxStrain`); `shift_p`/`shift_n` are the isotropic
    /// hardening shifts applied to the bounding lines; `strain`/`stress`
    /// are the last committed point (needed because the bounding formula
    /// `c = Cstress + E0*dStrain` references it directly, not just the
    /// "obvious" history fields — see the M7 stage-2 handoff notes).
    Steel01 {
        fy: f64,
        e0: f64,
        b: f64,
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
        min_strain: f64,
        max_strain: f64,
        shift_p: f64,
        shift_n: f64,
        loading: Steel01Loading,
        strain: f64,
        stress: f64,
    },
    /// Kent-Scott-Park compression envelope with Karsan-Jirsa degrading
    /// unload/reload and zero tensile strength, ported from OpenSees'
    /// `Concrete01`. Construct via [`Material::concrete01`] — `fpc`,
    /// `epsc0`, `fpcu`, `epscu` are normalized negative (compression) at
    /// construction, matching the OpenSees convention.
    ///
    /// `min_strain` is the most negative strain ever committed (drives the
    /// envelope); `end_strain` is where the current unload/reload line
    /// would cross zero stress (needed by `reload`'s branch to tell
    /// whether a trial strain is still on that line or must re-enter the
    /// envelope); `unload_slope` is that line's slope; `strain`/`stress`
    /// are the last committed point (same reason as `Steel01` above).
    Concrete01 {
        fpc: f64,
        epsc0: f64,
        fpcu: f64,
        epscu: f64,
        min_strain: f64,
        end_strain: f64,
        unload_slope: f64,
        strain: f64,
        stress: f64,
    },
    /// Composite: every child sees the same strain, stress/tangent are the
    /// factor-weighted sum. Ported from OpenSees' `ParallelMaterial` — the
    /// factor defaults to `1.0` (see [`Material::parallel`]) but is kept
    /// explicit per OpenSees' API.
    Parallel(Vec<(Material, f64)>),
    /// Composite: every child sees the same *stress*, total strain is the
    /// sum of the children's individual strains. Not directly invertible
    /// (unlike `Parallel`), so `evaluate` runs its own small internal
    /// flexibility-based Newton loop, entirely self-contained within one
    /// call — see the type's doc and the M7 stage-2 handoff notes for why
    /// this doesn't violate the trial/commit purity rule. `child_strains`
    /// is the last converged per-child strain split, kept explicitly
    /// (mirroring OpenSees' own `SeriesMaterial::strain[]` array) since not
    /// every leaf material tracks its own committed strain. Construct via
    /// [`Material::series`].
    Series {
        children: Vec<Material>,
        child_strains: Vec<f64>,
        stress: f64,
        max_iter: usize,
        tol: f64,
    },
    /// Wraps one material with strain bounds: once a committed strain
    /// exits `[min_strain, max_strain]`, `failed` becomes (and permanently
    /// stays) `true` — stress reads as `0.0`, tangent as `1e-8 *
    /// inner.initial_tangent()` (not exactly zero, to avoid handing back a
    /// structurally singular row/column — see OpenSees'
    /// `MinMaxMaterial::getTangent`). Once failed, the inner material's
    /// state is frozen (no further trial/commit reaches it), matching
    /// OpenSees' `Cfailed` short-circuit. Construct via
    /// [`Material::min_max`].
    MinMax {
        inner: Box<Material>,
        min_strain: f64,
        max_strain: f64,
        failed: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Steel01Loading {
    None,
    Loading,
    Unloading,
}

impl Material {
    pub fn elastic_pp(e: f64, eyp: f64) -> Self {
        Material::ElasticPP { e, eyp, ep: 0.0 }
    }

    pub fn steel01(fy: f64, e0: f64, b: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> Self {
        Material::Steel01 {
            fy,
            e0,
            b,
            a1,
            a2,
            a3,
            a4,
            min_strain: 0.0,
            max_strain: 0.0,
            shift_p: 1.0,
            shift_n: 1.0,
            loading: Steel01Loading::None,
            strain: 0.0,
            stress: 0.0,
        }
    }

    /// `fpc`/`epsc0`/`fpcu`/`epscu` are normalized to be negative
    /// (compression), matching OpenSees' `Concrete01` constructor.
    pub fn concrete01(fpc: f64, epsc0: f64, fpcu: f64, epscu: f64) -> Self {
        let fpc = -fpc.abs();
        let epsc0 = -epsc0.abs();
        let fpcu = -fpcu.abs();
        let epscu = -epscu.abs();
        let ec0 = 2.0 * fpc / epsc0;
        Material::Concrete01 {
            fpc,
            epsc0,
            fpcu,
            epscu,
            min_strain: 0.0,
            end_strain: 0.0,
            unload_slope: ec0,
            strain: 0.0,
            stress: 0.0,
        }
    }

    /// Equal-factor (`1.0` each) parallel composite. For custom factors,
    /// build the `Material::Parallel(Vec<(Material, f64)>)` variant
    /// directly.
    pub fn parallel(materials: Vec<Material>) -> Self {
        Material::Parallel(materials.into_iter().map(|m| (m, 1.0)).collect())
    }

    /// Series composite with OpenSees' default tolerance (`1e-10`) and a
    /// generous internal iteration cap (`30`) — since, unlike OpenSees'
    /// `SeriesMaterial`, `evaluate` can't lean on many external Newton
    /// calls to slowly converge a warm-started trial state (see the type's
    /// doc), the internal loop needs to reach convergence by itself.
    pub fn series(children: Vec<Material>) -> Self {
        let n = children.len();
        Material::Series {
            children,
            child_strains: vec![0.0; n],
            stress: 0.0,
            max_iter: 30,
            tol: 1e-10,
        }
    }

    pub fn min_max(inner: Material, min_strain: f64, max_strain: f64) -> Self {
        Material::MinMax {
            inner: Box::new(inner),
            min_strain,
            max_strain,
            failed: false,
        }
    }

    /// `(stress, tangent)` for a trial evaluation at `strain` — read-only,
    /// safe to call any number of times from any Newton iteration. See the
    /// type-level doc comment for why this must never mutate `self`.
    pub fn trial_stress_tangent(&self, strain: f64) -> (f64, f64) {
        let (stress, tangent, _next) = self.evaluate(strain);
        (stress, tangent)
    }

    /// The new committed `Material` if `strain` were accepted as the
    /// converged state. Called once per converged step, by `Domain::commit`
    /// — see the type-level doc comment.
    pub fn commit(&self, strain: f64) -> Material {
        let (_stress, _tangent, next) = self.evaluate(strain);
        next
    }

    /// The tangent modulus at zero history — OpenSees' `getInitialTangent`.
    /// Only needed by `MinMax` (for its post-failure tangent) and the
    /// composites (to propagate it to a wrapping `MinMax`).
    pub fn initial_tangent(&self) -> f64 {
        match self {
            Material::Elastic { e } => *e,
            Material::ElasticPP { e, .. } => *e,
            Material::Gap { e, .. } => *e,
            Material::Ent { e } => *e,
            Material::Steel01 { e0, .. } => *e0,
            Material::Concrete01 { fpc, epsc0, .. } => 2.0 * fpc / epsc0,
            Material::Parallel(children) => children.iter().map(|(m, f)| f * m.initial_tangent()).sum(),
            Material::Series { children, .. } => {
                let total_flex: f64 = children.iter().map(|m| 1.0 / m.initial_tangent()).sum();
                if total_flex.abs() > 1e-12 {
                    1.0 / total_flex
                } else {
                    0.0
                }
            }
            Material::MinMax { inner, .. } => inner.initial_tangent(),
        }
    }

    /// Shared core: `(stress, tangent, next_committed_state)`. Leaf
    /// materials without history (`Elastic`, `Gap`, `Ent`) return `self`'s
    /// own fields unchanged as `next_committed_state` — committing them is
    /// a no-op, correctly, since they have nothing to remember.
    fn evaluate(&self, strain: f64) -> (f64, f64, Material) {
        match self {
            Material::Elastic { e } => (e * strain, *e, self.clone()),
            Material::ElasticPP { e, eyp, ep } => {
                let (e, eyp, ep) = (*e, *eyp, *ep);
                let yield_stress = e * eyp;
                // Elastic predictor relative to the *committed* plastic
                // strain — standard 1D return mapping.
                let trial_stress = e * (strain - ep);
                if trial_stress.abs() <= yield_stress {
                    (trial_stress, e, Material::ElasticPP { e, eyp, ep })
                } else {
                    let stress = yield_stress.copysign(trial_stress);
                    // Return to the yield surface: the plastic strain
                    // absorbs whatever additional strain the elastic
                    // predictor overshot by.
                    let next_ep = strain - stress / e;
                    (
                        stress,
                        0.0,
                        Material::ElasticPP {
                            e,
                            eyp,
                            ep: next_ep,
                        },
                    )
                }
            }
            Material::Gap { e, gap } => {
                // `<=`, not `<`: at exactly zero closure the tangent must
                // stay nonzero (matching the just-closed regime) or a
                // system starting from zero displacement forms a singular
                // initial tangent — there is no elastic-vs-open distinction
                // to make at a single point anyway.
                let closure = strain + gap;
                if closure <= 0.0 {
                    (e * closure, *e, self.clone())
                } else {
                    (0.0, 0.0, self.clone())
                }
            }
            Material::Ent { e } => {
                // See `Gap` above for why this is `<=`, not `<`.
                if strain <= 0.0 {
                    (e * strain, *e, self.clone())
                } else {
                    (0.0, 0.0, self.clone())
                }
            }
            Material::Steel01 { .. } => self.evaluate_steel01(strain),
            Material::Concrete01 { .. } => self.evaluate_concrete01(strain),
            Material::Parallel(children) => {
                let mut stress = 0.0;
                let mut tangent = 0.0;
                let mut next_children = Vec::with_capacity(children.len());
                for (child, factor) in children {
                    let (s, t, next) = child.evaluate(strain);
                    stress += factor * s;
                    tangent += factor * t;
                    next_children.push((next, *factor));
                }
                (stress, tangent, Material::Parallel(next_children))
            }
            Material::Series { .. } => self.evaluate_series(strain),
            Material::MinMax {
                inner,
                min_strain,
                max_strain,
                failed,
            } => {
                if *failed || strain <= *min_strain || strain >= *max_strain {
                    let tangent = 1.0e-8 * inner.initial_tangent();
                    (
                        0.0,
                        tangent,
                        Material::MinMax {
                            inner: inner.clone(),
                            min_strain: *min_strain,
                            max_strain: *max_strain,
                            failed: true,
                        },
                    )
                } else {
                    let (s, t, next_inner) = inner.evaluate(strain);
                    (
                        s,
                        t,
                        Material::MinMax {
                            inner: Box::new(next_inner),
                            min_strain: *min_strain,
                            max_strain: *max_strain,
                            failed: false,
                        },
                    )
                }
            }
        }
    }

    /// Ported from `Steel01::determineTrialState` (`Steel01.cpp`). `self`
    /// stands in for the `C*` fields; `strain` is the trial (`Tstrain`).
    fn evaluate_steel01(&self, strain: f64) -> (f64, f64, Material) {
        let Material::Steel01 {
            fy,
            e0,
            b,
            a1,
            a2,
            a3,
            a4,
            min_strain,
            max_strain,
            shift_p,
            shift_n,
            loading,
            strain: cstrain,
            stress: cstress,
        } = self
        else {
            unreachable!()
        };
        let (fy, e0, b, a1, a2, a3, a4) = (*fy, *e0, *b, *a1, *a2, *a3, *a4);
        let (mut min_strain, mut max_strain, mut shift_p, mut shift_n) =
            (*min_strain, *max_strain, *shift_p, *shift_n);
        let mut loading = *loading;
        let (cstrain, cstress) = (*cstrain, *cstress);

        let dstrain = strain - cstrain;

        let fy_one_minus_b = fy * (1.0 - b);
        let esh = b * e0;
        let epsy = fy / e0;

        let c1 = esh * strain;
        let c2 = shift_n * fy_one_minus_b;
        let c3 = shift_p * fy_one_minus_b;
        let c = cstress + e0 * dstrain;

        // The three-way clamp `max(c1-c2, min(c1+c3, c))` — see the
        // commented-out simplification in `Steel01.cpp`, which the source
        // says only exists split into two `if`s as an optimizer
        // workaround; the single-expression form is used here.
        let stress = (c1 - c2).max((c1 + c3).min(c));
        let tangent = if (stress - c).abs() < f64::EPSILON { e0 } else { esh };

        if loading == Steel01Loading::None && dstrain != 0.0 {
            loading = if dstrain > 0.0 {
                Steel01Loading::Loading
            } else {
                Steel01Loading::Unloading
            };
        }

        // Transition from loading to unloading.
        if loading == Steel01Loading::Loading && dstrain < 0.0 {
            loading = Steel01Loading::Unloading;
            if cstrain > max_strain {
                max_strain = cstrain;
            }
            shift_n = 1.0 + a1 * ((max_strain - min_strain) / (2.0 * a2 * epsy)).powf(0.8);
        }

        // Transition from unloading to loading.
        if loading == Steel01Loading::Unloading && dstrain > 0.0 {
            loading = Steel01Loading::Loading;
            if cstrain < min_strain {
                min_strain = cstrain;
            }
            shift_p = 1.0 + a3 * ((max_strain - min_strain) / (2.0 * a4 * epsy)).powf(0.8);
        }

        (
            stress,
            tangent,
            Material::Steel01 {
                fy,
                e0,
                b,
                a1,
                a2,
                a3,
                a4,
                min_strain,
                max_strain,
                shift_p,
                shift_n,
                loading,
                strain,
                stress,
            },
        )
    }

    /// Ported from `Concrete01::setTrialStrain` (`Concrete01.cpp`), the
    /// version actually driving the reload/envelope/unload state machine
    /// (`determineTrialState` is dead code in the original — never called
    /// from `setTrialStrain`).
    fn evaluate_concrete01(&self, strain: f64) -> (f64, f64, Material) {
        let Material::Concrete01 {
            fpc,
            epsc0,
            fpcu,
            epscu,
            min_strain,
            end_strain,
            unload_slope,
            strain: cstrain,
            stress: cstress,
        } = self
        else {
            unreachable!()
        };
        let (fpc, epsc0, fpcu, epscu) = (*fpc, *epsc0, *fpcu, *epscu);
        let (min_strain, end_strain, unload_slope) = (*min_strain, *end_strain, *unload_slope);
        let (cstrain, cstress) = (*cstrain, *cstress);

        // Quick return: tension, zero stiffness, no compression-envelope
        // history update — but strain/stress still commit to the new
        // point (point 2 of the M7 stage-2 porting recipe).
        if strain > 0.0 {
            return (
                0.0,
                0.0,
                Material::Concrete01 {
                    fpc,
                    epsc0,
                    fpcu,
                    epscu,
                    min_strain,
                    end_strain,
                    unload_slope,
                    strain,
                    stress: 0.0,
                },
            );
        }

        let temp_stress = cstress + unload_slope * strain - unload_slope * cstrain;

        let (stress, tangent, min_strain, end_strain, unload_slope) = if strain < cstrain {
            // Further into compression: `reload()`.
            let (mut stress, mut tangent, min_strain, end_strain, unload_slope) = if strain <= min_strain {
                let new_min_strain = strain;
                let (stress, tangent) = concrete01_envelope(new_min_strain, fpc, epsc0, fpcu, epscu);
                let (end_strain, unload_slope) = concrete01_unload(new_min_strain, stress, fpc, epsc0, epscu);
                (stress, tangent, new_min_strain, end_strain, unload_slope)
            } else if strain <= end_strain {
                let tangent = unload_slope;
                let stress = tangent * (strain - end_strain);
                (stress, tangent, min_strain, end_strain, unload_slope)
            } else {
                (0.0, 0.0, min_strain, end_strain, unload_slope)
            };

            if temp_stress > stress {
                stress = temp_stress;
                tangent = unload_slope;
            }

            (stress, tangent, min_strain, end_strain, unload_slope)
        } else if temp_stress <= 0.0 {
            (temp_stress, unload_slope, min_strain, end_strain, unload_slope)
        } else {
            (0.0, 0.0, min_strain, end_strain, unload_slope)
        };

        (
            stress,
            tangent,
            Material::Concrete01 {
                fpc,
                epsc0,
                fpcu,
                epscu,
                min_strain,
                end_strain,
                unload_slope,
                strain,
                stress,
            },
        )
    }

    /// Ported from `SeriesMaterial::setTrialStrain` (`SeriesMaterial.cpp`)
    /// — flexibility-based iteration, but self-contained: applies the
    /// stress correction `dq` every inner iteration (a proper Newton
    /// update) rather than only once after the loop, since — unlike
    /// OpenSees, which carries `Tstress`/`Ttangent` across many external
    /// `setTrialStrain` calls to slowly converge — this must converge
    /// within one `evaluate` call. See the type's doc comment and the M7
    /// stage-2 handoff notes.
    fn evaluate_series(&self, total_strain: f64) -> (f64, f64, Material) {
        let Material::Series {
            children,
            child_strains,
            stress: cstress,
            max_iter,
            tol,
        } = self
        else {
            unreachable!()
        };

        let n = children.len();
        let mut e = child_strains.clone();
        let mut child_stress = vec![0.0; n];
        let mut flex = vec![0.0; n];
        let mut tstress = *cstress;
        let mut series_tangent = 0.0;

        let safe_inv = |t: f64| {
            if t.abs() > 1.0e-12 {
                1.0 / t
            } else if t < 0.0 {
                -1.0e12
            } else {
                1.0e12
            }
        };

        for _ in 0..*max_iter {
            let mut f = 0.0;
            let mut vr = 0.0;
            for i in 0..n {
                let ds = tstress - child_stress[i];
                e[i] += flex[i] * ds;
                let (s, t) = children[i].trial_stress_tangent(e[i]);
                child_stress[i] = s;
                flex[i] = safe_inv(t);
                let ds = tstress - child_stress[i];
                let de = flex[i] * ds;
                f += flex[i];
                vr += e[i] + de;
            }
            series_tangent = safe_inv(f);
            let dv = total_strain - vr;
            let dq = series_tangent * dv;
            tstress += dq;
            if (dq * dv).abs() < *tol {
                break;
            }
        }

        let next_children = children
            .iter()
            .zip(e.iter())
            .map(|(child, &strain)| child.evaluate(strain).2)
            .collect();

        (
            tstress,
            series_tangent,
            Material::Series {
                children: next_children,
                child_strains: e,
                stress: tstress,
                max_iter: *max_iter,
                tol: *tol,
            },
        )
    }
}

/// `Concrete01::envelope()`: the Kent-Scott-Park backbone, evaluated at a
/// strain at or past the previous minimum (i.e. loading further into
/// compression along the virgin curve).
fn concrete01_envelope(strain: f64, fpc: f64, epsc0: f64, fpcu: f64, epscu: f64) -> (f64, f64) {
    if strain > epsc0 {
        let eta = strain / epsc0;
        let stress = fpc * (2.0 * eta - eta * eta);
        let ec0 = 2.0 * fpc / epsc0;
        let tangent = ec0 * (1.0 - eta);
        (stress, tangent)
    } else if strain > epscu {
        let tangent = (fpc - fpcu) / (epsc0 - epscu);
        let stress = fpc + tangent * (strain - epsc0);
        (stress, tangent)
    } else {
        (fpcu, 0.0)
    }
}

/// `Concrete01::unload()`: the Karsan-Jirsa degrading unload-slope formula,
/// given the new `min_strain` and the envelope stress just computed there.
fn concrete01_unload(min_strain: f64, stress_at_min: f64, fpc: f64, epsc0: f64, epscu: f64) -> (f64, f64) {
    let mut temp_strain = min_strain;
    if temp_strain < epscu {
        temp_strain = epscu;
    }
    let eta = temp_strain / epsc0;

    let mut ratio = 0.707 * (eta - 2.0) + 0.834;
    if eta < 2.0 {
        ratio = 0.145 * eta * eta + 0.13 * eta;
    }

    let mut end_strain = ratio * epsc0;
    let temp1 = min_strain - end_strain;
    let ec0 = 2.0 * fpc / epsc0;
    let temp2 = stress_at_min / ec0;

    let unload_slope;
    if temp1 > -f64::EPSILON {
        // Should always be negative in practice.
        unload_slope = ec0;
    } else if temp1 <= temp2 {
        end_strain = min_strain - temp1;
        unload_slope = stress_at_min / temp1;
    } else {
        end_strain = min_strain - temp2;
        unload_slope = ec0;
    }

    (end_strain, unload_slope)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elastic_pp_clamps_symmetrically_on_first_loading() {
        // First loading from ep=0: identical to the old reversible-envelope
        // behavior, since nothing has been committed yet to make history
        // matter.
        let m = Material::elastic_pp(100.0, 0.01);
        assert_eq!(m.trial_stress_tangent(0.005), (0.5, 100.0));
        assert_eq!(m.trial_stress_tangent(0.02), (1.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.02), (-1.0, 0.0));
    }

    #[test]
    fn elastic_pp_remembers_permanent_set_after_commit() {
        let m = Material::elastic_pp(100.0, 0.01);

        // Load past yield (strain=0.02 well beyond eyp=0.01) and commit.
        let m = m.commit(0.02);
        let Material::ElasticPP { ep, .. } = m else { panic!() };
        assert!((ep - 0.01).abs() < 1e-12, "expected ep=0.01, got {ep}");

        // Partially unload to strain=0.015 (still positive, but now less
        // than the peak). A stateless envelope would still clamp this to
        // the yield force (|100*0.015|=1.5 > fy=1.0), but real plasticity
        // unloads *elastically* from the committed ep, not back through
        // the origin: trial_stress = e*(strain-ep) = 100*(0.015-0.01)=0.5,
        // well inside the yield surface.
        let (stress, tangent) = m.trial_stress_tangent(0.015);
        assert!((stress - 0.5).abs() < 1e-12, "expected elastic unload stress=0.5, got {stress}");
        assert_eq!(tangent, 100.0, "should be elastic (unloading), not the plastic tangent");
    }

    #[test]
    fn gap_only_engages_past_closure() {
        let m = Material::Gap { e: 100.0, gap: 0.01 };
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.005), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.02), (100.0 * (-0.02 + 0.01), 100.0));
    }

    #[test]
    fn ent_carries_no_tension() {
        let m = Material::Ent { e: 100.0 };
        assert_eq!(m.trial_stress_tangent(0.01), (0.0, 0.0));
        assert_eq!(m.trial_stress_tangent(-0.01), (-1.0, 100.0));
    }

    #[test]
    fn steel01_matches_elastic_on_first_loading_within_yield() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        let (stress, tangent) = m.trial_stress_tangent(0.001);
        assert!((stress - 29000.0 * 0.001).abs() < 1e-9, "got {stress}");
        assert_eq!(tangent, 29000.0);
    }

    #[test]
    fn steel01_bilinear_hardening_past_yield() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        // epsy = fy/E0 = 60/29000 ≈ 0.002069. Push well past it.
        let (stress, tangent) = m.trial_stress_tangent(0.02);
        let esh = 0.01 * 29000.0;
        // On first loading (shift_p=shift_n=1, min=max=0) the upper
        // bounding line is c1+c3 = Esh*strain + fy*(1-b).
        let expected = esh * 0.02 + 60.0 * (1.0 - 0.01);
        assert!((stress - expected).abs() < 1e-6, "expected {expected}, got {stress}");
        assert!((tangent - esh).abs() < 1e-9);
    }

    #[test]
    fn steel01_isotropic_hardening_shifts_the_bound_on_reversal() {
        let m = Material::steel01(60.0, 29000.0, 0.01, 0.9, 5.0, 0.9, 5.0);
        // Load well into the hardening range, commit, then reverse twice —
        // the first reversal (loading -> unloading) shifts shift_n, the
        // second (unloading -> loading) shifts shift_p.
        let m = m.commit(0.02);
        let m = m.commit(-0.02);
        let Material::Steel01 { shift_n, .. } = m else { panic!() };
        assert!(shift_n > 1.0, "shift_n should have grown after the first reversal: {shift_n}");

        let m = m.commit(0.02);
        let Material::Steel01 { shift_p, .. } = m else { panic!() };
        assert!(shift_p > 1.0, "shift_p should have grown after the second reversal: {shift_p}");
    }

    #[test]
    fn concrete01_matches_parabolic_envelope_on_first_loading() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        // Envelope is symmetric-normalized negative; strain -0.001 is
        // within the parabola (|strain| < |epsc0|).
        let strain = -0.001;
        let (stress, _tangent) = m.trial_stress_tangent(strain);
        let eta = strain / -0.002;
        let expected = -4.0 * (2.0 * eta - eta * eta);
        assert!((stress - expected).abs() < 1e-9, "expected {expected}, got {stress}");
    }

    #[test]
    fn concrete01_carries_no_tension() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        assert_eq!(m.trial_stress_tangent(0.001), (0.0, 0.0));
    }

    #[test]
    fn concrete01_unloads_with_degrading_stiffness_below_initial() {
        let m = Material::concrete01(4.0, 0.002, 3.0, 0.006);
        let ec0 = 2.0 * 4.0 / 0.002;

        // Push past peak stress (epsc0=-0.002) then commit, then partially
        // unload — but stay on the unload line (between min_strain and
        // end_strain), not far enough to cross back through zero stress.
        let m = m.commit(-0.003);
        let (_stress, unload_tangent) = m.trial_stress_tangent(-0.002);
        assert!(
            unload_tangent < ec0 && unload_tangent > 0.0,
            "Karsan-Jirsa unload slope should be positive but degraded below Ec0={ec0}, got {unload_tangent}"
        );
    }

    #[test]
    fn parallel_of_two_elastics_equals_one_elastic_with_summed_stiffness() {
        let combined = Material::parallel(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 50.0 }]);
        let single = Material::Elastic { e: 150.0 };
        assert_eq!(combined.trial_stress_tangent(0.01), single.trial_stress_tangent(0.01));
    }

    #[test]
    fn parallel_history_is_independent_per_child() {
        // Two ElasticPP with different yield points in parallel: combined
        // response is the sum of each independently-clamped branch.
        let combined = Material::parallel(vec![Material::elastic_pp(100.0, 0.01), Material::elastic_pp(50.0, 0.02)]);
        let (stress, _tangent) = combined.trial_stress_tangent(0.05);
        // Both branches are well past yield: 100*0.01 + 50*0.02 = 2.0.
        assert!((stress - 2.0).abs() < 1e-9, "got {stress}");
    }

    #[test]
    fn series_of_two_equal_elastics_halves_the_stiffness() {
        let series = Material::series(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 100.0 }]);
        let (stress, tangent) = series.trial_stress_tangent(0.02);
        // Two equal springs in series: combined stiffness = E/2, so for a
        // given total strain the stress equals a single spring of
        // stiffness 50 at that same total strain.
        assert!((tangent - 50.0).abs() < 1e-6, "expected tangent 50, got {tangent}");
        assert!((stress - 50.0 * 0.02).abs() < 1e-6, "got {stress}");
    }

    #[test]
    fn series_splits_strain_so_each_child_sees_equal_stress() {
        let series = Material::series(vec![Material::Elastic { e: 100.0 }, Material::Elastic { e: 300.0 }]);
        let Material::Series { child_strains, .. } = series.commit(0.04) else {
            panic!()
        };
        // Equilibrium requires equal stress: 100*e0 = 300*e1, and
        // compatibility requires e0+e1=0.04 => e0=0.03, e1=0.01.
        assert!((child_strains[0] - 0.03).abs() < 1e-6, "got {child_strains:?}");
        assert!((child_strains[1] - 0.01).abs() < 1e-6, "got {child_strains:?}");
    }

    #[test]
    fn min_max_fails_permanently_once_strain_exits_bounds() {
        let m = Material::min_max(Material::Elastic { e: 100.0 }, -0.01, 0.01);
        let m = m.commit(0.02); // exceeds max_strain
        let Material::MinMax { failed, .. } = m else { panic!() };
        assert!(failed);

        let (stress, tangent) = m.trial_stress_tangent(0.0); // back within bounds
        assert_eq!(stress, 0.0, "failed material must read zero stress even if strain re-enters bounds");
        assert!((tangent - 1.0e-8 * 100.0).abs() < 1e-12, "got {tangent}");
    }

    #[test]
    fn min_max_passes_through_inner_material_before_failure() {
        let m = Material::min_max(Material::Elastic { e: 100.0 }, -0.01, 0.01);
        assert_eq!(m.trial_stress_tangent(0.005), (0.5, 100.0));
    }
}
