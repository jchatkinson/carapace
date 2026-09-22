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
mod composite;
mod concrete01;
mod concrete02;
mod elastic_like;
mod hysteretic;
mod steel01;
mod steel02;

pub use hysteretic::HystereticFields;
pub use steel01::Steel01Loading;
pub use steel02::Steel02Kon;

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
    /// Menegotto-Pinto transition curve with Filippou isotropic hardening,
    /// ported from OpenSees' `Steel02`. Unlike `Steel01`'s piecewise-linear
    /// bilinear rule, the loading/unloading branches are smooth curved
    /// asymptotes to the elastic and strain-hardening lines (parameters
    /// `r0`/`cr1`/`cr2` control the transition sharpness). Construct via
    /// [`Material::steel02`] — `sigini` (OpenSees' optional initial-stress
    /// offset) is not supported, matching `Material::steel01`'s scope.
    ///
    /// `min_strain`/`max_strain` are the peak strain reached on each side so
    /// far (`Steel02`'s `epsmaxP`/`epsminP` — note the reversed initial
    /// values vs. `Steel01`: `Steel02` seeds them to `±fy/e0`, not `0.0`, so
    /// the very first reversal has a nonzero isotropic-hardening baseline);
    /// `pl_strain`/`asymptote_strain`/`asymptote_stress` are the strain at
    /// the last reversal's plastic-excursion bound and the strain/stress
    /// where the current curve's two asymptotes (elastic and
    /// strain-hardening) intersect (`Steel02`'s `epsplP`/`epss0P`/`sigs0P`);
    /// `reversal_strain`/`reversal_stress` are the point of the last load
    /// reversal (`epssrP`/`sigsrP`) — the curve's other anchor. `kon`
    /// mirrors `Steel02`'s `konP`, minus its `sigini`-only pinned state.
    Steel02 {
        fy: f64,
        e0: f64,
        b: f64,
        r0: f64,
        cr1: f64,
        cr2: f64,
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
        min_strain: f64,
        max_strain: f64,
        pl_strain: f64,
        asymptote_strain: f64,
        asymptote_stress: f64,
        reversal_strain: f64,
        reversal_stress: f64,
        kon: Steel02Kon,
        strain: f64,
        stress: f64,
    },
    /// `Concrete01`'s Kent-Scott-Park compression envelope plus linear
    /// tension softening and a stiffer reloading rule, ported from
    /// OpenSees' `Concrete02`. Construct via [`Material::concrete02`] —
    /// `fc`/`epsc0`/`fcu`/`epscu` are normalized negative, same as
    /// `Concrete01`.
    ///
    /// `min_strain` is the most negative strain ever committed (`Concrete02`'s
    /// `ecminP`, same role as `Concrete01`'s `min_strain`); `tension_strain`
    /// is the peak tensile-excursion strain beyond the unload line's
    /// zero-stress crossing (`Concrete02`'s `deptP` — the tension-side
    /// analogue of `min_strain`, with no `Concrete01` equivalent since
    /// `Concrete01` carries no tension branch at all).
    Concrete02 {
        fc: f64,
        epsc0: f64,
        fcu: f64,
        epscu: f64,
        rat: f64,
        ft: f64,
        ets: f64,
        min_strain: f64,
        tension_strain: f64,
        strain: f64,
        stress: f64,
    },
    /// Multi-linear (up to 3 points per side) backbone with pinching,
    /// deformation/energy-based damage, and ductility-degraded unloading
    /// stiffness, ported from OpenSees' `HystereticMaterial`. Construct via
    /// [`Material::hysteretic`].
    ///
    /// Boxed (`HystereticFields` bundles ~30 fields — see its own doc
    /// comment) rather than inlined into the enum like every other leaf
    /// variant: `Material`'s size is its largest variant's size, and
    /// without boxing this one, *every* `Material` — including the huge
    /// majority that are a plain `Elastic { e: f64 }` — would pay for the
    /// most field-heavy variant's stack space, everywhere a `Material`
    /// is stored (`Fiber`, `ZeroLength`'s per-DOF array, every composite's
    /// `Vec<Material>`). `Pinching4` (stage 4's other material) boxes the
    /// same way, for the same reason.
    Hysteretic(Box<HystereticFields>),
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

impl Material {
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
            Material::Steel02 { e0, .. } => *e0,
            Material::Concrete02 { fc, epsc0, .. } => 2.0 * fc / epsc0,
            Material::Hysteretic(h) => h.e1p,
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
    pub(crate) fn evaluate(&self, strain: f64) -> (f64, f64, Material) {
        match self {
            Material::Elastic { e } => elastic_like::evaluate_elastic(*e, strain),
            Material::ElasticPP { e, eyp, ep } => elastic_like::evaluate_elastic_pp(*e, *eyp, *ep, strain),
            Material::Gap { e, gap } => elastic_like::evaluate_gap(*e, *gap, strain),
            Material::Ent { e } => elastic_like::evaluate_ent(*e, strain),
            Material::Steel01 { .. } => steel01::evaluate(self, strain),
            Material::Concrete01 { .. } => concrete01::evaluate(self, strain),
            Material::Steel02 { .. } => steel02::evaluate(self, strain),
            Material::Concrete02 { .. } => concrete02::evaluate(self, strain),
            Material::Hysteretic(_) => hysteretic::evaluate(self, strain),
            Material::Parallel(children) => composite::evaluate_parallel(children, strain),
            Material::Series { .. } => composite::evaluate_series(self, strain),
            Material::MinMax {
                inner,
                min_strain,
                max_strain,
                failed,
            } => composite::evaluate_min_max(inner, *min_strain, *max_strain, *failed, strain),
        }
    }
}
