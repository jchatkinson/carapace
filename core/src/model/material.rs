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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Steel01Loading {
    None,
    Loading,
    Unloading,
}

/// `Steel02`'s `kon` history flag — `Initial` (0) until the first strain
/// increment, then `Positive`/`Negative` (1/2) for which asymptote the
/// current Menegotto-Pinto curve is heading toward. OpenSees also has a
/// `kon == 3` "pinned to `sigini`" state, not ported since `Material`
/// doesn't support `Steel02`'s optional initial-stress offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Steel02Kon {
    Initial,
    Positive,
    Negative,
}

/// `HystereticMaterial`'s `CloadIndicator`/`TloadIndicator` — which
/// direction the material last loaded in (`None` only before the first
/// nonzero strain increment ever).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HystereticLoad {
    None,
    Positive,
    Negative,
}

/// The full field set behind `Material::Hysteretic` (see that variant's
/// doc comment for why it's boxed rather than inlined). `Copy` despite the
/// field count — every field is itself `f64`/`HystereticLoad`, so this is
/// cheap to copy out of its `Box` for the duration of one `evaluate` call,
/// which is simpler than borrow-juggling a destructure through `&Box<_>`.
///
/// The backbone (`mom1p..mom3n`/`rot1p..rot3n`) and pinching/damage
/// parameters (`pinch_x`, `pinch_y`, `damfc1`, `damfc2`, `beta`) are fixed
/// at construction, along with the envelope slopes/limit (`e1p..e3n`,
/// `eup`, `eun`, `energy_a`) `HystereticMaterial::setEnvelope` derives from
/// them once — recomputing those every `evaluate` call would be wasted
/// work, since they never change. History: `rot_max`/`rot_min` are the
/// peak strain reached on each side (`CrotMax`/`CrotMin`); `rot_pu`/
/// `rot_nu` are the pinching reference points captured at the last
/// reversal into negative/positive loading respectively (`CrotPu`/
/// `CrotNu`); `energy_d` is cumulative dissipated energy, used only to
/// compute the damage factor on a reversal; `load` mirrors
/// `CloadIndicator`.
#[derive(Debug, Clone, Copy)]
pub struct HystereticFields {
    pub mom1p: f64,
    pub rot1p: f64,
    pub mom2p: f64,
    pub rot2p: f64,
    pub mom3p: f64,
    pub rot3p: f64,
    pub mom1n: f64,
    pub rot1n: f64,
    pub mom2n: f64,
    pub rot2n: f64,
    pub mom3n: f64,
    pub rot3n: f64,
    pub pinch_x: f64,
    pub pinch_y: f64,
    pub damfc1: f64,
    pub damfc2: f64,
    pub beta: f64,
    pub e1p: f64,
    pub e2p: f64,
    pub e3p: f64,
    pub e1n: f64,
    pub e2n: f64,
    pub e3n: f64,
    pub eup: f64,
    pub eun: f64,
    pub energy_a: f64,
    pub rot_max: f64,
    pub rot_min: f64,
    pub rot_pu: f64,
    pub rot_nu: f64,
    pub energy_d: f64,
    pub load: HystereticLoad,
    pub strain: f64,
    pub stress: f64,
}

/// `HystereticMaterial`'s fixed (non-history) backbone parameters, bundled
/// so `hysteretic_positive_increment`/`hysteretic_negative_increment` don't
/// need a 20-parameter signature each. Mirrors `posEnvlpStress`/
/// `negEnvlpStress`/`posEnvlpTangent`/`negEnvlpTangent`/`posEnvlpRotlim`/
/// `negEnvlpRotlim` from `HystereticMaterial.cpp` exactly.
struct HystereticEnvelope {
    mom1p: f64,
    rot1p: f64,
    mom2p: f64,
    rot2p: f64,
    mom3p: f64,
    rot3p: f64,
    e1p: f64,
    e2p: f64,
    e3p: f64,
    mom1n: f64,
    rot1n: f64,
    mom2n: f64,
    rot2n: f64,
    mom3n: f64,
    rot3n: f64,
    e1n: f64,
    e2n: f64,
    e3n: f64,
}

/// `HystereticMaterial.cpp`'s `POS_INF_STRAIN`/`NEG_INF_STRAIN` sentinels
/// (from `UniaxialMaterial.h`) — a strain limit this large means "no
/// limit," never an actual physical bound.
const HYSTERETIC_POS_INF_STRAIN: f64 = 1.0e16;
const HYSTERETIC_NEG_INF_STRAIN: f64 = -1.0e16;

impl HystereticEnvelope {
    fn pos_stress(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            0.0
        } else if strain <= self.rot1p {
            self.e1p * strain
        } else if strain <= self.rot2p {
            self.mom1p + self.e2p * (strain - self.rot1p)
        } else if strain <= self.rot3p || self.e3p > 0.0 {
            self.mom2p + self.e3p * (strain - self.rot2p)
        } else {
            self.mom3p
        }
    }

    fn neg_stress(&self, strain: f64) -> f64 {
        if strain >= 0.0 {
            0.0
        } else if strain >= self.rot1n {
            self.e1n * strain
        } else if strain >= self.rot2n {
            self.mom1n + self.e2n * (strain - self.rot1n)
        } else if strain >= self.rot3n || self.e3n > 0.0 {
            self.mom2n + self.e3n * (strain - self.rot2n)
        } else {
            self.mom3n
        }
    }

    fn pos_tangent(&self, strain: f64) -> f64 {
        if strain < 0.0 {
            self.e1p * 1.0e-9
        } else if strain <= self.rot1p {
            self.e1p
        } else if strain <= self.rot2p {
            self.e2p
        } else if strain <= self.rot3p || self.e3p > 0.0 {
            self.e3p
        } else {
            self.e1p * 1.0e-9
        }
    }

    fn neg_tangent(&self, strain: f64) -> f64 {
        if strain > 0.0 {
            self.e1n * 1.0e-9
        } else if strain >= self.rot1n {
            self.e1n
        } else if strain >= self.rot2n {
            self.e2n
        } else if strain >= self.rot3n || self.e3n > 0.0 {
            self.e3n
        } else {
            self.e1n * 1.0e-9
        }
    }

    /// The strain at which the positive backbone, extended past `strain`,
    /// would cross zero stress on a softening segment — `POS_INF_STRAIN`
    /// if the backbone never softens (or hasn't yet, from `strain`).
    ///
    /// The two branches below returning the same `HYSTERETIC_POS_INF_STRAIN`
    /// literal are checking genuinely different conditions from
    /// `HystereticMaterial.cpp`'s `posEnvlpRotlim` (no strain limit was
    /// ever found vs. the one found doesn't correspond to a real softening
    /// crossing) — collapsing them would lose that distinction.
    #[allow(clippy::if_same_then_else)]
    fn pos_rotlim(&self, strain: f64) -> f64 {
        let mut strain_limit = HYSTERETIC_POS_INF_STRAIN;

        if strain <= self.rot1p {
            return HYSTERETIC_POS_INF_STRAIN;
        }
        if strain > self.rot1p && strain <= self.rot2p && self.e2p < 0.0 {
            strain_limit = self.rot1p - self.mom1p / self.e2p;
        }
        if strain > self.rot2p && self.e3p < 0.0 {
            strain_limit = self.rot2p - self.mom2p / self.e3p;
        }

        if strain_limit == HYSTERETIC_POS_INF_STRAIN {
            HYSTERETIC_POS_INF_STRAIN
        } else if self.pos_stress(strain_limit) > 0.0 {
            HYSTERETIC_POS_INF_STRAIN
        } else {
            strain_limit
        }
    }

    /// See `pos_rotlim` above for why the same-looking branches aren't
    /// collapsed.
    #[allow(clippy::if_same_then_else)]
    fn neg_rotlim(&self, strain: f64) -> f64 {
        let mut strain_limit = HYSTERETIC_NEG_INF_STRAIN;

        if strain >= self.rot1n {
            return HYSTERETIC_NEG_INF_STRAIN;
        }
        if strain < self.rot1n && strain >= self.rot2n && self.e2n < 0.0 {
            strain_limit = self.rot1n - self.mom1n / self.e2n;
        }
        if strain < self.rot2n && self.e3n < 0.0 {
            strain_limit = self.rot2n - self.mom2n / self.e3n;
        }

        if strain_limit == HYSTERETIC_NEG_INF_STRAIN {
            HYSTERETIC_NEG_INF_STRAIN
        } else if self.neg_stress(strain_limit) < 0.0 {
            HYSTERETIC_NEG_INF_STRAIN
        } else {
            strain_limit
        }
    }
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

    /// `sigini` (OpenSees' optional initial-stress offset) is not
    /// supported — see the variant's doc comment. Ten parameters is
    /// `Steel02`'s actual OpenSees signature (minus `sigini`/`density`),
    /// not something to trim for its own sake.
    #[allow(clippy::too_many_arguments)]
    pub fn steel02(fy: f64, e0: f64, b: f64, r0: f64, cr1: f64, cr2: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> Self {
        let epsy = fy / e0;
        Material::Steel02 {
            fy,
            e0,
            b,
            r0,
            cr1,
            cr2,
            a1,
            a2,
            a3,
            a4,
            min_strain: -epsy,
            max_strain: epsy,
            pl_strain: 0.0,
            asymptote_strain: 0.0,
            asymptote_stress: 0.0,
            reversal_strain: 0.0,
            reversal_stress: 0.0,
            kon: Steel02Kon::Initial,
            strain: 0.0,
            stress: 0.0,
        }
    }

    /// `fc`/`epsc0`/`fcu`/`epscu` are normalized to be negative
    /// (compression), matching OpenSees' `Concrete02` constructor. `rat` is
    /// the ratio between the unloading slope at `epscu` and the initial
    /// slope; `ft`/`ets` are the tensile strength and tension-softening
    /// modulus.
    pub fn concrete02(fc: f64, epsc0: f64, fcu: f64, epscu: f64, rat: f64, ft: f64, ets: f64) -> Self {
        let fc = -fc.abs();
        let epsc0 = -epsc0.abs();
        let fcu = -fcu.abs();
        let epscu = -epscu.abs();
        Material::Concrete02 {
            fc,
            epsc0,
            fcu,
            epscu,
            rat,
            ft,
            ets,
            min_strain: 0.0,
            tension_strain: 0.0,
            strain: 0.0,
            stress: 0.0,
        }
    }

    /// The 3-point-per-side backbone constructor (OpenSees also offers a
    /// 2-point constructor that derives the midpoint as the average of the
    /// two given points — not ported, since a caller can just pass that
    /// average as `mom2p`/`rot2p`/`mom2n`/`rot2n` directly).
    ///
    /// Debug builds assert the backbone is monotonic and one-to-one
    /// (`0 < rot1p < rot2p < rot3p`, `rot1n < rot2n < rot3n < 0`) —
    /// OpenSees hard-exits the process on this same precondition failure,
    /// which isn't a library-appropriate response here.
    #[allow(clippy::too_many_arguments)]
    pub fn hysteretic(
        mom1p: f64,
        rot1p: f64,
        mom2p: f64,
        rot2p: f64,
        mom3p: f64,
        rot3p: f64,
        mom1n: f64,
        rot1n: f64,
        mom2n: f64,
        rot2n: f64,
        mom3n: f64,
        rot3n: f64,
        pinch_x: f64,
        pinch_y: f64,
        damfc1: f64,
        damfc2: f64,
        beta: f64,
    ) -> Self {
        debug_assert!(rot1p > 0.0 && rot2p > rot1p && rot3p > rot2p, "positive backbone must be increasing");
        debug_assert!(rot1n < 0.0 && rot2n < rot1n && rot3n < rot2n, "negative backbone must be decreasing");

        let energy_a = 0.5
            * (rot1p * mom1p
                + (rot2p - rot1p) * (mom2p + mom1p)
                + (rot3p - rot2p) * (mom3p + mom2p)
                + rot1n * mom1n
                + (rot2n - rot1n) * (mom2n + mom1n)
                + (rot3n - rot2n) * (mom3n + mom2n));

        let e1p = mom1p / rot1p;
        let e2p = (mom2p - mom1p) / (rot2p - rot1p);
        let e3p = (mom3p - mom2p) / (rot3p - rot2p);
        let e1n = mom1n / rot1n;
        let e2n = (mom2n - mom1n) / (rot2n - rot1n);
        let e3n = (mom3n - mom2n) / (rot3n - rot2n);
        let eup = e1p.max(e2p).max(e3p);
        let eun = e1n.max(e2n).max(e3n);

        Material::Hysteretic(Box::new(HystereticFields {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            pinch_x,
            pinch_y,
            damfc1,
            damfc2,
            beta,
            e1p,
            e2p,
            e3p,
            e1n,
            e2n,
            e3n,
            eup,
            eun,
            energy_a,
            rot_max: 0.0,
            rot_min: 0.0,
            rot_pu: 0.0,
            rot_nu: 0.0,
            energy_d: 0.0,
            load: HystereticLoad::None,
            strain: 0.0,
            stress: 0.0,
        }))
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
            Material::Steel02 { .. } => self.evaluate_steel02(strain),
            Material::Concrete02 { .. } => self.evaluate_concrete02(strain),
            Material::Hysteretic(_) => self.evaluate_hysteretic(strain),
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

    /// Ported from `Steel02::setTrialStrain` (`Steel02.cpp`). OpenSees'
    /// `sigini`-pinned quick return (`kon == 0 || kon == 3` with `|dStrain|`
    /// under a `DBL_EPSILON`-scale tolerance) is not ported — it exists to
    /// let the material sit at a nonzero initial stress before any strain
    /// increment, a feature `Material::steel02` doesn't expose. Letting a
    /// zero (or tiny) `dstrain` fall through the general formula instead is
    /// safe here for the same reason it is in `evaluate_steel01`/
    /// `evaluate_concrete01`: plugging `strain == self.strain` back in is
    /// idempotent, since the committed state was itself produced by this
    /// same formula.
    fn evaluate_steel02(&self, strain: f64) -> (f64, f64, Material) {
        let Material::Steel02 {
            fy,
            e0,
            b,
            r0,
            cr1,
            cr2,
            a1,
            a2,
            a3,
            a4,
            min_strain,
            max_strain,
            pl_strain,
            asymptote_strain,
            asymptote_stress,
            reversal_strain,
            reversal_stress,
            kon,
            strain: cstrain,
            stress: cstress,
        } = self
        else {
            unreachable!()
        };
        let (fy, e0, b, r0, cr1, cr2, a1, a2, a3, a4) = (*fy, *e0, *b, *r0, *cr1, *cr2, *a1, *a2, *a3, *a4);
        let (mut min_strain, mut max_strain, mut pl_strain) = (*min_strain, *max_strain, *pl_strain);
        let (mut asymptote_strain, mut asymptote_stress) = (*asymptote_strain, *asymptote_stress);
        let (mut reversal_strain, mut reversal_stress) = (*reversal_strain, *reversal_stress);
        let mut kon = *kon;
        let (cstrain, cstress) = (*cstrain, *cstress);

        let esh = b * e0;
        let epsy = fy / e0;
        let dstrain = strain - cstrain;

        if kon == Steel02Kon::Initial {
            max_strain = epsy;
            min_strain = -epsy;
            if dstrain < 0.0 {
                kon = Steel02Kon::Negative;
                asymptote_strain = min_strain;
                asymptote_stress = -fy;
                pl_strain = min_strain;
            } else {
                kon = Steel02Kon::Positive;
                asymptote_strain = max_strain;
                asymptote_stress = fy;
                pl_strain = max_strain;
            }
        }

        // Load reversal: recompute the intersection of the elastic and
        // strain-hardening asymptotes, shifted by the isotropic-hardening
        // factor (a1/a2 compression side, a3/a4 tension side) — see
        // `Steel02.cpp`'s comments on the corresponding branches.
        if kon == Steel02Kon::Negative && dstrain > 0.0 {
            kon = Steel02Kon::Positive;
            reversal_strain = cstrain;
            reversal_stress = cstress;
            if cstrain < min_strain {
                min_strain = cstrain;
            }
            let d1 = (max_strain - min_strain) / (2.0 * (a4 * epsy));
            let shift = 1.0 + a3 * d1.powf(0.8);
            asymptote_strain = (fy * shift - esh * epsy * shift - reversal_stress + e0 * reversal_strain) / (e0 - esh);
            asymptote_stress = fy * shift + esh * (asymptote_strain - epsy * shift);
            pl_strain = max_strain;
        } else if kon == Steel02Kon::Positive && dstrain < 0.0 {
            kon = Steel02Kon::Negative;
            reversal_strain = cstrain;
            reversal_stress = cstress;
            if cstrain > max_strain {
                max_strain = cstrain;
            }
            let d1 = (max_strain - min_strain) / (2.0 * (a2 * epsy));
            let shift = 1.0 + a1 * d1.powf(0.8);
            asymptote_strain = (-fy * shift + esh * epsy * shift - reversal_stress + e0 * reversal_strain) / (e0 - esh);
            asymptote_stress = -fy * shift + esh * (asymptote_strain + epsy * shift);
            pl_strain = min_strain;
        }

        // Menegotto-Pinto transition curve between the reversal point
        // `(reversal_strain, reversal_stress)` and the asymptote
        // intersection `(asymptote_strain, asymptote_stress)`; `xi` (the
        // normalized plastic excursion) softens the curvature `r` further
        // from `r0` the larger the excursion has been.
        let xi = ((pl_strain - asymptote_strain) / epsy).abs();
        let r = r0 * (1.0 - (cr1 * xi) / (cr2 + xi));
        let eps_ratio = (strain - reversal_strain) / (asymptote_strain - reversal_strain);
        let dum1 = 1.0 + eps_ratio.abs().powf(r);
        let dum2 = dum1.powf(1.0 / r);

        let stress = (b * eps_ratio + (1.0 - b) * eps_ratio / dum2) * (asymptote_stress - reversal_stress) + reversal_stress;
        let tangent = (b + (1.0 - b) / (dum1 * dum2)) * (asymptote_stress - reversal_stress) / (asymptote_strain - reversal_strain);

        (
            stress,
            tangent,
            Material::Steel02 {
                fy,
                e0,
                b,
                r0,
                cr1,
                cr2,
                a1,
                a2,
                a3,
                a4,
                min_strain,
                max_strain,
                pl_strain,
                asymptote_strain,
                asymptote_stress,
                reversal_strain,
                reversal_stress,
                kon,
                strain,
                stress,
            },
        )
    }

    /// Ported from `Concrete02::setTrialStrain` (`Concrete02.cpp`).
    fn evaluate_concrete02(&self, strain: f64) -> (f64, f64, Material) {
        let Material::Concrete02 {
            fc,
            epsc0,
            fcu,
            epscu,
            rat,
            ft,
            ets,
            min_strain,
            tension_strain,
            strain: cstrain,
            stress: cstress,
        } = self
        else {
            unreachable!()
        };
        let (fc, epsc0, fcu, epscu, rat, ft, ets) = (*fc, *epsc0, *fcu, *epscu, *rat, *ft, *ets);
        let (mut min_strain, mut tension_strain) = (*min_strain, *tension_strain);
        let (cstrain, cstress) = (*cstrain, *cstress);

        let ec0 = 2.0 * fc / epsc0;

        let (stress, tangent) = if strain < min_strain {
            // Further into compression than ever committed: monotonic
            // compression envelope, and a new low-water mark.
            let (stress, tangent) = concrete02_compr_envlp(strain, fc, epsc0, fcu, epscu);
            min_strain = strain;
            (stress, tangent)
        } else {
            // Unloading/reloading branch: point R (`reload_strain`,
            // `reload_stress`) is where the reloading slope, extended from
            // the previous minimum, would cross the *initial* tangent line
            // — the classic EERC-report construction (`Concrete02.cpp`'s
            // comments cite Eq. 2.31-2.36).
            let reload_strain = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
            let reload_stress = ec0 * reload_strain;
            let (min_stress, _) = concrete02_compr_envlp(min_strain, fc, epsc0, fcu, epscu);
            let reload_slope = (min_stress - reload_stress) / (min_strain - reload_strain);
            // Strain where the reloading line (through the previous
            // minimum) crosses zero stress.
            let zero_stress_strain = min_strain - min_stress / reload_slope;

            if strain <= zero_stress_strain {
                let stress_lower_bound = min_stress + reload_slope * (strain - min_strain);
                let stress_upper_bound = reload_slope * 0.5 * (strain - zero_stress_strain);
                let mut stress = cstress + ec0 * (strain - cstrain);
                let mut tangent = ec0;
                if stress <= stress_lower_bound {
                    stress = stress_lower_bound;
                    tangent = reload_slope;
                }
                if stress >= stress_upper_bound {
                    stress = stress_upper_bound;
                    tangent = 0.5 * reload_slope;
                }
                (stress, tangent)
            } else {
                // Tension side: `tension_peak_strain` is the strain at the
                // peak of the tensile stress-strain relation reached so
                // far, shifted by `zero_stress_strain` (Eq. 2.42-2.43).
                let tension_peak_strain = zero_stress_strain + tension_strain;
                if strain <= tension_peak_strain {
                    let (peak_stress, _) = concrete02_tens_envlp(tension_strain, ft, ec0, ets);
                    let tangent = if tension_strain != 0.0 { peak_stress / tension_strain } else { ec0 };
                    (tangent * (strain - zero_stress_strain), tangent)
                } else {
                    tension_strain = strain - zero_stress_strain;
                    concrete02_tens_envlp(tension_strain, ft, ec0, ets)
                }
            }
        };

        (
            stress,
            tangent,
            Material::Concrete02 {
                fc,
                epsc0,
                fcu,
                epscu,
                rat,
                ft,
                ets,
                min_strain,
                tension_strain,
                strain,
                stress,
            },
        )
    }

    /// Ported from `HystereticMaterial::setTrialStrain`
    /// (`HystereticMaterial.cpp`). Two of OpenSees' early-return guards
    /// aren't ported, for reasons consistent with the other materials
    /// above:
    ///
    /// - `TloadIndicator == 0 && strain == 0.0`: returns the untouched
    ///   initial `(Tstress, Ttangent) = (0.0, E1p)` — the general formula
    ///   below reproduces this exactly anyway (checked by
    ///   `hysteretic_reproduces_initial_tangent_at_zero_strain`).
    /// - `|dStrain| < DBL_EPSILON`: exists to skip re-deriving a stress/
    ///   tangent OpenSees would just recompute identically; the one place
    ///   this isn't quite idempotent (unlike every other material here) is
    ///   the interior branch's `dStrain < 0.0` / `dStrain > 0.0` split,
    ///   which has no case at all for `dStrain == 0.0` — see
    ///   `evaluate_hysteretic`'s own comment below for how that's resolved.
    fn evaluate_hysteretic(&self, strain: f64) -> (f64, f64, Material) {
        let Material::Hysteretic(fields) = self else {
            unreachable!()
        };
        let HystereticFields {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            pinch_x,
            pinch_y,
            damfc1,
            damfc2,
            beta,
            e1p,
            e2p,
            e3p,
            e1n,
            e2n,
            e3n,
            eup,
            eun,
            energy_a,
            mut rot_max,
            mut rot_min,
            mut rot_pu,
            mut rot_nu,
            energy_d: energy_d_committed,
            mut load,
            strain: cstrain,
            stress: cstress,
        } = **fields;

        let envelope = HystereticEnvelope {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            e1p,
            e2p,
            e3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            e1n,
            e2n,
            e3n,
        };

        let dstrain = strain - cstrain;

        if load == HystereticLoad::None {
            load = if dstrain < 0.0 { HystereticLoad::Negative } else { HystereticLoad::Positive };
        }

        let (stress, tangent) = if strain >= rot_max {
            rot_max = strain;
            load = HystereticLoad::Positive;
            (envelope.pos_stress(strain), envelope.pos_tangent(strain))
        } else if strain <= rot_min {
            rot_min = strain;
            load = HystereticLoad::Negative;
            (envelope.neg_stress(strain), envelope.neg_tangent(strain))
        } else if dstrain < 0.0 {
            let (stress, tangent, new_rot_min, new_rot_pu, new_load) = hysteretic_negative_increment(
                dstrain, strain, cstrain, cstress, rot_max, rot_min, rot_pu, load, rot1p, rot1n, beta, eup, eun,
                pinch_x, pinch_y, damfc1, damfc2, energy_a, energy_d_committed, &envelope,
            );
            rot_min = new_rot_min;
            rot_pu = new_rot_pu;
            load = new_load;
            (stress, tangent)
        } else {
            // `dstrain >= 0.0`: OpenSees only calls `positiveIncrement` for
            // `dStrain > 0.0` (the `dStrain == 0.0` case falls through with
            // whatever `Tstress`/`Ttangent` happened to hold from before —
            // meaningless for a pure function with no "before"). Treating
            // `dstrain == 0.0` as belonging to this branch instead gives a
            // well-defined answer that's consistent with the tie-break
            // this same file already uses for the *first ever* increment
            // (`dStrain < 0.0 ? Negative : Positive`, just above) — and
            // `hysteretic_positive_increment` handles `dstrain == 0.0`
            // gracefully (`tmpmo1` reduces to `cstress` exactly).
            let (stress, tangent, new_rot_max, new_rot_nu, new_load) = hysteretic_positive_increment(
                dstrain, strain, cstrain, cstress, rot_max, rot_min, rot_nu, load, rot1p, rot1n, beta, eup, eun,
                pinch_x, pinch_y, damfc1, damfc2, energy_a, energy_d_committed, &envelope,
            );
            rot_max = new_rot_max;
            rot_nu = new_rot_nu;
            load = new_load;
            (stress, tangent)
        };

        let energy_d = energy_d_committed + 0.5 * (cstress + stress) * dstrain;

        (
            stress,
            tangent,
            Material::Hysteretic(Box::new(HystereticFields {
                mom1p,
                rot1p,
                mom2p,
                rot2p,
                mom3p,
                rot3p,
                mom1n,
                rot1n,
                mom2n,
                rot2n,
                mom3n,
                rot3n,
                pinch_x,
                pinch_y,
                damfc1,
                damfc2,
                beta,
                e1p,
                e2p,
                e3p,
                e1n,
                e2n,
                e3n,
                eup,
                eun,
                energy_a,
                rot_max,
                rot_min,
                rot_pu,
                rot_nu,
                energy_d,
                load,
                strain,
                stress,
            })),
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

/// `Concrete02::Compr_Envlp`: same Kent-Scott-Park compression envelope as
/// `concrete01_envelope`, but note the flat crushed branch's tangent floor
/// here is `1e-10` (matching `Concrete02.cpp`), not the exact `0.0`
/// `concrete01_envelope` uses for the same branch — a genuine difference
/// between the two ports, not a simplification to reconcile away.
fn concrete02_compr_envlp(strain: f64, fc: f64, epsc0: f64, fcu: f64, epscu: f64) -> (f64, f64) {
    let ec0 = 2.0 * fc / epsc0;
    let ratio = strain / epsc0;
    if strain >= epsc0 {
        (fc * ratio * (2.0 - ratio), ec0 * (1.0 - ratio))
    } else if strain > epscu {
        let tangent = (fcu - fc) / (epscu - epsc0);
        (tangent * (strain - epsc0) + fc, tangent)
    } else {
        (fcu, 1.0e-10)
    }
}

/// `Concrete02::Tens_Envlp`: linear elastic up to `ft`, then linear
/// softening at slope `-ets` down to zero, then a flat zero branch (with
/// the same `1e-10` tangent floor as the compression envelope above, for
/// the same singular-tangent-avoidance reason).
fn concrete02_tens_envlp(strain: f64, ft: f64, ec0: f64, ets: f64) -> (f64, f64) {
    let eps0 = ft / ec0;
    let epsu = ft * (1.0 / ets + 1.0 / ec0);
    if strain <= eps0 {
        (strain * ec0, ec0)
    } else if strain <= epsu {
        (ft - ets * (strain - eps0), -ets)
    } else {
        (0.0, 1.0e-10)
    }
}

/// Ported from `HystereticMaterial::positiveIncrement`
/// (`HystereticMaterial.cpp`). Called only for the "interior" branch —
/// `rot_min < strain < rot_max` hasn't been hit yet by the caller's
/// envelope-extreme checks. `rot_max`/`rot_min`/`rot_nu`/`load` are the
/// values *entering* this call (`rot_min` and `rot_max` are read but only
/// `rot_max` is ever updated here); returns
/// `(stress, tangent, new_rot_max, new_rot_nu, new_load)`.
#[allow(clippy::too_many_arguments)]
fn hysteretic_positive_increment(
    dstrain: f64,
    strain: f64,
    cstrain: f64,
    cstress: f64,
    rot_max: f64,
    rot_min: f64,
    rot_nu: f64,
    mut load: HystereticLoad,
    rot1p: f64,
    rot1n: f64,
    beta: f64,
    eup: f64,
    eun: f64,
    pinch_x: f64,
    pinch_y: f64,
    damfc1: f64,
    damfc2: f64,
    energy_a: f64,
    energy_d: f64,
    envelope: &HystereticEnvelope,
) -> (f64, f64, f64, f64, HystereticLoad) {
    let kn_raw = (rot_min / rot1n).powf(beta);
    let kn = if kn_raw < 1.0 { 1.0 } else { 1.0 / kn_raw };
    let kp_raw = (rot_max / rot1p).powf(beta);
    let kp = if kp_raw < 1.0 { 1.0 } else { 1.0 / kp_raw };

    let mut rot_max = rot_max;
    let mut rot_nu = rot_nu;

    // Reversal from negative to positive loading: capture the pinching
    // reference point `rot_nu` and apply ductility-based damage to
    // `rot_max`, but only if the material had actually crossed into
    // tension-side stress (`cstress <= 0.0`) before reversing.
    if load == HystereticLoad::Negative && cstress <= 0.0 {
        rot_nu = cstrain - cstress / (eun * kn);
        let energy = energy_d - 0.5 * cstress / (eun * kn) * cstress;
        let mut damfc = 0.0;
        if rot_min < rot1n {
            damfc = damfc2 * energy / energy_a;
            damfc += damfc1 * (rot_min - rot1n) / rot1n;
        }
        rot_max *= 1.0 + damfc;
    }
    load = HystereticLoad::Positive;

    if rot_max > HYSTERETIC_POS_INF_STRAIN {
        rot_max = HYSTERETIC_POS_INF_STRAIN;
    }
    rot_max = if rot_max > rot1p { rot_max } else { rot1p };

    let maxmom = envelope.pos_stress(rot_max);
    let rotlim = envelope.neg_rotlim(rot_min);
    let rotrel = if rotlim > rot_nu { rotlim } else { rot_nu };

    let rotmp2 = rot_max - (1.0 - pinch_y) * maxmom / (eup * kp);
    let rotch = rotrel + (rotmp2 - rotrel) * pinch_x;

    let (stress, tangent) = if strain < rot_nu {
        let t = eun * kn;
        let s = cstress + t * dstrain;
        if s >= 0.0 {
            (0.0, eun * 1.0e-9)
        } else {
            (s, t)
        }
    } else if strain < rotch {
        if strain <= rotrel {
            (0.0, eup * 1.0e-9)
        } else {
            let t = maxmom * pinch_y / (rotch - rotrel);
            let tmpmo1 = cstress + eup * kp * dstrain;
            let tmpmo2 = (strain - rotrel) * t;
            if tmpmo1 < tmpmo2 {
                (tmpmo1, eup * kp)
            } else {
                (tmpmo2, t)
            }
        }
    } else {
        let t = (1.0 - pinch_y) * maxmom / (rot_max - rotch);
        let tmpmo1 = cstress + eup * kp * dstrain;
        let tmpmo2 = pinch_y * maxmom + (strain - rotch) * t;
        if tmpmo1 < tmpmo2 {
            (tmpmo1, eup * kp)
        } else {
            (tmpmo2, t)
        }
    };

    (stress, tangent, rot_max, rot_nu, load)
}

/// Ported from `HystereticMaterial::negativeIncrement`
/// (`HystereticMaterial.cpp`) — the mirror image of
/// `hysteretic_positive_increment` above; see its doc comment.
#[allow(clippy::too_many_arguments)]
fn hysteretic_negative_increment(
    dstrain: f64,
    strain: f64,
    cstrain: f64,
    cstress: f64,
    rot_max: f64,
    rot_min: f64,
    rot_pu: f64,
    mut load: HystereticLoad,
    rot1p: f64,
    rot1n: f64,
    beta: f64,
    eup: f64,
    eun: f64,
    pinch_x: f64,
    pinch_y: f64,
    damfc1: f64,
    damfc2: f64,
    energy_a: f64,
    energy_d: f64,
    envelope: &HystereticEnvelope,
) -> (f64, f64, f64, f64, HystereticLoad) {
    let kn_raw = (rot_min / rot1n).powf(beta);
    let kn = if kn_raw < 1.0 { 1.0 } else { 1.0 / kn_raw };
    let kp_raw = (rot_max / rot1p).powf(beta);
    let kp = if kp_raw < 1.0 { 1.0 } else { 1.0 / kp_raw };

    let mut rot_min = rot_min;
    let mut rot_pu = rot_pu;

    if load == HystereticLoad::Positive && cstress >= 0.0 {
        rot_pu = cstrain - cstress / (eup * kp);
        let energy = energy_d - 0.5 * cstress / (eup * kp) * cstress;
        let mut damfc = 0.0;
        if rot_max > rot1p {
            damfc = damfc2 * energy / energy_a;
            damfc += damfc1 * (rot_max - rot1p) / rot1p;
        }
        rot_min *= 1.0 + damfc;
    }
    load = HystereticLoad::Negative;

    if rot_min < HYSTERETIC_NEG_INF_STRAIN {
        rot_min = HYSTERETIC_NEG_INF_STRAIN;
    }
    rot_min = if rot_min < rot1n { rot_min } else { rot1n };

    let minmom = envelope.neg_stress(rot_min);
    let rotlim = envelope.pos_rotlim(rot_max);
    let rotrel = if rotlim < rot_pu { rotlim } else { rot_pu };

    let rotmp2 = rot_min - (1.0 - pinch_y) * minmom / (eun * kn);
    let rotch = rotrel + (rotmp2 - rotrel) * pinch_x;

    let (stress, tangent) = if strain > rot_pu {
        let t = eup * kp;
        let s = cstress + t * dstrain;
        if s <= 0.0 {
            (0.0, eup * 1.0e-9)
        } else {
            (s, t)
        }
    } else if strain > rotch {
        if strain >= rotrel {
            (0.0, eun * 1.0e-9)
        } else {
            let t = minmom * pinch_y / (rotch - rotrel);
            let tmpmo1 = cstress + eun * kn * dstrain;
            let tmpmo2 = (strain - rotrel) * t;
            if tmpmo1 > tmpmo2 {
                (tmpmo1, eun * kn)
            } else {
                (tmpmo2, t)
            }
        }
    } else {
        let t = (1.0 - pinch_y) * minmom / (rot_min - rotch);
        let tmpmo1 = cstress + eun * kn * dstrain;
        let tmpmo2 = pinch_y * minmom + (strain - rotch) * t;
        if tmpmo1 > tmpmo2 {
            (tmpmo1, eun * kn)
        } else {
            (tmpmo2, t)
        }
    };

    (stress, tangent, rot_min, rot_pu, load)
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
    fn steel02_gives_exact_initial_tangent_at_zero_strain() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        // At strain=reversal_strain=0 the Menegotto-Pinto curve's
        // eps_ratio is exactly 0 regardless of R, so this must be exact,
        // not just approximately E0.
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 29000.0));
    }

    #[test]
    fn steel02_tangent_approaches_hardening_slope_far_past_yield() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        let esh = 0.01 * 29000.0;
        // Far enough past the yield asymptote (epsy ≈ 0.00207) that the
        // Menegotto-Pinto curve has converged onto the strain-hardening
        // line.
        let (_stress, tangent) = m.trial_stress_tangent(0.05);
        assert!((tangent - esh).abs() < 1.0, "expected tangent near Esh={esh}, got {tangent}");
    }

    #[test]
    fn steel02_reversal_starts_out_near_elastic_stiffness() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        // Push well past yield and commit, then reverse by a small strain
        // step — right after a reversal the curve is close to its elastic
        // (E0) asymptote, unlike deep into either the loading or
        // hardening branch.
        let m = m.commit(0.02);
        let (_stress, tangent) = m.trial_stress_tangent(0.0199);
        assert!(tangent > 0.9 * 29000.0, "expected near-elastic unloading tangent, got {tangent}");
    }

    #[test]
    fn steel02_tracks_peak_strain_across_a_reversal() {
        let m = Material::steel02(60.0, 29000.0, 0.01, 18.5, 0.925, 0.15, 0.9, 5.0, 0.9, 5.0);
        let m = m.commit(0.02);
        let m = m.commit(-0.01);
        let Material::Steel02 { max_strain, .. } = m else { panic!() };
        // max_strain should have latched onto the 0.02 peak on reversal,
        // not stayed at the virgin epsy ≈ 0.00207.
        assert!((max_strain - 0.02).abs() < 1e-9, "got {max_strain}");
    }

    #[test]
    fn concrete02_matches_compression_envelope_on_first_loading() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let strain = -0.001;
        let (stress, tangent) = m.trial_stress_tangent(strain);
        let eta = strain / -0.002;
        let expected_stress = -4.0 * (2.0 * eta - eta * eta);
        assert!((stress - expected_stress).abs() < 1e-9, "expected {expected_stress}, got {stress}");
        assert!((tangent - 2000.0).abs() < 1e-9, "got {tangent}");
    }

    #[test]
    fn concrete02_is_linear_elastic_in_tension_below_ft() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let ec0 = 2.0 * 4.0 / 0.002;
        // ft/Ec0 = 0.4/4000 = 0.0001, so 0.00005 is within the linear
        // range — unlike Concrete01, Concrete02 carries real tension.
        let (stress, tangent) = m.trial_stress_tangent(0.00005);
        assert!((stress - ec0 * 0.00005).abs() < 1e-9, "got {stress}");
        assert!((tangent - ec0).abs() < 1e-9, "got {tangent}");
    }

    #[test]
    fn concrete02_softens_in_tension_past_ft() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        // Past eps0=ft/Ec0=0.0001: linear softening at slope -Ets from ft.
        let (stress, tangent) = m.trial_stress_tangent(0.0005);
        let expected_stress = 0.4 - 50.0 * (0.0005 - 0.0001);
        assert!((stress - expected_stress).abs() < 1e-9, "expected {expected_stress}, got {stress}");
        assert_eq!(tangent, -50.0);
    }

    #[test]
    fn concrete02_unloads_with_degrading_stiffness_below_initial() {
        let m = Material::concrete02(4.0, 0.002, 3.0, 0.006, 0.1, 0.4, 50.0);
        let ec0 = 2.0 * 4.0 / 0.002;

        let m = m.commit(-0.003);
        let (_stress, unload_tangent) = m.trial_stress_tangent(-0.002);
        assert!(
            unload_tangent < ec0 && unload_tangent > 0.0,
            "reload slope should be positive but degraded below Ec0={ec0}, got {unload_tangent}"
        );
    }

    /// A symmetric, 3-point-per-side, bilinear-then-flat backbone shared by
    /// every `Hysteretic` test below: `mom1p/rot1p=10/0.01` (E1p=1000),
    /// `mom2p/rot2p=15/0.02` (E2p=500), `mom3p/rot3p=15/0.03` (E3p=0, flat
    /// plateau) — negative side the exact mirror. `beta=0` makes the
    /// ductility-degradation factors `kn`/`kp` identically `1.0`
    /// throughout (`x^0 == 1`, and `1.0 < 1.0` is false so the `1.0/kn`
    /// branch never engages either) — deliberately excluded from these
    /// tests to keep the hand computation tractable; `damfc1=damfc2=0`
    /// disables the deformation/energy damage factor for the same reason.
    fn hysteretic_test_material() -> Material {
        Material::hysteretic(
            10.0, 0.01, 15.0, 0.02, 15.0, 0.03, -10.0, -0.01, -15.0, -0.02, -15.0, -0.03, 0.5, 0.5, 0.0, 0.0, 0.0,
        )
    }

    #[test]
    fn hysteretic_gives_exact_initial_tangent_at_zero_strain() {
        let m = hysteretic_test_material();
        assert_eq!(m.trial_stress_tangent(0.0), (0.0, 1000.0));
    }

    #[test]
    fn hysteretic_matches_bilinear_envelope_on_first_loading() {
        let m = hysteretic_test_material();
        // Within the first segment (strain <= rot1p=0.01): E1p=1000.
        assert_eq!(m.trial_stress_tangent(0.005), (5.0, 1000.0));
        // Within the second segment (rot1p < strain <= rot2p=0.02):
        // mom1p + E2p*(strain-rot1p) = 10 + 500*0.005 = 12.5.
        assert_eq!(m.trial_stress_tangent(0.015), (12.5, 500.0));
    }

    #[test]
    fn hysteretic_pinches_on_reload_after_a_full_reversal() {
        let m = hysteretic_test_material();
        // Push to the positive envelope's second breakpoint (rot_max
        // becomes exactly rot2p=0.02, stress=mom2p=15), then reverse all
        // the way past the negative breakpoint too (rot_min becomes
        // exactly rot2n=-0.02, stress=mom2n=-15) — both reach the
        // envelope-extreme branch directly, since rot_min/rot_max both
        // start at 0.0 (any first excursion in either direction is "new").
        let m = m.commit(0.02).commit(-0.02);

        // Now reload from -0.02 toward 0.02 — strictly interior to
        // [rot_min, rot_max], so this hits `hysteretic_positive_increment`
        // for the first time, which captures the pinching reference point
        // `rot_nu = cstrain - cstress/(Eun*kn) = -0.02 - (-15)/1000 =
        // -0.005` off the just-completed reversal.
        //
        // At strain=-0.01 (still below rot_nu=-0.005): pure elastic
        // reloading at slope Eun*kn=1000 from the reversal point:
        // stress = cstress + 1000*dstrain = -15 + 1000*0.01 = -5.
        let (stress, tangent) = m.trial_stress_tangent(-0.01);
        assert!((stress - -5.0).abs() < 1e-9, "got {stress}");
        assert!((tangent - 1000.0).abs() < 1e-9, "got {tangent}");

        // At strain=0.0 (past rot_nu, in the pinched-plateau region below
        // rotch): the pinched-target line `(strain-rotrel)*maxmom*pinchY/
        // (rotch-rotrel)` governs rather than the steep elastic-unload
        // line — the hallmark of pinching (a much softer path back through
        // zero than straight elastic reload would give).
        let (stress, tangent) = m.trial_stress_tangent(0.0);
        assert!((stress - 4.285714285714286).abs() < 1e-9, "got {stress}");
        assert!((tangent - 857.1428571428571).abs() < 1e-9, "got {tangent}");
        assert!(tangent < 1000.0, "pinched stiffness should be softer than the elastic Eun, got {tangent}");
    }

    #[test]
    fn hysteretic_energy_accumulates_as_half_trapezoid_of_stress_over_strain() {
        let m = hysteretic_test_material();
        let m = m.commit(0.01); // exactly rot1p: stress=mom1p=10.
        let Material::Hysteretic(fields) = m else { panic!() };
        let energy_d = fields.energy_d;
        // Triangular area under the elastic segment: 0.5*10*0.01 = 0.05.
        assert!((energy_d - 0.05).abs() < 1e-9, "got {energy_d}");
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
