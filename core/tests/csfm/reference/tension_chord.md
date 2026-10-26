# Independent reference: tension chord

Status: written in phase 0, before any Carapace `CsfmRebar` exists (plan
`docs/csfm-plan.md` S4 C). This is the a priori target the eventual
`CsfmRebar::trial_stress_tangent` inverse must reproduce on its TCM branch, to
`rel 1e-8` on the average strain, per the plan. It is derived from the Tension
Chord Model itself (Marti/Alvarez/Kaufmann/Sigrist), not read off calcs.app's
TypeScript.

## Setup

A single bar of diameter `d`, effective reinforcement ratio
`rho = As / Ac,eff`, embedded in a uniform concrete prism and pulled axially to
a uniform crack pattern (stabilized cracking: cracks have formed, cracked bond
applies everywhere, and the crack spacing is constant along the member). Bond
is rigid-plastic: shear stress `tau_b0` (elastic, before first cracking) and
`tau_b1` (cracked, stabilized), both proportional to `fctm`:
`tau_b0 = 2 fctm`, `tau_b1 = fctm` (the same proportionality CSFM's rebar law
uses — this is the standard TCM bond assumption, not a CSFM-specific choice).

## Equilibrium: steel stress profile between cracks

Let `x = 0` at a crack, where the concrete carries zero stress and the bar
carries the full force, stress `sigma0`. Axial equilibrium of the bar alone
(force lost to bond, transferred into the concrete) gives

```
d(As * sigma_s(x))/dx = -tau_b1 * (pi * d)        (bond shear over the bar's perimeter)
=> dsigma_s/dx = -tau_b1 * (pi * d) / As = -4 tau_b1 / d      (As = pi d^2/4)
=> sigma_s(x) = sigma0 - (4 tau_b1 / d) x
```

linear, until bond saturates at the midpoint between cracks, `x = sr/2`, where
by symmetry with the neighboring crack `dsigma_s/dx` reverses:

```
sigma_s_min = sigma0 - (4 tau_b1 / d)(sr/2) = sigma0 - 2 tau_b1 sr / d
```

## Mean strain

The bar stays elastic pointwise (bond only redistributes force along the bar,
it does not change the bar's own `sigma = Es eps` law), so
`eps_s(x) = sigma_s(x)/Es`, linear like the stress. Averaging over the half
spacing (by symmetry, the same as averaging over a full spacing):

```
eps_m(sigma0) = (1/Es) * (sigma0 + sigma_s_min)/2 = sigma0/Es - tau_b1 * sr / (Es * d)
```

This holds for `sigma0` in `[sigma_cr, fy]` — the *stabilized-cracking* branch.
`tension_chord.py::eps_m` implements exactly this, with `tau_b1` and `sr`
computed from `TcmParams`.

## Crack spacing

The theoretical minimum spacing comes from the *first*-crack bond condition
(uncracked bond `tau_b0`, since the member has only just cracked once): the
concrete regains `fctm` a distance `sr0/2` from the first crack, by the same
force-transfer argument used above but with the concrete's own stress instead
of the steel's:

```
sigma_c(x) = (rho/(1-rho)) * (sigma0 - sigma_s(x))     (section force balance: As*sigma0 = As*sigma_s(x) + Ac*sigma_c(x))
fctm = sigma_c(sr0/2) = (rho/(1-rho)) * (2 tau_b0 sr0 / d)
=> sr0 = d * fctm * (1 - rho) / (2 rho * tau_b0)
```

`lambda >= 1` is a calibration factor over this theoretical minimum (the CSFM
UI offers `lambda` in `0.5-1.0`, acknowledging real crack spacing scatters
above the theoretical value): `sr = lambda * sr0`. This matches the sr formula
in the plan exactly, because it **is** the standard TCM crack-spacing
derivation — not a CSFM-specific quantity.

## First-crack stress and the TCM/POM routing boundary

Just before the first crack, the bar and the concrete share one strain
`eps_cr = fctm/Ecm` (uncracked, bonded, both elastic). The bar there carries
`Es * fctm/Ecm` over its area, the concrete carries `fctm` over its net area
`(1-rho)`. Total axial force is continuous across the crack forming, but right
at the new crack the concrete there now carries zero, so everything the
concrete was carrying transfers to the bar:

```
rho * Es * (fctm/Ecm) + (1-rho) * fctm = rho * sigma_cr
sigma_cr = fctm/rho + fctm * (Es/Ecm - 1)
```

Setting `sigma_cr(rho) = fy` and solving for `rho` gives exactly

```
rho_crit = fctm / (fy - (Es/Ecm - 1) * fctm)
```

which is the plan's TCM/POM routing boundary (S1): below `rho_crit` the bar
would need to yield before a second crack could form, so TCM's own assumptions
(stabilized cracking under an elastic bar) break down and routing falls back
to the Pull-Out Model instead. `tension_chord.py::self_checks` verifies this
identity numerically (`rho_crit` from the formula above reproduces
`sigma_cr(rho_crit) == fy`), which is what caught an earlier, cruder
`sigma_cr = fctm*(1-rho)/rho` approximation (it ignored the bar's own elastic
strain before cracking) as inconsistent with the plan's `rho_crit` formula —
the fixed version above is the one checked in.

## What is checked, and the tolerance

- **Internal consistency** (`self_checks`, no external reference needed):
  the slope of `eps_m` in `sigma0` is exactly `1/Es` (the bond term does not
  depend on `sigma0`); tension stiffening strictly reduces the mean strain
  below the bare-steel value at the same crack stress; the forward-then-
  inverse round trip through `sigma_of_eps_m` (bisection) closes to residual
  `< rel 1e-8` (observed `~1e-12` at `tol = 1e-14`); and the `rho_crit`
  identity above.
- **Parity target.** The forward relation `eps_m(sigma0)` above is the
  closed-form average-strain relation `CsfmRebar`'s TCM branch must invert.
  Tolerance: `rel 1e-8` on the round trip `sigma -> eps_m -> sigma`, fixed by
  the plan (S4 C) before any Carapace result exists. This is tight because
  the relation is literally affine in `sigma0` (closed form, no iteration
  needed on the Carapace side either) — `rel 1e-8` only has to absorb
  floating-point rounding, not model uncertainty.

## Running it

```
python3 -I core/tests/csfm/reference/tension_chord.py
```

No dependencies beyond the standard library (numpy is not installed in this
environment; none of this needs it — everything is a scalar closed form).
