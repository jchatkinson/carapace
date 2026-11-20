# Independent reference: CMM panel

Status: written in phase 0, before any Carapace `CsfmConcrete`/`CsfmRebar`/
`EmbeddedBar` exists (plan `docs/csfm-plan.md` S4 C). This is the a priori
target for the combination of all three, once phases 1-2 land a single
uniform-strain element can be checked against. Derived from the Cracked
Membrane Model (Kaufmann & Marti 1998) equilibrium/compatibility equations and
the Tension Chord Model (`tension_chord.md`), not from calcs.app's TypeScript.

## Model

A single membrane element, orthogonal reinforcement in global `x` and `y`,
with displacements prescribed so the strain `[eps_x, eps_y, gamma_xy]`
(engineering shear) is uniform across it — a closed-form check on a strain
state, not a solved model.

**Equilibrium** (rotating-angle, the assumption CMM and CSFM share — one
isotropic concrete law acting at the strain's own principal direction):

```
sigma_x = rho_x * sigma_sx + sigma_c_x
sigma_y = rho_y * sigma_sy + sigma_c_y
tau_xy  =                    tau_c_xy
```

`(sigma_c_x, sigma_c_y, tau_c_xy)` is the ordinary Mohr's-circle rotation of `(sigma1, sigma2)`
back to `x`-`y`: `tau_c_xy = (sigma1 - sigma2)/2 * sin(2 theta)`, the tensor shear component
directly — stress has no engineering/tensor shear distinction (only *strain* does, through
`gamma_xy = 2 * eps_xy`). An earlier version of `cmm_panel.py` used one `from_principal` helper
for both strain and stress and doubled the shear term unconditionally (correct for recovering
`gamma_xy` from a strain decomposition, wrong for stress), so this file's own `tau_c_xy`/`tau_xy`
were exactly double what they should be. Caught when Carapace's `CsfmConcrete::trial_stress_tangent`
— an ordinary tensor stress rotation — returned exactly half of this file's original number; now
split into `from_principal_stress` (undoubled) and `from_principal_strain` (doubled, used only by
`assert_round_trip`'s strain check).

**Compatibility.** Bar strain in each direction equals the membrane's own
normal strain there (smeared reinforcement, shared by CMM and CSFM). Bar
stress is the Tension Chord Model inverse (`tension_chord.py`): both cases
below are chosen with `eps_x`, `eps_y` on the TCM branch
(`eps_m(sigma_cr) <= eps <= eps_m(fy)`), away from the first-crack transition
(the plan's "startup bridge" discussion is a separate, non-closed-form
concern, deliberately avoided here).

**Concrete** (`cmm_panel.py`): zero tension on the major principal direction;
parabola-rectangle compression (peak `-eta_fc*fck` at `eps_c2 = -0.002`, down
to `eps_cu2 = -0.0035`) on the minor direction, scaled by Kaufmann's
compression-softening coefficient
`kc2 = min(1, 1/(0.8 + 140*eps1))` when `eps1 > 0`. All of this is stated
directly from the plan's S1/S2.1 formulas (Eurocode/Kaufmann 1998 constants),
not calcs.app internals — and `kc2` in particular **is** the original
Kaufmann & Marti formula, shared verbatim between CMM and CSFM.

## Cases

Both use `fck = 30` MPa (`Ecm`, `fctm` derived), `Es = 200000` MPa,
`fy = 500`, `ft = 550` MPa, `lambda = 0.67`.

**A — pure shear, symmetric reinforcement.** Equal reinforcement both ways
(`rho_x = rho_y = 0.01`, `d = 16` mm), strain `eps_x = eps_y = 0.0016`,
`gamma_xy = 0.0045`. By symmetry `theta = 45 deg` exactly and the two bar
directions are identical.

**B — biaxial, asymmetric reinforcement.** `rho_x = 0.012` (`d = 16` mm),
`rho_y = 0.008` (`d = 12` mm); strain `eps_x = 0.0020`, `eps_y = 0.0017`,
`gamma_xy = 0.0045`. (`eps_x` is `0.0020`, not `0.0022`: at `0.0022` it falls
outside bar_x's TCM branch under the corrected `sigma_min`-anchored formula
below — see that section's note.)

## Expected quantities (from `cmm_panel.py`, run 2026-10-26, after the
`tension_chord.py` TCM-intercept fix and the `from_principal` shear-doubling fix)

| Quantity | A | B |
|---|---|---|
| `eps1` | 3.850000e-3 | 4.104994e-3 |
| `eps2` | -6.500000e-4 | -4.049945e-4 |
| `theta` | 45.000 deg | 43.093 deg |
| `kc2` | 0.746826 | 0.727432 |
| `sigma1` | 0 MPa | 0 MPa |
| `sigma2` | -12.196602 MPa | -7.943323 MPa |
| `sigma_c_x, sigma_c_y, tau_c_xy` | -6.098301, -6.098301, 6.098301 | -3.707470, -4.235852, 3.962865 |
| `sigma_sx, sigma_sy` | 416.0614, 416.0614 | 479.8894, 460.3193 |
| `sigma_x, sigma_y, tau_xy` | -1.937687, -1.937687, 6.098301 | 2.051203, -0.553298, 3.962865 |

Re-derive, don't hand-copy: run `python3 -I core/tests/csfm/reference/cmm_panel.py`.

## What must match to rounding, and what only within a band

CMM and CSFM differ by design (plan S4 C): CSFM has no concrete tension, and
its crack spacing and `kc2` are simplified relative to a full 2D treatment.
This reference separates the two explicitly.

**To rounding** (pure identities/shared formulas, independent of any
material-law uncertainty — tolerance `abs 1e-9 MPa`, `rel 1e-10`, the plan's
stress default sharpened since there is no modeling slack here at all):

- The Mohr round trip itself (`assert_round_trip`, checked at `1e-12`
  relative in the script): decomposing a strain into principal values and
  angle, then recombining, must reproduce the input exactly.
- `kc2(eps1)`: the original Kaufmann formula, shared verbatim by both models.
- Zero tension on the major principal direction, when `eps1 > 0`: both models
  make exactly this choice for CSFM's concrete law (the plan says CMM and
  CSFM "differ by design" specifically because *classical* CMM does carry
  some tension there, but *this* comparison is explicitly against the CSFM
  concrete law the plan specifies, not the original Kaufmann paper's tension
  branch — so there is no slack to budget here either).
- The equilibrium sum itself: `sigma_x`, `sigma_y`, `tau_xy` as computed above
  from `(sigma_c, sigma_s)` **is** the CMM solution by construction, so a
  Carapace model built to represent the same prescribed-strain state (once
  `CsfmConcrete`/`CsfmRebar`/`EmbeddedBar` exist) that independently assembles
  concrete Gauss stress and bar stress at that point must sum to the same
  panel stress to rounding — this is an equilibrium identity, not a
  material-model comparison.

**Within a band** (the modeling simplification the plan calls out — CSFM
evaluates each bar's tension stiffening from its own 1D crack spacing, not a
crack-angle-corrected 2D spacing):

- **Bar stress**, `sigma_sx` and `sigma_sy`. `crack_spacing_band` in
  `cmm_panel.py` quantifies this: it recomputes each bar's TCM stress using
  the angle-weighted combined spacing
  `1/sr_full(theta) = |cos(theta)|/sr_x + |sin(theta)|/sr_y` (the standard
  Vecchio-Collins-style combination for a crack oblique to both reinforcement
  layers) in place of that bar's own 1D `sr`, and reports the resulting shift.
  Observed: **6.76%** (case A, both directions, by symmetry) and **4.18%** /
  **8.81%** (case B, x / y). The largest observed shift is just under 9%.
- **Tolerance: `rel 10%`** on bar stress, fixed here (before any Carapace
  result exists) with a modest margin over the largest measured shift
  (8.81%), since the two cases above do not exhaustively cover every crack
  angle and reinforcement ratio a real test might use. If a future test needs
  a case further from these two (e.g. a much more oblique crack with very
  unequal reinforcement), recompute `crack_spacing_band` for it before relying
  on this 10% figure, since this bound is case-specific, not universal.
- Concrete principal stresses (`sigma1`, `sigma2`) are NOT in this band —
  they are checked to rounding above, since both models share `kc2` and zero
  tension exactly. Only the *panel-level* `sigma_x`/`sigma_y`/`tau_xy`
  inherit the bar-stress band through the `rho * sigma_s` term; the
  concrete-only contribution to each is still exact.

## Running it

```
python3 -I core/tests/csfm/reference/cmm_panel.py
```

Depends only on `tension_chord.py` (same directory) and the standard library.
