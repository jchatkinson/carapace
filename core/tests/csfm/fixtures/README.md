# CSFM golden fixtures

Generated from the CalcsApp TypeScript CSFM implementation, as the parity reference for the Carapace port
([docs/csfm-plan.md](../../../../docs/csfm-plan.md), phase 0). These fixtures show that the port is faithful. They do
not show that the method is right, which is what the independent references in §4 C of the plan are for.

**Never hand-edit these files.** To regenerate them, change the generator and run it from a CalcsApp checkout that
sits next to this repository:

```
cd ../calcs.app && npm run export:csfm-golden     # writes ../carapace/core/tests/csfm/fixtures
```

Every JSON file carries `provenance`: the CalcsApp commit and whether its `lib/` or `test/` had uncommitted changes.
Commit the CalcsApp source before regenerating, so that `calcsAppDirty` is `false`.

## Conventions

- Units are mm, N and MPa.
- Strain is `[eps_x, eps_y, gamma_xy]` with engineering shear. Stress is `[sigma_x, sigma_y, tau_xy]`. Tangents are
  3x3, row-major.
- Displacements are `[ux0, uy0, ux1, uy1, ...]`.
- Gauss-point states list concrete first, element by element, then rod segments, path by path. Q4 points follow the
  node order of the 2x2 rule, as Carapace's `Quad4` does. T3 has one point. `locations` in each model file names
  every entry.
- Supports in TS are penalty constraints (1e10 × the largest tangent diagonal), and Carapace uses exact supports.
  The difference is measured below.

## Files

| File | Contents |
|---|---|
| `concrete-points.json` | `makeConcreteContinuumLaw` stress and TS finite-difference tangent on a principal-strain grid × angles. It covers both laws, softening on and off, and regularization off, 0.01 and 0.001 at fck 30, plus fck 20 and 50. `tangentComparable` is false at the origin seed (`policy: "origin-seed"`) and wherever the TS stencil can straddle a branch kink. |
| `rebar-points.json` | Rod-law stress and tangent per strain for TCM at λ = 0.5, 0.67 and 1, POM, bare steel, and the deep-beam #4 bar. Each point is tagged with its `region`, `clampActive` and `nearBoundary`. `forward` holds TS forward-relation (stress, strain) pairs for inverse round trips; it excludes the TCM `[0, sigma_min]` interval and the startup bridge. |
| `embedding.json` | `resolveRebarPaths` output: nodes with host element, `xi`, `eta`, `N`, and rods. Cases: crossing, along and on shared edges, a T3 diagonal and through-node path, a distorted Q4, the deep beam, the span beam, and an outside-mesh error case. TS binds each end on its own, while Carapace binds a segment to the host containing its midpoint. |
| `models/*.json` | One TS load path each, under the parity settings: the model, the schedule, the solver settings, and for every increment the displacements, Gauss strain and stress, kc2, status and retries. |
| `summary.json` | Parity outcomes, sensitivity, the γ sweep, and outcomes under the UI's solver settings. |
| `critical-load.json` | `findCriticalLoad` bracketing (plan §4 G) on the deep-beam model: the load fraction at which the concrete tensile-strain criterion is first reached, bisected to `loadFractionTolerance`. Carapace's own bracketing on the same model and stop criteria should reproduce `outcome.loadFraction` within that tolerance. |
| `cracked-tangent.json` | A free-DOF tangent stiffness matrix and right-hand side (plan §4 H), dumped from a converged, cracked state of `span-beam-5kip-gamma-0.001` at full load. Boundary conditions are eliminated exactly (no TS penalty), so it is directly comparable to what Carapace's exact supports produce. `maxAsymmetry` confirms the matrix is genuinely nonsymmetric (kc2 coupling); faer's LU should solve it to `rel` 1e-10 residual. |
| `linear-elastic.json` | `solveLinearElastic` on the deep-beam topology (plan §4 E), at full load: isotropic concrete at the deep-beam's own `Ecm(fck)` (`nu = 0`, the same origin-seed tangent `CsfmConcrete`/`CsfmRebar` fall back to there) and plain elastic rods through the same `rebarPaths`/embedding. TS imposes boundary conditions by penalty; the deep-beam penalty sensitivity measured below bounds how closely this can match Carapace's exact supports. |

Models: `q4-startup`, `patch-2in-q4`, `patch-2in-t3` (first 5%), `patch-8in-{1x1,2x2}-*` (the six option
combinations), `deep-beam` (12x6, 7.5 kip, γ 0.01) and `span-beam-5kip-gamma-0.001` (15x7, 10 kip). All reach their
full requested load.

## Measured at generation (calcs.app `2f521f3`)

**Sensitivity of the parity paths.** This is the largest componentwise `|du| / max|u|` over all 20 increments:

| Model | Penalty 1e10 → 1e11 | Residual tolerance 1e-10 → 1e-12 |
|---|---|---|
| `deep-beam` | 4.2e-12 | 1.4e-11 |
| `span-beam-5kip-gamma-0.001` | 1.5e-12 | 6.2e-12 |

Both are about five orders of magnitude below the plan's parity tolerance (displacements to `rel` 1e-6). Neither
the penalty nor the convergence test can account for a parity failure.

**Span beam γ sweep** (10 kip; hand check about 67 MPa). The maximum bar stress is 6.34 MPa at γ = 0.1, 43.7 at 0.01,
50.5 at 1e-3 and 51.2 at 1e-4. At γ = 0.1 the concrete carries the tie.

**UI solver settings**, recorded as historical outcomes and not as parity:

| Load per point | Outcome | Total iterations | Maximum bar stress (MPa) |
|---|---|---|---|
| 5 kip | 100% | 65 | |
| 15 kip | 100% | 70 | |
| 30 kip | 100% | 75 | |
| 50 kip | aborts at 97.5% | | 473 |

The 50 kip abort is the load-control limit point. It is not a located peak and not a material limit (plan §2.5).

The independent reference calculations for §4 C (the CMM panel, the tension chord) are not golden fixtures - they
don't come from calcs.app - and live instead under
[core/tests/csfm/reference](../reference/README.md).
