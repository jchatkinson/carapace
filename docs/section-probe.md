# Section probe

`createSectionProbe(config)` (wasm-bridge `section_probe`) evaluates a fiber section at prescribed strain states
and solves the equilibrium problems built on that, without a `Domain`, nodes, elements or load history. It is the
section analogue of the material probe, built for section previews (PySees) and for CalcsApp's PMM and
moment-curvature analysis.

Evaluation is **stateless**: nothing commits, so the same strain state always returns the same resultants. For the
monotone-from-origin states a section analysis uses, that is the materials' envelopes, and it avoids the unloading
history that a Domain run with an axial preload picks up (Concrete01 unloads plastically, so it disagrees with a
reversible stress-strain law).

## Conventions

- Tension positive. Fiber `(y, z)`; strain `eps = eps0 - y kz + z ky` (as `FiberSection3`).
- `n = sum(f)`, `mz = -sum(y f)`, `my = sum(z f)`, `f = stress * area`.
- Also reported: `m2 = sum(y f) = -mz` and `m3 = sum(z f) = my` (CalcsApp's `M2`/`M3` when its `x` is the probe's `z`).
- `theta` defines `r = y cos(theta) + z sin(theta)`. A strain profile `(et, eb)` has `et` at the largest `r` and
  `eb` at the smallest, linear between. Moment-curvature curvature `kappa` gives `eps = eps0 - kappa r`.

## Config

```ts
createSectionProbe({
  coreMaterials?: MaterialSpec[],     // Carapace arena (OpenSees-compatible), may be empty
  materials: ProbeMaterialSpec[],     // what fibers refer to
  fibers: { y, z, area, material, group? },   // parallel arrays; material indexes `materials`
})
```

`ProbeMaterialSpec`:

| `kind` | Meaning |
|---|---|
| `core` `{ index }` | A `coreMaterials` entry, evaluated from its virgin state |
| `combine` `{ terms: [[materialIndex, factor], ...] }` | Weighted sum of earlier probe materials; a bar fiber is steel minus the host concrete: `[[steel, 1], [concrete, -1]]` |
| `concrete` | CalcsApp laws: `elastic`, `elasticNoTension`, `elasticPlastic`, `parabolic`, `whitney`, `hognestad`, `userCurve`. Fields `fc, ec, ecu, eu?, fu, alpha, beta1?, tensionCapacity, ft, ets, curvePoints, dropAfterLastPoint` (`whitney` requires `beta1`) |
| `steel` | `elastic`, `elasticNoCompression`, `elasticPlastic` (with `hardeningRatio`), `userCurve`; `epsUC`/`epsUT` strain cutoffs |
| `general` | `elastic`, `userCurve`; cutoffs |

`concrete`, `steel` and `general` have **no OpenSees equivalent and exist only in the probe**: they are not
`MaterialSpec` variants, so Carapace's core, arena and wire format stay OpenSees-compatible. Where a Carapace
material matches a law exactly it can be used through `core` instead: Hognestad is `concrete01`, Parabolic is
`concrete01` with `fpcu = fpc`, Hognestad with a tension branch is `concrete02`, and elastic-plastic steel with
hardening is `steel01`.

## Methods

| Method | |
|---|---|
| `evaluate(eps0, kz, ky)` | `{ n, mz, my, m2, m3 }` |
| `evaluateBatch(states, withGroups)` | `(eps0, kz, ky)` triples in; `(n, mz, my)` out, with `withGroups` followed by each group's `(n, mz, my)` (stride `3 * (1 + groups)`) |
| `fiberResponses(eps0, kz, ky)` | `[strain, stress]` per fiber, for contour previews |
| `profileStates(angles, profiles)` | `(eps0, kz, ky)` for every (angle, `(et, eb)` profile) pair |
| `axialExtremes(maxStrain)` | Maximum compressive and tensile axial force over uniform strain, with the strain and resultants there. The compressive one is the concentric capacity `P0`; it is interior to the crushing-limited PMM boundary and its moment is non-zero for a section that is not symmetric about its centroid |
| `buildGrid(angles, profiles)` | Samples the PMM boundary for `solveIsolines`; returns the number of states |
| `gridAxialRange()` | `[min, max]` axial force over the grid |
| `solveIsolines({ levels, momentAngles, options? })` | For each axial force and moment-vector angle, the crushing-limited state with `N = P` and the moment along that angle (the outer root) |
| `momentCurvature({ axial, angle, curvatures, tolerance? })` | `eps0` solved at each curvature for constant axial force |

### Isolines

`profiles` is the closed loop of `(et, eb)` strain profiles around the boundary, in order, and must be symmetric under
`(et, eb) -> (eb, et)` with reversed order (so that `(theta + pi, s)` and `(theta, last - s)` are the same field).
`angles` are the sampled profile angles in `[0, pi)`. `levels` may be in any order; points come back as
`levels.length x momentAngles.length` in request order, `found: false` where no root exists (for example above the
boundary's axial apex). Seeding, continuation and a known limitation near the axial apex are described in
`wasm-bridge/src/section_probe/isolines.rs`.

### Errors

Rejections are `{ kind, ... }`: `invalidFibers`, `invalidMaterial`, `degenerateSection`, `invalidRequest`, `noGrid`.

## Tests

`wasm-bridge/tests/wire/section_probe.rs` checks closed-form cases. `section_probe_parity.rs` checks the probe against
CalcsApp's JavaScript engine using generated fixtures (`wasm-bridge/tests/fixtures/section_probe/README.md`).
