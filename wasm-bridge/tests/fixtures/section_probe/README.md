# Section probe golden fixtures

Generated from CalcsApp's JavaScript sectional-analysis engine (`app/concrete/sectional-analysis/analysis`), as the
parity reference for `wasm-bridge/src/section_probe`. They show that the probe reproduces the engine, not that the
engine is right.

**Never hand-edit these files.** To regenerate them, change the generator and run it from a CalcsApp checkout that sits
next to this repository (the isolines take about a minute):

```
cd ../calcs.app && npm run export:section-probe-golden     # writes ../carapace/wasm-bridge/tests/fixtures/section_probe
```

Every JSON file carries `provenance`: the CalcsApp commit and whether the engine or generator had uncommitted changes.
Commit the CalcsApp source before regenerating, so that `calcsAppDirty` is `false`.

## Conventions

- Units are lb and in. The section is CalcsApp's `test_05` (8 x 8 in column, four bars, Hognestad concrete,
  elastic-plastic steel with 2% hardening) at mesh size 0.25 in: 418 fibers.
- Probe `z` is the engine's `x`. Engine `M2 = sum(f y)` and `M3 = sum(f x)` are the probe's `m2` and `m3`.
- Tension is positive. The engine's `P` is `sum(f)`, so it is already tension positive.
- `axialForce` in `mphi.json` is `-axial * 1000` (the engine's input is kip, compression positive).

## Files

| File | Content | Checked by |
|---|---|---|
| `section.json` | Probe config: concrete, steel and "bar minus displaced concrete" materials; fibers with concrete group 0 and bars group 1 | all tests |
| `states.json` | Every 13th PMM boundary state: angle, top and bottom strain, `n`, `m2`, `m3` | `pmm_states_match_the_engine` |
| `materials.json` | Stress at 85+ strains for 24 material laws (every concrete, steel and general model, with tension, cutoffs and user curves) | `material_laws_match_the_engine` |
| `mphi.json` | Moment-curvature, 150 steps, four axial loads and angles | `moment_curvature_matches_the_engine` |
| `isolines.json` | The PMM grid (72 angles, 53 closed-loop profiles), 20 axial levels and 36 moment angles, with the engine's outer-branch root for each | `isolines_match_the_engine` |
