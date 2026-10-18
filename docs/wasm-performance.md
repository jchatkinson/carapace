# Wasm build and performance

How the shipped wasm is built and what was measured while choosing the flags. The build is
`scripts/build-wasm.sh`; PySees (`../pysees/scripts/build-carapace.sh`), CI and the comparison harness all use it,
so every one of them runs the same binary. Only speed matters here: bundle size is not a concern (all variants
below are 1.7-2.2 MB).

## Decisions

- **No `wasm-opt`.** Measured under Node/V8 (the Chrome engine), `wasm-opt -O2` and above make the solve 30-45%
  slower than the unoptimized module. `-O1` matches it. The likely cause is that binaryen's heavier passes give V8's
  JIT code it optimizes worse; that was not investigated. Do not add `wasm-opt -O3` back on the grounds that "release
  wasm should be optimized".
- **`+simd128` stays on, but is not known to help.** It made no measurable difference on the frame models. The Quad4
  panel and continuum models were not measured. faer's runtime SIMD dispatch covers x86 and ARM only, so on wasm
  its kernels are probably scalar whatever the flags.
- **Binaryen older than 120 breaks the module.** The apt package (v108) mis-emits the table section: the
  `__wbindgen_externrefs` export points at the 345-entry funcref table, which has a fixed maximum, so the glue's
  `table.grow(4)` throws on load, in Node and in browsers. The build script refuses old binaryen when `wasm-opt` is
  enabled and always instantiates the result afterwards.
- **Release profile:** fat LTO, `codegen-units = 1`, `panic = "abort"`, plus `opt-level` 3 (the default).

## Measurements

Frame benchmark (`comparison/benchmark.py`: 10 stories, 3 bays, 50 steps), Node 20.19, WSL2, median of interleaved
runs. "Warm" is the solve after 5 untimed in-process runs; "cold" is a fresh Node process, so it includes JIT
compilation.

Optimization level (with simd128), warm solve in ms:

| `wasm-opt` | Size (MB) | Elastic frame | Fiber frame |
|---|---|---|---|
| none | 2.17 | 29.0 | 47.8 |
| `-O1` | 1.76 | 28.3 | 46.7 |
| `-O2` | 1.74 | 40.6 | 57.4 |
| `-O3` | 1.73 | 42.9 | 57.6 |
| `-Os` | 1.72 | 41.4 | 56.1 |
| `-Oz` | 1.71 | 43.9 | 56.9 |

SIMD (no `wasm-opt`): elastic 28.7 with, 27.7 without; fiber 49.1 with, 47.3 without. Within noise.

Against native and OpenSees, shipped build (no `wasm-opt`; `comparison/benchmark.py --runs 9 --wasm-warm-runs 5`, ms,
10 x 3 frame; native and OpenSees do not depend on the wasm build). OpenSees runs `SparseGeneral` (SuperLU) and
`BandSPD` + RCM, the fastest OpenSees configuration in the published benchmarks:

| Case | Native | Wasm warm | OpenSees SuperLU | OpenSees BandSPD |
|---|---|---|---|---|
| elastic, load | 16.6 | 30.7 | 45.3 | 11.3 |
| elastic, displacement | 16.5 | 29.0 | 44.6 | 10.0 |
| fiber, load | 25.9 | 42.8 | 63.7 | 25.7 |
| fiber, displacement | 27.7 | 47.2 | 71.1 | 26.0 |

Warm wasm is 1.6-1.8x native and faster than SuperLU, but 2.7-2.9x slower than BandSPD on the elastic frame. Native
Carapace is also slower than BandSPD there (about 1.5x) and level with it on the fiber frame: this frame is small
(120 free DOFs) and narrowly banded, which suits a band solver better than a sparse Cholesky.

Large frame, 80 x 8 (2,160 free DOFs, 1,360 elements, 20 steps; `benchmark.py --runs 3 --wasm-warm-runs 2`, ms).
The fiber frame uses 400 N per story index and 110 mm per displacement step: the small frame's 25,000 N overloads
80 stories and the native engine fails to converge above about 1,000 N. Roof drift is about 0.9%, so the fiber
sections are only mildly nonlinear here, unlike the published nonlinear benchmark, which drives to 2% drift.

| Case | Native | Wasm warm | OpenSees SuperLU | OpenSees BandSPD |
|---|---|---|---|---|
| elastic, load | 181.9 | 255.1 | 496.8 | 111.5 |
| elastic, displacement | 165.7 | 253.4 | 486.4 | 112.2 |
| fiber, load | 242.1 | 360.5 | 709.9 | 238.2 |
| fiber, displacement | 248.3 | 365.2 | 702.7 | 274.0 |

On the large frame warm wasm is about 1.4-1.5x native (the small frame's ratio is 1.6-2.3x; it shrinks as the solve
dominates), 1.4-2x faster than SuperLU, and 1.3-2.3x slower than BandSPD, which wins the elastic cases by 1.5x over
native and ties the fiber cases. The small-frame figures in this run (3 timed runs) are noisier than the 9-run table
above.

Cold wasm (a fresh Node process, so it includes JIT compilation; earlier run) is 112-154 ms, 5-6.6x native; a cold
first solve costs roughly 80-110 ms more than a warm one, so a warm-up solve on worker start may be worth testing.
"Warm" means 5 untimed in-process solves for wasm only; every engine gets one discarded repetition, OpenSees runs in
one long-lived Python process, and native Carapace is a fresh process per run (no JIT to warm).

Caveats: only V8 was measured; Firefox (SpiderMonkey) and Safari (JavaScriptCore) may behave differently and are
untested. Re-measure before changing flags, and do not infer wasm speed from native speed: wasm has no FMA or AVX, no
runtime CPU-feature dispatch, a slower allocator (dlmalloc), and a 1 MB default stack.

## Reproducing

```
scripts/build-wasm.sh                                   # the shipped build into pkg/
WASM_OPT_LEVEL=-O1 WASM_SIMD=0 scripts/build-wasm.sh /tmp/variant   # a variant (wasm-opt needs binaryen >= 120)
CARAPACE_PKG=/tmp/variant node comparison/wasm_frame.mjs --model fiber --displacement-increment 12 --warm-runs 5
comparison/.venv/bin/python comparison/benchmark.py --wasm-warm-runs 5   # native vs wasm vs OpenSees; small 10 x 3 and
                                                                         # large 80 x 8 frames (--no-large for small only)
```

`comparison/load_wasm.mjs` loads a build under Node (the web-target glue plus `initSync`). Interleave variants when
comparing; run-to-run noise is several percent.
