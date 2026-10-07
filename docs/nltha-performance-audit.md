# NLTHA performance audit

Benchmark: `core/examples/benchmark_nltha.rs` (`cargo run --release -p carapace-core --example benchmark_nltha -- stories bays steps pga_g points tangent`).
Multi-story RC frame, `DispBeamColumn` (5-pt Lobatto), 26-fiber RC sections (Concrete02 cover/core + Steel02 rebar), Newmark + Rayleigh 5%, uniform ground motion, `Algorithm::Newton`. Env knobs: `TEST=du TOL=1e-3` (convergence test), `CONC=01|02|el`, `STEEL=02|01|el`.
Reference run: `12 4 400 0.1 5 current` (108 elements, 180 DOF, ~2.2 iterations/step, 83 mm peak roof).

## Where the time goes (callgrind, 6x3 frame, inclusive % of `TransientAnalysis::step`)

| Piece | Share |
|---|---|
| Element state determination (`assemble_stiffness_triplets`, all fiber/material evals) | 66% |
| ...of which `assemble_newmark_system` | 40% |
| ...of which `multiply_stiffness` (a second, redundant assembly) | 33% |
| Material `evaluate` (Concrete02 + Steel02), of which memcpy/drop/clone of the discarded next-state ≈ 23% of total | 38% |
| Sparse LU numeric factorization | 10% |
| `Domain::commit` | 9.5% |
| Triplet to CSC build | 6.6% |
| Per-step `Domain::clone` snapshot | 4.6% |

Factorization is not the bottleneck at these sizes; element/material state determination is.

## Enhancement 1: assemble once per Newton iteration (measured 1.55x)

`transient.rs` `form_system` calls `assemble_newmark_system` and then `multiply_stiffness(&v_trial)`, and both run `assemble_stiffness_triplets`, so every element/fiber is evaluated twice per iteration. The `beta_k * K * v` term can be formed from the same triplets (or skipped when `beta_k == 0`). Same fix applies to `Linear` (`multiply_stiffness` on `damp_vec`) and `recompute_initial_acceleration`.
Follow-ups: a fixed-pattern scatter map (triplet index to CSC slot, built once) removes the per-iteration sort/alloc in `try_new_from_triplets`; the last-iteration `K` is assembled and discarded.

## Enhancement 2: trial evaluations must not build the next material state (measured a further ~1.15x on top of 1)

`Material::trial_stress_tangent` calls `evaluate`, which constructs a ~170-byte next-state `Material` (and `.clone()`s it on early-out paths) only to drop it. That is ~23% of runtime in memcpy/drop glue. Making the leaf evaluators generic over a `NEXT: bool` const (trial returns `None`) removes it with identical numerics. Do this for every leaf material; the same pattern fits Hysteretic/Pinching4, which heap-allocate a `Box` per trial today.

## Combined prototype

`docs/nltha-perf-prototype.patch` (both changes, Steel02/Concrete02 only for #2). Reference run 2.05 s to 1.25 s (1.65x), peak roof displacement identical to 4 digits, `cargo test -p carapace-core` passes.

## Not recommended / lower value
- `ReuseAtStepStart`/`Initial` tangents: saves factorizations but adds iterations; neutral to slightly negative at these sizes (see `tangent` arg).
- Per-step `Domain::clone` snapshot (4.6%) and `commit` (9.5%): worth revisiting after 1 and 2.

## Side finding (not investigated)
With realistic RC sections, Concrete01 and Concrete02 cause Newton two-cycle limit cycles (residual alternating between two values with tiny `du`) at moderate demand (e.g. 0.2g), failing with `FailedToConverge`. Elastic concrete converges fine. Worth a separate look at material tangent consistency at branch switches.
