"""Compare native release Carapace, the wasm bundle (under Node) and OpenSees frame solves: SuperLU
(`SparseGeneral`) and `BandSPD` + RCM, the fastest OpenSees configuration in the published benchmarks.

Run: comparison/.venv/bin/python comparison/benchmark.py
Compilation and imports precede all timings; solve timings exclude setup.
"""

import argparse
import json
from pathlib import Path
import statistics
import subprocess

from opensees_frame import run as run_opensees

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stories", type=int, default=10)
    parser.add_argument("--bays", type=int, default=3)
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--large-stories", type=int, default=80,
                        help="stories of the large frame (80 x 8 is the published nonlinear benchmark's Large case)")
    parser.add_argument("--large-bays", type=int, default=8)
    parser.add_argument("--large-steps", type=int, default=20)
    parser.add_argument("--large-load-base", type=float, default=400.0,
                        help="lateral load per story index (N) of the large fiber frame; the small frame's 25000 N "
                             "overloads 80 stories (the native engine fails to converge above ~1000)")
    parser.add_argument("--large-displacement-increment", type=float, default=110.0,
                        help="roof displacement per step (mm) of the large frame; ~0.9%% drift over the steps")
    parser.add_argument("--no-large", action="store_true", help="run only the small frame")
    parser.add_argument("--runs", type=int, default=21)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--wasm-warm-runs", type=int, default=0,
                        help="untimed in-process wasm runs before the timed one (0 = cold JIT, as in a fresh worker)")
    parser.add_argument("--displacement-increment", type=float, default=12.0)
    parser.add_argument("--output", type=Path, default=ROOT / "comparison/out/frame_benchmark.json")
    args = parser.parse_args()
    if min(args.stories, args.bays, args.steps, args.large_stories, args.large_bays, args.large_steps, args.runs) < 1 or args.warmups < 0:
        parser.error("dimensions, steps and runs must be positive; warmups must be nonnegative")
    subprocess.run(["cargo", "build", "--release", "-p", "carapace-core",
                    "--example", "benchmark_frame"], cwd=ROOT, check=True)
    subprocess.run([str(ROOT / "scripts/build-wasm.sh"), str(ROOT / "pkg")],
                   cwd=ROOT, check=True)
    report = {"configuration": vars(args).copy(), "cases": []}
    report["configuration"]["output"] = str(args.output)
    sizes = [("small", args.stories, args.bays, args.steps, args.displacement_increment, None)]
    if not args.no_large:
        sizes.append(("large", args.large_stories, args.large_bays, args.large_steps,
                      args.large_displacement_increment, args.large_load_base))
    for size, stories, bays, steps, displacement_increment, fiber_load_base in sizes:
        for model in ["elastic", "fiber"]:
            for control in ["load", "displacement"]:
                increment = displacement_increment if control == "displacement" else None
                load_base = fiber_load_base if model == "fiber" else None
                command = [str(ROOT / "target/release/examples/benchmark_frame"),
                           "--model", model, "--stories", str(stories),
                           "--bays", str(bays), "--steps", str(steps), "--json"]
                if increment is not None:
                    command += ["--displacement-increment", str(increment)]
                if load_base is not None:
                    command += ["--load-base", str(load_base)]
                wasm_command = ["node", str(ROOT / "comparison/wasm_frame.mjs"), "--model", model,
                                "--stories", str(stories), "--bays", str(bays),
                                "--steps", str(steps), "--warm-runs", str(args.wasm_warm_runs)]
                if increment is not None:
                    wasm_command += ["--displacement-increment", str(increment)]
                if load_base is not None:
                    wasm_command += ["--load-base", str(load_base)]
                samples = {engine: [] for engine in ["carapace", "wasm", "SuperLU", "BandSPD"]}
                for repetition in range(args.warmups + args.runs):
                    # Rotate execution order to reduce systematic timing bias.
                    engines = list(samples)
                    rotation = repetition % len(engines)
                    for engine in engines[rotation:] + engines[:rotation]:
                        if engine == "carapace":
                            data = json.loads(subprocess.check_output(command, cwd=ROOT, text=True))
                        elif engine == "wasm":
                            data = json.loads(subprocess.check_output(wasm_command, cwd=ROOT, text=True))
                        else:
                            system, numberer = ("SparseGeneral", "Plain") if engine == "SuperLU" else ("BandSPD", "RCM")
                            data = run_opensees(stories, bays, steps,
                                                model, system, increment, numberer, load_base)
                        if repetition >= args.warmups:
                            samples[engine].append(data)
                baseline = samples["carapace"][-1]
                summary = {"size": size, "stories": stories, "bays": bays, "steps": steps,
                           "model": model, "control": control, "engines": {}}
                for engine, runs in samples.items():
                    last = runs[-1]
                    differences = {}
                    for field in ["roof_disp_history", "load_factor_history"]:
                        errors = [abs(a - b) for a, b in zip(baseline[field], last[field])]
                        differences[field] = max(errors)
                        scale = max(1.0, max(abs(value) for value in baseline[field]))
                        if max(errors) > 1e-6 + 1e-8 * scale:
                            raise RuntimeError(f"{size}/{model}/{control}/{engine}: {field} mismatch ({max(errors)})")
                    times = [sample["solve_ms"] for sample in runs]
                    summary["engines"][engine] = {
                        "median_ms": statistics.median(times), "min_ms": min(times),
                        "max_ms": max(times), "samples_ms": times,
                        "iterations": last["iterations"],
                        "final_roof_disp": last["final_roof_disp"],
                        "max_absolute_errors": differences,
                        "reference_run": last,
                    }
                report["cases"].append(summary)
                timings = summary["engines"]
                print(f"{size} {stories}x{bays} {model}/{control}: " + ", ".join(
                    f"{engine} {data['median_ms']:.3f} ms"
                    + (f" ({data['iterations']} iterations)" if data["iterations"] is not None else "")
                    for engine, data in timings.items()), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
