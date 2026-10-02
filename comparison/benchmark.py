"""Compare native release Carapace with OpenSees SuperLU frame solves.

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
    parser.add_argument("--runs", type=int, default=21)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--displacement-increment", type=float, default=12.0)
    parser.add_argument("--output", type=Path, default=ROOT / "comparison/out/frame_benchmark.json")
    args = parser.parse_args()
    if min(args.stories, args.bays, args.steps, args.runs) < 1 or args.warmups < 0:
        parser.error("dimensions, steps and runs must be positive; warmups must be nonnegative")
    subprocess.run(["cargo", "build", "--release", "-p", "carapace-core",
                    "--example", "benchmark_frame"], cwd=ROOT, check=True)
    report = {"configuration": vars(args).copy(), "cases": []}
    report["configuration"]["output"] = str(args.output)
    for model in ["elastic", "fiber"]:
        for control in ["load", "displacement"]:
            increment = args.displacement_increment if control == "displacement" else None
            command = [str(ROOT / "target/release/examples/benchmark_frame"),
                       "--model", model, "--stories", str(args.stories),
                       "--bays", str(args.bays), "--steps", str(args.steps), "--json"]
            if increment is not None:
                command += ["--displacement-increment", str(increment)]
            samples = {engine: [] for engine in ["carapace", "SuperLU"]}
            for repetition in range(args.warmups + args.runs):
                # Rotate execution order to reduce systematic timing bias.
                engines = list(samples)
                rotation = repetition % len(engines)
                for engine in engines[rotation:] + engines[:rotation]:
                    if engine == "carapace":
                        data = json.loads(subprocess.check_output(command, cwd=ROOT, text=True))
                    else:
                        data = run_opensees(args.stories, args.bays, args.steps,
                                            model, "SparseGeneral", increment)
                    if repetition >= args.warmups:
                        samples[engine].append(data)
            baseline = samples["carapace"][-1]
            summary = {"model": model, "control": control, "engines": {}}
            for engine, runs in samples.items():
                last = runs[-1]
                differences = {}
                for field in ["roof_disp_history", "load_factor_history"]:
                    errors = [abs(a - b) for a, b in zip(baseline[field], last[field])]
                    differences[field] = max(errors)
                    scale = max(1.0, max(abs(value) for value in baseline[field]))
                    if max(errors) > 1e-6 + 1e-8 * scale:
                        raise RuntimeError(f"{model}/{control}/{engine}: {field} mismatch ({max(errors)})")
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
            print(f"{model}/{control}: " + ", ".join(
                f"{engine} {data['median_ms']:.3f} ms ({data['iterations']} iterations)"
                for engine, data in timings.items()), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
