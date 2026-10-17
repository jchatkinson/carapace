"""Native Carapace vs OpenSeesPy benchmarks for the docs site.

Writes one JSON file per benchmark to docs-site/benchmarks/ (rendered by docs-site/gen_benchmarks.py).
Run from the repo root: comparison/.venv/bin/python comparison/run_benchmarks.py [--only ID ...] [--quick]

Both engines solve the same model with the same convergence test; a benchmark fails if their results
disagree. Timings are the median of --runs timed repetitions after --warmups discarded ones, with the
engine order alternating between repetitions. Only the solve is timed: process start-up, imports and model
building are excluded.
"""

import argparse
import datetime
import json
import os
import platform
import statistics
import subprocess
from pathlib import Path

import openseespy.opensees as ops

import opensees_frame
import opensees_panel
import opensees_scaling

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs-site" / "benchmarks"
EXAMPLES = ["benchmark_frame", "benchmark_scaling", "benchmark_panel"]
STORY_HEIGHT = 3000.0
TOLERANCE = 1e-6  # max relative difference between the engines' results


def example(name, **flags):
    """Runs a core example with `--flag value` arguments and returns its JSON output."""
    command = [str(ROOT / "target/release/examples" / name), "--json"]
    for flag, value in flags.items():
        command += [f"--{flag.replace('_', '-')}", repr(value) if isinstance(value, float) else str(value)]
    return json.loads(subprocess.check_output(command, cwd=ROOT, text=True))


def machine():
    cpu = next((line.split(":", 1)[1].strip() for line in open("/proc/cpuinfo") if line.startswith("model name")),
               platform.processor())
    def tool(*command):
        return subprocess.run(command, cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return {"cpu": cpu, "threads": os.cpu_count(), "os": platform.platform(),
            "rustc": tool("rustc", "--version"), "openseespy": ops.version(),
            "commit": tool("git", "rev-parse", "--short", "HEAD"),
            "dirty": bool(tool("git", "status", "--porcelain", "--", "core", "wasm-bridge")),
            "date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")}


def relative_difference(a, b):
    a, b = (a, b) if isinstance(a, list) else ([a], [b])
    scale = max(1.0, max(abs(v) for v in a))
    return max(abs(x - y) for x, y in zip(a, b)) / scale


def measure(case, fields, runs, warmups):
    """Times `case["carapace"]` and `case["opensees"]` (callables returning result dicts)."""
    samples = {"carapace": [], "opensees": []}
    for repetition in range(warmups + runs):
        engines = list(samples)
        rotation = repetition % len(engines)
        for engine in engines[rotation:] + engines[:rotation]:
            data = case[engine]()
            if repetition >= warmups:
                samples[engine].append(data)
    reference = samples["carapace"][-1]
    error = max(relative_difference(reference[f], samples["opensees"][-1][f]) for f in fields)
    if error > TOLERANCE:
        raise RuntimeError(f"{case['label']}: engines disagree (relative difference {error:.3g})")
    engines = {}
    for engine, results in samples.items():
        times = [r["solve_ms"] for r in results]
        engines[engine] = {"median_ms": statistics.median(times), "min_ms": min(times), "max_ms": max(times),
                           "iterations": results[-1].get("iterations")}
    return {"label": case["label"], "free_dofs": reference["free_dofs"], "elements": reference["elements"],
            "max_relative_difference": error, "engines": engines}


def frame_nonlinear(steps=50, drift=0.02):
    sizes = [("Small", 20, 4), ("Medium", 40, 6), ("Large", 80, 8)]
    cases = []
    for name, stories, bays in sizes:
        increment = drift * stories * STORY_HEIGHT / steps
        # Load-control target: the load Carapace reaches under displacement control at the same drift.
        target = example("benchmark_scaling", stories=stories, bays=bays, steps=steps,
                         displacement_increment=increment)["load_factor_history"][-1]
        for control in ["displacement", "load"]:
            flags = {"displacement_increment": increment} if control == "displacement" else {"load_base": target}
            cases.append({
                "label": f"{name}, {control} control",
                "carapace": lambda s=stories, b=bays, f=flags: example("benchmark_scaling", stories=s, bays=b, steps=steps, **f),
                "opensees": lambda s=stories, b=bays, f=flags: opensees_scaling.run(s, b, steps, **f)})
    return {
        "title": "Nonlinear frame, increasing size",
        "description": (
            "A fixed-base 2D moment frame of displacement-based beam-columns (5 Legendre points) with 20-fibre "
            "rectangular Steel02 sections (fy 355 MPa, 2% hardening), under inverted-triangle lateral loads on "
            "the left column line. Frames are 20 x 4, 40 x 6 and 80 x 8 (stories x bays) with 3 m storeys and "
            "6 m bays. Each is driven to 2% roof drift in 50 steps, once by displacement control on the roof "
            "and once by load control to the load that displacement-control run reached, so both end deep in "
            "yield. Newton with a displacement-increment test (1e-6, 30 iterations)."),
        "note": "Cost is dominated by fibre state updates, so it grows with element count more than with DOF count.",
        "solver": "OpenSees SparseGeneral -piv (SuperLU), Plain numberer",
        "fields": ["roof_disp_history", "load_factor_history"]}, cases


def frame_elastic(steps=20):
    cases = []
    for name, stories, bays in [("Small", 30, 4), ("Medium", 60, 6), ("Large", 100, 10)]:
        cases.append({
            "label": f"{name}, {stories} x {bays}",
            "carapace": lambda s=stories, b=bays: example("benchmark_frame", model="elastic", stories=s, bays=b, steps=steps),
            "opensees": lambda s=stories, b=bays: opensees_frame.run(s, b, steps, "elastic", "SparseGeneral", None)})
    return {
        "title": "Elastic frame, solver-dominated",
        "description": (
            "Elastic 2D moment frames of 30 x 4, 60 x 6 and 100 x 10 (stories x bays) under inverted-triangle "
            "lateral loads, in 20 load-control steps. There is no material nonlinearity, so element work is "
            "small and the sparse factorisation and solve dominate as the DOF count grows. Newton with a "
            "displacement-increment test (1e-6, 30 iterations)."),
        "note": "Both engines refactorise the tangent on every iteration.",
        "solver": "OpenSees SparseGeneral -piv (SuperLU), Plain numberer",
        "fields": ["roof_disp_history", "load_factor_history"]}, cases


def panel_quad4():
    cases = []
    for n in [40, 80, 140]:
        cases.append({
            "label": f"{n} x {n} mesh",
            "carapace": lambda n=n: example("benchmark_panel", n=n),
            "opensees": lambda n=n: opensees_panel.run(n)})
    return {
        "title": "Quad4 plane-stress panel",
        "description": (
            "An N x N mesh of unit square plane-stress Quad4 elements (E 200 GPa, nu 0.3, t 1) clamped on the "
            "left edge, with a unit downward load on each right-edge node, solved in one Newton load step "
            "(2 iterations). The timed region is assembly, factorisation and solve. Meshes are 40 x 40, "
            "80 x 80 and 140 x 140."),
        "note": "OpenSees uses its fastest sparse direct solver here, with RCM numbering.",
        "solver": "OpenSees UmfPack, RCM numberer",
        "fields": ["tip_uy", "disp_norm"]}, cases


BENCHMARKS = {"frame-nonlinear": frame_nonlinear, "frame-elastic": frame_elastic, "panel-quad4": panel_quad4}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", action="append", choices=list(BENCHMARKS), help="benchmark id (repeatable)")
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--quick", action="store_true", help="2 runs, no warmup (development only)")
    args = parser.parse_args()
    runs, warmups = (2, 0) if args.quick else (args.runs, args.warmups)
    if runs < 1 or warmups < 0:
        parser.error("runs must be positive; warmups must be nonnegative")
    cargo = ["cargo", "build", "--release", "-p", "carapace-core"] + [a for e in EXAMPLES for a in ("--example", e)]
    subprocess.run(cargo, cwd=ROOT, check=True)
    OUT.mkdir(parents=True, exist_ok=True)
    meta = machine()
    for benchmark_id in args.only or list(BENCHMARKS):
        info, cases = BENCHMARKS[benchmark_id]()
        fields = info.pop("fields")
        results = []
        for case in cases:
            result = measure(case, fields, runs, warmups)
            results.append(result)
            c, o = result["engines"]["carapace"], result["engines"]["opensees"]
            print(f"{benchmark_id} / {result['label']} ({result['free_dofs']} DOFs): carapace {c['median_ms']:.1f} ms, "
                  f"opensees {o['median_ms']:.1f} ms, ratio {c['median_ms'] / o['median_ms']:.2f}", flush=True)
        report = {"id": benchmark_id, **info, "runs": runs, "warmups": warmups, "machine": meta, "cases": results}
        (OUT / f"{benchmark_id}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved results to {OUT}")


if __name__ == "__main__":
    main()
