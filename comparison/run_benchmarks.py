"""Native Carapace vs OpenSeesPy benchmarks for the docs site.

Writes one JSON file per benchmark to docs-site/benchmarks/ (rendered by docs-site/gen_benchmarks.py).
Run from the repo root (it takes a long time):
    comparison/.venv/bin/python -u comparison/run_benchmarks.py [--only ID ...] [--quick]

Carapace and OpenSees solve the same model with the same convergence test. Carapace is timed once; OpenSees
is timed with every solver and numberer combination in CONFIGURATIONS, and all of them are published,
so no OpenSees configuration is chosen for the comparison. Only the solve is timed: process start-up,
imports and model building are excluded.

Each OpenSees combination runs in its own subprocess (opensees_worker.py) and is killed once a run exceeds
TIMEOUT_FACTOR x Carapace's median time (at least MIN_TIMEOUT_S), then recorded as timed out. A combination
whose results differ from Carapace's by more than TOLERANCE is recorded as "differs" instead of a time.
"""

import argparse
import datetime
import json
import os
import platform
import select
import statistics
import subprocess
import sys
from pathlib import Path

import openseespy.opensees as ops

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs-site" / "benchmarks"
EXAMPLES = ["benchmark_frame", "benchmark_scaling", "benchmark_panel"]
STORY_HEIGHT = 3000.0
TOLERANCE = 1e-6  # max relative difference between the engines' results
# Numberers each solver is run with. Band and profile (skyline) solvers store a band set by the numbering, so
# only the bandwidth-reducing numberers make sense for them; AMD minimises fill for sparse direct solvers and
# gives them an almost dense band (minutes to hours per run). This rule is about how each solver stores the
# matrix, not about measured speed.
SPARSE_NUMBERERS, BANDED_NUMBERERS = ["Plain", "RCM", "AMD"], ["Plain", "RCM"]
CONFIGURATIONS = [(solver, numberer)
                  for solver, numberers in [("SparseGeneral -piv", SPARSE_NUMBERERS), ("UmfPack", SPARSE_NUMBERERS),
                                            ("BandGeneral", BANDED_NUMBERERS), ("ProfileSPD", BANDED_NUMBERERS),
                                            ("SparseSPD", SPARSE_NUMBERERS), ("BandSPD", BANDED_NUMBERERS),
                                            ("SparseSYM", SPARSE_NUMBERERS)]
                  for numberer in numberers]
TIMEOUT_FACTOR, MIN_TIMEOUT_S = 25, 2.0


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


def time_carapace(run, runs, warmups):
    results = [run() for _ in range(warmups + runs)]
    times = [r["solve_ms"] for r in results[warmups:]]
    return results[-1], {"median_ms": statistics.median(times), "min_ms": min(times), "max_ms": max(times),
                         "iterations": results[-1].get("iterations")}


def time_opensees(job, reference, fields, limit_s):
    """Runs one OpenSees configuration in a subprocess and returns its published record.

    The worker prints a line per run; the process is killed as soon as a run takes longer than the limit.
    """
    record = {"solver": job["system"], "numberer": job["numberer"]}
    command = [sys.executable, str(ROOT / "comparison/opensees_worker.py"), json.dumps(job)]
    process = subprocess.Popen(command, cwd=ROOT / "comparison", stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    runs, wait = [], limit_s + 30  # the first run also covers interpreter start-up and the OpenSees import
    timed_out = False
    while True:
        if not select.select([process.stdout], [], [], wait)[0]:
            timed_out = True
            break
        line = process.stdout.readline()
        if not line:
            break
        if line.startswith("{"):
            runs.append(json.loads(line))
            if runs[-1]["solve_ms"] > limit_s * 1000.0:
                timed_out = True
                break
        wait = limit_s + 5
    if timed_out:
        process.kill()
    process.wait()
    if timed_out:
        return {**record, "status": "timeout", "limit_ms": limit_s * 1000.0}
    if process.returncode != 0 or len(runs) != job["warmups"] + job["runs"]:
        return {**record, "status": "error", "detail": f"exit code {process.returncode} after {len(runs)} runs"}
    times = [r["solve_ms"] for r in runs[job["warmups"]:]]
    result = runs[-1]["result"]
    difference = max(relative_difference(reference[f], result[f]) for f in fields)
    return {**record, "status": "ok" if difference <= TOLERANCE else "differs", "median_ms": statistics.median(times),
            "min_ms": min(times), "max_ms": max(times), "iterations": result.get("iterations"),
            "max_relative_difference": difference}


def measure(case, fields, runs, warmups, opensees_runs, opensees_warmups):
    reference, carapace = time_carapace(case["carapace"], runs, warmups)
    limit_s = max(MIN_TIMEOUT_S, TIMEOUT_FACTOR * carapace["median_ms"] / 1000.0)
    print(f"  {case['label']}: carapace {carapace['median_ms']:.1f} ms (OpenSees limit {limit_s:.1f} s per run)", flush=True)
    configurations = []
    for system, numberer in CONFIGURATIONS:
        job = {"model": case["model"], "kwargs": case["kwargs"], "system": system, "numberer": numberer,
               "runs": opensees_runs, "warmups": opensees_warmups}
        record = time_opensees(job, reference, fields, limit_s)
        configurations.append(record)
        shown = f"{record['median_ms']:.1f} ms" if "median_ms" in record else record["status"]
        print(f"    {system:20s} {numberer:6s} {shown}" + (" (differs)" if record["status"] == "differs" else ""),
              flush=True)
    return {"label": case["label"], "free_dofs": reference["free_dofs"], "elements": reference["elements"],
            "carapace": carapace, "opensees": configurations}


def frame_nonlinear(steps=50, drift=0.02):
    cases = []
    for name, stories, bays in [("Small", 20, 4), ("Medium", 40, 6), ("Large", 80, 8)]:
        increment = drift * stories * STORY_HEIGHT / steps
        # Load-control target: the load Carapace reaches under displacement control at the same drift.
        target = example("benchmark_scaling", stories=stories, bays=bays, steps=steps,
                         displacement_increment=increment)["load_factor_history"][-1]
        for control in ["displacement", "load"]:
            flags = {"displacement_increment": increment} if control == "displacement" else {"load_base": target}
            cases.append({
                "label": f"{name}, {control} control", "model": "scaling",
                "kwargs": {"stories": stories, "bays": bays, "steps": steps, **flags},
                "carapace": lambda s=stories, b=bays, f=flags: example(
                    "benchmark_scaling", stories=s, bays=b, steps=steps, **f)})
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
        "fields": ["roof_disp_history", "load_factor_history"]}, cases


def frame_elastic(steps=20):
    cases = []
    for name, stories, bays in [("Small", 30, 4), ("Medium", 60, 6), ("Large", 100, 10)]:
        cases.append({
            "label": f"{name}, {stories} x {bays}", "model": "frame",
            "kwargs": {"stories": stories, "bays": bays, "steps": steps},
            "carapace": lambda s=stories, b=bays: example("benchmark_frame", model="elastic", stories=s, bays=b,
                                                          steps=steps)})
    return {
        "title": "Elastic frame, solver-dominated",
        "description": (
            "Elastic 2D moment frames of 30 x 4, 60 x 6 and 100 x 10 (stories x bays) under inverted-triangle "
            "lateral loads, in 20 load-control steps. There is no material nonlinearity, so element work is "
            "small and the sparse factorisation and solve dominate as the DOF count grows. Newton with a "
            "displacement-increment test (1e-6, 30 iterations)."),
        "note": "Both engines refactorise the tangent on every iteration. A frame's stiffness matrix is narrowly "
                "banded, which suits OpenSees' band solvers.",
        "fields": ["roof_disp_history", "load_factor_history"]}, cases


def panel_quad4():
    cases = [{"label": f"{n} x {n} mesh", "model": "panel", "kwargs": {"n": n},
              "carapace": lambda n=n: example("benchmark_panel", n=n)} for n in [40, 80, 140]]
    return {
        "title": "Quad4 plane-stress panel",
        "description": (
            "An N x N mesh of unit square plane-stress Quad4 elements (E 200 GPa, nu 0.3, t 1) clamped on the "
            "left edge, with a unit downward load on each right-edge node, solved in one Newton load step "
            "(2 iterations). The timed region is assembly, factorisation and solve. Meshes are 40 x 40, "
            "80 x 80 and 140 x 140."),
        "note": "A 2D mesh is not narrowly banded, so band and profile solvers do poorly here.",
        "fields": ["tip_uy", "disp_norm"]}, cases


BENCHMARKS = {"frame-nonlinear": frame_nonlinear, "frame-elastic": frame_elastic, "panel-quad4": panel_quad4}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", action="append", choices=list(BENCHMARKS), help="benchmark id (repeatable)")
    parser.add_argument("--runs", type=int, default=5, help="timed Carapace runs per case")
    parser.add_argument("--warmups", type=int, default=1, help="discarded Carapace runs per case")
    parser.add_argument("--opensees-runs", type=int, default=3, help="timed OpenSees runs per configuration")
    parser.add_argument("--opensees-warmups", type=int, default=1, help="discarded OpenSees runs per configuration")
    parser.add_argument("--quick", action="store_true", help="1 timed run and no warmup everywhere (development only)")
    args = parser.parse_args()
    if args.quick:
        args.runs = args.opensees_runs = 1
        args.warmups = args.opensees_warmups = 0
    if min(args.runs, args.opensees_runs) < 1 or min(args.warmups, args.opensees_warmups) < 0:
        parser.error("runs must be positive; warmups must be nonnegative")
    cargo = ["cargo", "build", "--release", "-p", "carapace-core"] + [a for e in EXAMPLES for a in ("--example", e)]
    subprocess.run(cargo, cwd=ROOT, check=True)
    OUT.mkdir(parents=True, exist_ok=True)
    meta = machine()
    for benchmark_id in args.only or list(BENCHMARKS):
        print(benchmark_id, flush=True)
        info, cases = BENCHMARKS[benchmark_id]()
        fields = info.pop("fields")
        results = [measure(case, fields, args.runs, args.warmups, args.opensees_runs, args.opensees_warmups)
                   for case in cases]
        report = {"id": benchmark_id, **info, "runs": args.runs, "warmups": args.warmups,
                  "opensees_runs": args.opensees_runs, "opensees_warmups": args.opensees_warmups,
                  "timeout_factor": TIMEOUT_FACTOR, "min_timeout_s": MIN_TIMEOUT_S,
                  "machine": meta, "cases": results}
        (OUT / f"{benchmark_id}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved results to {OUT}")


if __name__ == "__main__":
    main()
