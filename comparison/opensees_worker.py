"""Runs one OpenSees benchmark configuration in this process and prints one JSON line per run (its time and results).

run_benchmarks.py launches this in a subprocess per (case, solver, numberer) so that a configuration that
stalls or crashes OpenSees can be stopped without losing the rest of the sweep.

Usage: opensees_worker.py '{"model": "panel", "kwargs": {...}, "system": "...", "numberer": "...", "runs": 3, "warmups": 1}'
"""

import json
import sys

import opensees_frame
import opensees_panel
import opensees_scaling


def main():
    job = json.loads(sys.argv[1])
    kwargs, system, numberer = job["kwargs"], job["system"], job["numberer"]
    run = {
        "frame": lambda: opensees_frame.run(**kwargs, model="elastic", system=system, numberer=numberer),
        "scaling": lambda: opensees_scaling.run(**kwargs, system=system, numberer=numberer),
        "panel": lambda: opensees_panel.run(**kwargs, system=system, numberer=numberer),
    }[job["model"]]
    for _ in range(job["warmups"] + job["runs"]):
        result = run()
        # One line per run, so the parent can stop a run that exceeds its time limit.
        print(json.dumps({"solve_ms": result["solve_ms"], "result": result}), flush=True)


if __name__ == "__main__":
    main()
