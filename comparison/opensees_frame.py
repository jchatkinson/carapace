"""OpenSees counterpart of core/examples/benchmark_frame.rs.

Timings exclude Python import/process startup. Both engines time one
analyze/step call and history recording per load step.
"""

import argparse
import json
import math
from time import perf_counter

import openseespy.opensees as ops


def run(stories=10, bays=3, steps=50, model="elastic", system="SparseGeneral",
        displacement_increment=None):
    start = perf_counter()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nodes = {}
    for story in range(stories + 1):
        for bay in range(bays + 1):
            tag = story * (bays + 1) + bay + 1
            nodes[story, bay] = tag
            ops.node(tag, bay * 6000.0, story * 3000.0)
            if story == 0:
                ops.fix(tag, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    if model == "fiber":
        ops.uniaxialMaterial("Steel01", 1, 355.0, 200000.0, 0.02,
                            0.0, 1.0, 0.0, 1.0)
        for section, area, inertia in [(1, 8000.0, 2.0e8), (2, 6000.0, 3.0e8)]:
            height = math.sqrt(inertia / area)
            ops.section("Fiber", section)
            ops.fiber(height, 0.0, area / 2.0, 1)
            ops.fiber(-height, 0.0, area / 2.0, 1)
            ops.beamIntegration("Legendre", section, section, 3)

    tag = 0
    for story in range(stories):
        for bay in range(bays + 1):
            tag += 1
            ni, nj = nodes[story, bay], nodes[story + 1, bay]
            if model == "fiber":
                ops.element("dispBeamColumn", tag, ni, nj, 1, 1)
            else:
                ops.element("elasticBeamColumn", tag, ni, nj,
                            8000.0, 200000.0, 2.0e8, 1)
    for story in range(1, stories + 1):
        for bay in range(bays):
            tag += 1
            ni, nj = nodes[story, bay], nodes[story, bay + 1]
            if model == "fiber":
                ops.element("dispBeamColumn", tag, ni, nj, 1, 2)
            else:
                ops.element("elasticBeamColumn", tag, ni, nj,
                            6000.0, 200000.0, 3.0e8, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    base_load = 25000.0 if model == "fiber" else 10000.0
    for story in range(1, stories + 1):
        ops.load(nodes[story, 0], base_load * story, 0.0, 0.0)
    roof = nodes[stories, 0]
    ops.constraints("Plain")
    ops.numberer("Plain")
    if system != "SparseGeneral":
        raise ValueError("this comparison uses SuperLU via SparseGeneral -piv")
    ops.system("SparseGeneral", "-piv")
    ops.test("NormDispIncr", 1e-6, 30)
    ops.algorithm("Newton")
    if displacement_increment is None:
        ops.integrator("LoadControl", 1.0 / steps)
    else:
        ops.integrator("DisplacementControl", roof, 1, displacement_increment)
    ops.analysis("Static")
    setup_ms = (perf_counter() - start) * 1000.0

    roof_history, load_history = [], []
    iterations = 0
    start = perf_counter()
    for step in range(steps):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"OpenSees {model}/{system} failed at step {step + 1}")
        iterations += ops.testIter()
        roof_history.append(ops.nodeDisp(roof, 1))
        load_history.append(ops.getLoadFactor(1))
    solve_ms = (perf_counter() - start) * 1000.0
    return {
        "engine": "opensees", "version": ops.version(), "system": "SparseGeneral -piv (SuperLU)",
        "model": model, "stories": stories, "bays": bays, "steps": steps,
        "free_dofs": stories * (bays + 1) * 3,
        "setup_ms": setup_ms, "solve_ms": solve_ms,
        "total_ms": setup_ms + solve_ms, "iterations": iterations,
        "final_roof_disp": roof_history[-1], "roof_disp_history": roof_history,
        "load_factor_history": load_history,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stories", type=int, default=10)
    parser.add_argument("--bays", type=int, default=3)
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--model", choices=["elastic", "fiber"], default="elastic")
    parser.add_argument("--displacement-increment", type=float)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
