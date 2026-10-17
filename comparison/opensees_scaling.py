"""OpenSees counterpart of core/examples/benchmark_scaling.rs (same model, loads and solver settings).

Timings cover only the analyze/step loop with history recording, excluding Python startup and model build.
"""

import math
from time import perf_counter

import openseespy.opensees as ops

STORY_HEIGHT, BAY_WIDTH, E = 3000.0, 6000.0, 200000.0
COLUMN, BEAM = (8000.0, 2.0e8), (6000.0, 3.0e8)



def run(stories, bays, steps, fibers=20, points=5, load_base=1.0, displacement_increment=None,
        system="SparseGeneral -piv", numberer="Plain"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nodes = {}
    for s in range(stories + 1):
        for b in range(bays + 1):
            tag = s * (bays + 1) + b + 1
            nodes[s, b] = tag
            ops.node(tag, b * BAY_WIDTH, s * STORY_HEIGHT)
            if s == 0:
                ops.fix(tag, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Steel02", 1, 355.0, E, 0.02, 18.5, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0)
    for tag, (area, inertia) in [(1, COLUMN), (2, BEAM)]:
        ops.section("Fiber", tag)
        depth = math.sqrt(12.0 * inertia / area)
        for i in range(fibers):
            ops.fiber(-depth / 2.0 + (i + 0.5) * depth / fibers, 0.0, area / fibers, 1)
        ops.beamIntegration("Legendre", tag, tag, points)
    tag = 0
    for s in range(stories):
        for b in range(bays + 1):
            tag += 1
            ops.element("dispBeamColumn", tag, nodes[s, b], nodes[s + 1, b], 1, 1)
    for s in range(1, stories + 1):
        for b in range(bays):
            tag += 1
            ops.element("dispBeamColumn", tag, nodes[s, b], nodes[s, b + 1], 1, 2)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in range(1, stories + 1):
        ops.load(nodes[s, 0], load_base * s, 0.0, 0.0)
    roof = nodes[stories, 0]
    ops.constraints("Plain")
    ops.numberer(numberer)
    ops.system(*system.split())
    ops.test("NormDispIncr", 1e-6, 30)
    ops.algorithm("Newton")
    if displacement_increment is None:
        ops.integrator("LoadControl", 1.0 / steps)
    else:
        ops.integrator("DisplacementControl", roof, 1, displacement_increment)
    ops.analysis("Static")

    roof_history, load_history, iterations = [], [], 0
    start = perf_counter()
    for step in range(steps):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"OpenSees failed at step {step + 1}")
        iterations += ops.testIter()
        roof_history.append(ops.nodeDisp(roof, 1))
        load_history.append(ops.getLoadFactor(1))
    return {"engine": "opensees", "solve_ms": (perf_counter() - start) * 1000.0, "iterations": iterations,
            "final_roof_disp": roof_history[-1], "roof_disp_history": roof_history,
            "load_factor_history": load_history}
