"""OpenSees counterpart of core/examples/benchmark_panel.rs (same mesh, material, loads and solver settings).

The timed region is the single analyze(1) call: assembly, factorisation and solve.
"""

import math
from time import perf_counter

import openseespy.opensees as ops


def run(n, system="UmfPack"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    at = lambda i, j: j * (n + 1) + i + 1
    for j in range(n + 1):
        for i in range(n + 1):
            ops.node(at(i, j), float(i), float(j))
            if i == 0:
                ops.fix(at(i, j), 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3)
    tag = 0
    for j in range(n):
        for i in range(n):
            tag += 1
            ops.element("quad", tag, at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), 1.0, "PlaneStress", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(n + 1):
        ops.load(at(n, j), 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system.split())
    ops.test("NormDispIncr", 1e-9, 10)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    start = perf_counter()
    if ops.analyze(1) != 0:
        raise RuntimeError("OpenSees panel analysis failed")
    solve_ms = (perf_counter() - start) * 1000.0
    norm = math.sqrt(sum(ops.nodeDisp(at(i, j), d) ** 2 for j in range(n + 1) for i in range(n + 1) for d in (1, 2)))
    return {"engine": "opensees", "solve_ms": solve_ms, "iterations": ops.testIter(),
            "tip_uy": ops.nodeDisp(at(n, n), 2), "disp_norm": norm}
