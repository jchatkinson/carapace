"""Cantilever-panel cross-check of Carapace's Tri3/Quad4 against OpenSees (`tri31`, `quad`,
`enhancedQuad`): same mesh, same plane-stress isotropic material, tip loads. Compares nodal
displacements, support reactions and Gauss-point stresses.

Usage (from the repo root, after `wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg`):
    comparison/.venv/bin/python comparison/panel_compare.py

The plain elements (`quad` vs Carapace `full`, `tri31` vs `tri3`) are the same formulation and
must agree to round-off. `enhancedQuad` is OpenSees' own enhanced-strain variant: agreement is
expected on rectangles and is reported, not assumed, on distorted meshes.
"""

import json
import math
import pathlib
import subprocess
import sys
import tempfile

import openseespy.opensees as ops

ROOT = pathlib.Path(__file__).resolve().parent.parent


def grid(nx, ny, length, height, shape=None):
    """Row-major node grid; `shape(i, j, x, y) -> (x, y)` distorts it."""
    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x, y = length * i / nx, height * j / ny
            nodes.append(list(shape(i, j, x, y)) if shape else [x, y])
    return nodes


def cells(nx, ny):
    at = lambda i, j: j * (nx + 1) + i
    return [[at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)] for j in range(ny) for i in range(nx)]


def split(quads):
    tris = []
    for k, (a, b, c, d) in enumerate(quads):
        tris += [[a, b, c], [a, c, d]] if k % 2 == 0 else [[a, b, d], [b, c, d]]
    return tris


def distort(nx, ny):
    def shape(i, j, x, y):
        if 0 < i < nx and 0 < j < ny:
            return x + 0.35 * math.sin(3 * i + 5 * j), y + 0.12 * math.cos(7 * i + 2 * j)
        return x, y
    return shape


def tip_loads(nodes, nx, ny, total):
    first = nx  # node index of the bottom tip node
    per = total / ny
    return [[first + j * (nx + 1), 1, -per * (0.5 if j in (0, ny) else 1.0)] for j in range(ny + 1)] + \
           [[first + j * (nx + 1), 0, 0.3 * per * (0.5 if j in (0, ny) else 1.0)] for j in range(ny + 1)]


def case(name, nodes, nx, ny, thickness, e, nu, loads, quads=None, triangles=None, formulation="full"):
    fixed = [j * (nx + 1) for j in range(ny + 1)]
    return dict(name=name, nodes=nodes, quads=quads or [], triangles=triangles or [], thickness=thickness,
                e=e, nu=nu, formulation=formulation, fixed=fixed, loads=loads)


def run_opensees(c, element):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in enumerate(c["nodes"], start=1):
        ops.node(tag, x, y)
    for n in c["fixed"]:
        ops.fix(n + 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, c["e"], c["nu"])
    tags = []
    for k, cell in enumerate(c["quads"] or c["triangles"], start=1):
        ops.element(element, k, *[n + 1 for n in cell], c["thickness"], "PlaneStress", 1)
        tags.append(k)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for node, dof, value in c["loads"]:
        ops.load(node + 1, *([value, 0.0] if dof == 0 else [0.0, value]))
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()
    disp = [ops.nodeDisp(n + 1, d + 1) for n in range(len(c["nodes"])) for d in (0, 1)]
    react = [ops.nodeReaction(n + 1, d + 1) for n in c["fixed"] for d in (0, 1)]
    stress = []
    for tag in tags:
        stress += list(ops.eleResponse(tag, "stress"))
    return dict(displacements=disp, reactions=react, stresses=stress)


def run_carapace(c):
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(c, f)
    out = subprocess.run(["node", str(ROOT / "comparison" / "wasm_panel.mjs"), f.name],
                         check=True, capture_output=True, text=True, cwd=ROOT)
    return json.loads(out.stdout)


def rel(a, b):
    scale = max(max(abs(v) for v in a), max(abs(v) for v in b), 1e-300)
    return max(abs(x - y) for x, y in zip(a, b)) / scale


def compare(c, element):
    o, w = run_opensees(c, element), run_carapace(c)
    if len(o["stresses"]) != len(w["stresses"]):
        return {"displacement": rel(o["displacements"], w["displacements"]),
                "reaction": rel(o["reactions"], w["reactions"]),
                "stress": f"layout differs ({len(o['stresses'])} vs {len(w['stresses'])} values)"}
    return {"displacement": rel(o["displacements"], w["displacements"]),
            "reaction": rel(o["reactions"], w["reactions"]),
            "stress": rel(o["stresses"], w["stresses"])}


def main():
    results = []
    # Cantilever strip, 6 x 2 elements: regular and distorted interior nodes.
    for label, shape in [("regular", None), ("distorted", distort(6, 2))]:
        nodes = grid(6, 2, 6.0, 1.0, shape)
        quads = cells(6, 2)
        base = dict(nodes=nodes, nx=6, ny=2, thickness=0.2, e=3000.0, nu=0.3, loads=tip_loads(nodes, 6, 2, 1.0))
        results.append((f"quad / full, {label}", "quad", compare(case("q", quads=quads, **base), "quad")))
        results.append((f"tri31 / tri3, {label}", "tri31",
                        compare(case("t", triangles=split(quads), **base), "tri31")))
        results.append((f"enhancedQuad / enhanced, {label}", "enhancedQuad",
                        compare(case("e", quads=quads, formulation="enhanced", **base), "enhancedQuad")))
    # MacNeal-Harder style: 6 elements in a row over a 0.2-deep strip; rectangles, parallelograms, trapezoids.
    for label, top, bottom in [("rectangles", lambda i: 0.0, lambda i: 0.0),
                               ("parallelograms", lambda i: 0.0 if i in (0, 6) else 0.3, lambda i: 0.0),
                               ("trapezoids", lambda i: 0.0 if i in (0, 6) else 0.2, lambda i: 0.0 if i in (0, 6) else -0.2)]:
        nodes = [[1.0 * i + bottom(i), 0.0] for i in range(7)] + [[1.0 * i + top(i), 0.2] for i in range(7)]
        quads = [[i, i + 1, 7 + i + 1, 7 + i] for i in range(6)]
        c = case("m", nodes, 6, 1, 0.1, 1e7, 0.3, [[6, 1, -0.5], [13, 1, -0.5]], quads=quads)
        for formulation, element in [("full", "quad"), ("enhanced", "enhancedQuad")]:
            c["formulation"] = formulation
            results.append((f"MacNeal-Harder {label}: {element}", element, compare(c, element)))

    print(f"{'case':44} {'u rel diff':>12} {'reaction':>12} {'Gauss stress':>14}")
    plain_ok = True
    for name, element, r in results:
        stress = r["stress"] if isinstance(r["stress"], str) else f"{r['stress']:.2e}"
        print(f"{name:44} {r['displacement']:12.2e} {r['reaction']:12.2e} {stress:>14}")
        if element in ("quad", "tri31") and not (r["displacement"] < 1e-9 and r["reaction"] < 1e-9):
            plain_ok = False
    if not plain_ok:
        print("FAIL: a plain element differs from OpenSees beyond 1e-9")
        sys.exit(1)


if __name__ == "__main__":
    main()
