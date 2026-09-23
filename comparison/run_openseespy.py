"""Runs one material through the same reverse-cyclic DisplacementControl
protocol as core/examples/cyclic_material.rs, on the openseespy side, and
writes the recorded (strain, stress) trace to CSV.

Usage: run_openseespy.py <material> <out.csv>
"""

import csv
import sys

import openseespy.opensees as ops

import materials
from protocol import load_protocol, cyclic_targets

NODE_FIXED = 1
NODE_FREE = 2
MAT_TAG = 1
ELE_TAG = 1
PATTERN_TAG = 1
DOF = 1  # openseespy dofs are 1-indexed: 1 = x translation


def build_model(material_name):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(NODE_FIXED, 0.0, 0.0)
    ops.node(NODE_FREE, 0.0, 0.0)
    ops.fix(NODE_FIXED, 1, 1, 1)
    ops.fix(NODE_FREE, 0, 1, 1)  # free only in x translation

    safety_k = materials.define(material_name, MAT_TAG)
    ops.element("zeroLength", ELE_TAG, NODE_FIXED, NODE_FREE, "-mat", MAT_TAG, "-dir", DOF)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", PATTERN_TAG, 1)
    ops.load(NODE_FREE, 1.0, 0.0, 0.0)  # unit reference load: load factor reads as spring force

    ops.system("BandGeneral")
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.test("NormUnbalance", 1e-10, 30)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", NODE_FREE, DOF, 0.0)
    ops.analysis("Static")  # built once; only the integrator's increment changes per leg

    return safety_k


def leg(target):
    current = ops.nodeDisp(NODE_FREE, DOF)
    ops.integrator("DisplacementControl", NODE_FREE, DOF, target - current)
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError(f"analyze failed stepping toward {target}")
    return ops.nodeDisp(NODE_FREE, DOF), ops.getLoadFactor(PATTERN_TAG)


def ramp_to(target, step, points):
    while True:
        current = ops.nodeDisp(NODE_FREE, DOF)
        remaining = target - current
        if abs(remaining) < 1e-12:
            return
        next_target = current + max(-step, min(step, remaining))
        points.append(leg(next_target))


def main():
    material_name, out_path = sys.argv[1], sys.argv[2]
    tiers, step = load_protocol()

    safety_k = build_model(material_name)

    points = [(0.0, 0.0)]
    for target in cyclic_targets(tiers):
        ramp_to(target, step, points)

    # Subtract the safety spring's linear contribution back out, matching
    # core/examples/cyclic_material.rs.
    points = [(strain, total_stress - safety_k * strain) for strain, total_stress in points]

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["strain", "stress"])
        writer.writerows(points)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
