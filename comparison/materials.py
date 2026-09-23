"""openseespy uniaxialMaterial definitions — must mirror
core/examples/cyclic_material.rs's `material()` fn exactly (same params,
same units) for the comparison plots to mean anything."""

import openseespy.opensees as ops

MATERIALS = (
    "elastic",
    "elastic_pp",
    "steel01",
    "concrete01",
    "steel02",
    "concrete02",
    "hysteretic",
    "pinching4",
    "parallel",
    "series",
    "min_max",
)

# Must match core/examples/cyclic_material.rs's `initial_tangent()` result
# for each material (its `Material::initial_tangent()`), and the same
# SAFETY_FRACTION — see that file's `material_with_safety` doc comment for
# why a small parallel safety stiffness is needed at all (a bare spring
# with zero tangent, e.g. elastic_pp fully yielded or concrete01 in
# tension, is a singular system under DisplacementControl).
SAFETY_FRACTION = 0.005
INITIAL_TANGENT = {
    "elastic": 200.0,
    "elastic_pp": 200.0,
    "steel01": 200.0,
    "concrete01": 4000.0,  # 2*fpc/epsc0 = 2*4.0/0.002
    "steel02": 200.0,
    "concrete02": 4000.0,  # 2*fc/epsc0 = 2*4.0/0.002
    "hysteretic": 300.0,  # mom1p/rot1p = 1.5/0.005
    "pinching4": 200.0,  # stress1p/strain1p = 1.0/0.005
    "parallel": 200.0,  # 150 (steel01 e0) + 50 (elastic)
    "series": 75.0,  # 1/(1/300 + 1/100)
    "min_max": 200.0,  # inner steel01's e0 — bounds don't affect initial tangent
}


def _define_primary(name, tag):
    if name == "elastic":
        ops.uniaxialMaterial("Elastic", tag, 200.0)
    elif name == "elastic_pp":
        ops.uniaxialMaterial("ElasticPP", tag, 200.0, 0.01)
    elif name == "steel01":
        ops.uniaxialMaterial("Steel01", tag, 2.0, 200.0, 0.02, 0.0, 1.0, 0.0, 1.0)
    elif name == "concrete01":
        # OpenSees' Concrete01 takes compression-negative values directly
        # (carapace's Material::concrete01 takes positive magnitudes and
        # normalizes internally to the same negative convention).
        ops.uniaxialMaterial("Concrete01", tag, -4.0, -0.002, -3.4, -0.006)
    elif name == "steel02":
        ops.uniaxialMaterial("Steel02", tag, 2.0, 200.0, 0.02, 18.0, 0.925, 0.15, 0.0, 1.0, 0.0, 1.0)
    elif name == "concrete02":
        ops.uniaxialMaterial("Concrete02", tag, -4.0, -0.002, -3.4, -0.006, 0.1, 0.4, 200.0)
    elif name == "hysteretic":
        ops.uniaxialMaterial(
            "Hysteretic", tag,
            1.5, 0.005, 2.0, 0.02, 2.2, 0.04,
            -1.5, -0.005, -2.0, -0.02, -2.2, -0.04,
            0.5, 0.3, 0.0, 0.0, 0.0,
        )
    elif name == "pinching4":
        ops.uniaxialMaterial(
            "Pinching4", tag,
            1.0, 0.005, 1.8, 0.015, 2.0, 0.025, 1.0, 0.04,
            -1.0, -0.005, -1.8, -0.015, -2.0, -0.025, -1.0, -0.04,
            0.5, 0.25, 0.05, 0.5, 0.25, 0.05,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            10.0, "energy",
        )
    elif name == "parallel":
        ops.uniaxialMaterial("Steel01", tag * 10 + 1, 1.0, 150.0, 0.02, 0.0, 1.0, 0.0, 1.0)
        ops.uniaxialMaterial("Elastic", tag * 10 + 2, 50.0)
        ops.uniaxialMaterial("Parallel", tag, tag * 10 + 1, tag * 10 + 2)
    elif name == "series":
        ops.uniaxialMaterial("Elastic", tag * 10 + 1, 300.0)
        ops.uniaxialMaterial("Steel01", tag * 10 + 2, 1.0, 100.0, 0.02, 0.0, 1.0, 0.0, 1.0)
        ops.uniaxialMaterial("Series", tag, tag * 10 + 1, tag * 10 + 2)
    elif name == "min_max":
        ops.uniaxialMaterial("Steel01", tag * 10 + 1, 2.0, 200.0, 0.02, 0.0, 1.0, 0.0, 1.0)
        ops.uniaxialMaterial("MinMax", tag, tag * 10 + 1, "-min", -0.025, "-max", 0.025)
    else:
        raise ValueError(f"unknown material {name!r}")


def define(name, tag):
    """Defines uniaxialMaterial `tag` as `name` in parallel with a small
    safety stiffness, in the current openseespy model. Returns safety_k,
    the safety spring's stiffness, so the caller can subtract its linear
    contribution back out of the recorded force."""
    primary_tag, safety_tag = tag * 10 + 1, tag * 10 + 2
    _define_primary(name, primary_tag)
    safety_k = SAFETY_FRACTION * INITIAL_TANGENT[name]
    ops.uniaxialMaterial("Elastic", safety_tag, safety_k)
    ops.uniaxialMaterial("Parallel", tag, primary_tag, safety_tag)
    return safety_k
