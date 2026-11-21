"""The same OpenSeesPy script run through OpenSeesPy and through ``carapace.Model``; results must agree.

Skipped when ``openseespy`` is not installed (it is in ``comparison/.venv``). Each script takes the module-like
object (``openseespy.opensees`` or a ``carapace.Model``) so it is literally one script, not two ports of it.
"""

import numpy as np
import pytest

import carapace

ops_real = pytest.importorskip("openseespy.opensees")


def run_both(script):
    ops_real.wipe()
    want = script(ops_real)
    ops_real.wipe()
    got = script(carapace.Model())
    return np.asarray(got, dtype=float), np.asarray(want, dtype=float)


def analysis_setup(ops, dlambda=0.1, algorithm="Newton"):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 20)
    ops.algorithm(algorithm)
    ops.integrator("LoadControl", dlambda)
    ops.analysis("Static")


def frame2d(ops, transform="Linear", steps=10):
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for tag, xy in enumerate(((0, 0), (0, 3000), (4000, 3000), (4000, 0)), start=1):
        ops.node(tag, *map(float, xy))
    ops.fix(1, 1, 1, 1)
    ops.fix(4, 1, 1, 0)
    ops.geomTransf(transform, 1)
    for tag, (i, j) in enumerate(((1, 2), (2, 3), (4, 3)), start=1):
        ops.element("elasticBeamColumn", tag, i, j, 8000.0, 200000.0, 1e8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 2.0e5, -3.0e5, 0.0)
    ops.load(3, 0.0, -3.0e5, 1.0e8)
    analysis_setup(ops)
    out = []
    for _ in range(steps):
        assert ops.analyze(1) == 0
        out.append([ops.getTime(), *ops.nodeDisp(2), *ops.nodeDisp(3)])
    ops.reactions()
    out.append([*ops.nodeReaction(1), *ops.nodeReaction(4)])
    return np.concatenate([np.ravel(r) for r in out])


PDELTA_APPROXIMATION = pytest.mark.xfail(
    strict=True,
    reason="elasticBeamColumn PDelta is a linearized tangent only (resisting force stays elastic, axial force "
           "read from the previous step), so it differs from OpenSees; see elastic_beam_column.rs "
           "`geometric_stiffness` and docs/python-api.md",
)


@pytest.mark.parametrize("transform", ["Linear", "Corotational", pytest.param("PDelta", marks=PDELTA_APPROXIMATION)])
def test_2d_frame_stepwise(transform):
    got, want = run_both(lambda ops: frame2d(ops, transform))
    np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-8 * np.abs(want).max())


def test_truss_with_newton_and_one_big_analyze():
    def script(ops):
        ops.model("basic", "-ndm", 2, "-ndf", 2) if hasattr(ops, "ping") else ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 4000.0, 0.0)
        ops.node(3, 2000.0, 1500.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 1)
        ops.fix(3, 0, 0, 1)
        ops.uniaxialMaterial("Elastic", 1, 200000.0)
        ops.element("truss", 1, 1, 3, 500.0, 1)
        ops.element("truss", 2, 2, 3, 500.0, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 1.0e4, -5.0e4, 0.0)
        analysis_setup(ops, 0.25, "Linear")
        assert ops.analyze(4) == 0
        ops.reactions()
        return [*ops.nodeDisp(3), *ops.nodeReaction(1), *ops.nodeReaction(2)]

    got, want = run_both(script)
    np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-8 * np.abs(want).max())


def test_3d_cantilever_with_two_patterns():
    def script(ops):
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 2000.0, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.element("elasticBeamColumn", 1, 1, 2, 120.0, 200000.0, 80000.0, 3.0e5, 1.0e5, 2.0e5, 1)
        ops.timeSeries("Linear", 1)
        ops.timeSeries("Constant", 2, "-factor", 0.5)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.0, 100.0, 50.0, 0.0, 0.0, 0.0)
        ops.pattern("Plain", 2, 2, "-fact", 2.0)
        ops.load(2, 0.0, 0.0, 0.0, 1.0e3, 0.0, 2.0e3)
        analysis_setup(ops, 0.5, "Linear")
        assert ops.analyze(2) == 0
        ops.reactions()
        return [*ops.nodeDisp(2), *ops.nodeReaction(1)]

    got, want = run_both(script)
    np.testing.assert_allclose(got, want, rtol=1e-8, atol=1e-8 * np.abs(want).max())


def test_displacement_control():
    def script(ops):
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 0.0, 3000.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Corotational", 1)
        ops.element("elasticBeamColumn", 1, 1, 2, 8000.0, 200000.0, 1e8, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0e3, -2.0e4, 0.0)
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-10, 20)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 2, 1, 4.0)
        ops.analysis("Static")
        assert ops.analyze(5) == 0
        return [ops.getTime(), *ops.nodeDisp(2)]

    got, want = run_both(script)
    np.testing.assert_allclose(got, want, rtol=1e-7, atol=1e-7 * np.abs(want).max())


def test_eigenvalues_and_shapes():
    def script(ops):
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        for tag in range(1, 5):
            ops.node(tag, 0.0, 1000.0 * (tag - 1))
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        for tag in range(1, 4):
            ops.element("elasticBeamColumn", tag, tag, tag + 1, 8000.0, 200000.0, 1e8, 1, "-mass", 0.05)
        ops.mass(4, 2.0, 2.0, 0.0)
        values = ops.eigen(3)
        vec = [ops.nodeEigenvector(4, 1, d) for d in (1, 2, 3)]
        return [*values, *vec]

    got, want = run_both(script)
    np.testing.assert_allclose(got[:3], want[:3], rtol=1e-6)
    # Mode shapes are only defined up to sign and scale: compare directions.
    g, w = got[3:], want[3:]
    assert abs(np.dot(g, w)) / (np.linalg.norm(g) * np.linalg.norm(w)) == pytest.approx(1.0, abs=1e-6)
