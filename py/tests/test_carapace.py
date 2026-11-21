import json
import pathlib

import numpy as np
import pytest

import carapace

EXAMPLE = pathlib.Path(__file__).resolve().parents[2] / "examples" / "models" / "portal_frame_2d.json"


def portal() -> carapace.Model:
    m = carapace.Model(ndm=2)
    for tag, xy in enumerate(((0, 0), (0, 3000), (4000, 3000), (4000, 0)), start=1):
        m.node(tag, *map(float, xy))
    m.fix(1, 1, 1, 1)
    m.fix(4, 1, 1, 1)
    m.geomTransf("Linear", 1)
    for tag, (i, j) in enumerate(((1, 2), (2, 3), (4, 3)), start=1):
        m.element("elasticBeamColumn", tag, i, j, 8000.0, 200000.0, 1e8, 1)
    m.timeSeries("Linear", 1)
    m.pattern("Plain", 1, 1)
    m.load(2, 10000.0, 0.0, 0.0)
    m.integrator("LoadControl", 1.0)
    m.analysis("Static")
    return m


def drift(results: carapace.Results) -> float:
    """Node 2 (table index 1), dof ux, at the first step: recorder 0 is every node's displacement."""
    return results.bulk(0)[1][0, 1, 0, 0]


def test_builder_matches_the_bundled_json_example():
    expected = json.loads(EXAMPLE.read_text())
    got = drift(carapace.run(portal()))
    ref = carapace.run(expected).recorder(0)[0, 1]
    assert got == pytest.approx(ref, rel=1e-12)
    assert 0 < got < 100


def test_validate_reports_a_path():
    bad = portal().to_dict()
    bad["header"]["ndm"] = "two"
    out = carapace.validate(bad)
    assert out["valid"] is False and out["errors"][0]["path"] == "header.ndm"
    with pytest.raises(carapace.ModelError):
        carapace.run(bad)


def test_decode_error_is_structured():
    m = portal()
    m.raw["nodalLoads"] = {"pattern": [0], "node": [99], "dof": [0], "value": [1.0], "stage": [0]}
    out = carapace.validate(m)
    assert out["valid"] is False and out["errors"][0]["stage"] == "decode"


def test_3d_cantilever_matches_closed_form():
    m = carapace.Model(ndm=3)
    e, g, a, iy, iz, j, length, load = 200000.0, 80000.0, 100.0, 1e5, 1e5, 2e5, 1000.0, 100.0
    m.node(1, 0.0, 0.0, 0.0)
    m.node(2, length, 0.0, 0.0)
    m.fix(1, 1, 1, 1, 1, 1, 1)
    m.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    m.element("elasticBeamColumn", 1, 1, 2, a, e, g, j, iy, iz, 1)
    m.timeSeries("Linear", 1)
    m.pattern("Plain", 1, 1)
    m.load(2, 0.0, load, 0.0, 0.0, 0.0, 0.0)
    m.integrator("LoadControl", 1.0)
    m.analysis("Static")
    assert m.analyze(1) == 0
    assert m.nodeDisp(2, 2) == pytest.approx(load * length**3 / (3 * e * iz), rel=1e-9)


def test_modal():
    m = carapace.Model(ndm=2)
    m.node(1, 0.0, 0.0)
    m.node(2, 0.0, 3000.0, "-mass", 1.0, 1.0, 0.0)
    m.fix(1, 1, 1, 1)
    m.geomTransf("Linear", 1)
    m.element("elasticBeamColumn", 1, 1, 2, 8000.0, 200000.0, 1e8, 1)
    (lam,) = m.eigen(1)
    assert np.isfinite(lam) and lam > 0


def test_run_many_and_threads_do_not_change_results():
    models = []
    for k in range(1, 7):
        m = portal()
        m._loads["value"][0] *= k
        models.append(m)
    models.append({"header": {"schemaVersion": 1, "ndm": "x", "engineVersion": "t"}})
    out = carapace.run_many(models, threads=3)
    assert out[-1]["valid"] is False
    singles = [drift(carapace.run(m)) for m in models[:-1]]
    assert [drift(r) for r in out[:-1]] == singles
    assert singles[3] == pytest.approx(4 * singles[0], rel=1e-9)
    assert drift(carapace.run(models[0], threads=4)) == singles[0]


def test_unsupported_commands_raise_instead_of_being_ignored():
    m = carapace.Model()
    with pytest.raises(NotImplementedError, match="equalDOF"):
        m.equalDOF(1, 2, 1)
    with pytest.raises(NotImplementedError, match="Steel02"):
        m.uniaxialMaterial("Steel02", 1, 1.0)
    with pytest.raises(NotImplementedError, match="analysis 'Transient'"):
        m.analysis("Transient")
    with pytest.raises(AttributeError):
        m.not_a_command


def test_model_is_fixed_once_analysis_starts():
    m = portal()
    assert m.analyze(1) == 0
    with pytest.raises(RuntimeError, match="after analyze"):
        m.node(9, 0.0, 0.0)
    m.wipe()
    m.node(9, 0.0, 0.0)
