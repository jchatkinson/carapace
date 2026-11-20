import json
import pathlib

import numpy as np
import pytest

import carapace

EXAMPLE = pathlib.Path(__file__).resolve().parents[2] / "examples" / "models" / "portal_frame_2d.json"


def portal() -> carapace.Model:
    m = carapace.Model(ndm=2)
    a = m.add_node(0, 0, fixed="all")
    b = m.add_node(0, 3000)
    c = m.add_node(4000, 3000)
    d = m.add_node(4000, 0, fixed="all")
    for i, j in ((a, b), (b, c), (d, c)):
        m.add_elastic_beam(i, j, e=200000, a=8000, iz=1e8)
    m.add_static_stage()
    m.add_nodal_load(b, 0, 10000)
    m.record_node_disp(b, 0)
    return m


def test_builder_matches_the_bundled_json_example():
    expected = json.loads(EXAMPLE.read_text())
    got = carapace.run(portal()).recorder(0)[0, 1]
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
    m.add_nodal_load(99, 0, 1.0)
    out = carapace.validate(m)
    assert out["valid"] is False and out["errors"][0]["stage"] == "decode"


def test_3d_cantilever_matches_closed_form():
    m = carapace.Model(ndm=3)
    base = m.add_node(0, 0, 0, fixed="all")
    tip = m.add_node(1000, 0, 0)
    e, g, a, iy, iz, j, length, load = 200000.0, 80000.0, 100.0, 1e5, 1e5, 2e5, 1000.0, 100.0
    m.add_elastic_beam(base, tip, e=e, g=g, a=a, iy=iy, iz=iz, j_torsion=j)
    m.add_static_stage()
    m.add_nodal_load(tip, 1, load)
    r = m.record_node_disp(tip, 1)
    assert carapace.run(m).recorder(r)[0, 1] == pytest.approx(load * length**3 / (3 * e * iz), rel=1e-9)


def test_modal_stage():
    m = carapace.Model(ndm=2)
    base = m.add_node(0, 0, fixed="all")
    tip = m.add_node(0, 3000)
    m.add_elastic_beam(base, tip, e=200000, a=8000, iz=1e8, density=7.85e-9)
    m.add_modal_stage(modes=1)
    freq = carapace.run(m).modes()[0]["frequency"]
    assert np.isfinite(freq) and freq > 0
