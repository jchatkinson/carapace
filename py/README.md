# carapace-fe

Python interface to the Carapace finite-element engine (PyO3 + maturin). The commands follow OpenSeesPy;
see [`docs/python-api.md`](../docs/python-api.md) for what is supported and every deviation, and
[`docs/native-api.md`](../docs/native-api.md) for the CLI and threads.

```bash
python -m venv .venv && .venv/bin/pip install maturin numpy pytest openseespy
.venv/bin/maturin develop --release
.venv/bin/python -m pytest tests          # the OpenSeesPy comparison tests need openseespy
```

```python
import carapace.opensees as ops
ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 3000.0); ops.fix(1, 1, 1, 1)
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, 8000.0, 200000.0, 1e8, 1)
ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1); ops.load(2, 1e4, 0.0, 0.0)
ops.integrator('LoadControl', 1.0); ops.analysis('Static')
ops.analyze(1)
print(ops.nodeDisp(2, 1))
```
