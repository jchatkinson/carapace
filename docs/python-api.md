# Python API (`carapace-fe`)

The Python commands follow **OpenSeesPy**: same names, same argument order, 1-based DOFs, free tags.
A script written for OpenSeesPy should run by changing the import (or by creating a `Model`). Where Carapace
differs, it says so: an unsupported command or option raises `NotImplementedError` naming it (it is never
ignored), and every deliberate deviation is listed under [Deviations](#deviations) below.
`py/tests/test_opensees_compat.py` runs the same script through OpenSeesPy and through Carapace and compares
the results, so the "same" claim is tested (it needs `openseespy`; CI installs it).

```python
import carapace.opensees as ops          # module-level, like openseespy.opensees

ops.model('basic', '-ndm', 2, '-ndf', 3)
ops.node(1, 0.0, 0.0); ops.node(2, 0.0, 3000.0)
ops.fix(1, 1, 1, 1)
ops.geomTransf('Linear', 1)
ops.element('elasticBeamColumn', 1, 1, 2, 8000.0, 200000.0, 1e8, 1)
ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1); ops.load(2, 1e4, 0.0, 0.0)
ops.constraints('Plain'); ops.numberer('RCM'); ops.system('BandGeneral')
ops.test('NormDispIncr', 1e-8, 10); ops.algorithm('Newton')
ops.integrator('LoadControl', 0.5); ops.analysis('Static')
ops.analyze(2)
ops.nodeDisp(2, 1)
```

or, with no hidden global state (needed for several models, `run_many` and threads):

```python
import carapace
m = carapace.Model(ndm=2, ndf=3)         # same methods as `ops`, e.g. m.node(1, 0.0, 0.0)
```

## Supported commands

| Command | Supported forms |
|---|---|
| `model` | `model('basic', '-ndm', 2, '-ndf', 3)` and `-ndm 3 -ndf 6` |
| `wipe` | clears the model and analysis |
| `node` | `node(tag, x, y[, z])`, `-mass m1 ... m_ndf` |
| `mass`, `fix` | all `ndf` values / flags |
| `uniaxialMaterial` | `'Elastic', tag, E` |
| `geomTransf` | 2D `'Linear'`, `'PDelta'`, `'Corotational'`; 3D `'Linear'`, `'PDelta'` with `vecxz` |
| `element` | `'truss'` (`-rho`), `'elasticBeamColumn'` with `A, E, Iz` (2D) or `A, E, G, J, Iy, Iz` (3D) and `-mass` |
| `timeSeries` | `'Linear'`, `'Constant'`, `-factor` |
| `pattern`, `load` | `'Plain'`, `-fact`; nodal loads |
| `constraints`, `numberer`, `system` | accepted (see deviations) |
| `test` | `NormUnbalance`, `NormDispIncr`, `EnergyIncr` |
| `algorithm` | `Linear`, `Newton`, `ModifiedNewton`, `-initial` |
| `integrator` | `LoadControl dLambda`, `DisplacementControl node dof dU` |
| `analysis`, `analyze` | `'Static'`, `analyze(n)` |
| `nodeDisp`, `nodeReaction`, `reactions`, `getTime` | as OpenSeesPy |
| `eigen`, `nodeEigenvector` | `eigen(n)` returns eigenvalues (omega squared) |
| `nodeCoord`, `getNodeTags`, `getEleTags` | as OpenSeesPy |

Anything else (e.g. `equalDOF`, `section`, `nDMaterial`, `recorder`, `eleResponse`, `loadConst`, other elements,
materials, integrators, `analysis('Transient')`) raises `NotImplementedError`. The wire format supports much more
than the Python layer exposes yet; for those features put the table into `model.raw` (a dict merged over the
generated document) or pass a `CarapaceInputV1` dict/JSON to `carapace.run`.

## Deviations

| Where | OpenSeesPy | Carapace | Why |
|---|---|---|---|
| State | One global domain | `Model` objects; `carapace.opensees` is one module-level instance | `run_many` and threads need independent models |
| After `analyze` | Model and analysis can keep changing | `node`, `element`, `load`, ... and `test`/`algorithm`/`integrator` raise once `analyze` has run (`wipe()` starts over) | A stage's definition is fixed when it starts; staged analysis (`loadConst`, a second pattern and integrator) is not exposed yet |
| `analyze(n)` | Returns `0`, `-1`, `-2`, `-3`... | `0` or `-1`; the error is in `model.last_error` | The engine reports a tagged error, not OpenSees' codes |
| `analyze(n, dt)` | Transient | `NotImplementedError` | Static only for now |
| `constraints`, `numberer`, `system` | Select algorithms | Accepted, no effect (`constraints` takes `Plain` or `Transformation`; names are checked) | One equation numbering and one sparse direct solver (Cholesky, LU fallback) |
| Missing analysis commands | `analyze` fails without `algorithm`/`test` | Defaults to `Linear` and `NormUnbalance 1e-8, 20`; `integrator` and `analysis` are required | Convenience; set them explicitly for a nonlinear run |
| `-mass`, `-rho` | Mass per unit length | Same argument; stored as mass per volume (`mass / A`) | The wire format's `density` is volumetric |
| `nodeReaction` | Needs `reactions()` first | Reactions are recorded each step; `reactions()` is a no-op. `Model(reactions=False)` turns the recording off | Reactions come from a recorder; recording them costs one pass over the model per step |
| `eigen` | Uses the current state | Uses the undeformed model; solver flags are ignored | Modal analysis is its own stage |
| `nodeEigenvector` | OpenSees' normalization | Carapace's normalization (compare directions, not values) | |
| `getTime` | Pseudo-time | The load factor reached | |
| `geomTransf('PDelta')` on `elasticBeamColumn` | Chord-rotation (P-Delta) terms only, in force and tangent | The consistent geometric stiffness `(K + Kg(N)) d`, so it includes P-small-delta and matches the exact beam-column solution with one element. **Results differ from OpenSees by about 1%** and converge to the same answer as the mesh is refined | A deliberate choice: more accurate, not bit-comparable. `Corotational` (2D) matches OpenSees exactly |
| Tags | Any integers | Any integers; mapped to table indices internally | The wire format is columnar and 0-based |
| `load` | Zero values are loads | Zero components are skipped | The wire format rejects a load on a DOF nothing stiffens |
| Results to arrays | `ops.nodeDisp` etc. only | Also `carapace.run(model, steps=n)` returning `Results` with recorder 0 = every node displacement per step, recorder 1 = every reaction (`Results.bulk`) | Batch runs and numpy |

## Running many models

`carapace.run(model, steps=n)` and `carapace.run_many([...], threads=k, steps=n)` take `Model` objects (or wire
dicts). They run the static stage as configured for `n` steps and return `Results`; this is how you sweep load
cases in parallel. `carapace.validate(model)` and `ModelError` report problems with the JSON path of the bad value.
