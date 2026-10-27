// Run with Node 22.18+ after scripts/build-wasm.sh (tests the shipped build); also type-check against generated .d.ts.
import { loadCarapace } from "../../comparison/load_wasm.mjs";
import type {
  AlgorithmSpec, BulkBatch, CarapaceInputV1, ContinuationDetail, DecodeError, IntegratorSpec, LimitSample, ModelOutput, SequenceSpec, StepOutcome,
} from "../../pkg/carapace_wasm.js";

const { decodeInput } = await loadCarapace();

function assert(condition: boolean, message: string): asserts condition {
  if (!condition) throw new Error(message);
}

// Every table is optional on the wire, so an empty model is just a header.
function emptyInput(ndm: number): CarapaceInputV1 {
  return { header: { schemaVersion: 1, ndm, engineVersion: "smoke" } };
}

// Exercise both profiles through the generated JS glue and typed API.
for (const ndm of [2, 3]) {
  const input = emptyInput(ndm);
  const nodes = {
    coords: ndm === 2 ? [0, 0, 100, 0] : [0, 0, 0, 100, 0, 0],
    fixed: ndm === 2 ? [0b111, 0b110] : [0b111111, 0b111110],
    massNodeIndex: [], mass: [],
  };
  const trusses = { nodeI: [0], nodeJ: [1], area: [2], material: [0], density: [0] };
  const loads = { pattern: [0], node: [1], dof: [0], value: [50], stage: [0] };
  const stages: SequenceSpec["stages"] = [{
    kind: "static", id: "load", steps: 2,
    integrator: { kind: "loadControl", increment: 0.5 }, algorithm: { kind: "linear" },
    convergence: undefined, holdPatternsAfter: [],
  }];
  input.nodes = nodes; input.trusses = trusses; input.nodalLoads = loads;
  input.sequence = { stages, recorders: [{ response: "nodeDisp", node: 1, dof: 0 }] };
  input.materials = [{ kind: "elastic", e: 30000 }];
  input.loadPatterns = { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] };
  const session = decodeInput(input);
  try {
    assert(session.currentStageId() === "load", "current stage");
    const paused: StepOutcome = session.advance(0);
    assert(paused.stepsTaken === 0 && paused.recorderBatches.length === 0, "zero budget");
    for (let step = 0; step < 2; step++) {
      const result: StepOutcome = session.advance(1);
      assert(result.error === undefined, "successful analysis");
      assert(result.done === (step === 1) && result.stageComplete === (step === 1), "completion");
      const batch = result.recorderBatches[0];
      assert(batch.firstSample === step && batch.samples.length === 1, "incremental samples");
      assert(batch.recorderIndex === 0 && batch.stageIndex === 0, "batch indices");
      const [time, displacement] = batch.samples[0];
      assert(time === (step + 1) / 2, "sample time");
      assert(Math.abs(displacement - (50 * 100 / (2 * 30000)) * time) < 1e-12, "truss displacement");
    }
    assert(session.currentStageId() === undefined, "completed stage");
    assert(session.advance(1).recorderBatches.length === 0, "no repeated samples");
  } finally {
    session.free();
  }
}

function expectThrow(call: () => unknown): unknown {
  try { call(); } catch (error) { return error; }
  throw new Error("expected an exception");
}

const error = expectThrow(() => decodeInput(emptyInput(4))) as DecodeError;
assert(error.kind === "unsupportedNdm" && error.got === 4, "structured decode error");
const malformed = expectThrow(() => {
  // @ts-expect-error Missing fields must be rejected by TypeScript and at runtime.
  decodeInput({ header: {} });
});
assert(String(malformed).includes("malformed CarapaceInputV1"), "malformed input error");

function yieldingInput(ndm: number, transient: boolean, algorithm?: AlgorithmSpec): CarapaceInputV1 {
  const input = emptyInput(ndm);
  const nodes = {
    coords: ndm === 2 ? [0, 0, 0, 0] : [0, 0, 0, 0, 0, 0],
    fixed: ndm === 2 ? [0b111, 0b110] : [0b111111, 0b111110],
    massNodeIndex: [1], mass: ndm === 2 ? [1, 0, 0] : [1, 0, 0, 0, 0, 0],
  };
  const springs = { nodeI: [0, 0], nodeJ: [1, 1], materials: [[0, 0, 0], [1, 0, 1]] as [number, number, number][] };
  const loads = { pattern: [0], node: [1], dof: [0], value: [transient ? 1000 : 5], stage: [0] };
  input.materials = [{ kind: "elastic", e: 50 }, { kind: "elasticPp", e: 100, eyp: 0.01 }];
  input.loadPatterns = { series: [transient ? { kind: "linear", slope: 10 } : { kind: "constant" }], scaleFactor: [1] };
  const convergence = { kind: "normUnbalance" as const, tol: transient ? 1e-6 : 1e-9, maxIter: 120 };
  const sequence: SequenceSpec = {
    stages: [transient ? {
      kind: "transient", id: "yield", steps: 1, dt: 0.1,
      damping: { alphaM: 0, betaK: 0 }, groundMotions: [],
      // Omit both fields to check the default linear transient solver.
      ...(algorithm === undefined ? {} : { algorithm, convergence }),
    } : {
      kind: "static", id: "yield", steps: 1,
      integrator: { kind: "loadControl", increment: 1 },
      algorithm: algorithm ?? { kind: "linear" }, convergence, holdPatternsAfter: [],
    }], recorders: [{ response: "nodeDisp", node: 1, dof: 0 }],
  };
  input.nodes = nodes; input.zeroLengths = springs; input.nodalLoads = loads; input.sequence = sequence;
  return input;
}

const algorithms: AlgorithmSpec[] = [{ kind: "newton" }];
for (const tangent of ["current", "reuseAtStepStart", "initial"] as const) {
  algorithms.push({ kind: "newton", tangent }, { kind: "krylovNewton", tangent, maxDimension: 3 });
  for (const kind of ["bisection", "regulaFalsi"] as const) {
    algorithms.push({ kind: "newton", tangent, lineSearch: { kind, tol: 1e-10, maxIter: 30, maxEta: 16 } });
  }
}
for (const ndm of [2, 3]) {
  for (const transient of [false, true]) {
    for (const algorithm of [undefined, { kind: "linear" }, ...algorithms] as const) {
      const session = decodeInput(yieldingInput(ndm, transient, algorithm));
      try {
        const result = session.advance(1);
        assert(result.done && result.error === undefined, `solver failed: ${JSON.stringify(algorithm)}`);
        const linear = algorithm === undefined || algorithm.kind === "linear";
        const expected = transient ? (linear ? 1000 / 550 : 999 / 450) : (linear ? 5 / 150 : 4 / 50);
        assert(Math.abs(result.recorderBatches[0].samples[0][1] - expected) < 1e-8, "nonlinear equilibrium");
      } finally { session.free(); }
    }
  }
}

const invalidAlgorithms: [AlgorithmSpec, string][] = [
  [{ kind: "krylovNewton", tangent: "initial", maxDimension: 0 }, "algorithm.maxDimension"],
  ...["tol", "maxIter", "maxEta"].map((field): [AlgorithmSpec, string] => [{
    kind: "newton", lineSearch: { kind: "bisection", tol: 1e-10, maxIter: 30, maxEta: 16, [field]: 0 },
  }, `algorithm.lineSearch.${field}`]),
];
for (const [algorithm, field] of invalidAlgorithms) {
  for (const transient of [false, true]) {
    const error = expectThrow(() => decodeInput(yieldingInput(2, transient, algorithm))) as DecodeError;
    assert(error.kind === "invalidAnalysisOption" && error.stage === "yield" && error.field === field, "invalid solver option");
  }
}
for (const field of ["tol", "maxIter"] as const) {
  for (const value of field === "tol" ? [0, -1, NaN, Infinity] : [0]) {
    const input = yieldingInput(2, true, { kind: "newton" });
    const stage = input.sequence!.stages[0];
    assert(stage.kind === "transient", "validation stage");
    stage.convergence = { kind: "normUnbalance", tol: 1e-6, maxIter: 100, [field]: value };
    const error = expectThrow(() => decodeInput(input)) as DecodeError;
    assert(error.kind === "invalidAnalysisOption" && error.field === `convergence.${field}`, "invalid convergence setting");
  }
}
for (const dt of [0, -0.1, NaN, Infinity]) {
  const input = yieldingInput(2, true);
  const stage = input.sequence!.stages[0];
  assert(stage.kind === "transient", "validation stage");
  stage.dt = dt;
  const error = expectThrow(() => decodeInput(input)) as DecodeError;
  assert(error.kind === "invalidAnalysisOption" && error.field === "dt", "invalid timestep");
}
const failureInput = yieldingInput(2, true, { kind: "newton" });
const stage = failureInput.sequence!.stages[0];
assert(stage.kind === "transient", "transient test stage");
stage.convergence = { kind: "normUnbalance", tol: 1e-6, maxIter: 1 };
const failureSession = decodeInput(failureInput);
try {
  const result = failureSession.advance(1);
  assert(result.done && result.error?.kind === "failedToConverge", "transient convergence setting");
  assert(failureSession.advance(1).error?.kind === "failedToConverge", "sticky failure");
} finally { failureSession.free(); }

// Arc-length continuation past a degrading-strength peak in 2D and 3D:
// load factors fall after the peak while sample indices keep increasing,
// the stop criterion ends the stage before its step cap, and continuation
// diagnostics cross the boundary.
function softeningInput(ndm: number, integrator: IntegratorSpec): CarapaceInputV1 {
  const input = emptyInput(ndm);
  const nodes = {
    coords: ndm === 2 ? [0, 0, 0, 0] : [0, 0, 0, 0, 0, 0],
    fixed: ndm === 2 ? [0b111, 0b110] : [0b111111, 0b111110],
    massNodeIndex: [], mass: [],
  };
  const springs = { nodeI: [0], nodeJ: [1], materials: [[0, 0, 0]] as [number, number, number][] };
  const loads = { pattern: [0], node: [1], dof: [0], value: [1], stage: [0] };
  input.materials = [{
    kind: "hysteretic", mom1p: 10, rot1p: 0.01, mom2p: 6, rot2p: 0.02, mom3p: 2, rot3p: 0.03,
    mom1n: -10, rot1n: -0.01, mom2n: -6, rot2n: -0.02, mom3n: -2, rot3n: -0.03,
    pinchX: 1, pinchY: 1, damfc1: 0, damfc2: 0, beta: 0,
  }];
  input.loadPatterns = { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] };
  const sequence: SequenceSpec = {
    stages: [{
      kind: "static", id: "push", steps: 200, integrator,
      algorithm: { kind: "newton" }, convergence: undefined, holdPatternsAfter: [],
    }],
    recorders: [{ response: "nodeDisp", node: 1, dof: 0 }],
  };
  input.nodes = nodes; input.zeroLengths = springs; input.nodalLoads = loads; input.sequence = sequence;
  return input;
}

const arcLength: IntegratorSpec = {
  kind: "arcLength", initialRadius: 0.05,
  scales: { kind: "explicit", displacement: 0.01, load: 10 },
  arcTolerance: 1e-9,
  stop: { displacement: { node: 1, dof: 0, value: 0.02, exact: true } },
};
for (const ndm of [2, 3]) {
  const session = decodeInput(softeningInput(ndm, arcLength));
  try {
    const result: StepOutcome = session.advance(1000);
    assert(result.done && result.stageComplete && result.error === undefined, "arc length completes");
    assert(result.stepsTaken < 200, "stop criterion ends the stage early");
    const continuation: ContinuationDetail | undefined = result.continuation;
    assert(continuation !== undefined && continuation.stop?.reason === "displacementTarget", "stop reported");
    assert(continuation.stop.landedExactly && continuation.arcError <= 1e-9, "exact landing");
    const batch = result.recorderBatches[0];
    assert(batch.firstSample === 0 && batch.samples.length === result.stepsTaken, "one sample per step");
    const lambdas = batch.samples.map(([lambda]) => lambda);
    const peak = lambdas.indexOf(Math.max(...lambdas));
    assert(peak > 0 && peak < lambdas.length - 5, "interior peak");
    assert(lambdas.slice(peak).every((l, i, tail) => i === 0 || l < tail[i - 1]), "load falls after the peak");
    const [lambda, u] = batch.samples[batch.samples.length - 1];
    assert(Math.abs(u - 0.02) < 1e-9 && Math.abs(lambda - 6) < 1e-6, "lands on the backbone target");
  } finally { session.free(); }
}

const arcErrors: [IntegratorSpec, AlgorithmSpec, string][] = [
  [{ ...arcLength, minRadius: 0.1 }, { kind: "newton" }, "integrator.initialRadius"],
  [arcLength, { kind: "linear" }, "algorithm"],
  [arcLength, { kind: "newton", tangent: "initial" }, "algorithm"],
];
for (const [integrator, algorithm, field] of arcErrors) {
  const input = softeningInput(2, integrator);
  const stage = input.sequence!.stages[0];
  assert(stage.kind === "static", "static stage");
  stage.algorithm = algorithm;
  const error = expectThrow(() => decodeInput(input)) as DecodeError;
  assert(error.kind === "invalidAnalysisOption" && error.field === field, `arc length rejects ${field}`);
}
const typo = expectThrow(() => decodeInput(softeningInput(2, {
  ...arcLength,
  // @ts-expect-error Unknown arc-length fields must be rejected by TypeScript and at runtime.
  initialRadus: 0.05,
})));
assert(String(typo).includes("malformed CarapaceInputV1"), "unknown arc-length field");

// Unified tables: a profile mismatch is a structured error; rigid links and general constraints decode.
const wrongProfile = emptyInput(2);
wrongProfile.elasticBeamColumns3d = { nodeI: [0], nodeJ: [], e: [], g: [], a: [], j: [], iy: [], iz: [], transform: [], density: [] };
const profileError = expectThrow(() => decodeInput(wrongProfile)) as DecodeError;
assert(profileError.kind === "tableNotInProfile" && profileError.table === "elastic_beam_columns_3d", "table not in profile");
for (const ndm of [2, 3]) {
  const input = emptyInput(ndm);
  input.nodes = {
    coords: ndm === 2 ? [0, 0, 1, 0, 2, 0] : [0, 0, 0, 1, 0, 0, 2, 0, 0],
    fixed: ndm === 2 ? [0b111, 0, 0] : [0b111111, 0, 0], massNodeIndex: [], mass: [],
  };
  input.rigidLinks = { master: [0], slave: [1] };
  input.linearConstraints = { slaveNode: [2], slaveDof: [0], termOffsets: [0, 1], termNode: [1], termDof: [0], termCoeff: [1] };
  const session = decodeInput(input);
  session.free();
}
const badComponent = emptyInput(2);
badComponent.sequence = { stages: [], recorders: [{ response: "elementLoad", elementKind: "truss", elementIndex: 0, component: 16 }] };
badComponent.nodes = { coords: [0, 0, 1, 0], fixed: [7, 6], massNodeIndex: [], mass: [] };
badComponent.materials = [{ kind: "elastic", e: 1 }];
badComponent.trusses = { nodeI: [0], nodeJ: [1], area: [1], material: [0], density: [0] };
const componentError = expectThrow(() => decodeInput(badComponent)) as DecodeError;
assert(componentError.kind === "invalidRecorderComponent" && componentError.width === 16, "recorder component bound");

// Continuum elements: a two-quad patch under an edge traction, Gauss-point stress crossing the boundary.
const continuum = emptyInput(2);
continuum.nodes = { coords: [0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1], fixed: [3, 0, 0, 3, 0, 0], massNodeIndex: [], mass: [] };
continuum.planeMaterials = [{ kind: "isotropic", e: 1000, nu: 0, state: "planeStress" }];
continuum.quads = { nodeIds: [0, 1, 4, 3, 1, 2, 5, 4], thickness: [1, 1], material: [0, 0], density: [0, 0], formulation: ["full", "enhanced"] };
continuum.loadPatterns = { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] };
continuum.elementLoads = { pattern: [0], elementKind: ["quad4"], elementIndex: [1], load: [{ kind: "edgeTraction", edge: 1, tx: 1, ty: 0 }], stage: [0] };
continuum.sequence = {
  stages: [{ kind: "static", id: "pull", steps: 1, integrator: { kind: "loadControl", increment: 1 }, algorithm: { kind: "linear" }, holdPatternsAfter: [] }],
  recorders: [{ response: "gaussPoint", elementKind: "quad4", elementIndex: 0, point: 2, quantity: "stress", component: 0 }],
};
{
  const session = decodeInput(continuum);
  try {
    const outcome = session.advance(1);
    assert(outcome.error === undefined, "continuum run");
    // Uniform tension of 1 per unit length over unit height: sigma_x = 1 everywhere.
    assert(Math.abs(outcome.recorderBatches[0].samples[0][1] - 1) < 1e-9, "quad stress");
  } finally { session.free(); }
}
const gaussError = expectThrow(() => decodeInput({
  ...continuum, sequence: { ...continuum.sequence!, recorders: [{ response: "gaussPoint", elementKind: "quad4", elementIndex: 0, point: 4, quantity: "strain", component: 0 }] },
})) as DecodeError;
assert(gaussError.kind === "invalidGaussPoint" && gaussError.count === 4, "gauss point bound");
const planeError = expectThrow(() => decodeInput({ ...continuum, planeMaterials: [{ kind: "isotropic", e: -1, nu: 0, state: "planeStress" }] })) as DecodeError;
assert(planeError.kind === "invalidPlaneMaterial", "plane material validation");

// Shells: a clamped unit quad under pressure; the four vertical reactions sum to -p A.
const shell = emptyInput(3);
shell.nodes = { coords: [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0], fixed: [63, 63, 63, 63], massNodeIndex: [], mass: [] };
shell.shellSections = [{ kind: "elasticMembranePlate", e: 1000, nu: 0.3, h: 0.1, rho: 0 }];
shell.shell4s = { nodeIds: [0, 1, 2, 3], section: [0] };
shell.loadPatterns = { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] };
shell.elementLoads = { pattern: [0], elementKind: ["shell4"], elementIndex: [0], load: [{ kind: "shellPressure", pressure: 2 }], stage: [0] };
shell.sequence = {
  stages: [{ kind: "static", id: "press", steps: 1, integrator: { kind: "loadControl", increment: 1 }, algorithm: { kind: "linear" }, holdPatternsAfter: [] }],
  recorders: [0, 1, 2, 3].map((node) => ({ response: "reaction" as const, node, dof: 2 })),
};
{
  const session = decodeInput(shell);
  try {
    const outcome = session.advance(1);
    assert(outcome.error === undefined, "shell run");
    const total = outcome.recorderBatches.reduce((sum, batch) => sum + batch.samples[0][1], 0);
    assert(Math.abs(total + 2) < 1e-9, "shell pressure reactions");
  } finally { session.free(); }
}
const shellError = expectThrow(() => decodeInput({ ...shell, shellSections: [{ kind: "elasticMembranePlate", e: 1000, nu: 0.3, h: 0, rho: 0 }] })) as DecodeError;
assert(shellError.kind === "invalidShellSection", "shell section validation");

// CSFM: one cracking Q4 (csfmConcrete) with a bonded csfmRebar bar (embeddedBars), pulled in tension (mm, N, MPa).
// The bar spans the quad at mid-height, so its strain is the mean of the right-edge displacements over the length.
const csfm = emptyInput(2);
csfm.nodes = { coords: [0, 0, 200, 0, 200, 100, 0, 100], fixed: [3, 0, 0, 3], massNodeIndex: [], mass: [] };
// A tensile regularization is needed for a single cracking element to converge (the unregularized softening diverges).
csfm.planeMaterials = [{ kind: "csfmConcrete", fck: 30, law: "parabolaRectangle", softening: true, tensileRegularization: 0.01 }];
csfm.quads = { nodeIds: [0, 1, 2, 3], thickness: [100], material: [0], density: [0], formulation: ["full"] };
csfm.materials = [{
  kind: "csfmRebar", es: 200000, fy: 500, ft: 540, epsU: 0.05, fctm: 2.9, ecm: 32800,
  diameter: 16, rhoEff: 0.02, lambda: 1, tensionStiffening: true,
}];
csfm.embeddedBars = [{ polyline: [0, 50, 200, 50], area: 400, material: 0 }];
csfm.loadPatterns = { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] };
csfm.nodalLoads = { pattern: [0, 0], node: [1, 2], dof: [0, 0], value: [20000, 20000], stage: [0, 0] };
csfm.sequence = {
  stages: [{
    kind: "static", id: "pull", steps: 4, integrator: { kind: "loadControl", increment: 0.25 },
    algorithm: { kind: "newton" }, convergence: { kind: "normUnbalance", tol: 1e-3, maxIter: 50 }, holdPatternsAfter: [],
  }],
  recorders: [
    { response: "bulk", table: "embeddedBars", quantity: "axialForce" },
    { response: "embeddedBar", bar: 0, segment: 0, quantity: "stress" },
    { response: "gaussPoint", elementKind: "quad4", elementIndex: 0, point: 0, quantity: "principalStress", component: 0 },
    { response: "gaussPoint", elementKind: "quad4", elementIndex: 0, point: 0, quantity: "principalStress", component: 1 },
    { response: "nodeDisp", node: 1, dof: 0 },
    { response: "nodeDisp", node: 2, dof: 0 },
    { response: "embeddedBar", bar: 0, segment: 0, quantity: "strain" },
    { response: "bulk", table: "nodes", quantity: "displacement" },
    { response: "bulk", table: "nodes", quantity: "reaction" },
  ],
};
{
  const session = decodeInput(csfm);
  try {
    const model: ModelOutput = session.modelOutput();
    assert(model.embeddedBars.length === 1 && model.embeddedBars[0].firstSegment === 0, "model output groups");
    const group = model.embeddedBars[0];
    assert(group.segments.length === 1 && Math.abs(group.segments[0].length - 200) < 1e-9, "model output segment");
    assert(group.segments[0].hostKind === "quad4", "model output host");
    assert(group.routing === "tcm" || group.routing === "pom", "model output routing");
    assert(group.sr !== undefined && group.sr >= 0, "model output crack spacing");

    const outcome: StepOutcome = session.advance(100);
    assert(outcome.done && outcome.error === undefined, `csfm run: ${JSON.stringify(outcome.error)}`);
    assert(Math.abs(outcome.loadFactor - 1) < 1e-12, "csfm load factor");
    const series = (index: number) => outcome.recorderBatches.find((batch) => batch.recorderIndex === index)!.samples.map(([, v]) => v);
    const [stress, p1, p2, ux1, ux2, strain] = [1, 2, 3, 4, 5, 6].map(series);
    assert(stress.length === 4 && strain.length === 4, "bar recorder samples");

    // Bulk recorder: one frame of [axialForce] per accepted step, one segment, one point.
    const bulk: BulkBatch = outcome.bulkBatches.find((batch) => batch.recorderIndex === 0)!;
    assert(bulk.elements === 1 && bulk.points === 1 && bulk.components === 1 && bulk.progress.length === 4, "bulk batch shape");
    assert(bulk.values.length === 4, "bulk batch values");
    for (let k = 0; k < 4; k++) {
      assert(Math.abs(bulk.values[k] - stress[k] * 400) <= 1e-9 * Math.abs(bulk.values[k]) + 1e-9, "bulk bar force = stress x area");
      assert(Math.abs(strain[k] - (ux1[k] + ux2[k]) / 2 / 200) <= 1e-9 * Math.abs(strain[k]), "bar strain follows the host displacements");
      assert(p1[k] >= p2[k], "principal stresses are ordered");
    }
    assert(bulk.values[3] > bulk.values[0] && bulk.values[3] > 0, "bar picks up force as the load grows");
    assert(p1[3] > 0, "the host is in tension");

    // Bulk output arrives as typed arrays on their own buffers, so a worker can transfer them.
    assert(bulk.values instanceof Float64Array && bulk.progress instanceof Float64Array, "bulk arrays are Float64Array");
    // Nodes table: every node in input order, one point, ndf components, per step.
    const nodeDisp: BulkBatch = outcome.bulkBatches.find((batch) => batch.recorderIndex === 7)!;
    const nodeReaction: BulkBatch = outcome.bulkBatches.find((batch) => batch.recorderIndex === 8)!;
    for (const table of [nodeDisp, nodeReaction]) {
      assert(table.elements === 4 && table.points === 1 && table.components === 3 && table.progress.length === 4, "nodes bulk shape");
      assert(table.values.length === 4 * 4 * 3, "nodes bulk values");
    }
    const frame = (table: BulkBatch, step: number, node: number, dof: number) => table.values[(step * 4 + node) * 3 + dof];
    for (let k = 0; k < 4; k++) {
      assert(frame(nodeDisp, k, 1, 0) === ux1[k] && frame(nodeDisp, k, 2, 0) === ux2[k], "nodes displacement matches the scalar recorders");
      assert(frame(nodeDisp, k, 0, 0) === 0 && Math.abs(frame(nodeReaction, k, 1, 0)) < 1e-2, "supported displacement is zero and a free DOF carries only the convergence residual");
    }
    const rx = frame(nodeReaction, 3, 0, 0) + frame(nodeReaction, 3, 3, 0);
    assert(Math.abs(rx + 40000) < 1e-3, `x reactions balance the applied load: ${rx}`);
    const detached = bulk.values.buffer as ArrayBuffer;
    structuredClone(detached, { transfer: [detached] });
    assert(detached.byteLength === 0, "bulk buffers are transferable");

    // Per-step solver work, aligned with the bulk frames.
    assert(outcome.stepDiagnostics.length === 4, "one step diagnostic per step");
    assert(outcome.stepDiagnostics.every((d, k) => d.iterations >= 1 && d.factorizations >= 1 && d.cutback === undefined && Math.abs(d.progress - bulk.progress[k]) < 1e-12), "step diagnostics");

    // Controlling material-limit report: one entry per accepted step, located by table row.
    const limits: LimitSample[] = outcome.limitHistory;
    assert(limits.length === 4 && limits.every((l) => l.utilization >= 0 && l.utilization < 1), "limit history");
    assert(limits.every((l) => l.elementKind === "quad4" || l.elementKind === "embeddedBar"), "limit location kinds");
  } finally { session.free(); }
}
// Bad CSFM constants are structured errors naming the arena row; a bar outside the mesh names its segment.
const badRebar = expectThrow(() => decodeInput({
  ...csfm, materials: [{ kind: "csfmRebar", es: 200000, fy: 500, ft: 540, epsU: 0.05, fctm: 2.9, ecm: 32800, diameter: 16, rhoEff: 0, lambda: 1, tensionStiffening: true }],
})) as DecodeError;
assert(badRebar.kind === "invalidMaterial" && badRebar.index === 0, "csfmRebar validation");
const outsideBar = expectThrow(() => decodeInput({ ...csfm, embeddedBars: [{ polyline: [0, 500, 200, 500], area: 400, material: 0 }] })) as DecodeError;
assert(outsideBar.kind === "invalidModel" && outsideBar.error.kind === "embeddedBarOutsideMesh", "bar outside the mesh");

console.log("Wasm boundary smoke passed: 2D and 3D, unified tables, continuum elements, shells, CSFM concrete/rebar/embedded bars, configurable static/transient solvers, arc length, defaults, errors.");
