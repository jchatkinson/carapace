// Run with Node 22.18+ after wasm-pack build; also type-check against generated .d.ts.
import { decodeInput, axial_displacement } from "../../pkg/carapace_wasm.js";
import type { AlgorithmSpec, CarapaceInputV1, DecodeError, StepOutcome } from "../../pkg/carapace_wasm.js";

function assert(condition: boolean, message: string): asserts condition {
  if (!condition) throw new Error(message);
}

function emptyInput(space: number): CarapaceInputV1 {
  const nodes = { coords: [], fixed: [], massNodeIndex: [], mass: [] };
  const trusses = { nodeI: [], nodeJ: [], area: [], material: [], density: [] };
  const fibers = { sectionOffsets: [], y: [], area: [], material: [] };
  const fiberBeams = { nodeI: [], nodeJ: [], fiberSection: [], integration: [], corotational: [], density: [] };
  const zeroLengths = { nodeI: [], nodeJ: [], materials: [], friction: [] };
  const zeroLengthSections = { nodeI: [], nodeJ: [], fiberSection: [], materials: [] };
  const equalDofs = { retained: [], constrained: [], dofs: [] };
  const nodalLoads = { pattern: [], node: [], dof: [], value: [], stage: [] };
  const elementLoads = { pattern: [], elementKind: [], elementIndex: [], load: [], stage: [] };
  const spatialFiberBeams = { nodeI: [], nodeJ: [], g: [], j: [], vecXz: [], fiberSection: [], integration: [], density: [] };
  return {
    header: { schemaVersion: 1, space, engineVersion: "smoke" },
    nodes, materials: [], fibers, trusses,
    elasticBeamColumns: { nodeI: [], nodeJ: [], e: [], a: [], iz: [], transform: [], density: [] },
    dispBeamColumns: fiberBeams, forceBeamColumns: fiberBeams,
    zeroLengths, zeroLengthSections, equalDofs,
    rigidDiaphragms: { retained: [], constrained: [] },
    loadPatterns: { series: [], scaleFactor: [] }, nodalLoads, elementLoads,
    sequence: { stages: [], recorders: [] },
    nodes3: nodes, fibers3: { ...fibers, z: [] }, trusses3: trusses,
    elasticBeamColumns3: { nodeI: [], nodeJ: [], e: [], g: [], a: [], j: [], iy: [], iz: [], transform: [], density: [] },
    dispBeamColumns3: spatialFiberBeams, forceBeamColumns3: spatialFiberBeams,
    zeroLengths3: zeroLengths, zeroLengthSections3: zeroLengthSections, equalDofs3: equalDofs,
    rigidDiaphragms3: { retained: [], normal: [], constrained: [] },
    nodalLoads3: nodalLoads, elementLoads3: elementLoads,
    sequence3: { stages: [], recorders: [] },
  };
}

// Exercise both profile decoders through the generated JS glue and typed API.
for (const space of [2, 3]) {
  const input = emptyInput(space);
  const nodes = {
    coords: space === 2 ? [0, 0, 100, 0] : [0, 0, 0, 100, 0, 0],
    fixed: space === 2 ? [0b111, 0b110] : [0b111111, 0b111110],
    massNodeIndex: [], mass: [],
  };
  const trusses = { nodeI: [0], nodeJ: [1], area: [2], material: [0], density: [0] };
  const loads = { pattern: [0], node: [1], dof: [0], value: [50], stage: [0] };
  const stages: CarapaceInputV1["sequence"]["stages"] = [{
    kind: "static", id: "load", steps: 2,
    integrator: { kind: "loadControl", increment: 0.5 }, algorithm: "linear",
    convergence: undefined, holdPatternsAfter: [],
  }];
  if (space === 2) {
    input.nodes = nodes; input.trusses = trusses; input.nodalLoads = loads;
    input.sequence = { stages, recorders: [{ response: "nodeDisp", node: 1, dof: 0 }] };
  } else {
    input.nodes3 = nodes; input.trusses3 = trusses; input.nodalLoads3 = loads;
    input.sequence3 = { stages, recorders: [{ response: "nodeDisp", node: 1, dof: 0 }] };
  }
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
assert(error.kind === "unsupportedSpace" && error.got === 4, "structured decode error");
const malformed = expectThrow(() => {
  // @ts-expect-error Missing fields must be rejected by TypeScript and at runtime.
  decodeInput({ header: {} });
});
assert(String(malformed).includes("malformed CarapaceInputV1"), "malformed input error");
assert(Math.abs(axial_displacement(50, 100, 2, 30000) - 1 / 12) < 1e-12, "legacy export");

function yieldingInput(space: number, transient: boolean, algorithm?: AlgorithmSpec): CarapaceInputV1 {
  const input = emptyInput(space);
  const nodes = {
    coords: space === 2 ? [0, 0, 0, 0] : [0, 0, 0, 0, 0, 0],
    fixed: space === 2 ? [0b111, 0b110] : [0b111111, 0b111110],
    massNodeIndex: [1], mass: space === 2 ? [1, 0, 0] : [1, 0, 0, 0, 0, 0],
  };
  const springs = { nodeI: [0, 0], nodeJ: [1, 1], materials: [[0, 0, 0], [1, 0, 1]] as [number, number, number][], friction: [] };
  const loads = { pattern: [0], node: [1], dof: [0], value: [transient ? 1000 : 5], stage: [0] };
  input.materials = [{ kind: "elastic", e: 50 }, { kind: "elasticPp", e: 100, eyp: 0.01 }];
  input.loadPatterns = { series: [transient ? { kind: "linear", slope: 10 } : { kind: "constant" }], scaleFactor: [1] };
  const convergence = { kind: "normUnbalance" as const, tol: transient ? 1e-6 : 1e-9, maxIter: 120 };
  const sequence: CarapaceInputV1["sequence"] = {
    stages: [transient ? {
      kind: "transient", id: "yield", steps: 1, dt: 0.1,
      damping: { alphaM: 0, betaK: 0 }, groundMotions: [],
      // Omit both fields in the legacy case to check backward compatibility.
      ...(algorithm === undefined ? {} : { algorithm, convergence }),
    } : {
      kind: "static", id: "yield", steps: 1,
      integrator: { kind: "loadControl", increment: 1 },
      algorithm: algorithm ?? "linear", convergence, holdPatternsAfter: [],
    }], recorders: [{ response: "nodeDisp", node: 1, dof: 0 }],
  };
  if (space === 2) {
    input.nodes = nodes; input.zeroLengths = springs; input.nodalLoads = loads; input.sequence = sequence;
  } else {
    input.nodes3 = nodes; input.zeroLengths3 = springs; input.nodalLoads3 = loads; input.sequence3 = sequence;
  }
  return input;
}

const algorithms: AlgorithmSpec[] = ["newtonRaphson", { kind: "newton" }];
for (const tangent of ["current", "reuseAtStepStart", "initial"] as const) {
  algorithms.push({ kind: "newton", tangent }, { kind: "krylovNewton", tangent, maxDimension: 3 });
  for (const kind of ["bisection", "regulaFalsi"] as const) {
    algorithms.push({ kind: "newton", tangent, lineSearch: { kind, tol: 1e-10, maxIter: 30, maxEta: 16 } });
  }
}
for (const space of [2, 3]) {
  for (const transient of [false, true]) {
    for (const algorithm of [undefined, "linear", { kind: "linear" }, ...algorithms] as const) {
      const session = decodeInput(yieldingInput(space, transient, algorithm));
      try {
        const result = session.advance(1);
        assert(result.done && result.error === undefined, `solver failed: ${JSON.stringify(algorithm)}`);
        const linear = algorithm === undefined || algorithm === "linear" || (typeof algorithm === "object" && algorithm.kind === "linear");
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
    const input = yieldingInput(2, true, "newtonRaphson");
    const stage = input.sequence.stages[0];
    assert(stage.kind === "transient", "validation stage");
    stage.convergence = { kind: "normUnbalance", tol: 1e-6, maxIter: 100, [field]: value };
    const error = expectThrow(() => decodeInput(input)) as DecodeError;
    assert(error.kind === "invalidAnalysisOption" && error.field === `convergence.${field}`, "invalid convergence setting");
  }
}
for (const dt of [0, -0.1, NaN, Infinity]) {
  const input = yieldingInput(2, true);
  const stage = input.sequence.stages[0];
  assert(stage.kind === "transient", "validation stage");
  stage.dt = dt;
  const error = expectThrow(() => decodeInput(input)) as DecodeError;
  assert(error.kind === "invalidAnalysisOption" && error.field === "dt", "invalid timestep");
}
const failureInput = yieldingInput(2, true, "newtonRaphson");
const stage = failureInput.sequence.stages[0];
assert(stage.kind === "transient", "transient test stage");
stage.convergence = { kind: "normUnbalance", tol: 1e-6, maxIter: 1 };
const failureSession = decodeInput(failureInput);
try {
  const result = failureSession.advance(1);
  assert(result.done && result.error?.kind === "failedToConverge", "transient convergence setting");
  assert(failureSession.advance(1).error?.kind === "failedToConverge", "sticky failure");
} finally { failureSession.free(); }

console.log("Wasm boundary smoke passed: both profiles, configurable static/transient solvers, legacy inputs, errors.");
