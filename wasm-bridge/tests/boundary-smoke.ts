// Run with Node 22.18+ after wasm-pack build; also type-check against generated .d.ts.
import { decodeInput, axial_displacement } from "../../pkg/carapace_wasm.js";
import type { CarapaceInputV1, DecodeError, StepOutcome } from "../../pkg/carapace_wasm.js";

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
console.log("Wasm boundary smoke passed: planar/spatial sessions, batches, errors, legacy export.");
