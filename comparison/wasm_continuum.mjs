// Times an N x N Quad4 panel (clamped left edge, shear on the right edge) through the wasm build.
// Usage: node comparison/wasm_continuum.mjs [N] [reps] [full|enhanced]
// Run after: wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg
import { createRequire } from "node:module";

const { decodeInput } = createRequire(import.meta.url)("../pkg/carapace_wasm.js");
const n = Number(process.argv[2] ?? 50), reps = Number(process.argv[3] ?? 3), formulation = process.argv[4] ?? "full";
const id = (i, j) => j * (n + 1) + i;
const coords = [], fixed = [], quadNodes = [], loadNodes = [];
for (let j = 0; j <= n; j++) for (let i = 0; i <= n; i++) {
  coords.push(i, j);
  fixed.push(i === 0 ? 0b011 : 0);
  if (i === n) loadNodes.push(id(i, j));
}
for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) quadNodes.push(id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
const count = n * n;
const input = {
  header: { schemaVersion: 1, ndm: 2, engineVersion: "bench" },
  nodes: { coords, fixed, massNodeIndex: [], mass: [] },
  planeMaterials: [{ kind: "isotropic", e: 200000, nu: 0.3, state: "planeStress" }],
  quads: { nodeIds: quadNodes, thickness: Array(count).fill(1), material: Array(count).fill(0), density: Array(count).fill(0), formulation: Array(count).fill(formulation) },
  loadPatterns: { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] },
  nodalLoads: { pattern: loadNodes.map(() => 0), node: loadNodes, dof: loadNodes.map(() => 1), value: loadNodes.map(() => -1), stage: loadNodes.map(() => 0) },
  sequence: {
    stages: [{ kind: "static", id: "panel", steps: 1, integrator: { kind: "loadControl", increment: 1 }, algorithm: { kind: "newton", tangent: "current" }, holdPatternsAfter: [] }],
    recorders: [{ response: "nodeDisp", node: id(n, n), dof: 1 }],
  },
};
const decode = [], run = [];
let tip;
for (let r = 0; r < reps; r++) {
  let t = performance.now();
  const session = decodeInput(input);
  decode.push(performance.now() - t);
  t = performance.now();
  const outcome = session.advance(10);
  run.push(performance.now() - t);
  if (outcome.error) throw new Error(JSON.stringify(outcome.error));
  tip = outcome.recorderBatches[0].samples.at(-1)[1];
  session.free();
}
const median = (v) => [...v].sort((a, b) => a - b)[Math.floor(v.length / 2)];
console.log(`${n}x${n} Quad4 (${formulation}), ${2 * n * (n + 1)} DOF, tip uy ${tip.toExponential(6)}`);
console.log(`wasm decode ${median(decode).toFixed(1)} ms, advance ${median(run).toFixed(1)} ms (medians of ${reps})`);
