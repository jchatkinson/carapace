// Runs one continuum panel case on Carapace (wasm) and prints displacements, support
// reactions and Gauss-point stresses as JSON. The case is the JSON file given as argv[2]:
// { nodes: [[x, y]], quads: [[n0..n3]], triangles: [[n0..n2]], thickness, e, nu,
//   formulation: "full" | "enhanced", fixed: [node], loads: [[node, dof, value]] }
// Run after: wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg
import { createRequire } from "node:module";
import { readFileSync } from "node:fs";

const { decodeInput } = createRequire(import.meta.url)("../pkg/carapace_wasm.js");
const c = JSON.parse(readFileSync(process.argv[2], "utf8"));

const quads = c.quads ?? [], triangles = c.triangles ?? [];
const fixed = new Set(c.fixed);
const recorders = [];
c.nodes.forEach((_, node) => [0, 1].forEach((dof) => recorders.push({ response: "nodeDisp", node, dof })));
c.fixed.forEach((node) => [0, 1].forEach((dof) => recorders.push({ response: "reaction", node, dof })));
const points = [...quads.map(() => 4), ...triangles.map(() => 1)];
const kinds = [...quads.map(() => "quad4"), ...triangles.map(() => "tri3")];
const indices = [...quads.map((_, i) => i), ...triangles.map((_, i) => i)];
points.forEach((count, e) => {
  for (let point = 0; point < count; point++) for (let component = 0; component < 3; component++) {
    recorders.push({ response: "gaussPoint", elementKind: kinds[e], elementIndex: indices[e], point, quantity: "stress", component });
  }
});

const input = {
  header: { schemaVersion: 1, ndm: 2, engineVersion: "panel" },
  nodes: { coords: c.nodes.flat(), fixed: c.nodes.map((_, n) => (fixed.has(n) ? 0b011 : 0)), massNodeIndex: [], mass: [] },
  planeMaterials: [{ kind: "isotropic", e: c.e, nu: c.nu, state: "planeStress" }],
  quads: { nodeIds: quads.flat(), thickness: quads.map(() => c.thickness), material: quads.map(() => 0), density: quads.map(() => 0), formulation: quads.map(() => c.formulation) },
  triangles: { nodeIds: triangles.flat(), thickness: triangles.map(() => c.thickness), material: triangles.map(() => 0), density: triangles.map(() => 0) },
  loadPatterns: { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] },
  nodalLoads: { pattern: c.loads.map(() => 0), node: c.loads.map((l) => l[0]), dof: c.loads.map((l) => l[1]), value: c.loads.map((l) => l[2]), stage: c.loads.map(() => 0) },
  sequence: {
    stages: [{ kind: "static", id: "panel", steps: 1, integrator: { kind: "loadControl", increment: 1 }, algorithm: "linear", holdPatternsAfter: [] }],
    recorders,
  },
};
const session = decodeInput(input);
const outcome = session.advance(10);
if (outcome.error) throw new Error(JSON.stringify(outcome.error));
const value = new Array(recorders.length).fill(null);
for (const batch of outcome.recorderBatches) value[batch.recorderIndex] = batch.samples.at(-1)[1];
session.free();

const nNodes = c.nodes.length, nFixed = c.fixed.length;
console.log(JSON.stringify({
  displacements: value.slice(0, 2 * nNodes),
  reactions: value.slice(2 * nNodes, 2 * nNodes + 2 * nFixed),
  stresses: value.slice(2 * nNodes + 2 * nFixed),
}));
