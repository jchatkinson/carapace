// Wasm counterpart of core/examples/benchmark_frame.rs (same frame, same JSON shape).
// Timing covers the advance(1) loop only; decodeInput is reported as setup_ms.
// Run after: wasm-pack build wasm-bridge --target nodejs --out-dir ../pkg
import { createRequire } from "node:module";
import { performance } from "node:perf_hooks";
import { parseArgs } from "node:util";

const { decodeInput } = createRequire(import.meta.url)("../pkg/carapace_wasm.js");

const { values: opts } = parseArgs({ options: {
  stories: { type: "string", default: "10" }, bays: { type: "string", default: "3" },
  steps: { type: "string", default: "50" }, model: { type: "string", default: "elastic" },
  "displacement-increment": { type: "string" }, "warm-runs": { type: "string", default: "0" },
} });
const stories = +opts.stories, bays = +opts.bays, steps = +opts.steps;
const model = opts.model, dispInc = opts["displacement-increment"] === undefined ? undefined : +opts["displacement-increment"];

function buildInput() {
  const E = 200000, colA = 8000, colI = 2e8, beamA = 6000, beamI = 3e8;
  const coords = [], fixed = [], idx = (s, b) => s * (bays + 1) + b;
  for (let s = 0; s <= stories; s++) for (let b = 0; b <= bays; b++) {
    coords.push(b * 6000, s * 3000); fixed.push(s === 0 ? 0b111 : 0);
  }
  const columns = [], beams = [];
  for (let s = 0; s < stories; s++) for (let b = 0; b <= bays; b++) columns.push([idx(s, b), idx(s + 1, b)]);
  for (let s = 1; s <= stories; s++) for (let b = 0; b < bays; b++) beams.push([idx(s, b), idx(s, b + 1)]);
  const all = [...columns, ...beams];
  const empty = () => ({ nodeI: [], nodeJ: [], density: [] });
  const input = {
    header: { schemaVersion: 1, ndm: 2, engineVersion: "bench" },
    nodes: { coords, fixed, massNodeIndex: [], mass: [] },
    materials: [], fibers: { sectionOffsets: [], y: [], area: [], material: [] },
    elasticBeamColumns2d: { ...empty(), e: [], a: [], iz: [], transform: [] },
    dispBeamColumns2d: { ...empty(), fiberSection: [], integration: [], corotational: [] },
    loadPatterns: { series: [{ kind: "linear", slope: 1 }], scaleFactor: [1] },
    nodalLoads: { pattern: [], node: [], dof: [], value: [], stage: [] },
  };
  if (model === "fiber") {
    input.materials = [{ kind: "steel01", fy: 355, e0: E, b: 0.02, a1: 0, a2: 1, a3: 0, a4: 1 }];
    const sections = [[colA, colI], [beamA, beamI]];
    for (const [area, inertia] of sections) {
      const h = Math.sqrt(inertia / area);
      input.fibers.sectionOffsets.push(input.fibers.y.length);
      input.fibers.y.push(h, -h); input.fibers.area.push(area / 2, area / 2); input.fibers.material.push(0, 0);
    }
    input.fibers.sectionOffsets.push(input.fibers.y.length);
    all.forEach(([i, j], k) => {
      const t = input.dispBeamColumns2d;
      t.nodeI.push(i); t.nodeJ.push(j); t.density.push(0);
      t.fiberSection.push(k < columns.length ? 0 : 1);
      t.integration.push({ kind: "legendre", points: 3 }); t.corotational.push(false);
    });
  } else {
    all.forEach(([i, j], k) => {
      const t = input.elasticBeamColumns2d, col = k < columns.length;
      t.nodeI.push(i); t.nodeJ.push(j); t.density.push(0); t.e.push(E);
      t.a.push(col ? colA : beamA); t.iz.push(col ? colI : beamI); t.transform.push("linear");
    });
  }
  const pBase = model === "fiber" ? 25000 : 10000;
  for (let s = 1; s <= stories; s++) {
    const l = input.nodalLoads;
    l.pattern.push(0); l.node.push(idx(s, 0)); l.dof.push(0); l.value.push(pBase * s); l.stage.push(0);
  }
  const roof = idx(stories, 0);
  input.sequence = {
    recorders: [{ response: "nodeDisp", node: roof, dof: 0 }],
    stages: [{
      kind: "static", id: "bench", steps,
      integrator: dispInc === undefined ? { kind: "loadControl", increment: 1 / steps }
        : { kind: "displacementControl", node: roof, dof: 0, increment: dispInc },
      algorithm: { kind: "newton", tangent: "current" },
      convergence: { kind: "normDispIncr", tol: 1e-6, maxIter: 30 }, holdPatternsAfter: [],
    }],
  };
  return { input, elements: all.length };
}

function run() {
  const t0 = performance.now();
  const { input, elements } = buildInput();
  const session = decodeInput(input);
  const setupMs = performance.now() - t0;
  const roof = [], load = [];
  const t1 = performance.now();
  for (let i = 0; i < steps; i++) {
    const out = session.advance(1);
    if (out.error) throw new Error("analysis failed: " + JSON.stringify(out.error));
    load.push(out.loadFactor);
    for (const batch of out.recorderBatches) for (const [, v] of batch.samples) roof.push(v);
  }
  const solveMs = performance.now() - t1;
  session.free();
  return { elements, setupMs, solveMs, roof, load };
}

for (let i = 0; i < +opts["warm-runs"]; i++) run();
const r = run();
console.log(JSON.stringify({
  engine: "wasm", model, stories, bays, nodes: (stories + 1) * (bays + 1), elements: r.elements,
  free_dofs: stories * (bays + 1) * 3, steps, setup_ms: r.setupMs, solve_ms: r.solveMs,
  total_ms: r.setupMs + r.solveMs, final_roof_disp: r.roof.at(-1),
  roof_disp_history: r.roof, load_factor_history: r.load,
  iterations: null, // StepOutcome does not expose Newton iteration counts
  warm_runs: +opts["warm-runs"],
}));
