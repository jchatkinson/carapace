// Replays the CSFM span-beam golden model through the built wasm (pkg/) and compares every increment with the native run.
// Usage: node comparison/wasm_span_beam.mjs [budget]   (advance budget per call, default 1000; exits 1 on a mismatch)
// Run after: scripts/build-wasm.sh, then
//   cargo test -p carapace-wasm --test wire export_span_beam -- --ignored
// which writes comparison/out/span_beam_{wire,native}.json from the same conversion as wasm-bridge/tests/wire/span_beam.rs.
import { readFileSync } from "node:fs";
import { loadCarapace } from "./load_wasm.mjs";

const { decodeInput } = await loadCarapace();
const read = (name) => {
  try { return JSON.parse(readFileSync(new URL(`./out/${name}`, import.meta.url), "utf8")); }
  catch { throw new Error(`comparison/out/${name} is missing: run the export test named in this file's header`); }
};
const input = read("span_beam_wire.json"), native = read("span_beam_native.json");
const budget = Number(process.argv[2] ?? 1000);
// Same engine and arithmetic, so only rounding is allowed: the wasm build may fuse or reorder (simd128).
const REL = 1e-9, ABS = 1e-12;

const bulk = Array.from({ length: 7 }, () => []), disp = [], limits = [];
let loadFactor;
const session = decodeInput(input);
try {
  for (;;) {
    const outcome = session.advance(budget);
    if (outcome.error) throw new Error(JSON.stringify(outcome.error));
    for (const b of outcome.bulkBatches) {
      const width = b.elements * b.points * b.components;
      for (let k = 0; k < b.progress.length; k++) bulk[b.recorderIndex].push(b.values.slice(k * width, (k + 1) * width));
    }
    for (const b of outcome.recorderBatches) (disp[b.recorderIndex - 7] ??= []).push(...b.samples.map(([, v]) => v));
    limits.push(...outcome.limitHistory);
    if (outcome.done) { loadFactor = outcome.loadFactor; break; }
  }
} finally { session.free(); }

const fields = ["strain", "stress", "kc2", "barForce", "barStrain", "barStress", "barRuptured"];
let worst = 0, failures = 0, ties = 0;
function compare(label, a, b) {
  if (a.length !== b.length) { console.error(`${label}: length ${a.length} vs ${b.length}`); failures++; return; }
  for (let i = 0; i < a.length; i++) {
    const err = Math.abs(a[i] - b[i]) / Math.max(Math.abs(b[i]), ABS / REL);
    worst = Math.max(worst, err);
    if (!(Math.abs(a[i] - b[i]) <= ABS + REL * Math.abs(b[i]))) { console.error(`${label}[${i}]: wasm ${a[i]} vs native ${b[i]}`); if (++failures > 20) process.exit(1); }
  }
}
if (bulk[0].length !== native.length) throw new Error(`increments: wasm ${bulk[0].length} vs native ${native.length}`);
native.forEach((frame, k) => {
  compare(`increment ${k} displacements`, disp.map((series) => series[k]), frame.displacements);
  fields.forEach((name, r) => compare(`increment ${k} ${name}`, bulk[r][k], frame[name]));
  const w = limits[k], n = frame.limit;
  // The span beam is mirror-symmetric, so two elements can tie for the controlling utilization to ~1e-15 and rounding picks
  // the location. The kind and the utilization must match; the element may differ only when the two utilizations tie.
  if (!w !== !n || (n && w.kind !== n.kind)) { console.error(`increment ${k}: controlling limit ${JSON.stringify(w)} vs ${JSON.stringify(n)}`); failures++; }
  else if (n) {
    compare(`increment ${k} utilization`, [w.utilization], [n.utilization]);
    if (w.elementKind !== n.elementKind || w.elementIndex !== n.elementIndex) { ties++; console.log(`increment ${k}: tied controlling element ${w.elementKind} ${w.elementIndex} (wasm) vs ${n.elementKind} ${n.elementIndex} (native)`); }
  }
});
if (Math.abs(loadFactor - 1) > 1e-12) { console.error(`final load factor ${loadFactor}`); failures++; }
console.log(`span beam: ${native.length} increments, advance budget ${budget}, worst relative difference wasm vs native ${worst.toExponential(2)}, ${ties} tied-location swaps`);
if (failures) { console.error(`${failures} mismatches`); process.exit(1); }
