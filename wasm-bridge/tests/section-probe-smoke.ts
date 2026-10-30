// Run with Node 22.18+ after scripts/build-wasm.sh (tests the shipped build); also type-check against generated .d.ts.
// Drives the section probe through the generated JS glue with the CalcsApp-generated fixtures
// (wasm-bridge/tests/fixtures/section_probe/README.md): typed-array transfer, Ts<T> in and out, error objects, timings.
import { loadCarapace } from "../../comparison/load_wasm.mjs";
import sectionFixture from "./fixtures/section_probe/section.json" with { type: "json" };
import statesFixture from "./fixtures/section_probe/states.json" with { type: "json" };
import mphiFixture from "./fixtures/section_probe/mphi.json" with { type: "json" };
import isolinesFixture from "./fixtures/section_probe/isolines.json" with { type: "json" };
import type { IsolineRequest, MomentCurvatureRequest, SectionProbeConfig } from "../../pkg/carapace_wasm.js";

const { createSectionProbe } = await loadCarapace();

function assert(condition: boolean, message: string): asserts condition {
  if (!condition) throw new Error(message);
}

const time = <T>(label: string, f: () => T): T => {
  const t0 = performance.now();
  const out = f();
  console.log(`${label}: ${(performance.now() - t0).toFixed(1)} ms`);
  return out;
};

const config = sectionFixture as unknown as SectionProbeConfig;
const probe = time("createSectionProbe (418 fibers)", () => createSectionProbe(config));
assert(probe.fiberCount() === config.fibers.y.length, "fiber count");
assert(probe.groupCount() === 2, "concrete and bar groups");

// PMM states: every engine state through the batch call (typed arrays in and out)
const { states } = statesFixture;
const angles = Float64Array.from(states.map((s) => s.theta));
let worst = 0;
const m_max = Math.max(...states.map((s) => Math.hypot(s.m2, s.m3)));
const triples = new Float64Array(states.length * 3);
states.forEach((s, i) => {
  // one profile per angle: profileStates takes the cross product, so build each state on its own
  const [e0, kz, ky] = probe.profileStates(Float64Array.of(s.theta), Float64Array.of(s.et, s.eb));
  triples.set([e0, kz, ky], 3 * i);
});
const resultants = time(`evaluateBatch (${states.length} states)`, () => probe.evaluateBatch(triples, false));
states.forEach((s, i) => {
  worst = Math.max(worst, Math.abs(resultants[3 * i] - s.n) / m_max, Math.abs(-resultants[3 * i + 1] - s.m2) / m_max, Math.abs(resultants[3 * i + 2] - s.m3) / m_max);
});
assert(worst < 1e-9, `batch parity ${worst}`);
assert(angles.length === states.length, "angles");

// grouped sums add up
const grouped = probe.evaluateBatch(triples.slice(0, 3), true);
assert(grouped.length === 9, "stride 3 * (1 + 2)");
assert(Math.abs(grouped[0] - grouped[3] - grouped[6]) < 1e-6 * Math.max(1, Math.abs(grouped[0])), "group sums");

// axial extremes: the concentric point
const extremes = probe.axialExtremes(0.01);
assert(extremes.compression.response.n < 0 && extremes.tension.response.n > 0, "axial extremes");
console.log(`concentric: strain ${extremes.compression.strain.toExponential(3)}, n ${extremes.compression.response.n.toFixed(0)}, m2 ${extremes.compression.response.m2.toFixed(0)}`);

// moment-curvature against the engine
for (const c of mphiFixture.cases) {
  const request: MomentCurvatureRequest = { axial: c.axialForce, angle: c.angleRad, curvatures: c.points.map((p) => p.curvature) };
  const result = time(`momentCurvature, ${c.name} (${request.curvatures.length} steps)`, () => probe.momentCurvature(request));
  assert(result.nonConverged === 0, `${c.name}: non-converged`);
  const peak = Math.max(...c.points.map((p) => Math.hypot(p.m2, p.m3)));
  result.points.forEach((p, i) => {
    if (i > 0 && c.points[i].converged) assert(Math.abs(p.m2 - c.points[i].m2) / peak < 1e-5, `${c.name}: step ${i}`);
  });
}

// isolines against the engine
const iso = isolinesFixture as unknown as { angles: number[]; profiles: number[]; levels: number[]; momentAnglesRad: number[]; roots: { momentAngleIndex: number; levelIndex: number; moment?: number }[] };
time("buildGrid (72 angles x 53 profiles)", () => probe.buildGrid(Float64Array.from(iso.angles), Float64Array.from(iso.profiles)));
const request: IsolineRequest = { levels: iso.levels, momentAngles: iso.momentAnglesRad };
const solved = time(`solveIsolines (${iso.levels.length} levels x ${iso.momentAnglesRad.length} angles)`, () => probe.solveIsolines(request));
let missing = 0, off = 0;
const mmax = Math.max(...iso.roots.map((r) => r.moment ?? 0));
for (const r of iso.roots) {
  if (r.moment === undefined) continue;
  const p = solved.points[r.momentAngleIndex * iso.levels.length + r.levelIndex];
  if (!p.found) missing++;
  else if (Math.abs(p.moment - r.moment) / mmax > 1e-4) off++;
}
assert(missing === 0, `${missing} engine roots not found`);
assert(off <= 3, `${off} points on another root`);
console.log(`isolines: 0 missing, ${off} on another root`);

// errors surface as { kind } objects
const fresh = createSectionProbe(config);
try { fresh.solveIsolines({ levels: [0], momentAngles: [0] }); assert(false, "solveIsolines before buildGrid must be rejected"); } catch (e) { assert((e as { kind?: string }).kind === "noGrid", `expected noGrid, got ${JSON.stringify(e)}`); }
try { createSectionProbe({ ...config, fibers: { ...config.fibers, area: [1] } }); assert(false, "mismatched fibers must be rejected"); } catch (e) { assert((e as { kind?: string }).kind === "invalidFibers", `got ${JSON.stringify(e)}`); }
fresh.free();
probe.free();
console.log("section-probe-smoke ok");
