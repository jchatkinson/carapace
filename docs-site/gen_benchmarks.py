#!/usr/bin/env python3
"""Builds the docs-site benchmark page from the committed results in docs-site/benchmarks/*.json.

Usage: gen_benchmarks.py <out.html>
  python3 docs-site/gen_benchmarks.py target/doc/benchmarks/index.html

The results are produced locally by comparison/run_benchmarks.py (they depend on the machine and need
OpenSeesPy), so CI only renders them. Standard library only.
"""

import html
import json
import math
import sys
from pathlib import Path

from gen_tests import CSS as BASE_CSS

HERE = Path(__file__).resolve().parent
ORDER = ["frame-nonlinear", "frame-elastic", "panel-quad4"]

CSS = BASE_CSS + """
.machine{background:var(--card);border:1px solid var(--border);border-radius:.6rem;padding:.6rem 1rem;margin:0 0 1rem;font-size:.92rem}
.machine dl{display:grid;grid-template-columns:max-content 1fr;gap:.15rem 1rem;margin:0}.machine dt{color:var(--muted)}.machine dd{margin:0}
.method{color:var(--muted);font-size:.92rem}.method li{margin:.2rem 0}
section.bench{margin:1.5rem 0 2.5rem}section.bench h2{margin-top:0}
.chartbox{background:var(--card);border:1px solid var(--border);border-radius:.6rem;padding:.6rem .8rem;margin:.75rem 0}
.chartbox{overflow-x:auto}.chartbox svg{display:block;width:100%;min-width:32rem;max-width:42rem;margin:0 auto;height:auto}.chartbox figcaption{color:var(--muted);font-size:.85rem;margin-top:.3rem}
svg text{fill:var(--fg);font:11px system-ui,sans-serif}svg .case{font-weight:600}svg .note,svg .tick{fill:var(--muted);font-size:10px}
svg .grid{stroke:var(--border);stroke-width:1;opacity:.6}svg .sep{stroke:var(--border);stroke-width:1}svg .zone{fill:var(--accent);fill-opacity:.06}svg .one{stroke:var(--fg);stroke-width:1.5;stroke-dasharray:4 3}
svg .fast{fill:var(--ok)}svg .slow{fill:var(--accent)}
table{border-collapse:collapse;width:100%;font-size:.9rem;margin:.5rem 0;overflow-x:auto;display:block}
th,td{padding:.35rem .6rem;border-bottom:1px solid var(--border);text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600}
td.faster{color:var(--accent);font-weight:600}td.failed{color:var(--muted);font-style:italic}tr.ours td{font-weight:700;background:var(--card)}
tr.count td{color:var(--muted);border-top:2px solid var(--border)}th small{font-weight:400}
svg .dot{fill:var(--muted);fill-opacity:.7;stroke:var(--card);stroke-width:1}svg .dot.best{fill:var(--accent);fill-opacity:1;stroke:var(--fg);stroke-width:1.5}
td.best{outline:2px solid var(--accent);outline-offset:-2px;font-weight:700}
"""


def fmt_ms(ms):
    return f"{ms / 1000:.2f} s" if ms >= 1000 else f"{ms:.1f} ms"


def cell(record, carapace_ms, carapace_iterations):
    """Table cell (class, text) for one OpenSees configuration at one case."""
    status = record["status"]
    if status == "ok":
        text = fmt_ms(record["median_ms"])
        if record["iterations"] != carapace_iterations:  # rounding can change how many Newton iterations a run needs
            text += f" <small>{record['iterations']} it</small>"
        return ("faster" if record["median_ms"] < carapace_ms else "slower"), text
    if status == "timeout":
        return "failed", f"&gt; {fmt_ms(record['limit_ms'])}"
    if status == "differs":
        return "failed", f"differs ({record['max_relative_difference']:.0e})"
    return "failed", "error"


def fastest(case):
    """Index of the fastest OpenSees configuration in a case, or None if none finished."""
    done = [k for k, c in enumerate(case["opensees"]) if c["status"] == "ok"]
    return min(done, key=lambda k: case["opensees"][k]["median_ms"], default=None)


WIDE_TICKS = [2.0 ** k for k in range(-6, 7)]
NARROW_TICKS = [round(0.5 + 0.1 * k, 1) for k in range(0, 26)]  # 0.5, 0.6, ... 3.0


def axis_ticks(low, high, position, gap=46):
    """Tick values within [low, high]: powers of two over a wide range, steps of 0.1 over a narrow one.

    Starts from 1 and works outward, keeping each tick at least `gap` units (of `position`) from the last.
    """
    values = WIDE_TICKS if high / low > 6 else NARROW_TICKS
    kept = [1.0]
    for direction in ([v for v in values if v > 1], [v for v in values if v < 1][::-1]):
        last = 1.0
        for value in direction:
            if low <= value <= high and abs(position(value) - position(last)) >= gap:
                kept.append(value)
                last = value
    return sorted(kept)


def swarm(dots, radius):
    """Vertical offsets that keep dots (x, is_best) from overlapping; the fastest is placed first, at zero."""
    placed, offsets = [], {}
    order = sorted(range(len(dots)), key=lambda k: (not dots[k][1], dots[k][0]))
    for k in order:
        x = dots[k][0]
        candidates = [0] + [sign * step for step in range(4, 21, 4) for sign in (-1, 1)]
        gap = lambda dy: min((math.hypot(x - px, dy - py) for px, py in placed), default=99.0)
        choice = next((dy for dy in candidates if gap(dy) >= 2 * radius + 1), max(candidates, key=gap))
        placed.append((x, choice))
        offsets[k] = choice
    return offsets


def strip_plot(cases):
    """One row per case: a dot for every OpenSees configuration on a log axis of its time relative to Carapace's.

    The axis fits the data of the whole benchmark so rows stay comparable. Each row's label sits above its dots
    so it can never clip, and overlapping dots are spread vertically. The fastest configuration is highlighted.
    """
    ratios = [c["median_ms"] / case["carapace"]["median_ms"] for case in cases for c in case["opensees"] if c["status"] == "ok"]
    low, high = min(ratios + [1.0]) / 1.2, max(ratios + [1.0]) * 1.2
    left, width, row, pad, radius = 14, 560, 72, 30, 3.5
    height = pad + row * len(cases) + 40
    span = math.log2(high / low)
    x = lambda ratio: left + width * math.log2(ratio / low) / span
    bottom = pad + row * len(cases)
    out = [f"<svg viewBox='0 0 {left + width + 14} {height}' role=img "
           "aria-label='OpenSees time relative to Carapace for every solver configuration, per case'>",
           f"<rect class=zone x={left} y={pad - 10} width={x(1) - left:.1f} height={bottom - pad + 10} />"]
    for tick in axis_ticks(low, high, x):
        out.append(f"<line class=grid x1={x(tick):.1f} x2={x(tick):.1f} y1={pad - 10} y2={bottom} />"
                   f"<text class=tick x={x(tick):.1f} y={bottom + 14} text-anchor=middle>{tick:g}&times;</text>")
    for k, case in enumerate(cases):
        top = pad + row * k
        centre = top + 46
        carapace, best = case["carapace"]["median_ms"], fastest(case)
        if k:
            out.append(f"<line class=sep x1={left} x2={left + width} y1={top} y2={top} />")
        stopped = sum(1 for c in case["opensees"] if c["status"] != "ok")
        note = []
        if best is not None:
            b = case["opensees"][best]
            note.append(f"fastest OpenSees: {html.escape(b['solver'])}, {html.escape(b['numberer'])} ({b['median_ms'] / carapace:.2f}&times;)")
        if stopped:
            note.append(f"{stopped} stopped")
        out.append(f"<text class=case x={left} y={top + 16}>{html.escape(case['label'])}</text>"
                   f"<text class=note x={left + width} y={top + 16} text-anchor=end>{' &middot; '.join(note)}</text>")
        ok = [(index, c) for index, c in enumerate(case["opensees"]) if c["status"] == "ok"]
        points = [(x(c["median_ms"] / carapace), index == best) for index, c in ok]
        offsets = swarm(points, radius)
        for n, (index, c) in enumerate(ok):
            ratio = c["median_ms"] / carapace
            title = html.escape(f"{c['solver']}, {c['numberer']}: {fmt_ms(c['median_ms'])} ({ratio:.2f}x Carapace)"
                                + (" - fastest OpenSees" if index == best else ""))
            out.append(f"<circle class='dot{' best' if index == best else ''}' cx={x(ratio):.1f} cy={centre + offsets[n]} "
                       f"r={5 if index == best else radius}><title>{title}</title></circle>")
    out.append(f"<line class=one x1={x(1):.1f} x2={x(1):.1f} y1={pad - 10} y2={bottom} />"
               f"<text class=tick x={x(1):.1f} y=14 text-anchor=middle>Carapace</text>"
               f"<text class=tick x={left} y={height - 6}>&larr; OpenSees faster</text>"
               f"<text class=tick x={left + width} y={height - 6} text-anchor=end>Carapace faster &rarr;</text></svg>")
    return "".join(out)


def bench_section(b):
    cases = b["cases"]
    head = "".join(f"<th>{html.escape(c['label'])}<br><small>{c['free_dofs']:,} DOFs</small></th>" for c in cases)
    rows = ["<tr class=ours><td>Carapace</td>" + "".join(
        f"<td>{fmt_ms(c['carapace']['median_ms'])} <small>{c['carapace']['iterations']} it</small></td>" for c in cases) + "</tr>"]
    best = [fastest(c) for c in cases]
    for k, first in enumerate(cases[0]["opensees"]):
        label = f"OpenSees {html.escape(first['solver'])}, {html.escape(first['numberer'])}"
        cells = []
        for c in cases:
            klass, text = cell(c["opensees"][k], c["carapace"]["median_ms"], c["carapace"]["iterations"])
            if k == best[cases.index(c)]:
                klass += " best"
            cells.append(f"<td class='{klass}'{' title=\"Fastest OpenSees configuration for this case\"' if 'best' in klass else ''}>{text}</td>")
        rows.append(f"<tr><td>{label}</td>{''.join(cells)}</tr>")
    counts = []
    for c in cases:
        faster = sum(1 for r in c["opensees"] if r["status"] == "ok" and r["median_ms"] < c["carapace"]["median_ms"])
        counts.append(f"<td>{faster} of {len(c['opensees'])}</td>")
    rows.append(f"<tr class=count><td>OpenSees configurations faster than Carapace</td>{''.join(counts)}</tr>")
    return (f"<section class=bench><h2>{html.escape(b['title'])}</h2><p>{html.escape(b['description'])}</p>"
            f"<p class=method>{html.escape(b.get('note', ''))}</p>"
            f"<figure class=chartbox>{strip_plot(cases)}<figcaption>Each dot is one OpenSees solver and numberer "
            "combination: its solve time divided by Carapace&rsquo;s (log scale). Dots left of the line are "
            "faster than Carapace. The highlighted dot is the fastest OpenSees configuration for that case; overlapping dots are spread vertically. Hover a dot "
            "for its configuration. Combinations stopped for exceeding the time "
            "limit, or that disagreed with Carapace, have no dot.</figcaption></figure>"
            f"<table><thead><tr><th>Median solve time</th>{head}</tr></thead><tbody>"
            + "".join(rows) + "</tbody></table></section>")


def machine_block(benchmarks):
    metas = {json.dumps(b["machine"], sort_keys=True): b["machine"] for b in benchmarks}
    out = []
    for m in metas.values():
        dirty = " (with uncommitted changes)" if m.get("dirty") else ""
        out.append("<div class=machine><dl>"
                   f"<dt>CPU</dt><dd>{html.escape(m['cpu'])}, {m['threads']} threads</dd>"
                   f"<dt>OS</dt><dd>{html.escape(m['os'])}</dd>"
                   f"<dt>Rust</dt><dd>{html.escape(m['rustc'])}</dd>"
                   f"<dt>OpenSees</dt><dd>OpenSeesPy {html.escape(m['openseespy'])}</dd>"
                   f"<dt>Carapace</dt><dd>commit <code>{html.escape(m['commit'])}</code>{dirty}, run {html.escape(m['date'])}</dd>"
                   "</dl></div>")
    return "".join(out)


def render(benchmarks):
    runs = benchmarks[0]["runs"]
    return ("<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
            f"<title>Carapace benchmarks</title><link rel=icon href=../favicon.png><style>{CSS}</style><main>"
            "<p><a href=../>&larr; Carapace</a></p><h1>Benchmarks</h1>"
            "<p class=sub>Carapace&rsquo;s native Rust engine against OpenSees on the same models, with OpenSees run "
            "under every solver configuration.</p>" + machine_block(benchmarks) +
            "<ul class=method>"
            "<li>Each benchmark builds the same model in both engines with the same convergence test. A configuration "
            "whose results differ from Carapace&rsquo;s by more than 1e-6 (relative) is shown as &ldquo;differs&rdquo; "
            "instead of a time.</li>"
            f"<li>Carapace times are the median of {runs} timed runs after {benchmarks[0]['warmups']} discarded warm-up. "
            f"OpenSees was timed with <b>every</b> solver and numberer combination, each the median of "
            f"{benchmarks[0]['opensees_runs']} runs after {benchmarks[0]['opensees_warmups']} warm-up, and all are published: "
            "none was chosen for the comparison. Only the solve is timed; start-up and model building are excluded.</li>"
            "<li>OpenSees&rsquo; band and profile solvers (BandGeneral, BandSPD, ProfileSPD) were run with the Plain and RCM "
            "numberers only. They store a band set by the node numbering, and AMD, which minimises fill for sparse "
            "solvers, gives them an almost dense band. The sparse solvers (SuperLU, UmfPack, SparseSPD, SparseSYM) "
            "were run with all three numberers.</li>"
            f"<li>A combination whose run exceeded {benchmarks[0]['timeout_factor']}&times; Carapace&rsquo;s time (at least "
            f"{benchmarks[0]['min_timeout_s']:g} s) was stopped and shown as &ldquo;&gt;&rdquo; its limit.</li>"
            "<li>These are native timings, not the WebAssembly build, which runs slower than native. They were measured on one "
            "machine and are not run in CI, so absolute times will differ elsewhere; the ratios are the useful part.</li>"
            "<li>Orange cells are OpenSees configurations faster than Carapace; the outlined cell in each column is the fastest "
            "OpenSees configuration for that case. <i>it</i> is the total number of Newton "
            "iterations; it is shown for Carapace and wherever an OpenSees run needed a different number. Round-off in "
            "the linear solver can change the count (the convergence test sits near round-off on the elastic "
            "frames), and a different count changes the time without changing solver speed.</li>"
            "<li>"
            "Reproduce with <code>comparison/run_benchmarks.py</code>.</li></ul>"
            + "".join(bench_section(b) for b in benchmarks) +
            "<footer>Generated from <code>docs-site/benchmarks/*.json</code> (<code>docs-site/gen_benchmarks.py</code>). "
            "<a href=https://github.com/jchatkinson/carapace>Source on GitHub</a></footer></main></html>")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    found = {json.loads(p.read_text())["id"]: json.loads(p.read_text()) for p in (HERE / "benchmarks").glob("*.json")}
    if not found:
        sys.exit("no benchmark results in docs-site/benchmarks/")
    stale = [i for i, b in found.items() if not all("opensees" in c and "carapace" in c for c in b["cases"])]
    if stale:
        sys.exit(f"results for {', '.join(stale)} use an old format: rerun comparison/run_benchmarks.py")
    benchmarks = [found[i] for i in ORDER if i in found] + [b for i, b in sorted(found.items()) if i not in ORDER]
    dest = Path(sys.argv[1])
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(render(benchmarks))
    print(f"wrote {dest}: {len(benchmarks)} benchmarks")


if __name__ == "__main__":
    main()
