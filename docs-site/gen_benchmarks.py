#!/usr/bin/env python3
"""Builds the docs-site benchmark page from the committed results in docs-site/benchmarks/*.json.

Usage: gen_benchmarks.py <out.html>
  python3 docs-site/gen_benchmarks.py target/doc/benchmarks/index.html

The results are produced locally by comparison/run_benchmarks.py (they depend on the machine and need
OpenSeesPy), so CI only renders them. Standard library only.
"""

import html
import json
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
.chartbox svg{display:block;width:100%;height:auto}.chartbox figcaption{color:var(--muted);font-size:.85rem;margin-top:.3rem}
svg text{fill:var(--fg);font:12px system-ui,sans-serif}svg .dim{fill:var(--muted)}
svg .grid{stroke:var(--border);stroke-width:1}svg .one{stroke:var(--fg);stroke-width:1.5;stroke-dasharray:4 3}
svg .fast{fill:var(--ok)}svg .slow{fill:var(--accent)}
table{border-collapse:collapse;width:100%;font-size:.9rem;margin:.5rem 0;overflow-x:auto;display:block}
th,td{padding:.35rem .6rem;border-bottom:1px solid var(--border);text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}th{color:var(--muted);font-weight:600}
td.ratio{font-weight:600}.fastc{color:var(--ok)}.slowc{color:var(--accent)}
"""


def fmt_ms(ms):
    return f"{ms / 1000:.2f} s" if ms >= 1000 else f"{ms:.1f} ms"


def ratio_chart(cases):
    """Horizontal bars of Carapace time / OpenSees time per case, with a dashed line at parity."""
    ratios = [c["engines"]["carapace"]["median_ms"] / c["engines"]["opensees"]["median_ms"] for c in cases]
    top = max(1.2, max(ratios) * 1.1)
    left, width, row, pad = 210, 460, 26, 24
    height = pad + row * len(cases) + 22
    x = lambda value: left + width * value / top
    out = [f"<svg viewBox='0 0 {left + width + 60} {height}' role=img aria-label='Carapace time relative to OpenSees, per case'>"]
    step = 0.25 if top <= 2 else 0.5
    tick = 0.0
    while tick <= top:
        out.append(f"<line class=grid x1={x(tick):.1f} x2={x(tick):.1f} y1={pad - 6} y2={pad + row * len(cases)} />"
                   f"<text class=dim x={x(tick):.1f} y={height - 6} text-anchor=middle>{tick:g}</text>")
        tick += step
    for k, (case, r) in enumerate(zip(cases, ratios)):
        y = pad + row * k
        out.append(f"<text x={left - 8} y={y + row / 2 + 4} text-anchor=end>{html.escape(case['label'])}</text>"
                   f"<rect class='{'fast' if r <= 1 else 'slow'}' x={left} y={y + 4} width={x(r) - left:.1f} height={row - 8} rx=3 />"
                   f"<text x={x(r) + 6:.1f} y={y + row / 2 + 4}>{r:.2f}</text>")
    out.append(f"<line class=one x1={x(1):.1f} x2={x(1):.1f} y1={pad - 8} y2={pad + row * len(cases)} />"
               f"<text class=dim x={x(1):.1f} y=12 text-anchor=middle>OpenSees = 1</text></svg>")
    return "".join(out)


def bench_section(b):
    cases = b["cases"]
    rows = []
    for c in cases:
        ca, os_ = c["engines"]["carapace"], c["engines"]["opensees"]
        r = ca["median_ms"] / os_["median_ms"]
        iters = ca["iterations"] if ca["iterations"] == os_["iterations"] else f"{ca['iterations']} / {os_['iterations']}"
        rows.append(f"<tr><td>{html.escape(c['label'])}</td><td>{c['free_dofs']:,}</td><td>{c['elements']:,}</td>"
                    f"<td>{fmt_ms(ca['median_ms'])}</td><td>{fmt_ms(os_['median_ms'])}</td>"
                    f"<td class='ratio {'fastc' if r <= 1 else 'slowc'}'>{r:.2f}</td><td>{iters}</td>"
                    f"<td>{c['max_relative_difference']:.1e}</td></tr>")
    return (f"<section class=bench><h2>{html.escape(b['title'])}</h2><p>{html.escape(b['description'])}</p>"
            f"<p class=method>OpenSees: {html.escape(b['solver'])}. {html.escape(b.get('note', ''))}</p>"
            f"<figure class=chartbox>{ratio_chart(cases)}<figcaption>Carapace solve time divided by OpenSees solve time "
            "(median). Below the dashed line, Carapace is faster.</figcaption></figure>"
            "<table><thead><tr><th>Case</th><th>Free DOFs</th><th>Elements</th><th>Carapace</th><th>OpenSees</th>"
            "<th>Ratio</th><th>Newton iterations</th><th>Result difference</th></tr></thead><tbody>"
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
            "<p class=sub>Carapace&rsquo;s native Rust engine against OpenSees on the same models, with results checked "
            "against each other.</p>" + machine_block(benchmarks) +
            "<ul class=method>"
            "<li>Each benchmark builds the same model in both engines with the same convergence test, and fails if their "
            "results differ by more than 1e-6 (relative).</li>"
            f"<li>Times are the median of {runs} timed runs after a discarded warm-up. Engine order alternates between runs. "
            "Only the solve is timed; start-up and model building are excluded.</li>"
            "<li>These are native timings, not the WebAssembly build, which runs slower than native. They were measured on one "
            "machine and are not run in CI, so absolute times will differ elsewhere; the ratios are the useful part.</li>"
            "<li>The ratio is Carapace time divided by OpenSees time, so below 1 means Carapace is faster. "
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
    benchmarks = [found[i] for i in ORDER if i in found] + [b for i, b in sorted(found.items()) if i not in ORDER]
    dest = Path(sys.argv[1])
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(render(benchmarks))
    print(f"wrote {dest}: {len(benchmarks)} benchmarks")


if __name__ == "__main__":
    main()
