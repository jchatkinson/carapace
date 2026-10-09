#!/usr/bin/env python3
"""Builds the docs-site test report from `cargo test` output and the test sources.

Usage: gen_tests.py <cargo-test.log> <out.html>
  cargo test --workspace --no-fail-fast 2>&1 | tee target/test.log
  python3 docs-site/gen_tests.py target/test.log target/doc/tests/index.html

Results come from the log. Descriptions come from the sources: a test's `///` block, the `//!`
header of the file it lives in, and (for the area overview) AREAS below. Standard library only.
"""

import datetime
import html
import os
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PACKAGES = {"carapace_core": ROOT / "core", "carapace_wasm": ROOT / "wasm-bridge"}

AREAS = {
    "elements": "Element formulations (truss, beams, zero-length) against closed forms and each other.",
    "analysis": "Static analysis: solvers, algorithms, load patterns, displacement control, arc-length.",
    "dynamics": "Modal and transient analysis against SDOF closed forms.",
    "constraints": "Equal-DOF, rigid diaphragm and general linear constraints.",
    "continuum": "Plane-stress/strain Tri3 and Quad4 elements: patch tests and beam-theory benchmarks.",
    "materials": "Uniaxial material state handling (trial/commit semantics).",
    "wire": "The CarapaceInputV1 wire format: decoding, validation and end-to-end runs.",
}

TEST_RE = re.compile(r"^test (\S+)(?: - should panic)? \.\.\. (ok|FAILED|ignored)")
RUN_RE = re.compile(r"^\s*Running (unittests )?(\S+) \((?:.*/)?([^/)]+?)(?:-[0-9a-f]{8,})?\)")


ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def parse_log(path):
    """Yields (package_dir, kind, area, test_path, status) for every test line."""
    pkg, kind, area = None, None, None
    out = []
    for line in open(path, errors="replace"):
        line = ANSI_RE.sub("", line)  # CI sets CARGO_TERM_COLOR=always
        m = RUN_RE.match(line)
        if m:
            unit, src, binary = m.group(1), m.group(2), m.group(3)
            if unit:
                pkg, kind, area = PACKAGES.get(binary), "unit", binary
            else:
                parts = Path(src).parts  # tests/<area>/main.rs
                area = parts[1] if len(parts) > 2 else Path(src).stem
                pkg = next((p for p in PACKAGES.values() if (p / src).exists()), None)
                kind = "integration"
            continue
        m = TEST_RE.match(line)
        if m and pkg:
            out.append((pkg, kind, area, m.group(1), m.group(2)))
    return out


_cache = {}


def lines_of(path):
    if path not in _cache:
        _cache[path] = path.read_text().split("\n") if path.exists() else []
    return _cache[path]


def find_file(pkg, kind, area, mods):
    base = (pkg / "tests" / area) if kind == "integration" else (pkg / "src")
    mods = [m for m in mods if not (kind == "unit" and m == "tests")]
    for k in range(len(mods), 0, -1):
        stem = base.joinpath(*mods[:k])
        for cand in (stem.with_suffix(".rs"), stem / "mod.rs"):
            if cand.exists():
                return cand
    return None


def file_header(path):
    out = []
    for line in lines_of(path):
        if line.startswith("//!"):
            out.append(line[3:].removeprefix(" "))
        elif out or line.strip():
            break
    return "\n".join(out).strip()


def test_doc(path, name):
    lines = lines_of(path)
    pat = re.compile(rf"^\s*(pub(\(\w+\))? )?fn {re.escape(name)}\b")
    for i, line in enumerate(lines):
        if pat.match(line):
            doc, j = [], i - 1
            while j >= 0 and (lines[j].strip().startswith("#[") or not lines[j].strip()) and not doc:
                if not lines[j].strip():
                    return ""
                j -= 1
            while j >= 0 and re.match(r"\s*//[/!]?(?!/)", lines[j]) and not lines[j].strip().startswith("//!"):
                doc.append(re.sub(r"^\s*///? ?", "", lines[j]))
                j -= 1
                while j >= 0 and lines[j].strip().startswith("#["):
                    j -= 1
            return "\n".join(reversed(doc)).strip()
    return ""


def inline(text):
    parts = re.split(r"(`[^`]*`)", text)
    out = []
    for i, p in enumerate(parts):
        if i % 2:
            out.append(f"<code>{html.escape(p[1:-1])}</code>")
        else:
            p = html.escape(p)
            p = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", p)
            p = re.sub(r"(?<![\w*])\*(?!\s)([^*\n]+?)(?<!\s)\*(?![\w*])", r"<i>\1</i>", p)
            out.append(p)
    return "".join(out)


def markdown(text):
    """Paragraphs, bullet/numbered lists and inline code/bold/italic — what the doc comments use."""
    blocks, cur, kind = [], [], None

    def flush():
        nonlocal cur, kind
        if cur:
            blocks.append((kind, cur))
        cur, kind = [], None

    for line in text.split("\n"):
        if not line.strip():
            flush()
            continue
        m = re.match(r"\s*([-*]|\d+\.) (.*)", line)
        if m:
            k = "ol" if m.group(1)[0].isdigit() else "ul"
            if kind != k:
                flush()
                kind = k
            cur.append(m.group(2))
        elif kind in ("ul", "ol") and line.startswith("  "):
            cur[-1] += " " + line.strip()
        else:
            if kind != "p":
                flush()
                kind = "p"
            cur.append(line.strip())
    flush()
    out = []
    for kind, items in blocks:
        if kind == "p":
            out.append(f"<p>{inline(' '.join(items))}</p>")
        else:
            out.append(f"<{kind}>" + "".join(f"<li>{inline(i)}</li>" for i in items) + f"</{kind}>")
    return "".join(out)


def split_summary(doc):
    paras = re.split(r"\n\s*\n", doc, maxsplit=1)
    return paras[0], (paras[1] if len(paras) > 1 else "")


def title(name):
    s = name.replace("_", " ").strip()
    return (s[:1].upper() + s[1:]) if s else name


def collect(entries):
    """area -> module -> [(name, status, doc)] with module headers."""
    groups = OrderedDict()
    for pkg, kind, area, path, status in entries:
        parts = path.split("::")
        name, mods = parts[-1], parts[:-1]
        if kind == "integration":
            key = ("Integration tests", area)
        else:
            key = ("Unit tests", area)
        src = find_file(pkg, kind, area, mods)
        doc = test_doc(src, name) if src else ""
        header = file_header(src) if src else ""
        modname = "::".join(mods if kind == "integration" else [m for m in mods if m != "tests"]) or "(root)"
        mod = groups.setdefault(key, OrderedDict()).setdefault(modname, {"header": header, "tests": []})
        mod["tests"].append((name, status, doc))
    return groups


CSS = """
:root{--bg:#e8eefc;--fg:#1b2333;--muted:#4b5670;--accent:#c8321f;--card:#f6f8fe;--border:#cfd8ee;--ok:#1d7a3b;--bad:#c8321f;--skip:#8a6d00}
@media (prefers-color-scheme:dark){:root{--bg:#131826;--fg:#e6ebf7;--muted:#9aa6c2;--accent:#ff7a63;--card:#1b2234;--border:#2c3650;--ok:#58c27d;--bad:#ff7a63;--skip:#d9b84a}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:56rem;margin:0 auto;padding:2rem 1rem 4rem}a{color:var(--accent)}
h1{font-size:2rem;margin:0 0 .25rem;letter-spacing:-.02em}.sub{color:var(--muted);margin:0 0 1.25rem}
h2{font-size:1.1rem;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:2.5rem 0 .5rem}
.area-desc{color:var(--muted);margin:0 0 .75rem}
.mod-doc{background:var(--card);border:1px solid var(--border);border-radius:.5rem;padding:.5rem .9rem;margin:.25rem 0 .5rem;font-size:.92rem;color:var(--muted)}
.mod-doc p{margin:.3rem 0}.stats{display:flex;gap:1rem;flex-wrap:wrap;margin:0 0 1rem}
.stat{background:var(--card);border:1px solid var(--border);border-radius:.6rem;padding:.4rem .9rem}.stat b{font-size:1.3rem;display:block}
.test{background:var(--card);border:1px solid var(--border);border-radius:.6rem;margin:.4rem 0;padding:.5rem .9rem}
.test[hidden]{display:none}.head{display:flex;gap:.6rem;align-items:baseline}
.chip{font-size:.72rem;font-weight:600;border-radius:1rem;padding:.05rem .55rem;border:1px solid currentColor;flex:none}
.ok{color:var(--ok)}.FAILED{color:var(--bad)}.ignored{color:var(--skip)}
.name{font-weight:600}.fn{display:block;font:.78rem ui-monospace,Menlo,monospace;color:var(--muted)}
.desc{color:var(--fg);margin:.35rem 0 0;font-size:.93rem}.desc p{margin:.3rem 0}.desc ul,.desc ol{margin:.3rem 0;padding-left:1.3rem}
code{font:.85em ui-monospace,Menlo,monospace;background:rgba(127,127,127,.15);border-radius:.25rem;padding:0 .25em}
details.more summary{cursor:pointer;color:var(--muted);font-size:.85rem;margin-top:.3rem}
.controls{display:flex;gap:.75rem;align-items:center;margin:1rem 0;flex-wrap:wrap}
input[type=search]{flex:1;min-width:12rem;padding:.45rem .7rem;border-radius:.5rem;border:1px solid var(--border);background:var(--card);color:var(--fg);font:inherit}
details.area{margin:.6rem 0;border:1px solid var(--border);border-radius:.7rem;background:var(--card)}
details.area>summary{padding:.7rem 1rem;font-size:1.15rem;font-weight:600;display:flex;gap:.6rem;align-items:baseline;flex-wrap:wrap}
details.area>summary small{font-weight:400;color:var(--muted);font-size:.85rem}
details.area[open]>summary{border-bottom:1px solid var(--border)}
details.area>.inner{padding:.25rem 1rem 1rem}
details.area .test{background:var(--bg)}
details.module>summary{margin:1rem 0 .25rem;font-family:ui-monospace,Menlo,monospace;font-size:.9rem;color:var(--accent);font-weight:400}
details.area summary,details.module>summary{cursor:pointer;list-style:none}
details.area summary::-webkit-details-marker,details.module>summary::-webkit-details-marker{display:none}
details.area>summary::before,details.module>summary::before{content:"\\25B8";display:inline-block;width:1.1em;color:var(--muted);transition:transform .15s}
details.area[open]>summary::before,details.module[open]>summary::before{transform:rotate(90deg)}
.btns button{font:inherit;font-size:.85rem;padding:.3rem .7rem;border-radius:.5rem;border:1px solid var(--border);background:var(--card);color:var(--fg);cursor:pointer}
footer{margin-top:3rem;color:var(--muted);font-size:.9rem}
"""

JS = """
const q=document.getElementById('q'),f=document.getElementById('failed');
const areas=[...document.querySelectorAll('details.area')],mods=[...document.querySelectorAll('details.module')];
const dflt=()=>{areas.forEach(a=>a.open=a.dataset.bad>0);mods.forEach(m=>m.open=true)};
function apply(){const s=q.value.toLowerCase(),on=!!s||f.checked;
document.querySelectorAll('.test').forEach(t=>{t.hidden=!(t.dataset.text.includes(s)&&(!f.checked||t.dataset.status!=='ok'))});
mods.forEach(m=>{const any=m.querySelector('.test:not([hidden])');m.hidden=!any;if(on&&any)m.open=true});
areas.forEach(a=>{const any=a.querySelector('.test:not([hidden])');a.hidden=!any;if(on&&any)a.open=true});
if(!on)dflt()}
q.oninput=f.onchange=apply;
document.getElementById('expand').onclick=()=>{areas.concat(mods).forEach(d=>d.open=true)};
document.getElementById('collapse').onclick=()=>{areas.forEach(d=>d.open=false)};
dflt();
"""


def render(groups, meta):
    allt = [t for g in groups.values() for m in g.values() for t in m["tests"]]
    n = lambda s: sum(1 for t in allt if t[1] == s)
    out = [f"<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
           f"<title>Carapace tests</title><link rel=icon href=../favicon.png><style>{CSS}</style><main>"
           "<p><a href=../>&larr; Carapace</a></p><h1>Test report</h1>"
           f"<p class=sub>What the test suite verifies, and whether it passes. Commit <code>{html.escape(meta['commit'])}</code>, {html.escape(meta['date'])}.</p>"
           f"<div class=stats><div class=stat><b>{len(allt)}</b>tests</div><div class=stat><b class=ok>{n('ok')}</b>passing</div>"
           f"<div class=stat><b class=FAILED>{n('FAILED')}</b>failing</div><div class=stat><b class=ignored>{n('ignored')}</b>ignored</div></div>"
           "<div class=controls><input id=q type=search placeholder='Filter tests by name or description'>"
           "<label><input id=failed type=checkbox> not passing only</label>"
           "<span class=btns><button id=expand>Expand all</button> <button id=collapse>Collapse all</button></span></div>"]
    last_section = None
    for (section, area), mods in groups.items():
        count = sum(len(m["tests"]) for m in mods.values())
        if section != last_section:
            out.append(f"<h2>{section}</h2>")
            last_section = section
        tests = [x for m in mods.values() for x in m["tests"]]
        bad = sum(1 for x in tests if x[1] == "FAILED")
        badge = f"<small class=FAILED>{bad} failing</small>" if bad else ""
        out.append(f"<details class=area data-bad={bad}><summary>{html.escape(area)}<small>{count} tests</small>{badge}</summary><div class=inner>")
        if section == "Integration tests" and area in AREAS:
            out.append(f"<p class=area-desc>{html.escape(AREAS[area])}</p>")
        for modname, mod in mods.items():
            out.append(f"<details class=module open><summary>{html.escape(modname)}<small> &middot; {len(mod['tests'])}</small></summary>")
            if mod["header"]:
                out.append(f"<div class=mod-doc>{markdown(mod['header'])}</div>")
            for name, status, doc in mod["tests"]:
                summary, rest = split_summary(doc)
                body = f"<div class=desc>{markdown(summary)}</div>" if summary else ""
                if rest:
                    body += f"<details class=more><summary>More</summary><div class=desc>{markdown(rest)}</div></details>"
                text = html.escape((name + " " + doc).lower(), quote=True)
                out.append(f"<div class=test data-status={status} data-text=\"{text}\"><div class=head>"
                           f"<span class='chip {status}'>{status.lower() if status != 'FAILED' else 'failed'}</span>"
                           f"<div><span class=name>{html.escape(title(name))}</span><span class=fn>{html.escape(name)}</span></div></div>{body}</div>")
            out.append("</details>")
        out.append("</div></details>")
    out.append("<footer>Generated from <code>cargo test</code> output and the test sources "
               "(<code>docs-site/gen_tests.py</code>). <a href=https://github.com/jchatkinson/carapace>Source on GitHub</a></footer></main>"
               f"<script>{JS}</script></html>")
    return "".join(out)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    log, dest = sys.argv[1], Path(sys.argv[2])
    entries = parse_log(log)
    if not entries:
        sys.exit(f"no test results found in {log}")
    commit = os.environ.get("GITHUB_SHA", "")[:7] or subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip() or "unknown"
    meta = {"commit": commit, "date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")}
    groups = collect(entries)
    # Integration tests first, then unit tests.
    groups = OrderedDict(sorted(groups.items(), key=lambda kv: (kv[0][0] != "Integration tests",)))
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(render(groups, meta))
    print(f"wrote {dest}: {len(entries)} tests")


if __name__ == "__main__":
    main()
