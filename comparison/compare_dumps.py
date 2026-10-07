#!/usr/bin/env python3
"""Diff two `baseline_dump` JSON files (core/examples/baseline_dump.rs).

Usage: compare_dumps.py BEFORE.json AFTER.json [--rtol 1e-9] [--atol 1e-12]

Keys must match exactly. Values match when |a-b| <= atol + rtol*max(|a|,|b|).
Iteration and factorization counters are part of the dump and must therefore
match as well (they are small integers, so the tolerance never hides a change).
Exits 0 when everything matches, 1 otherwise, and prints the worst offenders.
"""
import argparse
import json
import math
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("before")
    parser.add_argument("after")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()

    with open(args.before) as f:
        before = json.load(f)
    with open(args.after) as f:
        after = json.load(f)

    problems = []
    for key in sorted(set(before) - set(after)):
        problems.append((math.inf, f"missing in AFTER: {key}"))
    for key in sorted(set(after) - set(before)):
        problems.append((math.inf, f"extra in AFTER: {key}"))

    worst = 0.0
    compared = 0
    for key in sorted(set(before) & set(after)):
        a, b = before[key], after[key]
        if len(a) != len(b):
            problems.append((math.inf, f"{key}: length {len(a)} -> {len(b)}"))
            continue
        for i, (x, y) in enumerate(zip(a, b)):
            compared += 1
            diff = abs(x - y)
            scale = max(abs(x), abs(y))
            if math.isnan(x) != math.isnan(y) or diff > args.atol + args.rtol * scale:
                rel = diff / scale if scale else diff
                problems.append((rel, f"{key}[{i}]: {x!r} -> {y!r} (rel {rel:.3e})"))
            else:
                worst = max(worst, diff / scale if scale else 0.0)

    print(f"compared {compared} values in {len(set(before) & set(after))} entries; "
          f"largest accepted relative difference {worst:.3e}")
    if problems:
        problems.sort(key=lambda p: -p[0])
        print(f"{len(problems)} MISMATCHES (worst first):")
        for _, text in problems[:40]:
            print("  " + text)
        return 1
    print("OK: dumps match")
    return 0


if __name__ == "__main__":
    sys.exit(main())
