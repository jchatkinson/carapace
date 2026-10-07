#!/usr/bin/env bash
# Before/after check for the domain refactor (docs/domain-refactor-plan.md).
# Compares the current tree against results recorded at tag pre-domain-refactor:
#   1. baseline_dump JSON (comparison/baseline/dump_pre-domain-refactor.json)
#   2. cyclic material CSVs (comparison/baseline/csv/*_carapace.csv), byte-identical
# Run from anywhere: comparison/check_refactor.sh [--rtol 1e-9] [--atol 1e-12]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "== baseline_dump =="
cargo build -q --release -p carapace-core --example baseline_dump
"$ROOT/target/release/examples/baseline_dump" > "$TMP/dump.json"
python3 comparison/compare_dumps.py comparison/baseline/dump_pre-domain-refactor.json "$TMP/dump.json" "$@"

echo "== cyclic material CSVs =="
cargo build -q --release -p carapace-core --example cyclic_material
status=0
for expected in comparison/baseline/csv/*_carapace.csv; do
  name="$(basename "$expected" _carapace.csv)"
  "$ROOT/target/release/examples/cyclic_material" "$name" "$TMP/$name.csv" >/dev/null
  if cmp -s "$expected" "$TMP/$name.csv"; then echo "  identical: $name"; else echo "  DIFFERS:   $name"; status=1; fi
done
exit $status
