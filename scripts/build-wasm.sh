#!/usr/bin/env bash
# The one wasm build: cargo (+simd128) -> wasm-bindgen --target web -> load check.
# Used by PySees (../pysees/scripts/build-carapace.sh), CI and the comparison harness, so all use one binary.
# Usage: scripts/build-wasm.sh [out-dir]   (default: pkg)
# wasm-opt is off by default: measured under Node/V8, -O2 and above made the solve 30-45% slower than the
# unoptimized module and -O1 matched it (size is not a concern), and SIMD made no measurable difference.
# Experiments: WASM_SIMD=0 drops simd128; WASM_OPT_LEVEL=-O1|-O2|... runs wasm-opt (binaryen >= 120, `npm i -g
# binaryen`; WASM_OPT=/path overrides the one on PATH).
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out="${1:-$root/pkg}"
wasm_opt="${WASM_OPT:-wasm-opt}"
simd="${WASM_SIMD:-1}"
opt_level="${WASM_OPT_LEVEL:-none}"

command -v wasm-bindgen >/dev/null || { echo "error: wasm-bindgen not found" >&2; exit 1; }
[ "$opt_level" = none ] || command -v "$wasm_opt" >/dev/null || { echo "error: wasm-opt not found (binaryen >= 120: npm i -g binaryen)" >&2; exit 1; }
# binaryen 108 (Ubuntu's apt package) mis-emits the module's tables: the externref table export points at the
# funcref table, so the module throws in `Table.grow` on load. Refuse anything old enough to do that.
version="$([ "$opt_level" = none ] && echo 999 || "$wasm_opt" --version | grep -oE '[0-9]+' | head -1)"
[ "${version:-0}" -ge 120 ] || { echo "error: $wasm_opt is version $version; need binaryen >= 120 (npm i -g binaryen)" >&2; exit 1; }

cd "$root"
features=""; [ "$simd" = 0 ] || features="-C target-feature=+simd128"
RUSTFLAGS="${RUSTFLAGS:-} $features" \
  cargo build -p carapace-wasm --target wasm32-unknown-unknown --release
wasm-bindgen target/wasm32-unknown-unknown/release/carapace_wasm.wasm --out-dir "$out" --target web
[ "$opt_level" = none ] || "$wasm_opt" "$opt_level" --enable-simd --enable-bulk-memory --enable-nontrapping-float-to-int \
  --enable-reference-types --enable-multivalue --enable-sign-ext --enable-mutable-globals \
  "$out/carapace_wasm_bg.wasm" -o "$out/carapace_wasm_bg.wasm"

# The module must still instantiate (wasm-opt 108 broke it this way).
node --input-type=module -e '
  import { pathToFileURL } from "node:url";
  import { loadCarapace } from "./comparison/load_wasm.mjs";
  await loadCarapace(pathToFileURL(process.argv[1] + "/"));
' "$(cd "$out" && pwd)"
