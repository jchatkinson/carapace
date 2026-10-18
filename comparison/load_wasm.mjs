// Loads the shipped wasm build (pkg/, from scripts/build-wasm.sh: the same binary PySees bundles) under Node.
// The glue is `--target web` ESM; it is imported from a data: URL so a missing "type": "module" in
// pkg/package.json does not matter, and the bytes are handed to initSync instead of fetch.
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

// CARAPACE_PKG=<dir> selects another build (e.g. a variant from scripts/build-wasm.sh); default is pkg/.
export async function loadCarapace(dir = process.env.CARAPACE_PKG ? pathToFileURL(process.env.CARAPACE_PKG + "/") : new URL("../pkg/", import.meta.url)) {
  const glue = readFileSync(new URL("carapace_wasm.js", dir), "utf8");
  const mod = await import(`data:text/javascript;base64,${Buffer.from(glue).toString("base64")}`);
  mod.initSync({ module: readFileSync(new URL("carapace_wasm_bg.wasm", dir)) });
  return mod;
}
