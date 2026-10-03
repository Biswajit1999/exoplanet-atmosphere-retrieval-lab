"use strict";

const crypto = require("node:crypto");
const fs = require("node:fs");
const path = require("node:path");

const failures = [];
const required = ["README.md", "index.html", "styles.css", "app.js", "inference.js", "physicsWorker.js", "data/wasp39b_benchmark.json"];
for (const file of required) if (!fs.existsSync(file)) failures.push(`${file} missing`);

const dataset = JSON.parse(fs.readFileSync("data/wasp39b_benchmark.json", "utf8"));
if (dataset.points.length !== 355) failures.push(`expected 355 included bins, found ${dataset.points.length}`);
if (dataset.excluded_saturated_prism_bins !== 35) failures.push("unexpected saturated-bin exclusion count");
if (dataset.products.length !== 6) failures.push("expected six source products");
for (const product of dataset.products) {
  const file = path.join("data", "source_products", product.file);
  if (!fs.existsSync(file)) { failures.push(`${file} missing`); continue; }
  const digest = crypto.createHash("sha256").update(fs.readFileSync(file)).digest("hex");
  if (digest !== product.sha256) failures.push(`${product.file} SHA-256 mismatch`);
}
for (const point of dataset.points) {
  if (![point.wavelength_um, point.transit_depth_ppm, point.uncertainty_ppm].every(Number.isFinite)) failures.push("non-finite spectral bin");
  if (point.uncertainty_ppm <= 0) failures.push("non-positive uncertainty");
}
const publicText = ["README.md", "index.html", "app.js"].map((file) => fs.readFileSync(file, "utf8")).join("\n").toLowerCase();
for (const phrase of ["research quality upgrade", "before and after", "maturity score"]) if (publicText.includes(phrase)) failures.push(`public self-scoring phrase: ${phrase}`);
if (!publicText.includes("10.5281/zenodo.10161743")) failures.push("dataset DOI missing from public documentation");
if (!publicText.includes("model inadequacy")) failures.push("model-inadequacy boundary missing");
if (failures.length) { console.error(failures.join("\n")); process.exit(1); }
console.log("Validated six pinned products, 355 finite bins, inference boundaries, and public presentation.");
