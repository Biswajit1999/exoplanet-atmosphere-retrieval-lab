"use strict";

const assert = require("node:assert/strict");
const test = require("node:test");
require("../inference.js");

function synthetic() {
  const points = [];
  for (const [instrument, offset] of [["A", 20_000], ["B", 20_250]]) {
    for (let i = 0; i < 90; i += 1) {
      const wavelength = 0.6 + 4.6 * i / 89;
      const t = global.Inference.templates(wavelength, 1);
      points.push({
        wavelength_um: wavelength,
        transit_depth_ppm: offset + 420 * t.h2o + 700 * t.co2 + 80 * t.so2 + 22 * Math.log(wavelength / 2.2),
        uncertainty_ppm: 35,
        instrument,
      });
    }
  }
  return points;
}

test("profile likelihood recovers injected template amplitudes", () => {
  const result = global.Inference.infer(synthetic(), { gridSize: 41, amplitudeMax: 800 });
  assert.equal(result.metrics.h2o_amplitude_ppm, 420);
  assert.equal(result.metrics.co2_amplitude_ppm, 700);
  assert.ok(result.metrics.reduced_chi2 < 1e-12);
});

test("instrument offsets are profiled independently", () => {
  const result = global.Inference.infer(synthetic(), { gridSize: 21, amplitudeMax: 800 });
  assert.ok(Math.abs((result.offsets.B - result.offsets.A) - 250) < 1e-6);
});

test("an additional error floor cannot increase chi-squared", () => {
  const data = synthetic();
  data[10].transit_depth_ppm += 300;
  const strict = global.Inference.infer(data, { gridSize: 21, amplitudeMax: 800, errorFloor: 0 });
  const inflated = global.Inference.infer(data, { gridSize: 21, amplitudeMax: 800, errorFloor: 100 });
  assert.ok(inflated.metrics.chi2 < strict.metrics.chi2);
});

test("short selections fail closed", () => {
  assert.throws(() => global.Inference.infer(synthetic().slice(0, 9)), /ten spectral bins/);
});
