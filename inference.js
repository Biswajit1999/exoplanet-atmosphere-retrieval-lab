(function (root) {
  "use strict";

  const gaussian = (x, center, width) => Math.exp(-0.5 * ((x - center) / width) ** 2);

  function templates(wavelength, widthScale = 1) {
    const g = (center, width) => gaussian(wavelength, center, width * widthScale);
    return {
      h2o: 0.62 * g(1.15, 0.08) + g(1.40, 0.12) + 0.42 * g(1.90, 0.14) + 0.48 * g(2.70, 0.18),
      co2: g(4.30, 0.16),
      so2: g(4.05, 0.07),
    };
  }

  function solve(matrix, vector) {
    const n = vector.length;
    const a = matrix.map((row, i) => [...row, vector[i]]);
    for (let col = 0; col < n; col += 1) {
      let pivot = col;
      for (let row = col + 1; row < n; row += 1) {
        if (Math.abs(a[row][col]) > Math.abs(a[pivot][col])) pivot = row;
      }
      if (Math.abs(a[pivot][col]) < 1e-12) return null;
      [a[col], a[pivot]] = [a[pivot], a[col]];
      const scale = a[col][col];
      for (let j = col; j <= n; j += 1) a[col][j] /= scale;
      for (let row = 0; row < n; row += 1) {
        if (row === col) continue;
        const factor = a[row][col];
        for (let j = col; j <= n; j += 1) a[row][j] -= factor * a[col][j];
      }
    }
    return a.map((row) => row[n]);
  }

  function profile(data, h2oAmplitude, co2Amplitude, options) {
    const instruments = [...new Set(data.map((point) => point.instrument))].sort();
    const columns = [...instruments.map((name) => ({ type: "offset", name }))];
    if (options.includeSlope) columns.push({ type: "slope", name: "spectral_slope" });
    columns.push({ type: "so2", name: "so2_template" });
    const normal = Array.from({ length: columns.length }, () => Array(columns.length).fill(0));
    const rhs = Array(columns.length).fill(0);
    const rows = [];

    for (const point of data) {
      const t = templates(point.wavelength_um, options.widthScale);
      const fixed = h2oAmplitude * t.h2o + co2Amplitude * t.co2;
      const y = point.transit_depth_ppm - fixed;
      const sigma = Math.hypot(point.uncertainty_ppm, options.errorFloor);
      const weight = 1 / sigma ** 2;
      const x = columns.map((column) => {
        if (column.type === "offset") return Number(column.name === point.instrument);
        if (column.type === "slope") return Math.log(point.wavelength_um / 2.2);
        return t.so2;
      });
      rows.push({ point, t, fixed, sigma, x });
      for (let i = 0; i < columns.length; i += 1) {
        rhs[i] += weight * x[i] * y;
        for (let j = 0; j < columns.length; j += 1) normal[i][j] += weight * x[i] * x[j];
      }
    }

    const coefficients = solve(normal, rhs);
    if (!coefficients) return null;
    let chi2 = 0;
    const predictions = [];
    const residuals = [];
    for (const row of rows) {
      const nuisance = row.x.reduce((sum, value, index) => sum + value * coefficients[index], 0);
      const model = row.fixed + nuisance;
      const residual = row.point.transit_depth_ppm - model;
      chi2 += (residual / row.sigma) ** 2;
      predictions.push(model);
      residuals.push(residual);
    }
    return {
      chi2,
      columns,
      coefficients,
      predictions,
      residuals,
      nuisanceCount: columns.length,
    };
  }

  function quantile(axis, weights, probability) {
    const total = weights.reduce((sum, value) => sum + value, 0);
    let cumulative = 0;
    for (let i = 0; i < axis.length; i += 1) {
      cumulative += weights[i] / total;
      if (cumulative >= probability) return axis[i];
    }
    return axis.at(-1);
  }

  function infer(input, options = {}) {
    const settings = {
      errorFloor: Number(options.errorFloor || 0),
      widthScale: Number(options.widthScale || 1),
      includeSlope: options.includeSlope !== false,
      gridSize: Number(options.gridSize || 51),
      amplitudeMax: Number(options.amplitudeMax || 1400),
      instrument: options.instrument || "all",
    };
    const data = input.filter((point) => settings.instrument === "all" || point.instrument === settings.instrument);
    if (data.length < 10) throw new Error("At least ten spectral bins are required for inference.");
    const axis = Array.from({ length: settings.gridSize }, (_, i) => settings.amplitudeMax * i / (settings.gridSize - 1));
    const fits = [];
    let minimum = Infinity;
    let best = null;
    for (const co2 of axis) {
      for (const h2o of axis) {
        const fit = profile(data, h2o, co2, settings);
        if (!fit) throw new Error("The nuisance design matrix is singular for this selection.");
        fits.push({ h2o, co2, fit });
        if (fit.chi2 < minimum) {
          minimum = fit.chi2;
          best = { h2o, co2, ...fit };
        }
      }
    }
    const posterior = fits.map((entry) => Math.exp(-0.5 * (entry.fit.chi2 - minimum)));
    const total = posterior.reduce((sum, value) => sum + value, 0);
    const normalized = posterior.map((value) => value / total);
    const hMarginal = axis.map((h2o) => fits.reduce((sum, entry, index) => sum + (entry.h2o === h2o ? normalized[index] : 0), 0));
    const cMarginal = axis.map((co2) => fits.reduce((sum, entry, index) => sum + (entry.co2 === co2 ? normalized[index] : 0), 0));
    const nullFit = profile(data, 0, 0, settings);
    const dof = data.length - best.nuisanceCount - 2;
    const residualRms = Math.sqrt(best.residuals.reduce((sum, value) => sum + value ** 2, 0) / data.length);
    const so2Index = best.columns.findIndex((column) => column.type === "so2");
    const offsets = Object.fromEntries(best.columns
      .map((column, index) => [column, best.coefficients[index]])
      .filter(([column]) => column.type === "offset")
      .map(([column, value]) => [column.name, value]));
    const points = data.map((point, index) => ({ ...point, model_ppm: best.predictions[index], residual_ppm: best.residuals[index] }));
    return {
      settings,
      points,
      instruments: [...new Set(data.map((point) => point.instrument))].sort(),
      posterior: { axis, values: normalized, n: axis.length, minimumChi2: minimum },
      intervals: {
        h2o: [quantile(axis, hMarginal, 0.16), quantile(axis, hMarginal, 0.5), quantile(axis, hMarginal, 0.84)],
        co2: [quantile(axis, cMarginal, 0.16), quantile(axis, cMarginal, 0.5), quantile(axis, cMarginal, 0.84)],
      },
      metrics: {
        bins: data.length,
        h2o_amplitude_ppm: best.h2o,
        co2_amplitude_ppm: best.co2,
        so2_nuisance_ppm: best.coefficients[so2Index],
        chi2: minimum,
        dof,
        reduced_chi2: minimum / dof,
        residual_rms_ppm: residualRms,
        delta_bic_vs_no_h2o_co2: nullFit.chi2 - minimum - 2 * Math.log(data.length),
      },
      offsets,
    };
  }

  root.Inference = { infer, profile, templates };
})(typeof self !== "undefined" ? self : globalThis);
