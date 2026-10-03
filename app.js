"use strict";

const COLORS = ["#0b665f", "#a65c24", "#46658a", "#7a536f", "#55733d", "#8d3f32"];
const state = { dataset: null, result: null, worker: null, requestId: 0, timer: null };
const $ = (id) => document.getElementById(id);
const format = (value, digits = 0) => Number(value).toLocaleString("en-GB", { maximumFractionDigits: digits });

function instrumentColours(names) {
  return Object.fromEntries([...names].sort().map((name, index) => [name, COLORS[index % COLORS.length]]));
}

function prepareCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();
  const ratio = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.max(320, Math.round(rect.width));
  const height = Math.max(300, Math.round(rect.height));
  if (canvas.width !== width * ratio || canvas.height !== height * ratio) {
    canvas.width = width * ratio;
    canvas.height = height * ratio;
  }
  const context = canvas.getContext("2d");
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  context.clearRect(0, 0, width, height);
  return { context, width, height };
}

function drawAxes(context, width, height, bounds, labels) {
  const margin = { left: 82, right: 22, top: 24, bottom: 56 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const x = (value) => margin.left + (value - bounds.xMin) / (bounds.xMax - bounds.xMin) * plotWidth;
  const y = (value) => margin.top + (bounds.yMax - value) / (bounds.yMax - bounds.yMin) * plotHeight;
  context.font = "11px IBM Plex Mono, Consolas, monospace";
  context.fillStyle = "#526164";
  context.strokeStyle = "#dfe3df";
  context.lineWidth = 1;
  for (let i = 0; i <= 5; i += 1) {
    const xValue = bounds.xMin + (bounds.xMax - bounds.xMin) * i / 5;
    const yValue = bounds.yMin + (bounds.yMax - bounds.yMin) * i / 5;
    const px = x(xValue);
    const py = y(yValue);
    context.beginPath(); context.moveTo(px, margin.top); context.lineTo(px, height - margin.bottom); context.stroke();
    context.beginPath(); context.moveTo(margin.left, py); context.lineTo(width - margin.right, py); context.stroke();
    context.textAlign = "center";
    context.fillText(format(xValue, 1), px, height - 33);
    context.textAlign = "right";
    context.fillText(format(yValue, 0), margin.left - 10, py + 4);
  }
  context.strokeStyle = "#71807d";
  context.beginPath();
  context.moveTo(margin.left, margin.top); context.lineTo(margin.left, height - margin.bottom); context.lineTo(width - margin.right, height - margin.bottom); context.stroke();
  context.fillStyle = "#182426";
  context.font = "600 12px Inter, Segoe UI, sans-serif";
  context.textAlign = "center";
  context.fillText(labels.x, margin.left + plotWidth / 2, height - 10);
  context.save(); context.translate(16, margin.top + plotHeight / 2); context.rotate(-Math.PI / 2); context.fillText(labels.y, 0, 0); context.restore();
  context.textAlign = "left";
  return { x, y, margin };
}

function drawSpectrum() {
  if (!state.result) return;
  const { context, width, height } = prepareCanvas($("spectrumCanvas"));
  const points = state.result.points;
  const xValues = points.map((point) => point.wavelength_um);
  const yLow = points.map((point) => point.transit_depth_ppm - point.uncertainty_ppm);
  const yHigh = points.map((point) => point.transit_depth_ppm + point.uncertainty_ppm);
  const padding = (Math.max(...yHigh) - Math.min(...yLow)) * 0.06;
  const axes = drawAxes(context, width, height, {
    xMin: Math.min(...xValues), xMax: Math.max(...xValues),
    yMin: Math.min(...yLow) - padding, yMax: Math.max(...yHigh) + padding,
  }, { x: "Wavelength (µm)", y: "Transit depth (ppm)" });
  const colours = instrumentColours(state.result.instruments);
  for (const instrument of state.result.instruments) {
    const selected = points.filter((point) => point.instrument === instrument).sort((a, b) => a.wavelength_um - b.wavelength_um);
    context.strokeStyle = colours[instrument];
    context.globalAlpha = 0.78;
    context.lineWidth = 1.5;
    context.beginPath();
    selected.forEach((point, index) => {
      const px = axes.x(point.wavelength_um), py = axes.y(point.model_ppm);
      if (index) context.lineTo(px, py); else context.moveTo(px, py);
    });
    context.stroke();
    context.globalAlpha = 0.48;
    context.lineWidth = 1;
    for (const point of selected) {
      const px = axes.x(point.wavelength_um), py = axes.y(point.transit_depth_ppm);
      context.beginPath(); context.moveTo(px, axes.y(point.transit_depth_ppm - point.uncertainty_ppm)); context.lineTo(px, axes.y(point.transit_depth_ppm + point.uncertainty_ppm)); context.stroke();
      context.globalAlpha = 0.9;
      context.fillStyle = colours[instrument];
      context.beginPath(); context.arc(px, py, 2.2, 0, Math.PI * 2); context.fill();
      context.globalAlpha = 0.48;
    }
    context.globalAlpha = 1;
  }
  $("spectrumSummary").textContent = `${points.length} spectral bins across ${state.result.instruments.length} products. The profiled model has residual RMS ${format(state.result.metrics.residual_rms_ppm, 1)} ppm and reduced chi-squared ${format(state.result.metrics.reduced_chi2, 2)}.`;
}

function posteriorColour(value, maximum) {
  const t = Math.max(0, Math.min(1, value / maximum));
  const stops = [[247, 246, 240], [163, 203, 193], [24, 111, 105], [22, 45, 53]];
  const scaled = t * (stops.length - 1);
  const index = Math.min(stops.length - 2, Math.floor(scaled));
  const fraction = scaled - index;
  return `rgb(${stops[index].map((channel, i) => Math.round(channel + (stops[index + 1][i] - channel) * fraction)).join(",")})`;
}

function drawPosterior() {
  if (!state.result) return;
  const { context, width, height } = prepareCanvas($("posteriorCanvas"));
  const posterior = state.result.posterior;
  const maxValue = Math.max(...posterior.values);
  const axes = drawAxes(context, width, height, { xMin: 0, xMax: posterior.axis.at(-1), yMin: 0, yMax: posterior.axis.at(-1) }, { x: "H₂O amplitude (ppm)", y: "CO₂ amplitude (ppm)" });
  const plotWidth = width - axes.margin.left - axes.margin.right;
  const plotHeight = height - axes.margin.top - axes.margin.bottom;
  const cellWidth = plotWidth / posterior.n;
  const cellHeight = plotHeight / posterior.n;
  for (let row = 0; row < posterior.n; row += 1) {
    for (let column = 0; column < posterior.n; column += 1) {
      const value = posterior.values[row * posterior.n + column];
      context.fillStyle = posteriorColour(value, maxValue);
      context.fillRect(axes.margin.left + column * cellWidth, axes.margin.top + (posterior.n - row - 1) * cellHeight, cellWidth + .5, cellHeight + .5);
    }
  }
  for (const [ratio, colour] of [[Math.exp(-1.15), "#182426"], [Math.exp(-3), "#a65c24"]]) {
    context.strokeStyle = colour; context.lineWidth = 1.2;
    for (let row = 1; row < posterior.n - 1; row += 1) for (let column = 1; column < posterior.n - 1; column += 1) {
      const here = posterior.values[row * posterior.n + column] / maxValue;
      const right = posterior.values[row * posterior.n + column + 1] / maxValue;
      const up = posterior.values[(row + 1) * posterior.n + column] / maxValue;
      if ((here - ratio) * (right - ratio) <= 0 || (here - ratio) * (up - ratio) <= 0) {
        context.strokeRect(axes.margin.left + column * cellWidth, axes.margin.top + (posterior.n - row - 1) * cellHeight, cellWidth, cellHeight);
      }
    }
  }
  const [hLow, hMid, hHigh] = state.result.intervals.h2o;
  const [cLow, cMid, cHigh] = state.result.intervals.co2;
  $("posteriorSummary").textContent = `Conditional H2O amplitude median ${format(hMid)} ppm with 16th to 84th grid interval ${format(hLow)} to ${format(hHigh)} ppm. CO2 amplitude median ${format(cMid)} ppm with interval ${format(cLow)} to ${format(cHigh)} ppm.`;
}

function renderMetrics() {
  const result = state.result;
  const m = result.metrics;
  const intervals = result.intervals;
  const entries = [
    ["H₂O template amplitude", `${format(intervals.h2o[1])} [${format(intervals.h2o[0])}, ${format(intervals.h2o[2])}] ppm`],
    ["CO₂ template amplitude", `${format(intervals.co2[1])} [${format(intervals.co2[0])}, ${format(intervals.co2[2])}] ppm`],
    ["SO₂ nuisance amplitude", `${format(m.so2_nuisance_ppm, 1)} ppm`],
    ["Reduced χ²", `${format(m.reduced_chi2, 2)} (${format(m.chi2, 1)} / ${m.dof})`],
    ["Residual RMS", `${format(m.residual_rms_ppm, 1)} ppm`],
    ["ΔBIC vs no H₂O/CO₂", format(m.delta_bic_vs_no_h2o_co2, 1)],
    ["Included bins", format(m.bins)],
    ["Profiled offsets", format(Object.keys(result.offsets).length)],
  ];
  $("metrics").innerHTML = entries.map(([label, value]) => `<div class="metric"><span>${label}</span><strong>${value}</strong></div>`).join("");
  const poor = m.reduced_chi2 > 1.5 || m.reduced_chi2 < 0.65;
  $("fitFlag").className = `fit-flag ${poor ? "bad" : "good"}`;
  $("fitFlag").textContent = poor ? "model inadequate" : "variance compatible";
  $("interpretationText").textContent = poor
    ? `The reduced χ² of ${format(m.reduced_chi2, 2)} rejects this template family as an adequate generative description under the selected errors. The narrow grid intervals and ΔBIC therefore describe only conditional feature contrast; they must not be reported as atmospheric constraints or molecular significances.`
    : `With the selected ${format(result.settings.errorFloor)} ppm error floor, residual variance is broadly compatible with the declared likelihood. This is a sensitivity result: the amplitudes remain phenomenological and do not map to mixing ratios.`;
}

function renderLegend() {
  const colours = instrumentColours(state.result.instruments);
  $("legend").innerHTML = state.result.instruments.map((name) => `<span><i style="--swatch:${colours[name]}"></i>${name.replaceAll("_", " ")}</span>`).join("");
}

function renderCoverage() {
  const groups = state.dataset.points.reduce((accumulator, point) => {
    (accumulator[point.instrument] ||= []).push(point);
    return accumulator;
  }, {});
  $("coverageRows").innerHTML = Object.entries(groups).sort(([a], [b]) => a.localeCompare(b)).map(([name, points]) => {
    const errors = points.map((point) => point.uncertainty_ppm).sort((a, b) => a - b);
    const median = errors[Math.floor(errors.length / 2)];
    const wavelengths = points.map((point) => point.wavelength_um);
    return `<tr><th scope="row">${name.replaceAll("_", " ")}</th><td>${points.length}</td><td>${format(Math.min(...wavelengths), 2)}–${format(Math.max(...wavelengths), 2)}</td><td>${format(median, 1)}</td></tr>`;
  }).join("");
  $("provenanceNote").textContent = `Archive MD5 ${state.dataset.archive.md5}. Six source CSV files are pinned by SHA-256; ${state.dataset.excluded_saturated_prism_bins} flagged PRISM bins are excluded. ${state.dataset.covariance}`;
}

function render() {
  renderLegend(); renderMetrics(); drawSpectrum(); drawPosterior();
}

function requestAnalysis() {
  if (!state.dataset) return;
  clearTimeout(state.timer);
  state.timer = setTimeout(() => {
    state.requestId += 1;
    $("workerStatus").textContent = "Evaluating likelihood grid";
    $("runId").textContent = `analysis ${String(state.requestId).padStart(2, "0")}`;
    state.worker.postMessage({
      requestId: state.requestId,
      spectrum: state.dataset.points,
      options: {
        instrument: $("instrument").value,
        errorFloor: Number($("errorFloor").value),
        widthScale: Number($("widthScale").value),
        includeSlope: $("includeSlope").checked,
        gridSize: 51,
        amplitudeMax: 1400,
      },
    });
  }, 90);
}

async function initialise() {
  const response = await fetch("data/wasp39b_benchmark.json", { cache: "no-cache" });
  if (!response.ok) throw new Error(`Dataset request failed (${response.status})`);
  state.dataset = await response.json();
  const instruments = [...new Set(state.dataset.points.map((point) => point.instrument))].sort();
  for (const instrument of instruments) {
    const option = document.createElement("option"); option.value = instrument; option.textContent = instrument.replaceAll("_", " "); $("instrument").appendChild(option);
  }
  renderCoverage();
  state.worker = new Worker("physicsWorker.js");
  state.worker.onmessage = (event) => {
    if (event.data.requestId !== state.requestId) return;
    if (!event.data.ok) { $("workerStatus").textContent = event.data.error; return; }
    state.result = event.data.result;
    $("workerStatus").textContent = `${state.result.metrics.bins} bins · likelihood complete`;
    render();
  };
  state.worker.onerror = (event) => { $("workerStatus").textContent = `Worker error: ${event.message}`; };
  requestAnalysis();
}

$("errorFloor").addEventListener("input", () => { $("errorFloorValue").textContent = `${$("errorFloor").value} ppm`; requestAnalysis(); });
$("widthScale").addEventListener("input", () => { $("widthScaleValue").textContent = `${Number($("widthScale").value).toFixed(2)}×`; requestAnalysis(); });
$("instrument").addEventListener("change", requestAnalysis);
$("includeSlope").addEventListener("change", requestAnalysis);
$("reset").addEventListener("click", () => {
  $("instrument").value = "all"; $("errorFloor").value = 0; $("errorFloorValue").textContent = "0 ppm";
  $("widthScale").value = 1; $("widthScaleValue").textContent = "1.00×"; $("includeSlope").checked = true; requestAnalysis();
});
window.addEventListener("resize", () => { if (state.result) { drawSpectrum(); drawPosterior(); } });

initialise().catch((error) => { $("workerStatus").textContent = error.message; });
