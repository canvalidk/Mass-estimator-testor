/* Offline results browser. Numerical work stays in the shared Python engine. */
"use strict";

const $ = (id) => document.getElementById(id);
const state = {runs: [], run: null, metric: "squared_log", method: "norm_ratio", path: "", axis: "a", row: null, token: null, job: null, loadSequence: 0, explicitMeasurement: false};
const colors = {ours: "#087f79", other: "#b47b16", ink: "#193238", muted: "#718083", grid: "#e5ebe6"};
const numberFormat = new Intl.NumberFormat("en", {maximumFractionDigits: 3});
const wholeFormat = new Intl.NumberFormat("en", {maximumFractionDigits: 0});
let pollTimer = null;
let toastTimer = null;

function element(tag, attrs = {}, children = []) {
  const node = document.createElement(tag);
  Object.entries(attrs).forEach(([key, value]) => {
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else if (value !== null && value !== undefined) node.setAttribute(key, value);
  });
  (Array.isArray(children) ? children : [children]).forEach(child => {
    if (child !== null && child !== undefined) node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  });
  return node;
}

function svgElement(tag, attrs = {}, content) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
  if (content !== undefined) node.textContent = content;
  return node;
}

function finite(value) { return typeof value === "number" && Number.isFinite(value); }
function format(value, significant = 5) {
  if (!finite(value)) return "Unavailable";
  if (value === 0) return "0";
  if (Math.abs(value) < 0.0001 || Math.abs(value) >= 1e6) return value.toExponential(significant - 1);
  return Number(value.toPrecision(significant)).toLocaleString("en", {maximumFractionDigits: 12});
}
function signed(value) { return finite(value) ? `${value > 0 ? "+" : ""}${format(value)}` : "Unavailable"; }
function errorText(error) { return error instanceof Error ? error.message : String(error); }
function metricInfo() { return (state.run?.metrics || []).find(metric => metric.id === state.metric) || {id: state.metric, label: state.metric, better: "lower"}; }
function methodInfo(id) { return (state.run?.methods || []).find(method => method.id === id) || {id, label: id.replaceAll("_", " ")}; }
function distributionMetric() { return ["log_crps", "log_density_score"].includes(state.metric) || /^(interval_score_|coverage_|log_width_|below_interval_|above_interval_)/.test(state.metric); }
function methodLabel(id) { return distributionMetric() ? ({flat_joint: "Flat posterior law (supplied noise)", tube_joint: "Tube posterior law", calibrated_joint: "Calibration-integrated law", oracle_joint: "Known-noise oracle law"}[id] || methodInfo(id).label) : methodInfo(id).label; }
function diagonalRun() { return state.run?.config?.noise_model === "diagonal_gaussian_3d" || state.run?.metadata?.noise_model === "diagonal_gaussian_3d"; }
function measurementInfo() { return state.run?.config?.measurement || state.run?.metadata?.measurement || {repeats: 1, calibration: {mode: "known"}}; }
function repeatCount() { return measurementInfo().repeats || 1; }
function calibrationLabel(calibration = measurementInfo().calibration) { return {known: "Exact noise SDs", scaled: "Miscalibrated noise SDs", estimated: `Coordinate noise SDs estimated from ${calibration?.samples ?? "—"} vectors`, pooled_isotropic: `Pooled isotropic calibration from ${calibration?.samples ?? "—"} vectors`}[calibration?.mode || "known"]; }
function snrName() { return `${repeatCount() > 1 ? "per-reading " : ""}${diagonalRun() ? "RMS-based SNR" : "SNR"}`; }
function massName() { return diagonalRun() ? "True physical mass" : "True mass"; }
function valueFor(row, method, metric = state.metric) { return row?.values?.[method]?.[metric] || null; }
function pairFor(row, method, metric = state.metric) { return row?.paired?.[method]?.[metric] || null; }
function sourceText(value) { return typeof value === "string" ? value : JSON.stringify(value); }

async function api(url, options = {}) {
  const response = await fetch(url, options);
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(typeof body.error === "string" ? body.error : body.message || `Request failed (${response.status}).`);
  return body;
}

function showError(message) { $("global-error").textContent = message; $("global-error").hidden = !message; }
function toast(message) {
  $("toast").textContent = message;
  $("toast").hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { $("toast").hidden = true; }, 4500);
}

function download(name, text, type = "application/json") {
  const url = URL.createObjectURL(new Blob([text], {type}));
  const link = element("a", {href: url, download: name});
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function fillSelect(node, items, selected) {
  node.replaceChildren(...items.map(item => element("option", {value: item.id, text: item.label})));
  if (items.some(item => item.id === selected)) node.value = selected;
}

async function refreshRuns(selectId = null) {
  const result = await api("/api/runs");
  state.runs = result.runs || [];
  renderRunList();
  if (Object.hasOwn(result, "active")) updateJob(result.active);
  if (selectId) await loadRun(selectId);
  else if (!state.run && state.runs.length) await loadRun(state.runs[0].id);
  $("empty-state").hidden = state.runs.length > 0;
  $("results-content").hidden = !state.run;
}

function renderRunList() {
  const terms = $("run-search").value.trim().toLowerCase().split(/\s+/).filter(Boolean);
  const matches = state.runs.filter(run => terms.every(term => `${run.title} ${run.kind || ""}`.toLowerCase().includes(term)));
  $("run-count").textContent = `${matches.length} of ${state.runs.length} experiments${state.run && !matches.some(run => run.id === state.run.id) ? " · selected result stays open" : ""}`;
  $("run-list").replaceChildren(...matches.map(run => {
    const button = element("button", {class: `run-item${state.run?.id === run.id ? " active" : ""}`, "aria-current": state.run?.id === run.id ? "true" : "false", onclick: () => loadRun(run.id).catch(error => showError(errorText(error)))}, [
      element("span", {class: "run-type", text: String(run.kind || "experiment").replaceAll("_", " ")}),
      element("strong", {text: run.title}),
      element("small", {text: `${finite(run.samples) ? wholeFormat.format(run.samples) + " trials / scenario" : "Reference calculation"} · ${run.scenario_count ?? "—"} scenarios`})
    ]);
    return button;
  }));
  if (!matches.length) $("run-list").append(element("p", {class: "muted", text: state.runs.length ? "No experiments match this search." : "No saved runs yet."}));
}

async function loadRun(id) {
  const sequence = ++state.loadSequence;
  $("results-content").classList.add("loading");
  showError("");
  try {
    const run = await api(`/api/runs/${encodeURIComponent(id)}`);
    if (sequence !== state.loadSequence) return;
    state.run = run;
    run.rows ||= [];
    run.methods ||= [];
    run.metrics ||= [];
    state.row = null;
    state.path = "";
    renderRun();
    $("results-content").hidden = false;
    $("empty-state").hidden = true;
  } finally {
    if (sequence === state.loadSequence) $("results-content").classList.remove("loading");
  }
}

function renderRun() {
  const run = state.run;
  renderRunList();
  $("run-title").textContent = run.title;
  $("run-kind").textContent = String(run.kind || "Saved experiment").replaceAll("_", " ");
  const meta = run.metadata || {};
  $("run-subtitle").textContent = diagonalRun()
    ? "Noise specified per coordinate · SNRs use each channel’s RMS standard deviation; mass is physical."
    : meta.antithetic
    ? "Antithetic trials · uncertainty is computed across independent trial pairs."
    : "Compare accuracy against known true mass across controlled signal conditions.";
  if (run.config?.measurement || meta.measurement) $("run-subtitle").textContent += ` ${repeatCount()} ${repeatCount() === 1 ? "reading" : "readings"} per trial at one fixed excitation. ${calibrationLabel()}.`;
  const facts = [
    ["Trials / scenario", finite(run.samples) ? wholeFormat.format(run.samples) : "Reference"],
    ["Scenarios", run.scenario_count ?? run.rows.length],
    ["Estimators", run.methods.length],
    ["Available scores", run.metrics.length]
  ];
  if (repeatCount() > 1) facts.push(["Readings / trial", repeatCount()]);
  $("run-facts").replaceChildren(...facts.map(([label, value]) => element("div", {class: "run-fact"}, [element("span", {text: label}), element("strong", {text: String(value)})])));
  if (!run.metrics.some(metric => metric.id === state.metric)) state.metric = run.metrics.find(metric => metric.id === "squared_log")?.id || run.metrics[0]?.id || "squared_log";
  const others = run.methods.filter(method => method.id !== "flat_joint");
  if (!others.some(method => method.id === state.method)) state.method = others.find(method => method.id === "norm_ratio")?.id || others[0]?.id || "norm_ratio";
  fillSelect($("metric-select"), run.metrics, state.metric);
  fillSelect($("method-select"), others, state.method);
  fillSelect($("axis-select"), ["a", "strength", "f", "mass"].map(id => ({id, label: axisLabel(id)})), state.axis);
  const paths = getPaths(run.rows);
  state.path = paths.find(path => path.id === "force:8")?.id || paths.find(path => path.id === "equal")?.id || paths.find(path => path.id !== "all")?.id || "all";
  fillSelect($("path-select"), paths, state.path);
  $("reuse-config").hidden = !run.config;
  renderMethods();
  renderProvenance();
  renderComparison();
}

function keyNumber(number) { return String(Number(number.toPrecision(10))); }
function closeNumber(left, right) { return Math.abs(left - right) <= 1e-8 * Math.max(Number.MIN_VALUE, Math.abs(left), Math.abs(right)); }
function getPaths(rows) {
  const forceGroups = new Map();
  const massGroups = new Map();
  rows.forEach(row => {
    if (finite(row.f)) { const key = keyNumber(row.f); forceGroups.set(key, (forceGroups.get(key) || 0) + 1); }
    if (finite(row.mass)) { const key = keyNumber(row.mass); massGroups.set(key, (massGroups.get(key) || 0) + 1); }
  });
  const paths = [{id: "all", label: "All scenarios · separate points"}];
  [...forceGroups].filter(([, count]) => count > 1).sort((a, b) => +a[0] - +b[0]).forEach(([number, count]) => paths.push({id: `force:${number}`, label: `Force ${snrName()} = ${format(+number)} · ${count} points`}));
  if (rows.filter(row => finite(row.f) && finite(row.a) && closeNumber(row.f, row.a)).length > 1) paths.push({id: "equal", label: `Equal force & acceleration ${snrName()}`});
  [...massGroups].filter(([, count]) => count > 1).sort((a, b) => +a[0] - +b[0]).forEach(([number, count]) => paths.push({id: `mass:${number}`, label: `${massName()} = ${format(+number)} · ${count} points`}));
  return paths;
}

function selectedRows() {
  let rows = [...(state.run?.rows || [])];
  if (state.path === "equal") rows = rows.filter(row => finite(row.f) && finite(row.a) && closeNumber(row.f, row.a));
  else if (state.path.startsWith("force:")) rows = rows.filter(row => finite(row.f) && closeNumber(row.f, Number(state.path.slice(6))));
  else if (state.path.startsWith("mass:")) rows = rows.filter(row => finite(row.mass) && closeNumber(row.mass, Number(state.path.slice(5))));
  return rows.sort((left, right) => (left[state.axis] ?? 0) - (right[state.axis] ?? 0));
}

function cellConcern(value) {
  if (!value) return "Individual score not saved";
  if (value.status === "not_applicable") return value.note ? `Not applicable: ${value.note}` : "Not applicable to this scenario";
  const text = `${value.status || ""} ${value.note || ""}`;
  if (/infinite.*population|population.*infinite|infinite_under_declared/i.test(text)) return "Infinite population risk; a finite sample mean does not estimate a finite risk";
  if (/infinite.*variance|no.*finite.variance|finite.variance.*not|mcse.*descriptive/i.test(text)) return "Empirical MCSE has no usual finite-variance justification";
  if (!finite(value.mean)) return /inf|unbounded|invalid/i.test(text) ? "Unbounded or undefined under this score’s failure policy" : "Unavailable in this run";
  return "";
}

function validWhisker(value) {
  if (!value || !finite(value.mean) || !finite(value.mcse) || value.mcse < 0) return false;
  const text = `${value.status || ""} ${value.note || ""}`;
  return !/infinite.*variance|no.*finite.variance|finite.variance.*not|infinite.*population|population.*infinite|infinite_under_declared|mcse.*descriptive/i.test(text);
}

function renderComparison() {
  if (!state.run) return;
  const metric = metricInfo();
  const rows = selectedRows();
  const ourLabel = distributionMetric() ? "Our flat law" : "Our equation";
  const otherLabel = methodLabel(state.method);
  fillSelect($("method-select"), state.run.methods.filter(method => method.id !== "flat_joint").map(method => ({id: method.id, label: methodLabel(method.id)})), state.method);
  $("metric-description").textContent = metric.description || "Expected score is estimated by averaging across simulated trials.";
  $("risk-heading").textContent = metric.label;
  $("risk-subheading").textContent = metric.better === "target" ? `Compare with the nominal target ${finite(metric.target) ? format(metric.target) : "shown"}; a higher value is not automatically better.` : metric.better === "context" ? "A diagnostic of estimator behavior." : "Smaller average loss indicates greater accuracy under this score.";
  $("paired-chart").previousElementSibling.querySelector("p").textContent = metric.better === "target"
    ? "Our equation − comparator · assess each against the nominal target"
    : metric.better === "context" ? "Our equation − comparator · signed diagnostic difference" : "Our equation − comparator · below zero favors ours";
  $("risk-legend").replaceChildren(...[[ourLabel, "ours"], [otherLabel, "other"]].map(([label, color]) => element("span", {}, [element("i", {class: `legend-${color}`}), label])));
  const warnings = new Set();
  rows.forEach(row => [valueFor(row, "flat_joint"), valueFor(row, state.method), pairFor(row, state.method)].forEach(value => { if (value && cellConcern(value)) warnings.add(cellConcern(value)); }));
  if (state.path === "all") warnings.add("All scenarios are shown as separate points. Choose a signal path to compare a continuous sequence.");
  $("score-warning").hidden = !warnings.size;
  $("score-warning").textContent = [...warnings].join(". ") + (warnings.size ? "." : "");
  if (!rows.some(row => row.id === state.row)) state.row = rows[Math.floor(rows.length / 2)]?.id || null;
  fillSelect($("scenario-select"), rows.map(row => ({id: row.id, label: `F ${format(row.f, 4)} · a ${format(row.a, 4)} · mass ${format(row.mass, 4)}`})), state.row);
  drawChart($("risk-chart"), rows, [
    {label: ourLabel, color: colors.ours, accessor: row => valueFor(row, "flat_joint")},
    {label: otherLabel, color: colors.other, accessor: row => valueFor(row, state.method)}
  ], {target: metric.better === "target" ? metric.target : null, empty: "Individual scores were not saved for this score.", subtext: "The paired difference may still be available in the adjacent chart."});
  drawChart($("paired-chart"), rows, [{label: "Our − comparator", color: colors.ours, accessor: row => pairFor(row, state.method)}], {paired: true, empty: "No paired difference is available.", subtext: "Paired uncertainty needs the score differences on the same trials; it cannot be reconstructed by subtracting separate standard errors."});
  renderScenario();
}

function drawChart(host, rows, series, options = {}) {
  host.replaceChildren();
  const points = series.flatMap(item => rows.map(row => ({row, value: item.accessor(row), series: item})).filter(point => finite(point.value?.mean) && finite(point.row[state.axis])));
  if (!points.length) {
    host.append(element("div", {class: "chart-empty"}, [element("strong", {text: options.empty}), element("span", {text: options.subtext || ""})]));
    return;
  }
  const width = 580, height = 330;
  const margin = {left: 68, right: 20, top: 19, bottom: 54};
  const xValues = rows.map(row => row[state.axis]).filter(finite);
  let xMin = Math.min(...xValues), xMax = Math.max(...xValues);
  if (xMin === xMax) { xMin -= Math.max(0.1, Math.abs(xMin) * 0.1); xMax += Math.max(0.1, Math.abs(xMax) * 0.1); }
  const yValues = points.flatMap(point => validWhisker(point.value) ? [point.value.mean - 1.96 * point.value.mcse, point.value.mean + 1.96 * point.value.mcse] : [point.value.mean]);
  if (options.paired) yValues.push(0);
  if (finite(options.target)) yValues.push(options.target);
  let yMin = Math.min(...yValues), yMax = Math.max(...yValues);
  if (yMin === yMax) { yMin -= Math.max(0.01, Math.abs(yMin) * 0.1); yMax += Math.max(0.01, Math.abs(yMax) * 0.1); }
  const pad = (yMax - yMin) * 0.12;
  yMin -= pad; yMax += pad;
  if (!options.paired && Math.min(...yValues) >= 0 && yMin < 0) yMin = 0;
  const sx = x => margin.left + (x - xMin) / (xMax - xMin) * (width - margin.left - margin.right);
  const sy = y => height - margin.bottom - (y - yMin) / (yMax - yMin) * (height - margin.top - margin.bottom);
  const svg = svgElement("svg", {viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": `${metricInfo().label}: ${options.paired ? "paired difference" : "individual scores"} by ${axisLabel()}`});
  svg.append(svgElement("title", {}, `${metricInfo().label}. ${points.length} recorded values. Select a point to inspect its scenario.`));
  for (let index = 0; index <= 4; index++) {
    const y = yMin + (yMax - yMin) * index / 4;
    const lineY = sy(y);
    svg.append(svgElement("line", {x1: margin.left, y1: lineY, x2: width - margin.right, y2: lineY, stroke: colors.grid, "stroke-width": 1}));
    svg.append(svgElement("text", {x: margin.left - 10, y: lineY + 3.5, fill: colors.muted, "font-size": 10, "text-anchor": "end"}, format(y, 3)));
  }
  for (let index = 0; index <= 4; index++) {
    const x = xMin + (xMax - xMin) * index / 4;
    svg.append(svgElement("text", {x: sx(x), y: height - margin.bottom + 21, fill: colors.muted, "font-size": 10, "text-anchor": "middle"}, format(x, 3)));
  }
  const reference = options.paired ? 0 : options.target;
  if (finite(reference)) {
    svg.append(svgElement("line", {x1: margin.left, y1: sy(reference), x2: width - margin.right, y2: sy(reference), stroke: colors.ink, "stroke-width": 1.1, "stroke-dasharray": "4 4", opacity: .6}));
    if (!options.paired) svg.append(svgElement("text", {x: width - margin.right - 2, y: sy(reference) - 5, fill: colors.muted, "font-size": 9, "text-anchor": "end"}, `Nominal ${format(reference)}`));
  }
  svg.append(svgElement("text", {x: (margin.left + width - margin.right) / 2, y: height - 7, fill: colors.muted, "font-size": 10, "text-anchor": "middle"}, axisLabel()));
  const tooltip = element("div", {class: "chart-tooltip", role: "tooltip", hidden: ""});
  function showTooltip(point, target) {
    const value = point.value;
    tooltip.textContent = `${point.series.label}\nForce ${snrName()} ${format(point.row.f)} · acceleration ${snrName()} ${format(point.row.a)}\n${massName()} ${format(point.row.mass)}\n${options.paired ? "Difference" : "Score"}: ${format(value.mean, 7)}${finite(value.mcse) ? `\nEmpirical MCSE: ${format(value.mcse, 5)}` : ""}${cellConcern(value) ? `\n${cellConcern(value)}` : ""}`;
    tooltip.hidden = false;
    const box = target.getBoundingClientRect(), parent = host.getBoundingClientRect();
    tooltip.style.left = `${Math.max(0, Math.min(box.left - parent.left + 12, parent.width - Math.min(245, tooltip.offsetWidth) - 4))}px`;
    tooltip.style.top = `${Math.max(0, box.top - parent.top - tooltip.offsetHeight - 12)}px`;
  }
  series.forEach(item => {
    // Missing values break a path; unrelated scenarios are displayed without a path.
    if (state.path !== "all") {
      let segment = [];
      const flush = () => {
        if (segment.length > 1) svg.append(svgElement("polyline", {points: segment.join(" "), fill: "none", stroke: item.color, "stroke-width": 1.8, opacity: .85, "stroke-linejoin": "round"}));
        segment = [];
      };
      rows.forEach(row => {
        const value = item.accessor(row);
        if (finite(value?.mean) && finite(row[state.axis])) segment.push(`${sx(row[state.axis])},${sy(value.mean)}`);
        else flush();
      });
      flush();
    }
    rows.forEach(row => {
      const value = item.accessor(row);
      if (!finite(value?.mean) || !finite(row[state.axis])) return;
      const x = sx(row[state.axis]), y = sy(value.mean);
      if (validWhisker(value)) {
        const upper = sy(value.mean + 1.96 * value.mcse), lower = sy(value.mean - 1.96 * value.mcse);
        svg.append(svgElement("path", {d: `M ${x} ${upper} V ${lower} M ${x - 3} ${upper} H ${x + 3} M ${x - 3} ${lower} H ${x + 3}`, fill: "none", stroke: item.color, "stroke-width": 1, opacity: .45}));
      }
      if (row.id === state.row) svg.append(svgElement("circle", {cx: x, cy: y, r: 8, fill: item.color, opacity: .12}));
      const point = svgElement("circle", {cx: x, cy: y, r: row.id === state.row ? 4.5 : 3.3, fill: item.color, stroke: "white", "stroke-width": 1.2, class: "chart-point", tabindex: "0", role: "button", "aria-label": `${item.label}, force ${format(row.f)}, acceleration ${format(row.a)}, score ${format(value.mean)}. Inspect scenario.`});
      point.append(svgElement("title", {}, `${item.label}: ${format(value.mean, 7)}`));
      const select = () => { state.row = row.id; renderComparison(); };
      point.addEventListener("click", select);
      point.addEventListener("keydown", event => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(); } });
      point.addEventListener("mouseenter", () => showTooltip({row, value, series: item}, point));
      point.addEventListener("focus", () => showTooltip({row, value, series: item}, point));
      point.addEventListener("mouseleave", () => { tooltip.hidden = true; });
      point.addEventListener("blur", () => { tooltip.hidden = true; });
      svg.append(point);
    });
  });
  host.append(svg, tooltip);
}

function axisLabel(axis = state.axis) {
  const label = diagonalRun() ? {a: "Acceleration / RMS noise SD", f: "Force / RMS noise SD", strength: "Combined RMS-based signal / noise", mass: "True physical mass"}[axis] : {a: "Acceleration signal / noise", f: "Force signal / noise", strength: "Total signal / noise", mass: "True mass ratio"}[axis];
  return axis !== "mass" && repeatCount() > 1 ? `${label} · per reading` : label;
}

function renderScenario() {
  const row = state.run.rows.find(item => item.id === state.row);
  $("table-score-label").textContent = metricInfo().label;
  $("scenario-table").replaceChildren();
  $("scenario-notes").replaceChildren();
  if (!row) {
    $("scenario-title").textContent = "No scenario selected";
    $("scenario-detail").textContent = "Choose a run with recorded trial results.";
    return;
  }
  $("scenario-title").textContent = `${massName()} ${format(row.mass)} · force ${snrName()} ${format(row.f)} · acceleration ${snrName()} ${format(row.a)}`;
  const magnitudes = diagonalRun() && finite(row.true_force_magnitude) && finite(row.true_acceleration_magnitude) ? `True magnitudes: force ${format(row.true_force_magnitude)}, acceleration ${format(row.true_acceleration_magnitude)}. ` : "";
  const effective = repeatCount() > 1 && finite(row.effective_force_snr) && finite(row.effective_acceleration_snr) ? `${repeatCount()} readings combined: mean force SNR ${format(row.effective_force_snr)}, mean acceleration SNR ${format(row.effective_acceleration_snr)} (using the true noise). ` : "";
  $("scenario-detail").textContent = `${diagonalRun() ? "Combined RMS-based" : "Total"} signal / noise ${format(row.strength)}${repeatCount() > 1 ? " per reading" : ""}. ${magnitudes}${effective}${metricInfo().better === "target" ? "Coverage is assessed against its nominal target." : "Values are recorded sample means, with unavailable or unbounded results retained."}`;
  const methods = [...state.run.methods].sort((left, right) => (left.id === "flat_joint" ? -1 : left.id === state.method ? 0 : 1) - (right.id === "flat_joint" ? -1 : right.id === state.method ? 0 : 1));
  const notes = new Set();
  methods.forEach(method => {
    const value = valueFor(row, method.id), paired = pairFor(row, method.id), diagnostic = row.diagnostics?.[method.id] || {};
    const concern = cellConcern(value);
    const status = [];
    if (finite(diagnostic.invalid_count)) status.push(`${wholeFormat.format(diagnostic.invalid_count)} invalid outputs`);
    else if (finite(diagnostic.finite_positive_fraction)) status.push(`${format(100 * diagnostic.finite_positive_fraction, 5)}% finite positive outputs`);
    if (concern) status.push(concern);
    else if (!status.length) status.push(value ? "Finite recorded score" : "Not recorded");
    if (value?.note) notes.add(`${method.label}: ${value.note}`);
    const first = element("td", {}, [methodLabel(method.id), method.id === "flat_joint" ? element("span", {class: "our-mark", text: "OURS"}) : null]);
    const score = element("td", {text: value && finite(value.mean) ? format(value.mean) : "—"});
    if (value && !finite(value.mean)) score.append(element("span", {class: "cell-subtext", text: concern || value.status || "Unavailable"}));
    const mcse = element("td", {text: value && finite(value.mcse) ? format(value.mcse, 3) : "—"});
    if (value && finite(value.mcse) && !validWhisker(value)) mcse.append(element("span", {class: "cell-subtext", text: "Descriptive only; whiskers omitted"}));
    const delta = element("td", {text: method.id === "flat_joint" ? "Reference" : finite(paired?.mean) ? signed(paired.mean) : "—", class: metricInfo().better === "lower" && finite(paired?.mean) ? (paired.mean < 0 ? "delta-negative" : "delta-positive") : ""});
    if (finite(paired?.mcse)) delta.append(element("span", {class: "cell-subtext", text: `MCSE ${format(paired.mcse, 3)}`}));
    if (paired?.note) notes.add(`Paired comparison with ${method.label}: ${paired.note}`);
    $("scenario-table").append(element("tr", {class: ["flat_joint", state.method].includes(method.id) ? "selected-method" : ""}, [first, score, mcse, delta, element("td", {class: `status-note${concern ? " caution" : ""}`, text: status.join(" · ")})]));
  });
  if (state.metric === "squared_log") {
    const diagnostics = ["flat_joint", state.method].map(id => ({id, data: row.diagnostics?.[id]})).filter(item => finite(item.data?.bias) && finite(item.data?.variance));
    if (diagnostics.length) $("scenario-notes").append(element("p", {}, [element("strong", {text: "Bias and spread: "}), ...diagnostics.map(({id, data}, index) => `${index ? " · " : ""}${id === "flat_joint" ? "Our equation" : methodInfo(id).label}: signed log bias ${format(data.bias, 4)}, log variance ${format(data.variance, 4)}`), ". Squared log error combines squared log bias with log variance."]));
  }
  if (notes.size) $("scenario-notes").append(element("details", {class: "interpretation-notes"}, [element("summary", {text: `Interpretation notes · ${notes.size}`}), ...[...notes].map(note => element("p", {text: note}))]));
}

function renderMethods() {
  $("method-cards").replaceChildren(...state.run.methods.map(method => {
    const card = element("article", {class: "method-card"}, [element("h3", {text: method.label})]);
    if (method.formula) card.append(element("p", {class: "formula", text: method.formula}));
    if (method.assumptions) card.append(element("p", {text: sourceText(method.assumptions)}));
    if (method.note || method.risk_note) card.append(element("p", {text: method.note || method.risk_note}));
    if (typeof method.source === "string" && /^https?:\/\//i.test(method.source)) card.append(element("a", {href: method.source, target: "_blank", rel: "noopener noreferrer", text: "Source paper ↗"}));
    else if (method.source) card.append(element("p", {text: sourceText(method.source)}));
    return card;
  }));
}

function renderProvenance() {
  const run = state.run, meta = run.metadata || {};
  const measurement = measurementInfo(), calibration = measurement.calibration || {mode: "known"};
  const list = element("dl", {class: "provenance-list"});
  const entries = [
    ["Source", run.source], ["Training seed", meta.training_seed ?? meta.seed ?? run.config?.training_seed],
    ["Held-out seed", meta.heldout_seed ?? run.config?.heldout_seed], ["Quadrature order", meta.quadrature_order ?? run.config?.order],
    ["Managed context", meta.managed_commit], ["Noise model", run.config?.noise_model ?? "See source assumptions"],
    ["True force noise SD / reading", run.config?.noise?.force_sd ?? meta.measurement?.true_force_sd], ["True acceleration noise SD / reading", run.config?.noise?.acceleration_sd ?? meta.measurement?.true_acceleration_sd],
    ["Truth direction", run.config?.direction ? `${JSON.stringify(run.config.direction)} · used to generate trials; unknown to the fitted estimators` : null],
    ["Readings / trial", repeatCount()], ["Noise information", calibrationLabel(calibration)],
    ["Force SD factors", calibration.force_scale], ["Acceleration SD factors", calibration.acceleration_scale],
    ["Calibration seed", calibration.seed], ["Calibration samples / coordinate", calibration.samples],
    ["True mean force noise SD", meta.measurement?.true_mean_force_sd], ["True mean acceleration noise SD", meta.measurement?.true_mean_acceleration_sd],
    ["Supplied mean force noise SD", meta.measurement?.supplied_mean_force_sd_summary], ["Supplied mean acceleration noise SD", meta.measurement?.supplied_mean_acceleration_sd_summary],
    ["Antithetic trials", typeof meta.antithetic === "boolean" ? (meta.antithetic ? "Yes; independent pairs are the uncertainty units" : "No") : null],
    ["Duration", finite(meta.elapsed_seconds) ? `${format(meta.elapsed_seconds, 4)} seconds` : finite(meta.seconds) ? `${format(meta.seconds, 4)} seconds` : null]
  ];
  entries.filter(([, value]) => value !== undefined && value !== null).forEach(([label, value]) => list.append(element("dt", {text: label}), element("dd", {text: sourceText(value)})));
  const content = [list];
  if (repeatCount() > 1) content.push(element("p", {text: "Each trial combines independent readings of the same latent force and acceleration. Charts use per-reading signal / noise; the scenario details also show the true signal / noise of the combined mean. Excitation does not vary within a trial."}));
  if (calibration.mode === "estimated") content.push(element("p", {text: "Noise calibration is independent of the mass readings and is repeated for every trial. The fitted methods treat each estimated noise SD as fixed, so their mass intervals omit calibration uncertainty. Calibration draws are shared across scenario cells; cells are not independent replications."}));
  if (calibration.mode === "pooled_isotropic") content.push(element("p", {text: "Each trial has fresh independent calibration, pooled over three coordinates within each channel. Compare the supplied-noise plug-in law, the law integrating calibration uncertainty, and the oracle given the true noise SDs. Inference uses reading means and external calibration; it does not use scatter among repeated readings. Calibration draws are reused across truth cells, which are not independent replications."}));
  ["assumptions", "failure_policy", "paired_difference_convention", "uncertainty_scope", "legacy_scope", "config_note"].forEach(key => { if (meta[key]) content.push(element("p", {text: sourceText(meta[key])})); });
  if (run.config) {
    content.push(element("div", {class: "detail-actions"}, [element("button", {class: "text-button", text: "Download configuration ↓", onclick: () => download(`${safeName(run.title)}-config.json`, JSON.stringify(run.config, null, 2))}), element("button", {class: "text-button", text: "Reuse these settings", onclick: () => { setConfig(run.config); openSetup(); }})]));
  } else content.push(element("p", {text: "This historical run predates the common configuration format. Its original results remain available for download."}));
  const details = element("details", {}, [element("summary", {text: "All saved metadata"}), element("pre", {text: JSON.stringify(meta, null, 2)})]);
  details.className = "metadata-details";
  content.push(details);
  $("provenance").replaceChildren(...content);
}

function safeName(text) { return String(text || "experiment").replace(/[^a-zA-Z0-9._-]+/g, "-").slice(0, 80).replace(/^-+|-+$/g, "") || "experiment"; }
function exportCSV() {
  if (!state.run) return;
  const header = ["run", "scenario", "force_snr", "acceleration_snr", "true_mass", "total_snr", "metric", "method", "sample_mean", "empirical_mcse", "status", "note", "paired_our_minus_method", "paired_mcse", "paired_status", "paired_note", "noise_model", "force_noise_sd", "acceleration_noise_sd", "truth_direction", "snr_definition", "readings_per_trial", "effective_force_snr", "effective_acceleration_snr", "effective_total_snr", "calibration"];
  const data = [header];
  selectedRows().forEach(row => state.run.methods.forEach(method => {
    const value = valueFor(row, method.id), paired = pairFor(row, method.id);
    data.push([state.run.id, row.id, row.f, row.a, row.mass, row.strength, state.metric, method.id, value?.mean, value?.mcse, value?.status || "not_recorded", value?.note, paired?.mean, paired?.mcse, paired?.status, paired?.note,
      state.run.config?.noise_model || state.run.metadata?.noise_model || "isotropic_unit_3d", JSON.stringify(state.run.config?.noise?.force_sd || [1, 1, 1]), JSON.stringify(state.run.config?.noise?.acceleration_sd || [1, 1, 1]), JSON.stringify(state.run.config?.direction || [1, 0, 0]), "per-reading magnitude / true channel RMS component noise SD", repeatCount(), row.effective_force_snr ?? row.f, row.effective_acceleration_snr ?? row.a, row.effective_total_signal_snr ?? row.strength, JSON.stringify(measurementInfo().calibration)]);
  }));
  const encode = value => {
    let text = value === undefined || value === null ? "" : String(value);
    // Spreadsheet formula characters in textual metadata should stay literal.
    if (typeof value === "string" && /^[=+@\-\t\r]/.test(text)) text = "'" + text;
    return /[",\n\r]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  };
  download(`${safeName(state.run.title)}-${state.metric}.csv`, data.map(row => row.map(encode).join(",")).join("\r\n"), "text/csv;charset=utf-8");
}

function parseNumbers(text, label, allowZero = false) {
  const tokens = text.trim().split(/[,;\s]+/), values = tokens.map(Number);
  if (!text.trim() || tokens.some(token => !token) || values.some(value => !finite(value) || (allowZero ? value < 0 : value <= 0))) throw new Error(`${label} must contain ${allowZero ? "nonnegative" : "positive"} finite numbers.`);
  return values;
}
function integerInput(id, label) {
  const value = Number($(id).value);
  if (!$(id).value.trim() || !Number.isSafeInteger(value) || value < (id.includes("seed") || ["config-training", "config-heldout"].includes(id) ? 0 : 1)) throw new Error(`${label} must be a valid whole number.`);
  return value;
}
function inputVector(prefix, label, positive, maxRatio = 16) {
  const values = ["x", "y", "z"].map(axis => {
    const raw = $(`${prefix}-${axis}`).value.trim();
    const value = Number(raw);
    if (!raw || !finite(value) || (positive && value <= 0)) throw new Error(`${label} needs three ${positive ? "positive " : ""}finite numbers.`);
    return value;
  });
  if (positive && Math.max(...values) / Math.min(...values) > maxRatio) throw new Error(`${label}: largest / smallest standard deviation must be at most ${maxRatio} for the supported coordinate-noise model.`);
  return values;
}
function rmsSD(values) {
  const scale = Math.max(...values);
  return scale * (Math.hypot(...values.map(value => value / scale)) / Math.sqrt(3));
}
function unitDirection(values) {
  const scale = Math.max(...values.map(Math.abs));
  if (!scale) throw new Error("True signal direction must be nonzero.");
  if (Math.abs(Math.hypot(...values) - 1) <= 5e-15) return values;
  const scaled = values.map(value => value / scale);
  const norm = Math.hypot(...scaled);
  return scaled.map(value => value / norm);
}
function getNoiseSettings() {
  return {
    noise: {force_sd: inputVector("noise-force", "Force noise SD", true), acceleration_sd: inputVector("noise-acceleration", "Acceleration noise SD", true)},
    direction: unitDirection(inputVector("direction", "True signal direction", false))
  };
}
function getMeasurementSettings() {
  const repeats = integerInput("config-repeats", "Readings per trial");
  if (repeats > 1000) throw new Error("Readings per trial must be between 1 and 1,000.");
  const mode = $("config-calibration").value;
  let calibration = {mode};
  if (mode === "scaled") {
    calibration = {mode, force_scale: inputVector("calibration-force", "Force calibration factors", true, Infinity), acceleration_scale: inputVector("calibration-acceleration", "Acceleration calibration factors", true, Infinity)};
    const noise = $("config-noise").value === "diagonal_gaussian_3d" ? getNoiseSettings().noise : {force_sd: [1, 1, 1], acceleration_sd: [1, 1, 1]};
    ["force", "acceleration"].forEach(channel => {
      const supplied = noise[`${channel}_sd`].map((value, index) => value * calibration[`${channel}_scale`][index]);
      if (supplied.some(value => !finite(value) || value <= 0) || Math.max(...supplied) / Math.min(...supplied) > 16) throw new Error(`Supplied ${channel} noise SDs must be positive and finite, with largest / smallest at most 16.`);
    });
  } else if (["estimated", "pooled_isotropic"].includes(mode)) {
    calibration = {mode, samples: integerInput("config-calibration-samples", "Calibration samples"), seed: integerInput("config-calibration-seed", "Calibration seed")};
    if (calibration.samples < 4 || calibration.samples > 100000) throw new Error("Calibration samples must be between 4 and 100,000.");
    if (calibration.seed > 4294967295 || [Number($("config-training").value), Number($("config-heldout").value)].includes(calibration.seed)) throw new Error("Calibration seed must be between 0 and 4,294,967,295 and differ from the training and held-out seeds.");
    if (mode === "pooled_isotropic") {
      const noise = $("config-noise").value === "diagonal_gaussian_3d" ? getNoiseSettings().noise : {force_sd: [1, 1, 1], acceleration_sd: [1, 1, 1]};
      calibration.calibration_order = integerInput("config-calibration-order", "Calibration quadrature order");
      if (calibration.calibration_order < 8 || calibration.calibration_order > 96) throw new Error("Calibration quadrature order must be between 8 and 96.");
      if (Object.values(noise).some(vector => vector.some(value => value !== vector[0]))) throw new Error("Pooled calibration requires equal coordinate noise SDs within each channel.");
    }
  } else if (mode !== "known") throw new Error("Choose a supported noise calibration mode.");
  return {repeats, calibration};
}
function getConfig() {
  let scenario;
  const type = $("config-type").value;
  if (type === "preset") scenario = {type, selection: $("config-preset").value};
  else if (type === "force_path") scenario = {type, force_snr: parseNumbers($("config-force").value, "Force SNR")[0], acceleration_snrs: parseNumbers($("config-accelerations").value, "Acceleration SNRs")};
  else if (type === "fixed_mass") scenario = {type, masses: parseNumbers($("config-masses").value, "Mass ratios"), signal_snrs: parseNumbers($("config-signals").value, "Total signal SNRs")};
  else if (type === "mass_excitation") scenario = {type, masses: parseNumbers($("config-excitation-masses").value, "Physical masses"), acceleration_snrs: parseNumbers($("config-excitation-accelerations").value, "Acceleration SNRs", true)};
  else {
    const pairs = $("config-pairs").value.trim().split(/\n+/).map((line, index) => {
      const pair = parseNumbers(line, `Pair ${index + 1}`);
      if (pair.length !== 2) throw new Error(`Pair ${index + 1} must contain exactly two numbers: force and acceleration SNR.`);
      return pair;
    });
    scenario = {type: "pairs", pairs};
  }
  const name = $("config-name").value.trim();
  if (!name) throw new Error("Give this experiment a name.");
  const config = {schema_version: 1, name, samples: integerInput("config-samples", "Trials"), order: integerInput("config-order", "Quadrature order"), batch_size: integerInput("config-batch", "Batch size"), training_seed: integerInput("config-training", "Training seed"), heldout_seed: integerInput("config-heldout", "Held-out seed"), save_trials: $("config-save").checked, noise_model: $("config-noise").value, scenario};
  if (config.noise_model === "diagonal_gaussian_3d") Object.assign(config, getNoiseSettings());
  const measurement = getMeasurementSettings();
  if (state.explicitMeasurement || measurement.repeats !== 1 || measurement.calibration.mode !== "known") config.measurement = measurement;
  return config;
}

function setConfig(config) {
  if (!config || typeof config !== "object" || config.schema_version !== 1) throw new Error("Expected a version 1 experiment configuration.");
  if (!["isotropic_unit_3d", "diagonal_gaussian_3d"].includes(config.noise_model)) throw new Error("Unsupported noise model; use isotropic_unit_3d or diagonal_gaussian_3d.");
  if (!["preset", "force_path", "fixed_mass", "mass_excitation", "pairs"].includes(config.scenario?.type)) throw new Error("The configuration has an unsupported signal design.");
  const measurement = config.measurement || {repeats: 1, calibration: {mode: "known"}};
  if (!["known", "scaled", "estimated", "pooled_isotropic"].includes(measurement.calibration?.mode)) throw new Error("The configuration has an unsupported calibration mode.");
  if (measurement.calibration.mode === "scaled" && ![measurement.calibration.force_scale, measurement.calibration.acceleration_scale].every(vector => Array.isArray(vector) && vector.length === 3)) throw new Error("Calibration factors each require three coordinates.");
  const mapping = {"config-name": "name", "config-samples": "samples", "config-order": "order", "config-batch": "batch_size", "config-training": "training_seed", "config-heldout": "heldout_seed", "config-noise": "noise_model"};
  Object.entries(mapping).forEach(([id, key]) => { if (config[key] !== undefined) $(id).value = config[key]; });
  $("config-save").checked = !!config.save_trials;
  state.explicitMeasurement = Object.hasOwn(config, "measurement");
  $("config-repeats").value = measurement.repeats;
  $("config-calibration").value = measurement.calibration.mode;
  $("config-calibration-samples").value = measurement.calibration.samples ?? 20;
  $("config-calibration-seed").value = measurement.calibration.seed ?? 2026091903;
  $("config-calibration-order").value = measurement.calibration.calibration_order ?? 24;
  ["x", "y", "z"].forEach((axis, index) => { $(`calibration-force-${axis}`).value = measurement.calibration.force_scale?.[index] ?? 1; $(`calibration-acceleration-${axis}`).value = measurement.calibration.acceleration_scale?.[index] ?? 1; });
  $("config-type").value = config.scenario.type;
  if (config.scenario.type === "preset") $("config-preset").value = config.scenario.selection;
  if (config.scenario.type === "force_path") { $("config-force").value = config.scenario.force_snr; $("config-accelerations").value = config.scenario.acceleration_snrs.join(", "); }
  if (config.scenario.type === "fixed_mass") { $("config-masses").value = config.scenario.masses.join(", "); $("config-signals").value = config.scenario.signal_snrs.join(", "); }
  if (config.scenario.type === "mass_excitation") { $("config-excitation-masses").value = config.scenario.masses.join(", "); $("config-excitation-accelerations").value = config.scenario.acceleration_snrs.join(", "); }
  if (config.scenario.type === "pairs") $("config-pairs").value = config.scenario.pairs.map(pair => pair.join(", ")).join("\n");
  if (config.noise_model === "diagonal_gaussian_3d") {
    const force = config.noise?.force_sd || [1, 1, 1], acceleration = config.noise?.acceleration_sd || [1, 1, 1], direction = config.direction || [1, 0, 0];
    if (![force, acceleration, direction].every(vector => Array.isArray(vector) && vector.length === 3)) throw new Error("Noise SDs and direction each require three coordinates.");
    ["x", "y", "z"].forEach((axis, index) => { $(`noise-force-${axis}`).value = force[index]; $(`noise-acceleration-${axis}`).value = acceleration[index]; $(`direction-${axis}`).value = direction[index]; });
  }
  updateSetup();
}

function updateSetup() {
  document.querySelectorAll(".conditional").forEach(node => { node.hidden = node.dataset.type !== $("config-type").value; });
  const diagonal = $("config-noise").value === "diagonal_gaussian_3d";
  $("diagonal-noise-settings").hidden = !diagonal;
  $("diagonal-noise-settings").querySelectorAll("input").forEach(input => { input.disabled = !diagonal; });
  const calibrationMode = $("config-calibration").value;
  ["scaled", "estimated"].forEach(mode => {
    const section = $(`${mode}-calibration-settings`);
    const active = calibrationMode === mode || (mode === "estimated" && calibrationMode === "pooled_isotropic");
    section.hidden = !active;
    section.querySelectorAll("input").forEach(input => { input.disabled = !active; });
  });
  $("calibration-order-setting").hidden = calibrationMode !== "pooled_isotropic";
  $("config-calibration-order").disabled = calibrationMode !== "pooled_isotropic";
  $("calibration-model-hint").textContent = calibrationMode === "pooled_isotropic" ? "Equal coordinate noise within each channel is assumed. Compare plug-in calibration, integrated calibration uncertainty and a true-noise oracle on identical observations. Only reading means and external calibration enter inference." : "Each trial receives independent coordinate SD estimates, treated as fixed in mass inference. These intervals omit calibration uncertainty.";
  $("config-samples").max = diagonal || calibrationMode !== "known" || Number($("config-repeats").value) > 1 ? "10000" : "1048576";
  $("noise-explanation").textContent = diagonal ? "Force and acceleration may have different noise on x, y, and z. Component errors remain independent and Gaussian." : "Unit standard deviation in every coordinate of both channels.";
  if (diagonal) {
    try {
      const settings = getNoiseSettings(), force = rmsSD(settings.noise.force_sd), acceleration = rmsSD(settings.noise.acceleration_sd);
      const massRule = $("config-type").value === "mass_excitation" ? "Acceleration = its SNR × acceleration RMS SD; force = declared mass × acceleration." : `Physical mass = SNR ratio × ${format(force / acceleration, 5)}.`;
      $("noise-summary").textContent = `RMS noise SD: force ${format(force, 5)} · acceleration ${format(acceleration, 5)}. ${massRule} Unit truth direction: [${settings.direction.map(value => format(value, 4)).join(", ")}].`;
      $("noise-summary").classList.remove("invalid");
    } catch (error) { $("noise-summary").textContent = errorText(error); $("noise-summary").classList.add("invalid"); }
  }
  try {
    const measurement = getMeasurementSettings();
    $("measurement-summary").textContent = `${measurement.repeats} ${measurement.repeats === 1 ? "reading" : "readings"} per trial · mean signal / noise is ${format(Math.sqrt(measurement.repeats), 4)}× the per-reading value. ${calibrationLabel(measurement.calibration)}.`;
    $("measurement-summary").classList.remove("invalid");
  } catch (error) { $("measurement-summary").textContent = errorText(error); $("measurement-summary").classList.add("invalid"); }
  try {
    const config = getConfig(), scenario = config.scenario;
    const count = scenario.type === "preset" ? {all: 46, original: 16, fixed_mass: 30}[scenario.selection] : scenario.type === "force_path" ? scenario.acceleration_snrs.length : scenario.type === "fixed_mass" ? scenario.masses.length * scenario.signal_snrs.length : scenario.type === "mass_excitation" ? scenario.masses.length * scenario.acceleration_snrs.length : scenario.pairs.length;
    const readings = config.measurement?.repeats || 1;
    $("setup-size").textContent = `${wholeFormat.format(count)} scenarios · ${wholeFormat.format(count * config.samples)} training trials${readings > 1 ? ` · ${readings} readings / trial` : ""}`;
  } catch { $("setup-size").textContent = "Complete the signal settings to calculate run size."; }
}

function applyNoiseProfile(profile) {
  const definitions = {equal: [[1, 1, 1], [1, 1, 1]], transverse: [[.5, 2, 2], [.5, 2, 2]], parallel: [[2, .5, .5], [2, .5, .5]], crossed: [[.5, 1, 2], [2, 1, .5]]};
  const values = definitions[profile];
  if (!values) return;
  const [force, acceleration] = values.map(vector => { const rms = rmsSD(vector); return vector.map(value => value / rms); });
  ["x", "y", "z"].forEach((axis, index) => { $(`noise-force-${axis}`).value = force[index]; $(`noise-acceleration-${axis}`).value = acceleration[index]; $(`direction-${axis}`).value = index === 0 ? "1" : "0"; });
  updateSetup();
  toast("Noise profile applied: RMS SD = 1 in each channel; true direction = x.");
}

function openSetup() {
  $("form-error").hidden = true;
  $("start-run").disabled = state.job?.status === "running";
  updateSetup();
  $("experiment-dialog").showModal();
}

async function startRun(event) {
  event.preventDefault();
  $("form-error").hidden = true;
  $("start-run").disabled = true;
  try {
    const config = getConfig();
    if (!state.token) state.token = (await api("/api/session")).token;
    const result = await api("/api/runs", {method: "POST", headers: {"Content-Type": "application/json", "X-Session-Token": state.token}, body: JSON.stringify(config)});
    updateJob(result.job);
    $("experiment-dialog").close();
    toast("Experiment started. You can keep exploring saved results.");
  } catch (error) {
    $("form-error").textContent = errorText(error);
    $("form-error").hidden = false;
    $("form-error").scrollIntoView({block: "nearest"});
  } finally { $("start-run").disabled = state.job?.status === "running"; }
}

function updateJob(job) {
  state.job = job;
  clearTimeout(pollTimer);
  const banner = $("job-banner");
  banner.hidden = !job;
  if (!job) return;
  banner.classList.toggle("failed", job.status === "failed");
  const title = job.config?.name || "Experiment";
  const label = job.status === "running" ? `Running ${title}` : job.status === "complete" ? `Completed ${title}` : `${title} needs attention`;
  const children = [element("strong", {text: label})];
  if (job.status === "running") {
    children.push(element("progress", {max: job.total || 1, value: job.completed || 0, "aria-label": "Completed scenarios"}));
    children.push(element("span", {text: `${job.completed || 0} / ${job.total || "…"} scenarios`}));
    pollTimer = setTimeout(pollJob, 1800);
  }
  if (job.message || job.error) children.push(element("span", {class: "job-message", text: job.error || job.message}));
  if (job.status === "complete" && job.run_id) children.push(element("button", {class: "text-button", text: "Open results →", onclick: () => loadRun(job.run_id).catch(error => showError(errorText(error)))}));
  if (job.status !== "running") children.push(element("button", {class: "text-button", text: "Dismiss", onclick: () => { banner.hidden = true; }}));
  banner.replaceChildren(...children);
  $("start-run").disabled = job.status === "running";
}

async function pollJob() {
  try {
    const wasRunning = state.job?.status === "running";
    const result = await api("/api/job");
    updateJob(result.job);
    if (wasRunning && result.job?.status === "complete") {
      await refreshRuns(result.job.run_id || null);
      toast("Your experiment is complete and saved in the library.");
    }
  } catch (error) {
    showError(`Could not update experiment progress: ${errorText(error)}`);
    pollTimer = setTimeout(pollJob, 5000);
  }
}

$("metric-select").addEventListener("change", event => {
  state.metric = event.target.value;
  if (!state.run.rows.some(row => valueFor(row, state.method) || pairFor(row, state.method))) {
    const available = state.run.methods.find(method => method.id !== "flat_joint" && state.run.rows.some(row => valueFor(row, method.id) || pairFor(row, method.id)));
    if (available) { state.method = available.id; $("method-select").value = state.method; }
  }
  renderComparison();
});
$("method-select").addEventListener("change", event => { state.method = event.target.value; renderComparison(); });
$("path-select").addEventListener("change", event => { state.path = event.target.value; state.row = null; renderComparison(); });
$("axis-select").addEventListener("change", event => { state.axis = event.target.value; renderComparison(); });
$("scenario-select").addEventListener("change", event => { state.row = event.target.value; renderComparison(); });
$("refresh-runs").addEventListener("click", () => refreshRuns().catch(error => showError(errorText(error))));
$("run-search").addEventListener("input", renderRunList);
$("export-csv").addEventListener("click", exportCSV);
$("download-raw").addEventListener("click", () => { if (state.run) { const link = element("a", {href: `/api/runs/${encodeURIComponent(state.run.id)}/raw`, download: `${safeName(state.run.title)}-results.json`}); document.body.append(link); link.click(); link.remove(); } });
$("new-run").addEventListener("click", openSetup);
$("empty-new-run").addEventListener("click", openSetup);
$("close-dialog").addEventListener("click", () => $("experiment-dialog").close());
$("experiment-form").addEventListener("submit", startRun);
$("experiment-form").addEventListener("input", updateSetup);
$("experiment-form").addEventListener("change", updateSetup);
document.querySelectorAll("[data-noise-profile]").forEach(button => button.addEventListener("click", () => applyNoiseProfile(button.dataset.noiseProfile)));
$("load-config").addEventListener("click", () => $("config-file").click());
$("config-file").addEventListener("change", async event => {
  const file = event.target.files[0];
  if (!file) return;
  try { setConfig(JSON.parse(await file.text())); $("form-error").hidden = true; toast("Configuration loaded. Review settings before running."); }
  catch (error) { $("form-error").textContent = errorText(error); $("form-error").hidden = false; }
  event.target.value = "";
});
$("save-config").addEventListener("click", () => {
  try { const config = getConfig(); download(`${safeName(config.name)}-config.json`, JSON.stringify(config, null, 2)); }
  catch (error) { $("form-error").textContent = errorText(error); $("form-error").hidden = false; }
});
$("reuse-config").addEventListener("click", () => { if (state.run?.config) { setConfig(state.run.config); toast("Selected run’s settings loaded."); } });

async function initialize() {
  try {
    await refreshRuns();
    state.token = (await api("/api/session")).token;
  } catch (error) { showError(errorText(error)); }
}
initialize();
