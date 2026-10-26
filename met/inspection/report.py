"""Presentation of an inspection result, independent of its scientific producers.

The only entry point is ``write_report(folder, result)``. Sections are optional,
unknown fields remain inspectable, and no inference or check is performed here.
Matplotlib is optional; the HTML and supplied check outcomes remain available
without it. All plots use supplied arrays, with display transformations labelled.
"""

from collections.abc import Mapping
from html import escape
from io import StringIO
import json
from pathlib import Path
import re

import numpy as np


_STYLE = """
body{font:15px/1.55 system-ui,sans-serif;color:#1e293b;background:#f1f5f9;margin:0}
main{max-width:1240px;margin:auto;padding:32px 24px 64px;background:white}
h1{font-size:30px;line-height:1.2}h2{margin-top:38px;border-top:1px solid #e2e8f0;padding-top:22px}
h3{margin-top:25px}p{max-width:1050px}a{color:#1d4ed8}code,pre{font-family:ui-monospace,monospace}
pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;background:#f8fafc;padding:12px;border-radius:5px}
.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #e2e8f0}
th{background:#f8fafc}td pre{margin:0;min-width:140px}details{margin:10px 0}summary{cursor:pointer}
.note{background:#eff6ff;padding:12px 16px;border-left:4px solid #3b82f6}
.warning{background:#fff7ed;padding:14px 18px;border-left:4px solid #ea580c;margin:20px 0}
.warning h2{border:0;margin:0;padding:0;font-size:20px}
.muted,figcaption{color:#475569;font-size:13px}.badge{display:inline-block;border-radius:4px;padding:2px 7px;font-weight:650}
.pass{background:#dcfce7;color:#166534}.fail,.error{background:#fee2e2;color:#991b1b}
.unresolved{background:#fef3c7;color:#92400e}.inapplicable{background:#e2e8f0;color:#475569}
figure{margin:22px 0}figure svg{display:block;width:100%;height:auto}figcaption{padding:7px 0}
ul{padding-left:22px}li{margin:3px 0}@media print{main{padding:0}details{break-inside:avoid}figure{break-inside:avoid}}
"""


def _plain(value):
    """JSON-friendly, compact display without imposing a scientific schema."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        if value.size <= 64:
            return _plain(value.tolist())
        return {"array_shape": list(value.shape), "dtype": str(value.dtype),
                "values": "See arrays.npz / run.json for the supplied array"}
    if isinstance(value, (list, tuple)):
        if len(value) > 64:
            return {"length": len(value), "first_items": [_plain(item) for item in value[:8]],
                    "values": "See run.json / arrays.npz for complete values"}
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _pre(value):
    return "<pre>" + escape(json.dumps(_plain(value), indent=2, default=str, ensure_ascii=False)) + "</pre>"


def _details(label, value, opened=False):
    return ("<details" + (" open" if opened else "") + "><summary>" + escape(str(label))
            + "</summary>" + _pre(value) + "</details>")


def _fmt(value):
    if value is None:
        return "—"
    if isinstance(value, (float, int, np.number)):
        return format(value, ".6g")
    return str(value)


def _table(headers, rows):
    return ('<div class="scroll"><table><thead><tr>'
            + "".join("<th>" + escape(str(item)) + "</th>" for item in headers)
            + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join("<td>" + item + "</td>" for item in row) + "</tr>" for row in rows)
            + "</tbody></table></div>")


def _mapping(value):
    return value if isinstance(value, Mapping) else {}


def _records(value):
    return value if isinstance(value, (list, tuple)) else []


def _array(value, ndim=None):
    if value is None:
        return None
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return None
    return array if ndim is None or array.ndim == ndim else None


def _log_mass(value):
    values = _array(value, 1)
    if values is None:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(values > 0, np.log(values), np.nan)


def _curve(ax, x, value, label, **kwargs):
    y = _array(value, 1)
    if x is None or y is None or y.shape != x.shape:
        return False
    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        return False
    # Mask, rather than drop, so missing interior values do not become bridges.
    ax.plot(x, np.ma.masked_where(~finite, y), label=label, **kwargs)
    return True


def _components(ax, x, values, label):
    array = _array(values, 2)
    if array is None or x is None or len(array) != len(x):
        return False
    plotted = False
    for index in range(array.shape[1]):
        plotted |= _curve(ax, x, array[:, index], f"{label} {index + 1}")
    return plotted


def _finish_axis(ax, title, ylabel, *, costs=False):
    ax.set(title=title, xlabel="ln(m) in declared mass units", ylabel=ylabel)
    if costs:
        ax.set_yscale("symlog", linthresh=1)
    ax.grid(True, alpha=.2)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, loc="best")
    else:
        ax.text(.5, .5, "No compatible supplied curves", transform=ax.transAxes, ha="center", color="#64748b")


def _save_figure(plt, folder, fig, stem):
    fig.savefig(folder / (stem + ".png"), dpi=150, facecolor="white")
    fig.savefig(folder / (stem + ".svg"), facecolor="white")
    stream = StringIO()
    fig.savefig(stream, format="svg", facecolor="white")
    plt.close(fig)
    svg = stream.getvalue()
    svg = svg[svg.index("<svg"):]
    # Embedded SVGs share an HTML namespace: make clipping and glyph IDs unique.
    ids = sorted(set(re.findall(r'\bid="([^"]+)"', svg)), key=len, reverse=True)
    for ident in ids:
        svg = svg.replace(f'id="{ident}"', f'id="{stem}_{ident}"')
        svg = svg.replace(f'#{ident})', f'#{stem}_{ident})')
        svg = svg.replace(f'="#{ident}"', f'="#{stem}_{ident}"')
    return svg


def _law_figure(plt, evaluation):
    laws = _mapping(evaluation.get("laws"))
    x = _log_mass(evaluation.get("mass"))
    prior = _array(evaluation.get("prior_log_density"), 1)
    if not laws and prior is None:
        return None
    fig, ax = plt.subplots(figsize=(12, 5.3), layout="constrained")
    if prior is not None:
        with np.errstate(under="ignore", over="ignore"):
            density = np.exp(prior)
        _curve(ax, x, np.where(density > 0, density, np.nan), "Prior", color="#64748b", linestyle="--")
    colors = {"unweighted": "#2563eb", "untempered": "#9333ea", "quoted": "#dc2626"}
    for name, law_value in laws.items():
        law = _mapping(law_value)
        native = _array(law.get("log_mass"), 1)
        density = _array(law.get("density"), 1)
        color = colors.get(name)
        if density is not None:
            _curve(ax, native, np.where(density > 0, density, np.nan), str(name), color=color)
        point = law.get("point")
        try:
            point = float(point)
            if np.isfinite(point) and point > 0:
                ax.axvline(np.log(point), color=color, linestyle=":", alpha=.8,
                           label=f"{name} point")
        except (TypeError, ValueError):
            pass
        interval = _array(law.get("interval"), 1)
        if interval is not None and interval.size == 2 and np.all(np.isfinite(interval) & (interval > 0)):
            ax.axvspan(*np.log(interval), color=color, alpha=.06)
    ax.set_yscale("log")
    _finish_axis(ax, "Mass laws on their complete supplied native grids", "Density per unit ln(m)")
    return fig


def _cost_figure(plt, diagnostics):
    x = _log_mass(diagnostics.get("mass"))
    if x is None:
        return None
    fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
    _curve(ax, x, diagnostics.get("total_residual_cost"), "Total Newton residual cost", color="#0f172a", linestyle="--")
    _curve(ax, x, diagnostics.get("total_declared_cost"), "Total declared-compatible correction cost", color="#0f172a")
    for index, row_value in enumerate(_records(diagnostics.get("readings"))):
        row = _mapping(row_value)
        name = str(row.get("id", index + 1))
        color = f"C{index % 10}"
        _curve(ax, x, row.get("residual_cost"), name + ": Newton residual", color=color, linestyle="--")
        _curve(ax, x, row.get("declared_cost"), name + ": declared fit", color=color)
        outside = row.get("outside_subspace_cost")
        if isinstance(outside, (float, int, np.number)) and np.isfinite(outside):
            ax.axhline(outside, color=color, linestyle=":", alpha=.65, label=name + ": outside declared span")
    _finish_axis(ax, "Independent geometry: required correction costs", "Squared correction in declared noise units", costs=True)
    return fig


def _evidence_figure(plt, evaluation):
    rows = _records(evaluation.get("per_reading"))
    x = _log_mass(evaluation.get("mass"))
    if not rows or x is None:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for index, row_value in enumerate(rows):
        row = _mapping(row_value)
        name = str(row.get("id", index + 1))
        evidence = _array(row.get("weighted_log_evidence"), 1)
        if evidence is not None and np.any(np.isfinite(evidence)):
            offset = float(np.max(evidence[np.isfinite(evidence)]))
            _curve(axes[0], x, evidence - offset, name + f" (offset {offset:.4g})")
        _curve(axes[1], x, row.get("explained_energy"), name)
    _finish_axis(axes[0], "Per-reading weighted log evidence", "Supplied weighted log evidence − finite maximum")
    _finish_axis(axes[1], "Per-reading explained energy", "Supplied explained energy")
    return fig


def _size_figure(plt, evaluation):
    x = _log_mass(evaluation.get("mass"))
    if x is None or _array(evaluation.get("size_weight"), 1) is None:
        return None
    fig, ax = plt.subplots(figsize=(12, 4), layout="constrained")
    _curve(ax, x, evaluation["size_weight"], "Global readout size weight")
    _finish_axis(ax, "Global size weight entering the readout", "Supplied conditional size weight")
    return fig


def _correction_figure(plt, x, reading):
    if x is None:
        return None
    fields = (("residual", "Newton residual r(m)", "Force units", "r"),
              ("force_correction", "Fitted force correction", "Force units", "δF"),
              ("acceleration_correction", "Fitted acceleration correction", "Acceleration units", "δa"))
    if not any(_array(reading.get(key), 2) is not None for key, *_ in fields):
        return None
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout="constrained")
    for ax, (key, title, units, label) in zip(axes, fields):
        _components(ax, x, reading.get(key), label)
        ax.axhline(0, color="#94a3b8", linewidth=.7)
        _finish_axis(ax, title, units)
    fig.suptitle("Reading " + str(reading.get("id", "")), fontweight="bold")
    return fig


def _ellipse(ax, center, covariance, color):
    from matplotlib.patches import Ellipse

    if center.shape != (2,) or covariance.shape != (2, 2) or not np.all(np.isfinite(covariance)):
        return
    values, vectors = np.linalg.eigh(covariance)
    if np.any(values < 0):
        return
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    ax.add_patch(Ellipse(center, 2 * np.sqrt(values[0]), 2 * np.sqrt(values[1]),
                         angle=angle, facecolor=color, edgecolor=color, alpha=.15))


def _geometry_figure(plt, observation, diagnostic, masses, evaluation):
    force = _array(observation.get("force"), 1)
    acceleration = _array(observation.get("acceleration"), 1)
    covariance = _array(observation.get("covariance"), 2)
    if (force is None or acceleration is None or covariance is None
            or force.shape != (2,) or acceleration.shape != (2,) or covariance.shape != (4, 4)):
        return None
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    laws = _mapping(evaluation.get("laws"))
    quoted = _mapping(laws.get("quoted"))
    point = quoted.get("point")
    selected = None
    if masses is not None and len(masses) and point is not None:
        try:
            point = float(point)
            valid = np.isfinite(masses) & (masses > 0)
            if point > 0 and np.isfinite(point) and np.any(valid):
                differences = np.where(valid, np.abs(np.log(masses) - np.log(point)), np.inf)
                selected = int(np.argmin(differences))
        except (TypeError, ValueError):
            pass
    for ax, measured, block, key, channel, color in (
            (axes[0], force, covariance[:2, :2], "candidate_force", "Force", "#2563eb"),
            (axes[1], acceleration, covariance[2:, 2:], "candidate_acceleration", "Acceleration", "#ea580c")):
        _ellipse(ax, measured, block, color)
        ax.plot(*measured, "o", color=color, label="Observed channel")
        ax.plot([], [], color=color, alpha=.3, linewidth=5, label="Marginal noise ellipse (radius 1)")
        candidates = _array(diagnostic.get(key), 2)
        if candidates is not None and candidates.shape[1] == 2:
            valid = np.all(np.isfinite(candidates), axis=1)
            masked = np.ma.masked_where(~np.repeat(valid[:, None], 2, axis=1), candidates)
            ax.plot(masked[:, 0], masked[:, 1], color="#64748b", linewidth=1, label="Compatible fit along mass grid")
            if selected is not None and selected < len(candidates) and valid[selected]:
                candidate = candidates[selected]
                ax.plot(*candidate, "x", color="#16a34a", markersize=8,
                        label="Fit at grid m = " + _fmt(masses[selected]))
                ax.annotate("", xy=measured, xytext=candidate,
                            arrowprops={"arrowstyle": "->", "color": color, "lw": 1.5})
        ax.set(title=channel + " channel", xlabel=channel + " coordinate 1", ylabel=channel + " coordinate 2")
        ax.set_aspect("equal", adjustable="datalim")
        ax.autoscale_view()
        ax.margins(.15)
        ax.grid(True, alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Observed vectors and declared marginal uncertainty: " + str(observation.get("id", "")), fontweight="bold")
    return fig


def _figures(folder, result):
    """Keep optional plotting failures local to a section, with visible diagnostics."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        return [], ["Figures unavailable: " + str(error)]
    evaluation = _mapping(result.get("evaluation"))
    diagnostics = _mapping(result.get("diagnostics"))
    case = _mapping(result.get("case"))
    jobs = [
        ("Mass laws", "laws", lambda: _law_figure(plt, evaluation),
         "Density is per unit ln(m). All finite positive native-grid values are retained; no peak-relative tail cutoff is applied. "
         "Dotted vertical lines mark supplied points; shading marks supplied intervals. Interval probability is recorded in the result metadata."),
        ("Correction costs", "costs", lambda: _cost_figure(plt, diagnostics),
         "The Newton cost uses rᵀV⁻¹r. The declared fit also enforces the supplied signal subspace. "
         "Each dotted outside-span baseline belongs to the named reading; the total curves retain every reading's outside-span cost. "
         "Costs are squared corrections in noise units, not p-values. The vertical scale is linear near zero and logarithmic farther away."),
        ("Reading contributions", "evidence", lambda: _evidence_figure(plt, evaluation),
         "Evidence curves have separate finite maxima subtracted for display. These shifts show relative shape; "
         "they do not establish absolute compatibility. Q, g, ν, projective φ and η are tabulated separately."),
        ("Readout size weight", "size_weight", lambda: _size_figure(plt, evaluation),
         "This is the supplied global size weight used by the readout. In ours29 it is E[‖stack(a*₁, …, a*ₙ)‖ | m, readings], "
         "the conditional expected norm of the stacked acceleration vector across readings. It is distinct from the mass-law density."),
    ]
    mass = _array(diagnostics.get("mass"), 1)
    x = _log_mass(diagnostics.get("mass"))
    observations = {str(row.get("id", index + 1)): row for index, row in enumerate(_records(case.get("readings")))
                    if isinstance(row, Mapping)}
    for index, reading_value in enumerate(_records(diagnostics.get("readings"))):
        reading = _mapping(reading_value)
        identifier = str(reading.get("id", index + 1))
        jobs.append(("Residual and corrections: " + identifier, f"reading_{index:03d}_components",
                     lambda row=reading: _correction_figure(plt, x, row),
                     "Each component remains in its channel's units. Corrections are observed minus the minimum-cost compatible fit. "
                     "Fitted corrections are not the injected instrument errors."))
        observation = _mapping(observations.get(identifier))
        jobs.append(("Uncertainty geometry: " + identifier, f"reading_{index:03d}_geometry",
                     lambda obs=observation, row=reading: _geometry_figure(plt, obs, row, mass, evaluation),
                     "The channels use separate axes and units. Ellipses show marginal covariance contours of Mahalanobis radius one, "
                     "not a stated joint confidence level. Cross-channel correlations remain in the full input covariance and correction costs. "
                     "The marked fit uses the supplied grid point nearest the quoted mass, when available."))
    fragments, errors = [], []
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "savefig.bbox": "tight"}):
        for title, stem, build, caption in jobs:
            # An optional or newly extended section must not suppress the report.
            before = set(plt.get_fignums())
            try:
                fig = build()
                if fig is not None:
                    svg = _save_figure(plt, folder, fig, stem)
                    fragments.append("<h3>" + escape(title) + "</h3><figure>" + svg
                                     + "<figcaption>" + escape(caption) + ' <a href="' + stem + '.svg">SVG</a> · '
                                     + '<a href="' + stem + '.png">PNG</a></figcaption></figure>')
            except Exception as error:  # Presentation failure, explicitly reported below.
                errors.append(title + ": " + type(error).__name__ + ": " + str(error))
            finally:
                for number in set(plt.get_fignums()) - before:
                    plt.close(number)
    return fragments, errors


def _check_section(checks):
    if not checks:
        return "<p>No checks supplied. The renderer assigns no scientific outcomes.</p>"
    rows = []
    for check_value in checks:
        check = _mapping(check_value)
        status = str(check.get("status", "unresolved"))
        css = status if status in {"pass", "fail", "unresolved", "inapplicable", "error"} else "unresolved"
        identity = escape(str(check.get("id", ""))) + "<br>" + escape(str(check.get("title", "")))
        source = "<b>" + escape(str(check.get("kind", "unspecified"))) + "</b><br>" + escape(str(check.get("source", "unspecified")))
        assumptions = check.get("assumptions")
        if assumptions:
            source += _details("Assumptions", assumptions)
        expected = escape(str(check.get("expected", "Not supplied")))
        observed = _pre(check.get("observed", {}))
        if check.get("tolerance"):
            observed += _details("Tolerance", check["tolerance"])
        if check.get("detail"):
            observed += _details("Detail", check["detail"])
        extras = {key: value for key, value in check.items() if key not in
                  {"id", "title", "status", "kind", "source", "assumptions", "expected", "observed", "tolerance", "detail"}}
        if extras:
            observed += _details("Additional check fields", extras)
        rows.append((identity, '<span class="badge ' + css + '">' + escape(status) + "</span>", source, expected, observed))
    return _table(("Check", "Supplied outcome", "Authority and assumptions", "Expected relationship", "Observed and tolerance"), rows)


def _support_warning(evaluation):
    """Present adapter-supplied support warnings without deciding grid adequacy."""
    rows = []
    for name, law_value in _mapping(evaluation.get("laws")).items():
        law = _mapping(law_value)
        if not law.get("support_edge_warning", False):
            continue
        limits = law.get("log_grid_bounds")
        if limits is None:
            grid = _array(law.get("log_mass"), 1)
            if grid is not None and grid.size:
                limits = [grid[0], grid[-1]]
        rows.append((escape(str(name)), escape(_fmt(law.get("support_edge_density_fraction"))),
                     _pre(limits if limits is not None else "Not supplied")))
    if not rows:
        return ""
    return ('<section class="warning"><h2>Native-grid support requires review</h2>'
            '<p>The adapter flags density remaining at one or both edges of these laws\' finite native grids. '
            'The normalized densities, points and intervals are supplied numerical outputs; the report has not established '
            'that their grids capture the full law. The edge fraction is the maximum endpoint density divided by the '
            'maximum grid density. It is a density ratio, not omitted probability.</p>'
            + _table(("Flagged law", "Supplied edge density fraction", "Native ln(m) grid limits"), rows)
            + '</section>')


def write_report(folder: Path, result: dict) -> bool:
    """Write self-contained HTML, optional PNG/SVG figures; return bool(figures).

    No section is required. The result producer owns all scientific quantities,
    decisions and statuses; this module only presents them. Unknown fields are
    retained in expandable views. ``run.json``, ``arrays.npz`` and ``case.yaml``
    are linked as companion artifacts and are written by the inspection runner.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    result = _mapping(result)
    case = _mapping(result.get("case"))
    evaluation = _mapping(result.get("evaluation"))
    diagnostics = _mapping(result.get("diagnostics"))
    title = str(case.get("name", "Case inspection"))
    parts = ["<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
             '<meta name="viewport" content="width=device-width,initial-scale=1">',
             "<title>" + escape(title) + " — case inspection</title><style>" + _STYLE + "</style></head><body><main>",
             "<h1>" + escape(title) + "</h1>"]
    if case.get("question"):
        parts.append("<p>" + escape(str(case["question"])) + "</p>")
    parts.append('<p><a href="case.yaml">Saved case</a> · <a href="run.json">Structured result</a> · '
                 '<a href="arrays.npz">Complete numeric arrays</a></p>')
    parts.append('<p class="note">Independent geometric diagnostics and estimator internals are shown separately. '
                 'The report assigns no pass/fail outcomes and introduces no required scientific behaviour. '
                 'Correction costs are not p-values. Minimum-cost fitted corrections are not recovered instrument errors.</p>')
    parts.append(_support_warning(evaluation))

    parts.append("<h2>Checks and their authority</h2>")
    parts.append(_check_section(_records(result.get("checks"))))

    parts.append("<h2>Permitted inputs and declarations</h2>")
    parts.append("<p>These are the readings and declarations supplied to the estimator. The joint covariance is ordered "
                 "as force coordinates followed by acceleration coordinates. The subspace is a supplied basis for the allowed signal directions. "
                 "Units, windows and other declarations are preserved below as provided.</p>")
    parts.append(_details("Prior", case.get("prior", {}), opened=True))
    for index, reading_value in enumerate(_records(case.get("readings"))):
        reading = _mapping(reading_value)
        parts.append(_details("Reading " + str(reading.get("id", index + 1)), reading, opened=True))
    case_extra = {key: value for key, value in case.items() if key not in {"name", "question", "prior", "readings"}}
    if case_extra:
        parts.append(_details("Additional case declarations", case_extra))

    if result.get("investigator_truth") is not None:
        parts.append("<h2>Investigator-only truth</h2><p>This separate record is available to the investigator. "
                     "It is not an estimator input. Any latent vectors and injected errors here must be distinguished from fitted corrections.</p>")
        parts.append(_details("Synthetic truth / injected errors", result["investigator_truth"], opened=True))

    parts.append("<h2>Estimator output and internal contributions</h2>")
    parts.append("<p>Adapter: <b>" + escape(str(evaluation.get("adapter", "unspecified"))) + "</b>. "
                 "For ours29, <b>unweighted</b> sets reading weights g and global tempering η to 1; "
                 "<b>untempered</b> uses the current g with η = 1; <b>quoted</b> uses the current g and η. "
                 "These are inspection stages of the supplied equation. The global size weight affects the point readout separately. "
                 "For ours29 this is E[‖stack(a*₁, …, a*ₙ)‖ | m, readings], the conditional expected norm of the "
                 "stacked acceleration vector across all readings.</p>")
    laws = _mapping(evaluation.get("laws"))
    law_rows = []
    for name, law_value in laws.items():
        law = _mapping(law_value)
        law_rows.append((escape(str(name)), escape(_fmt(law.get("point"))), escape(_fmt(law.get("median"))),
                         escape(str(_plain(law.get("interval", "—")))), escape(_fmt(law.get("eta"))),
                         _pre({key: law[key] for key in ("g", "Q") if key in law})))
    if law_rows:
        parts.append(_table(("Law", "Point", "Median", "Supplied interval", "η", "Reading weights / Q"), law_rows))
    contribution_rows = []
    for index, row_value in enumerate(_records(evaluation.get("per_reading"))):
        row = _mapping(row_value)
        contribution_rows.append(tuple(escape(_fmt(row.get(key, index + 1 if key == "id" else None)))
                                       for key in ("id", "Q", "g", "nu", "projective_phi")))
    if contribution_rows:
        parts.append(_table(("Reading", "Q", "g", "ν", "Projective φ"), contribution_rows))
    if evaluation:
        parts.append(_details("All supplied estimator quantities (including extra fields)", evaluation))

    parts.append("<h2>Independent geometry and curves</h2>")
    parts.append("<p>The independent geometry fits compatible candidates across the supplied positive-mass grid. "
                 "The residual r(m) = F̃ − m ã combines both channels' errors. Its covariance is "
                 "V(m) = ΣFF + m²Σaa − m(ΣFa + ΣaF). A declared-compatible candidate also respects the supplied subspace. "
                 "The diagnostic and estimator grids may have different ranges; each density is shown on its own native grid.</p>")
    figures, plot_errors = _figures(folder, result)
    parts.extend(figures)
    if not figures:
        parts.append("<p>No figures were produced. All supplied numeric results and checks remain available below and in the companion files.</p>")
    if plot_errors:
        parts.append(_details("Optional figure availability / rendering notes", plot_errors, opened=True))
    if diagnostics:
        parts.append(_details("All independent diagnostics (including conditional covariances and extra fields)", diagnostics))
    if result.get("provenance"):
        parts.append("<h2>Provenance</h2>" + _pre(result["provenance"]))
    extras = {key: value for key, value in result.items() if key not in
              {"case", "investigator_truth", "diagnostics", "evaluation", "checks", "provenance"}}
    if extras:
        parts.append("<h2>Additional supplied result fields</h2>" + _pre(extras))
    parts.append("</main></body></html>")
    (folder / "report.html").write_text("\n".join(parts), encoding="utf-8")
    return bool(figures)
