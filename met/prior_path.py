"""Prior-centre investigation composed from ordinary cases and inspection functions.

Only the declared half-Cauchy centre varies. The recipe's independent integration
reference is restricted to a single unrestricted isotropic 2D reading; that is not
a restriction on case storage or the general path builder.
"""
import csv
import hashlib
from html import escape
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import yaml

from .full import Reading
from .inspection.case import Case, parse_cases
from .inspection.paths import expand_paths, set_prior_center
from .misalignment import isotropic_minima


def validate(document):
    cases = parse_cases(document)
    centers = document.get("prior_centers")
    if not isinstance(centers, list):
        raise ValueError("prior_path needs prior_centers")
    if type(document.get("seed")) is not int or document["seed"] < 0:
        raise ValueError("prior_path needs a nonnegative seed (no random draws)")
    for case in cases:
        if len(case.inputs.readings) != 1 or case.inputs.readings[0].force.size != 2:
            raise ValueError("this prior-path reference requires one 2D reading")
        r = case.inputs.readings[0]
        isotropic_minima(Reading(r.force, r.acceleration, r.covariance, r.subspace))
    expanded = expand_paths(cases, centers, transform=set_prior_center, parameter="prior_center")
    detail = document.get("detail_centers", centers)
    if not isinstance(detail, list) or any(c not in centers for c in detail):
        raise ValueError("detail_centers must select prior_centers")
    band = document.get("mass_band", [0.5, 4.0])
    if (not isinstance(band, list) or len(band) != 2 or any(
            isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in band)
            or not 0 < band[0] < band[1]):
        raise ValueError("mass_band must be two increasing finite positive masses")
    grid_center = document.get("mass_grid_center", 1.0)
    if isinstance(grid_center, bool) or not isinstance(grid_center, (int, float)) or not math.isfinite(grid_center) or grid_center <= 0:
        raise ValueError("mass_grid_center must be finite and positive")
    numerics = dict(log_half_span=14.0, nodes=301, check_tolerance=1e-7)
    numerics.update(document.get("numerics", {}))
    nodes, half = numerics["nodes"], numerics["log_half_span"]
    if type(nodes) is not int or nodes < 3 or not math.isfinite(half) or half <= 0:
        raise ValueError("common grid needs at least three nodes and positive log_half_span")
    masses = np.asarray(numerics.get("mass_grid", np.exp(math.log(grid_center) + np.linspace(-half, half, nodes))))
    if masses.ndim != 1 or len(masses) < 3 or not np.all(np.isfinite(masses) & (masses > 0)) or not np.all(np.diff(masses) > 0):
        raise ValueError("mass_grid must be increasing finite positive masses")
    numerics["mass_grid"] = masses.tolist()
    display = document.get("plot_log_mass_bounds")
    if display is not None and (not isinstance(display, list) or len(display) != 2 or any(
            isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in display)
            or display[0] >= display[1]):
        raise ValueError("plot_log_mass_bounds must be two increasing finite log-masses")
    return dict(kind="prior_path", raw=document, cases=expanded, numerics=numerics, mass_band=band,
                cells=[dict(id=f"case_{i:03d}", name=c.name) for i, c in enumerate(expanded)])


def summarize(result, mass_band):
    from .inspection.references import angle_reference, cloud_readout_reference
    case = Case.from_dict(result["case"])
    r = case.inputs.readings[0]
    bounds = isotropic_minima(Reading(r.force, r.acceleration, r.covariance, r.subspace))
    state = result["evaluation"]["per_reading"][0]
    center = case.inputs.prior["center"]
    row = dict(case=case.name, path=case.metadata["path"]["source"], prior_center=center,
               Q=state["Q"], g=state["g"], **bounds, laws={})
    for name, law in result["evaluation"]["laws"].items():
        weight = law["eta"] * law["g"][0]
        reference = angle_reference(r, center, weight, np.array([*mass_band, 1.0, *law["interval"]]))
        point = cloud_readout_reference(r, center, weight)
        x, density = law["log_mass"], law["density"]
        native_cdf = np.r_[0, np.cumsum((density[1:] + density[:-1]) * np.diff(x) / 2)]
        native_cdf /= native_cdf[-1]
        interval_error = float(np.max(np.abs(reference["cdf"][-2:] - [0.025, 0.975])))
        log_point_error = abs(math.log(law["point"] / point["point"]))
        row["laws"][name] = dict(point=law["point"], median=law["median"], interval=law["interval"],
            eta=law["eta"], log_evidence_span=weight * bounds["positive_cost_span"] / 2,
            probability_below_one=float(reference["cdf"][2]),
            native_probability_below_one=float(np.interp(0, x, native_cdf)),
            band_probability=float(reference["cdf"][1] - reference["cdf"][0]),
            interval_cdf_error=interval_error, reference_point=point["point"],
            absolute_log_point_error=log_point_error,
            cdf_quadrature_error_estimate=reference["quadrature_error_bound"],
            point_quadrature_error_estimate=point["quadrature_error_estimate"],
            reference_point_numerator=point["numerator"], reference_point_denominator=point["denominator"],
            reference_numerator_error_estimate=point["numerator_error_estimate"],
            reference_denominator_error_estimate=point["denominator_error_estimate"],
            reference_warnings=list(dict.fromkeys(reference["quadrature_warnings"] + point["quadrature_warnings"])),
            support_edge_warning=law["support_edge_warning"],
            support_edge_density_fraction=law["support_edge_density_fraction"])
    return row


def compare_paths(rows, results):
    """Check data-only factors at identical physical masses across prior choices."""
    comparisons = []
    for group in dict.fromkeys(r["path"] for r in rows):
        pairs = [(row, result) for row, result in zip(rows, results) if row["path"] == group]
        base_row, base = pairs[0]
        base_eval = base["evaluation"]
        base_fit = base["diagnostics"]["readings"][0]
        errors = dict(data_only_log_evidence=0.0, Q=0.0, g=0.0, eta=0.0,
                      correction=0.0, conditional_covariance=0.0, cloud_size_weight=0.0)
        for row, result in pairs[1:]:
            evaluation = result["evaluation"]
            if not np.array_equal(result["diagnostics"]["mass"], base["diagnostics"]["mass"]):
                raise ValueError("prior comparisons require the same physical candidate grid")
            for name, law in evaluation["laws"].items():
                baseline = base_eval["laws"][name]
                difference = (law["log_kernel"] - evaluation["prior_log_density"]
                              - baseline["log_kernel"] + base_eval["prior_log_density"])
                difference -= difference[len(difference) // 2]
                errors["data_only_log_evidence"] = max(errors["data_only_log_evidence"], float(np.max(np.abs(difference))))
                errors["eta"] = max(errors["eta"], abs(law["eta"] - baseline["eta"]))
            errors["Q"] = max(errors["Q"], abs(row["Q"] - base_row["Q"]))
            errors["g"] = max(errors["g"], abs(row["g"] - base_row["g"]))
            fit = result["diagnostics"]["readings"][0]
            for key in ("force_correction", "acceleration_correction"):
                errors["correction"] = max(errors["correction"], float(np.max(np.abs(fit[key] - base_fit[key]))))
            errors["conditional_covariance"] = max(errors["conditional_covariance"], float(np.max(np.abs(
                fit["conditional_covariance"] - base_fit["conditional_covariance"]))))
            errors["cloud_size_weight"] = max(errors["cloud_size_weight"], float(np.max(np.abs(
                evaluation["size_weight"] - base_eval["size_weight"]))))
        ordered = sorted((r for r, _ in pairs), key=lambda r: r["prior_center"])
        sensitivity = (math.log(ordered[-1]["laws"]["quoted"]["point"] / ordered[0]["laws"]["quoted"]["point"])
                       / math.log(ordered[-1]["prior_center"] / ordered[0]["prior_center"])) if len(ordered) > 1 else None
        comparisons.append(dict(path=group, maximum_absolute_errors=errors,
            point_log_sensitivity=sensitivity,
            sensitivity_definition="log(point_at_largest_c/point_at_smallest_c) / log(largest_c/smallest_c); descriptive finite change, not a derivative or criterion",
            authority="Fixed declarations leave geometry and cloud weights unchanged. One-reading eta=1; a prior change only multiplies the mass kernel by the prior-density ratio. Curvature and the projective coordinate themselves need not be invariant.",
            point_range=[min(r["laws"]["quoted"]["point"] for r, _ in pairs), max(r["laws"]["quoted"]["point"] for r, _ in pairs)]))
    return comparisons


def write_comparison(out, rows, results, mass_band, display=None):
    relations = compare_paths(rows, results)
    flat = []
    for row in rows:
        for name, law in row["laws"].items():
            flat.append(dict({k: v for k, v in row.items() if k != "laws"}, law=name,
                **{k: v for k, v in law.items() if k != "interval"},
                interval_low=law["interval"][0], interval_high=law["interval"][1]))
    (out / "comparison.json").write_text(json.dumps(dict(mass_band=mass_band, rows=rows, path_relations=relations), indent=2, allow_nan=False), encoding="utf-8")
    with (out / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    figures = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        figures.append("<p>Matplotlib unavailable; all numerical results are saved.</p>")
    else:
        groups = list(dict.fromkeys(row["path"] for row in rows))
        for mode in ("laws", "points"):
            fig, axes = plt.subplots(math.ceil(len(groups) / 2), 2, figsize=(13, 3.6 * math.ceil(len(groups) / 2)), squeeze=False)
            for ax, group in zip(axes.flat, groups):
                selected = [(row, result) for row, result in zip(rows, results) if row["path"] == group]
                if mode == "laws":
                    for row, result in selected:
                        law, evaluation = result["evaluation"]["laws"]["quoted"], result["evaluation"]
                        line, = ax.plot(law["log_mass"], law["density"], label=f'law c={row["prior_center"]:g}')
                        ax.plot(np.log(evaluation["mass"]), np.exp(evaluation["prior_log_density"]), "--", color=line.get_color(), alpha=.6, label=f'prior c={row["prior_center"]:g}')
                    ax.set(xlabel="ln(mass)", ylabel="Density per ln(mass)")
                    if display is not None:
                        ax.set_xlim(display)
                else:
                    centers = [r["prior_center"] for r, _ in selected]
                    laws = [r["laws"]["quoted"] for r, _ in selected]
                    ax.fill_between(centers, [r["interval"][0] for r in laws], [r["interval"][1] for r in laws], alpha=.15, label="95% equal-tailed interval")
                    ax.plot(centers, [r["point"] for r in laws], "o-", label="cloud point")
                    ax.plot(centers, [r["median"] for r in laws], "x--", label="numerical median")
                    ax.plot(centers, centers, "k:", label="prior centre")
                    ax.set(xscale="log", yscale="log", xlabel="Declared prior centre", ylabel="Mass")
                ax.set_title(group)
                ax.grid(alpha=.15)
                ax.legend(fontsize=8, ncol=2)
            for ax in list(axes.flat)[len(groups):]:
                ax.set_visible(False)
            fig.tight_layout()
            for suffix in ("png", "svg"):
                fig.savefig(out / f"{mode}.{suffix}", dpi=140)
            plt.close(fig)
            figures.append("<svg" + (out / f"{mode}.svg").read_text(encoding="utf-8").split("<svg", 1)[1])
    header = ["Case", "Prior c", "Q", "g", "Point", "Median", "95% interval", "P(m < 1)", "Point / c", "Inspect"]
    table = ["<tr>" + "".join(f"<th>{escape(x)}</th>" for x in header) + "</tr>"]
    for index, row in enumerate(rows):
        law = row["laws"]["quoted"]
        values = [row["path"], f'{row["prior_center"]:g}', f'{row["Q"]:.6g}', f'{row["g"]:.6g}',
                  f'{law["point"]:.6g}', f'{law["median"]:.6g}', f'[{law["interval"][0]:.6g}, {law["interval"][1]:.6g}]',
                  f'{law["probability_below_one"]:.6g}', f'{law["point"]/row["prior_center"]:.6g}']
        table.append("<tr>" + "".join(f"<td>{escape(x)}</td>" for x in values) + f'<td><a href="case_{index:03d}/report.html">Open</a></td></tr>')
    laws = [law for row in rows for law in row["laws"].values()]
    html = ['<!doctype html><meta charset="utf-8"><title>Prior-centre paths</title>',
        '<style>body{font:16px/1.5 system-ui;max-width:1300px;margin:36px auto;padding:20px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:right}td:first-child{text-align:left}svg{width:100%;height:auto}.scroll{overflow:auto}</style>',
        '<h1>What does the declared prior determine?</h1>',
        '<p>Measured vectors, covariances and candidate subspaces are fixed. Only the half-Cauchy prior centre changes. '
        'Every case uses the same physical diagnostic mass grid; each law keeps its own native probability grid. Prior dependence is allowed by the declared input type.</p>',
        '<p>The equal-perpendicular law equals its supplied prior. The cloud point additionally uses a mass-dependent '
        'candidate-size weight, so equality of the law with the prior does not imply that point equals prior centre. '
        'The opposite case can have two separated modes; a central interval can bridge a nearly empty region. '
        'The known numerical median discrepancy at c=1 remains visible.</p>',
        f'<p>Independent angle integration: maximum interval-endpoint CDF error {max(x["interval_cdf_error"] for x in laws):.3g}; '
        f'maximum absolute log point error {max(x["absolute_log_point_error"] for x in laws):.3g}. '
        f'Native support warnings: {sum(x["support_edge_warning"] for x in laws)}; '
        f'reference warning stages: {sum(bool(x["reference_warnings"]) for x in laws)}. '
        'All three law stages and numerical details are in the comparison JSON. These are numerical checks, not newly adopted physical criteria.</p>',
        '<p>Plot window: ' + (escape(str(display)) + ' in ln(mass); tails outside this display are retained in native arrays.' if display else 'native grid limits.') + '</p>',
        '<p><a href="report.html">Case checks</a> · <a href="comparison.csv">CSV</a> · <a href="comparison.json">Full comparison</a> · '
        '<a href="recipe.yaml">Recipe</a> · <a href="preset.yaml">Expanded replay cases</a></p>',
        '<div class="scroll"><table>', *table, '</table></div>',
        '<details><summary>What stayed fixed across priors</summary><pre>' + escape(json.dumps(relations, indent=2)) + '</pre></details>', *figures]
    (out / "comparison.html").write_text("\n".join(html), encoding="utf-8")


def run(path, out_root="results", progress=print, *, reporter=None):
    from .inspection.report import write_report
    from .inspection.run import run as inspect_run
    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    document = dict(spec["raw"], kind="inspection", numerics=spec["numerics"],
                    cases=[case.to_dict(include_truth=True) for case in spec["cases"]])
    rows, results = [], []
    render = reporter or write_report

    def collect(folder, result):
        rows.append(summarize(result, spec["mass_band"]))
        results.append(result)
        if render is write_report:
            plots = result["case"]["prior"]["center"] in spec["raw"].get("detail_centers", spec["raw"]["prior_centers"])
            return render(folder, result, plots=plots)
        return render(folder, result)

    with TemporaryDirectory(prefix="met-prior-path-") as temp:
        expanded = Path(temp) / "expanded.yaml"
        expanded.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        out, bundle = inspect_run(expanded, out_root, progress=progress, reporter=collect)
    (out / "recipe.yaml").write_bytes(source)
    bundle["path_recipe"] = dict(sha256=hashlib.sha256(source).hexdigest(), file="recipe.yaml", transformation="met.inspection.paths.set_prior_center")
    (out / "run.json").write_text(json.dumps(bundle, indent=2, allow_nan=False), encoding="utf-8")
    write_comparison(out, rows, results, spec["mass_band"], spec["raw"].get("plot_log_mass_bounds"))
    if progress:
        progress(f"comparison: {out / 'comparison.html'}")
    return out, bundle
