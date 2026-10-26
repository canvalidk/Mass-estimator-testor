"""A replaceable investigation recipe composed from the general case inspector.

This first recipe asks what common noise scaling changes. Path construction is
general; this recipe's extra analytic summaries require one unrestricted isotropic
reading. Nothing here changes the estimator or declares physical adequacy.
"""
import csv
import hashlib
from html import escape
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import numpy as np
import yaml
from scipy.integrate import IntegrationWarning, quad

from .full import Reading
from .inspection.case import parse_cases
from .inspection.paths import expand_paths
from .misalignment import isotropic_minima


def _reading(observation):
    return Reading(observation.force, observation.acceleration, observation.covariance, observation.subspace)


def validate(data):
    cases = parse_cases(data)
    values = data.get("noise_scales")
    if not isinstance(values, list):
        raise ValueError("noise_path needs noise_scales")
    if type(data.get("seed")) is not int or data["seed"] < 0:
        raise ValueError("noise_path needs a nonnegative seed (no random draws)")
    band = data.get("mass_band", [0.5, 4.0])
    if (not isinstance(band, list) or len(band) != 2
            or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in band)
            or not 0 < band[0] < band[1]):
        raise ValueError("mass_band must be two increasing finite positive masses")
    for case in cases:
        if len(case.inputs.readings) != 1:
            raise ValueError("this noise-path summary requires one reading; the path builder supports multiple")
        isotropic_minima(_reading(case.inputs.readings[0]))
        if case.inputs.prior.get("family") != "half_cauchy":
            raise ValueError("this noise-path reference requires a half-Cauchy prior")
    expanded = expand_paths(cases, values)
    detail = data.get("detail_scales", values)
    if not isinstance(detail, list) or any(value not in values for value in detail):
        raise ValueError("detail_scales must select values from noise_scales")
    display = data.get("plot_log_mass_bounds")
    if display is not None and (not isinstance(display, list) or len(display) != 2
            or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in display)
            or display[0] >= display[1]):
        raise ValueError("plot_log_mass_bounds must be two increasing finite log-masses")
    return dict(kind="noise_path", raw=data, cases=expanded, mass_band=band,
                cells=[dict(id=f"case_{i:03d}", name=c.name) for i, c in enumerate(expanded)])


def angle_reference(observation, prior_center, weight, masses):
    """Independent adaptive integrals on the compact physical angle domain.

    m=s*tan(phi), s=sigma_F/sigma_a. The direct residual cost is
    ||cos(phi) F/sigma_F - sin(phi) a/sigma_a||². The half-Cauchy
    Jacobian is c*s/(c² cos²(phi)+s² sin²(phi)). Subtracting the
    positive-domain infimum prevents underflow of the shared constant.
    Q/g/eta are supplied by the inspected equation, not justified by this oracle.
    """
    reading = _reading(observation)
    bounds = isotropic_minima(reading)
    if not math.isfinite(weight) or weight < 0 or not math.isfinite(prior_center) or prior_center <= 0:
        raise ValueError("reference needs a nonnegative weight and positive prior center")
    masses = np.asarray(masses, dtype=float)
    if masses.ndim != 1 or not np.all(np.isfinite(masses) & (masses > 0)):
        raise ValueError("reference masses must be a finite positive vector")
    if weight == 0 or bounds["positive_cost_span"] == 0:
        # Constant evidence has a closed-form CDF. This also avoids quadrature of
        # an arbitrarily narrow prior Jacobian on a null reading.
        return dict(cdf=2 / math.pi * np.arctan2(masses, prior_center), normalizer=math.pi / 2,
                    quadrature_error_bound=0.0, quadrature_warnings=[],
                    note="Analytic half-Cauchy CDF: the weighted correction cost is constant in mass.")
    d = reading.dimension
    sf, sa = math.sqrt(reading.cov[0, 0]), math.sqrt(reading.cov[d, d])
    ratio, c = sf / sa, prior_center
    f, a = reading.force / sf, reading.acceleration / sa

    def kernel(phi):
        co, si = math.cos(phi), math.sin(phi)
        residual = co * f - si * a
        excess = max(0.0, float(residual @ residual) - bounds["positive_min_cost"])
        jacobian = c * ratio / (c * c * co * co + ratio * ratio * si * si)
        return jacobian * math.exp(-0.5 * weight * excess)

    # Supply narrow peak and boundary scales to the adaptive integrator. Repeating
    # integrals on this fixed partition avoids missing a narrow peak in a long CDF.
    energy = float(f @ f + a @ a)
    width = min(0.1, 1 / math.sqrt(1 + weight * energy))
    edge_width = min(0.1, 1 / (1 + weight * energy))
    peak = math.atan(bounds["positive_min_mass"] / ratio) if bounds["positive_min_mass"] else math.pi / 4
    cuts = [0.0, math.pi / 2, math.pi / 4, peak]
    for multiplier in (1, 4, 16):
        cuts.extend([edge_width * multiplier, math.pi / 2 - edge_width * multiplier,
                     peak - width * multiplier, peak + width * multiplier])
        # Resolve the prior Jacobian too, including weak/null signals whose
        # likelihood supplies no scale. atan2 avoids overflow in c/s ratios.
        cuts.extend([math.atan2(multiplier * c, ratio),
                     math.atan2(c, multiplier * ratio)])
    cuts = sorted(set(x for x in cuts if 0 <= x <= math.pi / 2))

    integration_warnings = []

    def integral(lo, hi):
        points = [lo, *[x for x in cuts if lo < x < hi], hi]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", IntegrationWarning)
            values = [quad(kernel, x, y, epsabs=1e-13, epsrel=2e-11, limit=100)
                      for x, y in zip(points[:-1], points[1:]) if x < y]
        integration_warnings.extend(str(w.message) for w in caught)
        return sum(x[0] for x in values), sum(x[1] for x in values)

    z, z_error = integral(0, math.pi / 2)
    if z <= 0:
        raise ValueError("angle reference has no resolved normalizer")
    probabilities, errors = [], []
    for mass in masses:
        phi = math.atan(float(mass) / ratio)
        # Integrate the shorter probability tail to avoid integrating across an
        # almost empty central region in a boundary-concentrated law.
        if phi <= math.pi / 4:
            value, error = integral(0, phi)
            probabilities.append(value / z)
        else:
            value, error = integral(phi, math.pi / 2)
            probabilities.append(1 - value / z)
        errors.append(error / z + value * z_error / z ** 2)
    return dict(cdf=np.array(probabilities), normalizer=z,
                quadrature_error_bound=max(errors, default=0.0),
                quadrature_warnings=list(dict.fromkeys(integration_warnings)),
                note="QUADPACK estimated error, not a rigorous bound; independent of the native log-mass grid.")


def summarize(result, mass_band):
    from .inspection.case import Case
    case = Case.from_dict(result["case"])
    observation = case.inputs.readings[0]
    bounds = isotropic_minima(_reading(observation))
    path = case.metadata["path"]
    row = dict(case=case.name, path=path["source"], noise_scale=path["value"], **bounds)
    row["noise_sd_force"] = math.sqrt(observation.covariance[0, 0])
    f, a = observation.force, observation.acceleration
    norm = np.linalg.norm(f) * np.linalg.norm(a)
    row["angle_degrees"] = math.degrees(math.acos(float(np.clip(f @ a / norm, -1, 1)))) if norm > 0 else None
    row["Q"] = result["evaluation"]["per_reading"][0]["Q"]
    row["g"] = result["evaluation"]["per_reading"][0]["g"]
    fixed = case.inputs.prior["center"]
    from .inspection.geometry import profile
    fit = profile(observation, [fixed])
    row.update(fixed_candidate_mass=fixed, fixed_candidate_cost=float(fit["declared_cost"][0]),
               fixed_force_correction_norm=float(np.linalg.norm(fit["force_correction"][0])),
               fixed_acceleration_correction_norm=float(np.linalg.norm(fit["acceleration_correction"][0])),
               fixed_conditional_covariance_trace=float(np.trace(fit["conditional_covariance"][0])))
    row["laws"] = {}
    for name, law in result["evaluation"]["laws"].items():
        x, density = law["log_mass"], law["density"]
        cdf = np.r_[0, np.cumsum((density[:-1] + density[1:]) * np.diff(x) / 2)]
        cdf /= cdf[-1]
        weight = law["eta"] * law["g"][0]
        queries = [*mass_band, *law["interval"]]
        reference = angle_reference(observation, case.inputs.prior["center"], weight, queries)
        native_band = float(np.diff(np.interp(np.log(mass_band), x, cdf))[0])
        oracle_band = float(reference["cdf"][1] - reference["cdf"][0])
        interval_error = float(np.max(np.abs(reference["cdf"][2:] - [0.025, 0.975])))
        row["laws"][name] = dict(point=law["point"], median=law["median"], interval=law["interval"],
            eta=law["eta"], curvature=law["curvature"] if math.isfinite(law["curvature"]) else None,
            J=law["J"] if math.isfinite(law["J"]) else None,
            log_evidence_span=weight * bounds["positive_cost_span"] / 2,
            band_probability=native_band, reference_band_probability=oracle_band,
            band_absolute_error=abs(native_band - oracle_band), interval_cdf_error=interval_error,
            quadrature_error_bound=reference["quadrature_error_bound"],
            quadrature_warnings=reference["quadrature_warnings"],
            support_edge_warning=law["support_edge_warning"],
            support_edge_density_fraction=law["support_edge_density_fraction"])
    return row


def compare_paths(rows, results):
    """Numerical path identities, separate from any required physical response."""
    from .inspection.checks import _law_errors
    comparisons = []
    for group in dict.fromkeys(row["path"] for row in rows):
        pairs = [(row, result) for row, result in zip(rows, results) if row["path"] == group]
        base_row, base = pairs[0]
        baseline = base["diagnostics"]["readings"][0]
        correction_error, covariance_error, cost_error = 0.0, 0.0, 0.0
        for row, result in pairs[1:]:
            ratio = (row["noise_scale"] / base_row["noise_scale"]) ** 2
            fit = result["diagnostics"]["readings"][0]
            for key in ("force_correction", "acceleration_correction"):
                correction_error = max(correction_error, float(np.max(np.abs(fit[key] - baseline[key]))))
            covariance_error = max(covariance_error, float(np.max(np.abs(
                fit["conditional_covariance"] / ratio - baseline["conditional_covariance"]))))
            cost_error = max(cost_error, float(np.max(np.abs(fit["declared_cost"] * ratio - baseline["declared_cost"])))
                            / (1 + float(np.max(baseline["declared_cost"]))))
        saturated = [(r, x) for r, x in pairs
                     if x["evaluation"]["per_reading"][0]["nu"] > 0
                     and r["Q"] >= x["evaluation"]["per_reading"][0]["nu"]
                     and r["laws"]["quoted"]["eta"] == 1]
        frozen = None
        if len(saturated) >= 2:
            errors = [_law_errors(saturated[0][1]["evaluation"]["laws"]["quoted"],
                                  other["evaluation"]["laws"]["quoted"]) for _, other in saturated[1:]]
            frozen = dict(scales=[r["noise_scale"] for r, _ in saturated],
                          maximum_log_kernel_error=max(e["sampled_log_kernel_error_modulo_constant"] for e in errors),
                          maximum_native_cdf_error=max(e["native_cdf_maximum_error"] for e in errors),
                          point_range=[min(r["laws"]["quoted"]["point"] for r, _ in saturated),
                                       max(r["laws"]["quoted"]["point"] for r, _ in saturated)])
        comparisons.append(dict(path=group, maximum_correction_change=correction_error,
            maximum_rescaled_covariance_error=covariance_error, maximum_rescaled_cost_error=cost_error,
            saturated_law_comparison=frozen,
            authority="Common covariance scaling leaves the least-squares minimizer fixed; K scales inversely. For one reading with nu>0, saturated g=nu/Q cancels that scale and eta=1. No physical adequacy verdict."))
    return comparisons


def _write_comparison(out, rows, results, band, plot_log_mass_bounds=None):
    """This recipe's presentation; the underlying case/report format stays generic."""
    flat = []
    for row in rows:
        for name, law in row["laws"].items():
            flat.append(dict({k: v for k, v in row.items() if k != "laws"}, law=name,
                             **{k: v for k, v in law.items() if k != "interval"},
                             interval_low=law["interval"][0], interval_high=law["interval"][1]))
    with (out / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    relations = compare_paths(rows, results)
    (out / "comparison.json").write_text(json.dumps(dict(mass_band=band, rows=rows, path_relations=relations), indent=2, allow_nan=False), encoding="utf-8")
    figure_html = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        figure_html.append("<p>Matplotlib unavailable; all comparisons are in CSV/JSON.</p>")
    else:
        groups = list(dict.fromkeys(row["path"] for row in rows))
        for mode in ("laws", "intervals"):
            fig, axes = plt.subplots(math.ceil(len(groups) / 2), 2, figsize=(13, 3.2 * math.ceil(len(groups) / 2)), squeeze=False)
            for ax, group in zip(axes.flat, groups):
                selected = [(row, result) for row, result in zip(rows, results) if row["path"] == group]
                if mode == "laws":
                    for row, result in selected:
                        law = result["evaluation"]["laws"]["quoted"]
                        ax.plot(law["log_mass"], law["density"], label=f'SD {row["noise_sd_force"]:g}', linewidth=1.3)
                    result = selected[0][1]
                    ax.plot(np.log(result["evaluation"]["mass"]), np.exp(result["evaluation"]["prior_log_density"]), "k--", label="prior", alpha=.7)
                    limits = plot_log_mass_bounds or [
                        min(float(x[1]["evaluation"]["laws"]["quoted"]["log_mass"][0]) for x in selected),
                        max(float(x[1]["evaluation"]["laws"]["quoted"]["log_mass"][-1]) for x in selected)]
                    ax.set(xlabel="ln(mass)", ylabel="Density per ln(mass)", xlim=limits)
                else:
                    sd = [r["noise_sd_force"] for r, _ in selected]
                    laws = [r["laws"]["quoted"] for r, _ in selected]
                    ax.fill_between(sd, [l["interval"][0] for l in laws], [l["interval"][1] for l in laws], alpha=.18, label="95% equal-tailed interval")
                    ax.plot(sd, [l["point"] for l in laws], "o-", label="point")
                    ax.plot(sd, [l["median"] for l in laws], "x--", label="numerical median")
                    ax.set(xscale="log", yscale="log", xlabel="Declared common SD (decreasing →)", ylabel="Mass")
                    ax.invert_xaxis()
                ax.set_title(group)
                ax.grid(alpha=.15)
                ax.legend(fontsize=7, ncol=2)
            for ax in list(axes.flat)[len(groups):]:
                ax.set_visible(False)
            fig.tight_layout()
            for suffix in ("png", "svg"):
                fig.savefig(out / f"{mode}.{suffix}", dpi=140)
            plt.close(fig)
            figure_html.append("<svg" + (out / f"{mode}.svg").read_text(encoding="utf-8").split("<svg", 1)[1])
        # Embed SVG rather than external assets so the comparison can be shared.
    header = ("Path", "SD", "Q+", "Q", "g", "η", "log-evidence span", "Point", "95% interval", f"P({band[0]:g} < m < {band[1]:g})", "Detail")
    table = ["<tr>" + "".join(f"<th>{escape(x)}</th>" for x in header) + "</tr>"]
    for index, row in enumerate(rows):
        law = row["laws"]["quoted"]
        values = [row["path"], f'{row["noise_sd_force"]:g}', f'{row["positive_min_cost"]:.6g}',
                  f'{row["Q"]:.6g}', f'{row["g"]:.6g}', f'{law["eta"]:.6g}',
                  f'{law["log_evidence_span"]:.6g}', f'{law["point"]:.6g}',
                  f'[{law["interval"][0]:.6g}, {law["interval"][1]:.6g}]', f'{law["band_probability"]:.6g}']
        table.append("<tr>" + "".join(f"<td>{escape(x)}</td>" for x in values)
                     + f'<td><a href="case_{index:03d}/report.html">Inspect</a></td></tr>')
    max_band = max(law["band_absolute_error"] for row in rows for law in row["laws"].values())
    max_interval = max(law["interval_cdf_error"] for row in rows for law in row["laws"].values())
    warnings = sum(law["support_edge_warning"] for row in rows for law in row["laws"].values())
    reference_warnings = sum(bool(law["quadrature_warnings"]) for row in rows for law in row["laws"].values())
    html = ['<!doctype html><meta charset="utf-8"><title>Shrinking-noise paths</title>',
            '<style>body{font:16px/1.5 system-ui;margin:36px auto;max-width:1300px;padding:20px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:right}td:first-child{text-align:left}svg{width:100%;height:auto}.scroll{overflow:auto}</style>',
            '<h1>What changes as declared noise shrinks?</h1>',
            '<p>Fixed measured vectors and prior; one reading per case. Every SD scales both instruments together. '
            'This is a deterministic pressure experiment, without a planted mass or repeated-experiment coverage claim.</p>',
            '<p>Q+ is the least correction cost over positive masses, including unattained boundary infima. '
            'Q searches the whole projective line. Log-evidence span measures relative preference after g and η; '
            'it excludes the prior and the mass-independent cost. Costs are not p-values.</p>',
            '<p>Inspect the full laws: stable points can coexist with frozen, narrowing, or escaping probability. '
            'Equal-tailed intervals can bridge a low-probability gap. Numerical medians on CDF plateaus are retained as supplied.</p>',
            '<p>Law plot display: ' + (f'explicit ln(mass) bounds {escape(str(plot_log_mass_bounds))}; '
            'tails outside this display remain in the saved native arrays.' if plot_log_mass_bounds else
            'limits span the saved native grids for each path.') + '</p>',
            f'<p>Independent compact-angle integration: largest band-probability difference {max_band:.3g}; '
            f'largest interval-endpoint CDF error {max_interval:.3g}. Native-grid support warnings: {warnings}; '
            f'reference quadrature warning stages: {reference_warnings} (details in comparison.json). '
            'These checks verify numerical integration for these cases, not the physical adequacy of the equation.</p>',
            '<p><a href="report.html">All case checks</a> · <a href="comparison.csv">CSV</a> · '
            '<a href="comparison.json">Full comparison</a> · <a href="recipe.yaml">Path recipe</a> · '
            '<a href="preset.yaml">Expanded replay cases</a></p>',
            '<div class="scroll"><table>', *table, '</table></div>',
            '<details><summary>Path identities and numerical differences</summary><pre>'
            + escape(json.dumps(relations, indent=2)) + '</pre></details>', *figure_html]
    (out / "comparison.html").write_text("\n".join(html), encoding="utf-8")


def run(path, out_root="results", progress=print, *, reporter=None):
    from .inspection.report import write_report
    from .inspection.run import run as inspect_run
    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    document = dict(spec["raw"], kind="inspection", cases=[c.to_dict(include_truth=True) for c in spec["cases"]])
    rows, results = [], []
    reporter = reporter or write_report

    def collect(folder, result):
        rows.append(summarize(result, spec["mass_band"]))
        results.append(result)
        if reporter is write_report:
            plots = result["case"]["metadata"]["path"]["value"] in spec["raw"].get("detail_scales", spec["raw"]["noise_scales"])
            return reporter(folder, result, plots=plots)
        return reporter(folder, result)

    # The inspector saves its exact expanded input before this temporary input is
    # removed. The original recipe and its hash travel beside that explicit replay.
    with TemporaryDirectory(prefix="met-noise-path-") as temp:
        expanded = Path(temp) / "expanded.yaml"
        expanded.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        out, bundle = inspect_run(expanded, out_root, progress=progress, reporter=collect)
    (out / "recipe.yaml").write_bytes(source)
    bundle["path_recipe"] = dict(sha256=hashlib.sha256(source).hexdigest(), file="recipe.yaml", transformation="met.inspection.paths.scale_noise")
    (out / "run.json").write_text(json.dumps(bundle, indent=2, allow_nan=False), encoding="utf-8")
    _write_comparison(out, rows, results, spec["mass_band"], spec["raw"].get("plot_log_mass_bounds"))
    if progress:
        progress(f"comparison: {out / 'comparison.html'}")
    return out, bundle
