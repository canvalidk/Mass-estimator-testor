"""Deterministic pressure study of one measured pair, without changing the estimator.

The diagnostic projection is derived directly from the declared covariance, separately
from met.full's candidate construction. Corrections mean observed minus compatible;
there is no planted mass or known actual noise in this fixed-reading experiment.
"""
import csv
import datetime as dt
import hashlib
from importlib import metadata
import json
import math
from pathlib import Path

import numpy as np
import yaml
from scipy.special import i0e, i1e

from .full import Reading, isotropic, ours29


def residual_profile(reading, masses):
    """Closest pair satisfying F*=m a*, at each finite positive mass.

    With R=[I,-mI], r=R y and V=R Sigma R.T, the required correction is
    delta=Sigma R.T V^-1 r. Its Mahalanobis cost is r.T V^-1 r. This profiles
    over an unrestricted latent acceleration; declared subspaces need extra constraints.
    No density normalization, tempering or goodness-of-fit calibration is applied.
    """
    if reading.k != reading.dimension:
        raise ValueError("residual_profile requires the full candidate subspace")
    masses = np.asarray(masses, dtype=float)
    if masses.ndim != 1 or not masses.size or not np.all(np.isfinite(masses) & (masses > 0)):
        raise ValueError("masses must be a nonempty vector of finite positive values")
    d = reading.dimension
    eye = np.broadcast_to(np.eye(d), (len(masses), d, d))
    constraint = np.concatenate((eye, -masses[:, None, None] * eye), axis=2)
    residual = reading.force[None] - masses[:, None] * reading.acceleration[None]
    cov = constraint @ reading.cov @ constraint.transpose(0, 2, 1)
    score = np.linalg.solve(cov, residual[..., None])[..., 0]
    correction = np.einsum("nij,nj->ni", reading.cov @ constraint.transpose(0, 2, 1), score)
    force_correction, acceleration_correction = correction[:, :d], correction[:, d:]
    return dict(residual=residual, residual_covariance=cov,
                force_correction=force_correction, acceleration_correction=acceleration_correction,
                candidate_force=reading.force[None] - force_correction,
                candidate_acceleration=reading.acceleration[None] - acceleration_correction,
                correction_cost=np.einsum("ni,ni->n", residual, score))


def isotropic_minima(reading):
    """Exact cost bounds on the projective line and on the positive half-line.

    The domain includes both limiting endpoints, without assigning either a finite
    fitted mass. Projective cost is the squared smallest singular value of [F/sF,a/sa].
    On m>0 an interior minimum exists when the whitened dot product is positive.
    The positive-domain penalty is the extra infimum cost caused by forbidding
    negative slopes. It is a geometric diagnostic, not a replacement for Q in g.
    """
    d, cov = reading.dimension, reading.cov
    if reading.k != d:
        raise ValueError("isotropic_minima requires the full candidate subspace")
    sf2, sa2 = float(cov[0, 0]), float(cov[d, d])
    expected = np.diag(np.r_[np.full(d, sf2), np.full(d, sa2)])
    if not np.allclose(cov, expected, rtol=1e-12, atol=0):
        raise ValueError("isotropic_minima requires independent isotropic channel errors")
    x, y = reading.force / math.sqrt(sf2), reading.acceleration / math.sqrt(sa2)
    p, q, dot = float(x @ x), float(y @ y), float(x @ y)
    singular = np.linalg.svd(np.column_stack((x, y)), compute_uv=False)
    projective = float(singular[-1] ** 2) if d >= 2 else 0.0
    mass = None
    if dot > 0:
        gap = math.hypot(p - q, 2 * dot)
        slope = ((p - q + gap) / (2 * dot) if p >= q else 2 * dot / (q - p + gap))
        mass = math.sqrt(sf2 / sa2) * slope
        cost = float(residual_profile(reading, [mass])["correction_cost"][0])
        location = "interior"
    elif p < q:
        cost, location = p, "mass_to_zero"
    elif q < p:
        cost, location = q, "mass_to_infinity"
    else:
        cost = p
        location = "all_masses" if dot == 0 else "both_boundaries"
    supremum = float(singular[0] ** 2) if dot < 0 else max(p, q)
    # Keep small contrast/penalty even when the shared baseline dwarfs it.
    contrast = math.hypot(p - q, 2 * dot) + abs(p - q)
    penalty = abs(dot) * (2 * abs(dot) / contrast) if dot < 0 else 0.0
    return dict(projective_min_cost=projective, positive_min_cost=cost,
                positive_min_location=location, positive_min_mass=mass,
                positive_cost_supremum=supremum,
                positive_cost_span=0.5 * contrast,
                positive_domain_penalty=penalty)


def candidate_mean_length(reading, masses):
    """Untruncated E[|a*| | m] for this study's isotropic 2D cloud (Rice mean).

    The estimator's A array contains computational cutoff zeros outside retained
    density support. Those zeros are not cloud sizes, so diagnostic plots use this
    independent closed form on the entire exported grid instead.
    """
    if reading.dimension != 2:
        raise ValueError("candidate_mean_length requires a 2D reading")
    isotropic_minima(reading)  # also checks the full subspace and isotropic covariance
    masses = np.asarray(masses, dtype=float)
    if masses.ndim != 1 or not masses.size or not np.all(np.isfinite(masses) & (masses > 0)):
        raise ValueError("masses must be a nonempty vector of finite positive values")
    sf2, sa2 = reading.cov[0, 0], reading.cov[2, 2]
    denominator = sf2 + masses ** 2 * sa2
    mu = (masses[:, None] * sa2 * reading.force + sf2 * reading.acceleration) / denominator[:, None]
    tau = np.sqrt(sf2 / denominator) * math.sqrt(sa2)
    z = np.sum(mu ** 2, axis=1) / (4 * tau ** 2)
    return tau * math.sqrt(math.pi / 2) * ((1 + 2 * z) * i0e(z) + 2 * z * i1e(z))


def _positive(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return float(value)


def validate(data):
    """A separate preset kind: fixed observations, not draws from a Newtonian world."""
    required = {"kind", "name", "question", "seed", "cases", "noise", "prior"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"misalignment preset needs {sorted(required)}")
    if set(data) - required - {"notes", "numerics", "plots"}:
        raise ValueError("unknown misalignment preset setting")
    if data["kind"] != "misalignment":
        raise ValueError("expected kind: misalignment")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for ch in name)):
        raise ValueError("preset.name must be a simple identifier, not a path")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("preset.question must be nonempty")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("preset.seed must be a nonnegative integer (no random draws are used)")
    if type(data.get("plots", True)) is not bool:
        raise ValueError("plots must be true or false")
    prior = data["prior"]
    if not isinstance(prior, dict) or set(prior) != {"family", "center"} or prior["family"] != "half_cauchy":
        raise ValueError("prior must declare family: half_cauchy and center")
    _positive(prior["center"], "prior.center")
    noise = data["noise"]
    if not isinstance(noise, dict) or set(noise) != {"force_sd", "acceleration_sd", "scales"}:
        raise ValueError("noise needs force_sd, acceleration_sd and scales")
    sf = _positive(noise["force_sd"], "noise.force_sd")
    sa = _positive(noise["acceleration_sd"], "noise.acceleration_sd")
    scales = noise["scales"]
    if not isinstance(scales, list) or not scales:
        raise ValueError("noise.scales must be a nonempty list")
    scales = [_positive(v, "noise scale") for v in scales]
    if len(set(scales)) != len(scales):
        raise ValueError("noise scales must be distinct")
    numerics = dict(profile_log_half_span=16.0, profile_nodes=1601)
    if not isinstance(data.get("numerics", {}), dict) or set(data.get("numerics", {})) - numerics.keys():
        raise ValueError("numerics accepts profile_log_half_span and profile_nodes")
    numerics.update(data.get("numerics", {}))
    _positive(numerics["profile_log_half_span"], "profile_log_half_span")
    if numerics["profile_log_half_span"] > 40:
        raise ValueError("profile_log_half_span must be <= 40")
    if type(numerics["profile_nodes"]) is not int or numerics["profile_nodes"] < 101:
        raise ValueError("profile_nodes must be an integer >= 101")
    if not isinstance(data["cases"], list) or not data["cases"]:
        raise ValueError("cases must be a nonempty list")
    names, cells = set(), []
    for case in data["cases"]:
        if not isinstance(case, dict) or set(case) != {"name", "force", "acceleration"}:
            raise ValueError("each case needs name, force and acceleration")
        if not isinstance(case["name"], str) or not case["name"].strip() or case["name"] in names:
            raise ValueError("case names must be nonempty and distinct")
        names.add(case["name"])
        force, acceleration = np.asarray(case["force"], dtype=float), np.asarray(case["acceleration"], dtype=float)
        if (force.shape != (2,) or acceleration.shape != (2,)
                or not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration))):
            raise ValueError("fixed force and acceleration must be finite 2D vectors")
        for scale in scales:
            for sd in (sf * scale, sa * scale):
                if not np.isfinite(sd) or sd > math.sqrt(np.finfo(float).max) or sd < math.sqrt(np.finfo(float).tiny):
                    raise ValueError("scaled noise variance must remain finite and positive in float64")
            cells.append(dict(id=f"c{len(cells):03d}", case=case["name"], force=force.tolist(),
                              acceleration=acceleration.tolist(), noise_scale=scale,
                              force_sd=sf * scale, acceleration_sd=sa * scale))
    return dict(kind="misalignment", raw=data, cells=cells, numerics=numerics)


def _density(result):
    """Convert the estimator's discrete log-grid weights to density per unit ln m."""
    return result["law"] / np.trapezoid(result["law"], result["lm"])


def _point_summary(reading, result):
    projection = residual_profile(reading, [result["readout"]])
    density = _density(result)
    # met.full currently searches only ln(c)+/-25. Preserve evidence of truncated
    # tails rather than silently trusting a finite endpoint as a physical boundary.
    edge_fraction = float(max(density[0], density[-1]) / density.max())
    return dict(point=result["readout"], median=result["median"], interval=list(result["interval"]),
                cost_at_point=float(projection["correction_cost"][0]),
                force_correction_at_point=projection["force_correction"][0].tolist(),
                acceleration_correction_at_point=projection["acceleration_correction"][0].tolist(),
                support_edge_density_fraction=edge_fraction,
                support_edge_warning=edge_fraction > 1e-7,
                log_grid_min=float(result["lm"][0]), log_grid_max=float(result["lm"][-1]))


def investigate_cell(cell, c, numerics):
    reading = isotropic(cell["force"], cell["acceleration"], cell["force_sd"], cell["acceleration_sd"])
    raw = ours29([reading], c, misfit=False, calibrate=False)
    quoted = ours29([reading], c)
    lm = math.log(c) + np.linspace(-numerics["profile_log_half_span"],
                                  numerics["profile_log_half_span"], numerics["profile_nodes"])
    arrays = dict(profile_lm=lm, profile_mass=np.exp(lm), **residual_profile(reading, np.exp(lm)))
    bounds = isotropic_minima(reading)
    arrays["positive_excess_cost"] = arrays["correction_cost"] - bounds["positive_min_cost"]
    arrays["prior_density"] = (2 / math.pi) * np.exp(math.log(c) + lm - np.logaddexp(2 * math.log(c), 2 * lm))
    for label, result in (("raw", raw), ("quoted", quoted)):
        arrays.update({f"{label}_lm": result["lm"], f"{label}_density": _density(result),
                       f"{label}_A": candidate_mean_length(reading, np.exp(result["lm"]))})
    f, a = reading.force, reading.acceleration
    norm = np.linalg.norm(f) * np.linalg.norm(a)
    angle = math.degrees(math.acos(float(np.clip(f @ a / norm, -1, 1)))) if norm > 0 else None
    summary = dict(cell, **bounds, angle_degrees=angle,
                   equation_Q=quoted["Q"][0], misfit_weight=quoted["g"][0], bartlett_eta=quoted["eta"],
                   raw=_point_summary(reading, raw), quoted=_point_summary(reading, quoted))
    for label, result in (("raw", raw), ("quoted", quoted)):
        summary[label]["log_evidence_span"] = (
            0.5 * result["eta"] * result["g"][0] * bounds["positive_cost_span"])
    return summary, arrays


def _summary(data):
    lines = [f"# {data['name']}", "", data["question"], "",
             "Deterministic fixed observations; no planted mass, random noise draws or coverage claim.",
             f"Half-Cauchy prior centre fixed at {data['prior']['center']:g}. Seed {data['seed']} (unused).",
             "Raw = neither misfit weighting nor Bartlett power. Quoted = current ours29.",
             "Corrections are observed minus the closest compatible candidate, not recovered actual errors.",
             "Costs are squared Mahalanobis distances, not p-values or marginal evidence.", "",
             "| case | noise scale | positive cost infimum | location | projective Q | g | raw point | quoted point | quoted median | quoted 95% interval |",
             "|---|---:|---:|---|---:|---:|---:|---:|---:|---|"]
    for cell in data["cells"]:
        q = cell["quoted"]
        lines.append(f"| {cell['case']} | {cell['noise_scale']:g} | {cell['positive_min_cost']:.6g} | "
                     f"{cell['positive_min_location']} | {cell['projective_min_cost']:.6g} | "
                     f"{cell['misfit_weight']:.6g} | {cell['raw']['point']:.6g} | {q['point']:.6g} | "
                     f"{q['median']:.6g} | [{q['interval'][0]:.6g}, {q['interval'][1]:.6g}] |")
    lines += ["", "`mass_to_zero` and `mass_to_infinity` are limiting infima, never finite fitted masses.",
              "Normalized laws show relative preference between masses; absolute correction costs are kept separately.",
              "The estimator searches the whole projective line for Q, including negative slopes.",
              "", "| case | noise scale | extra cost from positive mass | quoted log-evidence span |",
              "|---|---:|---:|---:|"]
    for cell in data["cells"]:
        lines.append(f"| {cell['case']} | {cell['noise_scale']:g} | "
                     f"{cell['positive_domain_penalty']:.6g} | {cell['quoted']['log_evidence_span']:.6g} |")
    lines += ["", "The extra cost is Q_positive - Q_projective; it is not used to alter g.",
              "Log-evidence span is sup(log(p/prior)) - inf(log(p/prior)) on m > 0, in nats.",
              "Normalization cancels the mass-independent positive minimum; its excess curve is in values.npz.",
              "No pass/fail standard for the equation's response is inferred from these plots.",
              "Full-precision values, support-edge warnings and code/environment provenance are in run.json;",
              "all curves and correction vectors are in values.npz. Plots use density per unit ln(m).", ""]
    return "\n".join(lines)


def _write_csv(path, cells):
    rows = []
    for cell in cells:
        row = {k: cell[k] for k in ("id", "case", "angle_degrees", "noise_scale", "force_sd", "acceleration_sd",
                                   "positive_min_cost", "positive_min_location", "positive_min_mass",
                                   "projective_min_cost", "positive_cost_supremum", "positive_cost_span",
                                   "positive_domain_penalty", "equation_Q", "misfit_weight", "bartlett_eta")}
        for name in ("raw", "quoted"):
            for key in ("point", "median", "cost_at_point", "support_edge_warning", "log_evidence_span"):
                row[f"{name}_{key}"] = cell[name][key]
            row[f"{name}_interval_low"], row[f"{name}_interval_high"] = cell[name]["interval"]
        rows.append(row)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(path, out_root="results", progress=print):
    """Run through the same `python -m met run PRESET` entry point as other studies."""
    from .runner import code_version
    with Path(path).open(encoding="utf-8") as handle:
        spec = validate(yaml.safe_load(handle))
    raw = spec["raw"]
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d-%H%M%S")
    parent = Path(out_root) / raw["name"]
    out = parent / stamp
    suffix = 1
    while out.exists():
        suffix += 1
        out = parent / f"{stamp}-{suffix}"
    cells, curves = [], {}
    for cell in spec["cells"]:
        summary, values = investigate_cell(cell, float(raw["prior"]["center"]), spec["numerics"])
        cells.append(summary)
        curves.update({f"{cell['id']}/{k}": v for k, v in values.items()})
        if progress:
            progress(f"{cell['id']} {cell['case']}, noise {cell['noise_scale']:g}: "
                     f"positive cost {summary['positive_min_cost']:.5g}, "
                     f"projective Q {summary['projective_min_cost']:.5g}, g {summary['misfit_weight']:.5g}")
    data = dict(schema_version=1, kind="misalignment", name=raw["name"], question=raw["question"],
                seed=raw["seed"], deterministic=True, prior=raw["prior"], numerics=spec["numerics"],
                code_version=code_version(), started_utc=stamp, cells=cells,
                law_definitions={"raw": "Unweighted baseline (ours27): misfit=False, calibrate=False",
                                 "quoted": "Current ours29: misfit=True, calibrate=True"})
    versions = {}
    for package in ("numpy", "scipy", "PyYAML", "matplotlib", "contourpy", "cycler", "fonttools",
                    "kiwisolver", "pillow", "pyparsing", "python-dateutil", "six", "tzdata"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    data["code_version"]["dependency_versions"] = versions
    data["code_version"]["misalignment_preset_sha256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    out.mkdir(parents=True)
    (out / "preset.yaml").write_bytes(Path(path).read_bytes())
    np.savez_compressed(out / "values.npz", **curves)
    _write_csv(out / "summary.csv", cells)
    (out / "summary.md").write_text(_summary(data), encoding="utf-8")
    data["plots_generated"] = False
    if raw.get("plots", True):
        from .investigation_report import write_report
        data["plots_generated"] = bool(write_report(out, data, curves))
        if not data["plots_generated"] and progress:
            progress("Plots unavailable: install requirements-plots.txt; all numeric outputs were saved.")
    (out / "run.json").write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")
    if progress:
        progress(f"wrote {out}")
    return out, data
