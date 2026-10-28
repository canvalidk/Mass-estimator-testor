"""O17's declared tensor worlds, independent scalar check and unchanged ours29.

Tensor and latent vectors are investigator-only. The estimator boundary accepts a
ReadingSet and an explicitly supplied prior centre. No scalar accuracy target is
invented for worlds whose excited directions admit no common scalar.
"""
from datetime import datetime
from functools import lru_cache
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.stats import chi2
import yaml

from .analyst import bartlett_weight, joint_log_density, misfit_weights, posterior_summaries, prior_log_density
from .inspection.sources import from_reading_set
from .law import ReadingSet
from .world import series_rngs

GEOMETRIES = ("principal_line", "principal_mix", "oblique_line", "broad")
LAW = {"reference": "cartesian2", "nuisance": "radius"}


def _positive(value, label):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value <= 0):
        raise ValueError(f"{label} must be finite and positive")
    return float(value)


def _tensor(value):
    tensor = np.asarray(value, dtype=float)
    if (tensor.shape != (3, 3) or not np.all(np.isfinite(tensor))
            or not np.allclose(tensor, tensor.T, rtol=1e-12, atol=0)):
        raise ValueError("tensor must be a finite symmetric 3 by 3 matrix")
    if np.linalg.eigvalsh(tensor).min() <= 0:
        raise ValueError("tensor must be positive definite")
    return tensor


def validate(data):
    required = {"kind", "name", "question", "tensors", "geometries", "readings", "noise_sd",
                "replicates", "seed", "acceleration_length", "prior_center"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"anisotropic_world needs {sorted(required)}")
    if set(data) - required - {"notes", "numerics"}:
        raise ValueError("unknown anisotropic_world setting")
    if data["kind"] != "anisotropic_world":
        raise ValueError("expected kind: anisotropic_world")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in name)):
        raise ValueError("anisotropic_world name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("anisotropic_world needs a nonempty question")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("anisotropic_world seed must be a nonnegative integer")
    if type(data["replicates"]) is not int or data["replicates"] < 2:
        raise ValueError("anisotropic_world replicates must be an integer >= 2")
    for key in ("tensors", "geometries", "readings", "noise_sd"):
        if not isinstance(data[key], list) or not data[key]:
            raise ValueError(f"anisotropic_world {key} must be a nonempty list")
    names, tensors = set(), []
    for item in data["tensors"]:
        if not isinstance(item, dict) or set(item) != {"name", "matrix"}:
            raise ValueError("each tensor needs a name and matrix")
        if not isinstance(item["name"], str) or not item["name"].strip() or item["name"] in names:
            raise ValueError("tensor names must be nonempty and unique")
        names.add(item["name"])
        tensor = _tensor(item["matrix"])
        if not np.allclose(tensor, np.diag(np.diag(tensor)), rtol=1e-12, atol=0):
            raise ValueError("this recipe requires tensors diagonal in the common principal frame")
        tensors.append(dict(name=item["name"], matrix=tensor.tolist()))
    if any(not isinstance(g, str) or g not in GEOMETRIES for g in data["geometries"]):
        raise ValueError(f"geometries must be drawn from {GEOMETRIES}")
    if len(set(data["geometries"])) != len(data["geometries"]):
        raise ValueError("geometries must be unique")
    if any(type(n) is not int or n < 3 for n in data["readings"]):
        raise ValueError("reading counts must be integers >= 3")
    if len(set(data["readings"])) != len(data["readings"]):
        raise ValueError("reading counts must be unique")
    noise = [_positive(value, "noise_sd") for value in data["noise_sd"]]
    if len(set(noise)) != len(noise):
        raise ValueError("noise_sd values must be unique")
    numerics = dict(representative_log_half_span=8.0, representative_nodes=4001)
    supplied = data.get("numerics", {})
    if not isinstance(supplied, dict) or set(supplied) - numerics.keys():
        raise ValueError("anisotropic_world numerics accepts representative_log_half_span and representative_nodes")
    numerics.update(supplied)
    numerics["representative_log_half_span"] = _positive(numerics["representative_log_half_span"], "representative_log_half_span")
    if type(numerics["representative_nodes"]) is not int or numerics["representative_nodes"] < 101:
        raise ValueError("representative_nodes must be an integer >= 101")
    settings = dict(acceleration_length=_positive(data["acceleration_length"], "acceleration_length"),
                    prior_center=_positive(data["prior_center"], "prior_center"), **numerics)
    cells = [dict(id=f"c{i:03d}", tensor_name=tensor["name"], tensor=tensor["matrix"],
                  geometry=geometry, readings=n, noise_sd=sd)
             for i, (tensor, geometry, n, sd) in enumerate(itertools.product(
                 tensors, data["geometries"], data["readings"], noise))]
    return dict(kind="anisotropic_world", raw=data, settings=settings, cells=cells)


@lru_cache(maxsize=32)
def _combine(prior_center, stage):
    combine = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": prior_center}}]}
    if stage in ("untempered", "quoted"):
        combine["misfit"] = "kink"
    if stage == "quoted":
        combine["calibrate"] = "bartlett"
    return combine


def generate_cell(cell, settings, replicates, seed):
    """Pair latent directions and errors across tensors; G is investigator-only.

    The seed key deliberately omits G and prior centre: equal geometry/noise/count
    uses common random numbers. Principal-line worlds with G e2=2 e2 therefore
    generate identical measured inputs for the two shipped tensors.
    """
    n, sd = cell["readings"], cell["noise_sd"]
    tensor = _tensor(cell["tensor"])
    key = {name: cell[name] for name in ("geometry", "readings", "noise_sd")}
    key["acceleration_length"] = settings["acceleration_length"]
    latent, force, acceleration = (np.empty((replicates, n, 3)) for _ in range(3))
    true_force = np.empty_like(latent)
    for i, rng in enumerate(series_rngs(seed, key, replicates)):
        geometry = cell["geometry"]
        if geometry == "principal_line":
            directions = np.broadcast_to([0., 1., 0.], (n, 3))
        elif geometry == "principal_mix":
            directions = np.eye(3)[np.arange(n) % 3]
        elif geometry == "oblique_line":
            directions = np.broadcast_to(np.array([1., 1., 0.]) / math.sqrt(2), (n, 3))
        else:
            directions = rng.normal(size=(n, 3))
            directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        latent[i] = settings["acceleration_length"] * directions
        true_force[i] = latent[i] @ tensor.T
        force[i] = true_force[i] + sd * rng.normal(size=(n, 3))
        acceleration[i] = latent[i] + sd * rng.normal(size=(n, 3))
    observed = ReadingSet(force, acceleration, np.full((replicates, n), sd), np.full((replicates, n), sd))
    return observed, dict(tensor=tensor, true_acceleration=latent, true_force=true_force)


def measured_scalar_diagnostics(readings):
    """Common positive scalar infimum from measured totals, without a tensor input.

    For the recipe's common equal channel SD, sum K is the total-Gram Rayleigh
    quotient. Its minimum is not treated as chi-square with fitted degrees of
    freedom. Comparing it to chi-square(3N) is conservative under scalar H0:
    C_min <= K(m_true), and K(m_true) has that exact distribution.
    """
    if not readings.isotropic or not np.all(readings.s == 1):
        raise ValueError("tensor-world scalar diagnostic requires equal isotropic supplied channel SDs")
    p, q, dot = readings.P, readings.Q, readings.D
    per_q = np.maximum(0.0, 0.5 * (p + q - np.hypot(p - q, 2 * dot)))
    a, b, c = p.sum(axis=1), q.sum(axis=1), dot.sum(axis=1)
    projective = np.maximum(0.0, 0.5 * (a + b - np.hypot(a - b, 2 * c)))
    positive = np.where(c > 0, projective, np.minimum(a, b))
    return dict(common_scalar_cost=positive, common_scalar_projective_cost=projective,
                sum_Q=per_q.sum(axis=1), per_reading_Q=per_q,
                extra_common_vs_individual_projective=np.maximum(0.0, projective - per_q.sum(axis=1)),
                scalar_h0_conservative_p=chi2.sf(positive, 3 * readings.readings))


def evaluate(readings, prior_center):
    """The whole estimator-facing boundary: observed readings and declared prior."""
    quoted, raw = _combine(float(prior_center), "quoted"), _combine(float(prior_center), "untempered")
    result = posterior_summaries(readings, LAW, quoted, ["ratio_of_means_rss", "median", "interval_95"])
    untempered = posterior_summaries(readings, LAW, raw, ["interval_95"])
    weights = misfit_weights(readings, quoted)
    low, high = result["interval_95"].T
    rlow, rhigh = untempered["interval_95"].T
    return dict(point=result["ratio_of_means_rss"], median=result["median"],
                interval_low=low, interval_high=high, raw_interval_low=rlow, raw_interval_high=rhigh,
                quoted_log_width=np.log(high / low), raw_log_width=np.log(rhigh / rlow),
                eta=bartlett_weight(readings, quoted).copy(), g=weights.copy(),
                min_g=weights.min(axis=1), mean_g=weights.mean(axis=1),
                unresolved=result["unresolved"].copy())


def known_tensor_cost(force, acceleration, tensor, noise_sd):
    """Unfitted investigator reference: exactly chi-square(3N) under the declared G world."""
    tensor, sd = _tensor(tensor), _positive(noise_sd, "noise_sd")
    force, acceleration = np.asarray(force, dtype=float), np.asarray(acceleration, dtype=float)
    if (force.ndim != 3 or force.shape[-1] != 3 or acceleration.shape != force.shape
            or not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration))):
        raise ValueError("known-tensor cost needs matching finite B by N by 3 observations")
    residual = force - acceleration @ tensor.T
    covariance = sd ** 2 * (np.eye(3) + tensor @ tensor.T)
    return np.einsum("bni,ij,bnj->b", residual, np.linalg.inv(covariance), residual)


def investigator_scalar_compatibility(true_force, true_acceleration, tolerance=1e-10):
    """Exact-latent representability; this is never passed to evaluate."""
    force, acceleration = np.asarray(true_force, dtype=float), np.asarray(true_acceleration, dtype=float)
    if (force.ndim != 3 or force.shape[-1] != 3 or acceleration.shape != force.shape
            or not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration))):
        raise ValueError("scalar compatibility needs matching finite B by N by 3 latent vectors")
    tolerance = _positive(tolerance, "scalar tolerance")
    denominator = np.sum(acceleration ** 2, axis=(1, 2))
    numerator = np.sum(force * acceleration, axis=(1, 2))
    scalar = np.divide(numerator, denominator, out=np.zeros_like(denominator), where=denominator > 0)
    residual = np.linalg.norm(force - scalar[:, None, None] * acceleration, axis=(1, 2))
    scale = np.linalg.norm(force, axis=(1, 2))
    every = (denominator == 0) & (scale == 0)
    one = (denominator > 0) & (scalar > 0) & (residual <= tolerance * scale)
    return dict(scalar_exact=one | every, scalar_status=np.where(every, "every", np.where(one, "one", "none")),
                compatible_scalar=np.where(one, scalar, np.nan),
                least_squares_scalar=np.where(denominator > 0, scalar, np.nan), scalar_residual_norm=residual)


def _summary(value):
    value = np.asarray(value)
    finite = value[np.isfinite(value)]
    return dict(mean=float(finite.mean()) if finite.size else None,
                sd=float(finite.std(ddof=1)) if finite.size > 1 else None,
                minimum=float(finite.min()) if finite.size else None,
                maximum=float(finite.max()) if finite.size else None,
                finite_count=int(finite.size), total_count=int(value.size))


def _representative(folder, readings, truth, cell, settings, arrays):
    """Keep native summaries separate from laws sampled on an explicit display grid."""
    folder.mkdir()
    center = settings["prior_center"]
    one = readings.subset(slice(0, 1))
    coarse = math.log(center) + np.linspace(-settings["representative_log_half_span"],
                                           settings["representative_log_half_span"], settings["representative_nodes"])
    grids = [coarse]
    for stage in ("unweighted", "untempered", "quoted"):
        native = posterior_summaries(one, LAW, _combine(center, stage), ["interval_95"])["interval_95"][0]
        if np.all(np.isfinite(native)) and 0 < native[0] < native[1]:
            lo, hi = np.log(native)
            width, middle = hi - lo, (lo + hi) / 2
            grids.append(np.linspace(middle - 4 * width, middle + 4 * width, 401))
    lm = np.unique(np.concatenate(grids))
    prefix = cell["id"] + "/representative"
    arrays[prefix + "/log_mass"] = lm
    prior_log = prior_log_density(_combine(center, "unweighted")["prior"], lm[None], one)[0]
    arrays[prefix + "/prior_density_per_log_mass"] = (2 / math.pi) * np.exp(prior_log)
    stages = {}
    for stage in ("unweighted", "untempered", "quoted"):
        lp = joint_log_density(one, lm[None], LAW, _combine(center, stage), alpha=False)[0][0]
        scaled = np.exp(lp - lp.max())
        ratio = float(max(scaled[0], scaled[-1]))
        density = scaled / np.trapezoid(scaled, lm)
        arrays[prefix + f"/{stage}_log_kernel"] = lp
        arrays[prefix + f"/{stage}_display_density"] = density
        stages[stage] = dict(endpoint_density_ratio=ratio, support_edge_warning=bool(ratio > 1e-7))
    arrays[prefix + "/cloud_size"] = one.stacked_cond_alpha(lm[None], reference=LAW["reference"])[0]
    case = from_reading_set(readings, prior={"family": "half_cauchy", "center": center},
                            name=cell["id"], question="Inspect the measured scalar analysis of a tensor world.",
                            truth={"true_force": truth["true_force"], "true_acceleration": truth["true_acceleration"]},
                            metadata={"recipe_cell": cell["id"], "investigator_tensor_name": cell["tensor_name"]})
    case.investigator_truth["tensor"] = truth["tensor"].tolist()
    replay = dict(case.to_dict(include_truth=True), kind="inspection",
                  numerics={"mass_grid": np.exp(lm).tolist()})
    (folder / "case.yaml").write_text(yaml.safe_dump(replay, sort_keys=False), encoding="utf-8")
    return dict(case_file=f'{cell["id"]}/case.yaml', arrays_prefix=prefix, stages=stages,
                prior_role="Exact half-Cauchy density per unit log mass, sampled on the diagnostic grid.",
                grid_role="Display/diagnostic union of a coarse global window and fine windows around each stage's native interval; not a native probability grid.",
                tail_role="Endpoint ratios are warnings, not omitted-probability bounds. Measured inputs plus g/eta reconstruct the full kernel.")


def _markdown(record):
    lines = ["# Tensor world: scalar compatibility and ours29", "", record["question"], "",
        "Tensor and latent vectors belong to the investigator. Ours29 receives observations, supplied SDs and an explicit prior only.",
        "A tensor world need not have a true scalar mass. No scalar error score is assigned when none exists on the excited subspace.",
        "The independent scalar H0 p-value uses chi-square(3N) conservatively, without subtracting fitted degrees of freedom; alpha=.01 is per series, not familywise.",
        "Known-tensor cost is an unfitted chi-square(3N) reference. Neither diagnostic is a production rejection gate.", "",
        "| Cell | Tensor | Excitation | N | SD | Latent scalar | Point mean | Quoted log width | eta | Mean g | Common scalar cost | Known G cost / df | Fraction scalar p<.01 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in record["cells"]:
        s = row["summary"]
        lines.append(f'| {row["id"]} | {row["tensor_name"]} | {row["geometry"]} | {row["readings"]} | '
                     f'{row["noise_sd"]:g} | {row["investigator"]["scalar_status"]} | {s["point"]["mean"]:.6g} | '
                     f'{s["quoted_log_width"]["mean"]:.6g} | {s["eta"]["mean"]:.6g} | {s["mean_g"]["mean"]:.6g} | '
                     f'{s["common_scalar_cost"]["mean"]:.6g} | {s["known_tensor_cost"]["mean"] / (3 * row["readings"]):.6g} | '
                     f'{row["investigator"]["scalar_small_p_fraction"]:.4g} |')
    lines += ["", "Each cell has a first-series `case.yaml` replay. Complete arrays and sampled three-stage law curves are in values.npz.",
              "Native point/quantile outputs come from the production adaptive pipeline; displayed densities are separately normalized on their stated finite grid.",
              "Tight laws, Q and eta are recorded separately. What the complete physical result must do with incompatibility remains a premise question."]
    return "\n".join(lines) + "\n"


def run(path, out_root="results", replicates=None, progress=print):
    """Run the seeded declared worlds serially; the production equation is unchanged."""
    from .runner import code_version

    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    reps = spec["raw"]["replicates"] if replicates is None else replicates
    if type(reps) is not int or reps < 2:
        raise ValueError("anisotropic_world replicates must be an integer >= 2")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base, suffix = Path(out_root) / spec["raw"]["name"], 1
    out = base / stamp
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    if replicates is None:
        (out / "preset.yaml").write_bytes(source)
    else:
        (out / "source_preset.yaml").write_bytes(source)
        effective = dict(spec["raw"], replicates=reps)
        (out / "preset.yaml").write_text(yaml.safe_dump(effective, sort_keys=False), encoding="utf-8")
    start, arrays, rows = time.perf_counter(), {}, []
    for cell in spec["cells"]:
        if progress:
            progress(f'cell {cell["id"]}: {cell["tensor_name"]}, {cell["geometry"]}, N={cell["readings"]}, SD={cell["noise_sd"]:g}')
        readings, truth = generate_cell(cell, spec["settings"], reps, spec["raw"]["seed"])
        values = evaluate(readings, spec["settings"]["prior_center"])
        values.update(measured_scalar_diagnostics(readings))
        investigator = investigator_scalar_compatibility(truth["true_force"], truth["true_acceleration"])
        cost = known_tensor_cost(readings.force, readings.acceleration, truth["tensor"], cell["noise_sd"])
        for key, value in values.items():
            arrays[f'{cell["id"]}/{key}'] = value
        arrays[f'observations/{cell["id"]}/force'] = readings.force
        arrays[f'observations/{cell["id"]}/acceleration'] = readings.acceleration
        for key in ("true_force", "true_acceleration"):
            arrays[f'investigator/{cell["id"]}/{key}'] = truth[key]
        for key, value in investigator.items():
            arrays[f'investigator/{cell["id"]}/{key}'] = value
        arrays[f'investigator/{cell["id"]}/known_tensor_cost'] = cost
        summaries = {key: _summary(value) for key, value in values.items() if np.asarray(value).ndim == 1}
        summaries["known_tensor_cost"] = _summary(cost)
        exact = np.isfinite(investigator["compatible_scalar"])
        coverage = ((values["interval_low"] <= investigator["compatible_scalar"])
                    & (investigator["compatible_scalar"] <= values["interval_high"]))
        statuses = investigator["scalar_status"]
        status = str(statuses[0]) if np.all(statuses == statuses[0]) else "mixed"
        representative = _representative(out / cell["id"], readings, truth, cell, spec["settings"], arrays)
        rows.append({**cell, "summary": summaries, "representative": representative,
                     "investigator": dict(scalar_status=status, compatible_scalar_count=int(exact.sum()),
                         compatible_scalar=_summary(investigator["compatible_scalar"]),
                         scalar_coverage=float(coverage[exact].mean()) if np.any(exact) else None,
                         scalar_small_p_fraction=float(np.mean(values["scalar_h0_conservative_p"] < .01)),
                         known_tensor_df=3 * cell["readings"]),
                     "operational_errors": int(np.sum(values["unresolved"]))})
    code = code_version()
    code.update(preset_sha256=hashlib.sha256(source).hexdigest(),
                effective_preset_sha256=hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest())
    record = dict(kind="anisotropic_world", preset=spec["raw"]["name"], question=spec["raw"]["question"],
                  seed=spec["raw"]["seed"], replicates=reps, workers=1, settings=spec["settings"],
                  code=code, started=stamp, seconds=time.perf_counter() - start, cells=rows,
                  assumptions=["Calibrated independent isotropic Gaussian instruments in this declared tensor world.",
                               "Paired random streams across tensors omit tensor and prior from the generator key intentionally.",
                               "Scalar H0 diagnostic is conservative under its stated model and is not a production gate.",
                               "The known-G diagnostic is unfitted and uses investigator truth only.",
                               "A small scalar p-value or a tight law does not establish a newly adopted physical response criterion."])
    np.savez_compressed(out / "values.npz", **arrays)
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    (out / "summary.md").write_text(_markdown(record), encoding="utf-8")
    with (out / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["cell", "tensor", "geometry", "readings", "noise_sd", "scalar_status", "point_mean",
                  "quoted_log_width", "eta_mean", "mean_g", "common_scalar_cost", "known_tensor_cost_per_df", "scalar_small_p_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            s = row["summary"]
            writer.writerow(dict(cell=row["id"], tensor=row["tensor_name"], geometry=row["geometry"],
                readings=row["readings"], noise_sd=row["noise_sd"], scalar_status=row["investigator"]["scalar_status"],
                point_mean=s["point"]["mean"], quoted_log_width=s["quoted_log_width"]["mean"],
                eta_mean=s["eta"]["mean"], mean_g=s["mean_g"]["mean"], common_scalar_cost=s["common_scalar_cost"]["mean"],
                known_tensor_cost_per_df=s["known_tensor_cost"]["mean"] / (3 * row["readings"]),
                scalar_small_p_fraction=row["investigator"]["scalar_small_p_fraction"]))
    if progress:
        progress(f"wrote {out}")
    return out, record
