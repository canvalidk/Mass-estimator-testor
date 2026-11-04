"""Seeded T54/T55 sanity study; run python -m met.joint_sanity PRESET.

Full positive joint profile-law quadrature is used in the small-noise arm.
The replication arm tests local Wald intervals, not an invented globally
calibrated joint law. Truth enters generation/scoring only. Every failure
counts against coverage; numerical refinement is reported separately.
"""
import argparse
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re

import numpy as np
from scipy import optimize, special, stats
import yaml


QUANTITIES = ("mass1", "mass2", "sum", "parallel")


def validate(config):
    config = deepcopy(config)
    if set(config) != {"kind", "name", "question", "seed", "replicates", "probability", "small_noise", "repeated"}:
        raise ValueError("joint_sanity requires its explicit study settings")
    if config["kind"] != "joint_sanity" or not 0 < config["probability"] < 1:
        raise ValueError("invalid study kind or probability")
    if type(config["seed"]) is not int or config["seed"] < 0:
        raise ValueError("seed must be a nonnegative integer")
    if type(config["replicates"]) is not int or config["replicates"] < 2:
        raise ValueError("replicates must be an integer >=2")
    if not isinstance(config["name"], str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", config["name"]):
        raise ValueError("study name must be one path component")
    small, repeated = config["small_noise"], config["repeated"]
    if set(small) != {"masses", "acceleration", "scales", "angle_nodes", "refinement_nodes", "refinement_replicates", "prior"}:
        raise ValueError("small_noise needs all its numerical and physical declarations")
    if set(repeated) != {"masses", "accelerations", "variance_diagonals", "blocks"}:
        raise ValueError("repeated needs all its physical declarations")
    for value in (small["masses"], repeated["masses"]):
        if np.shape(value) != (2,) or np.any(~np.isfinite(value)) or np.any(np.asarray(value) <= 0):
            raise ValueError("two positive finite masses are required")
    for values in (small["scales"], repeated["accelerations"], [small["acceleration"]]):
        if not values or np.any(~np.isfinite(values)) or np.any(np.asarray(values) <= 0):
            raise ValueError("this benchmark requires positive finite scales/signals")
    v = np.asarray(repeated["variance_diagonals"])
    if v.shape != (len(repeated["accelerations"]), 3) or np.any(~np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("each event needs three positive finite variances")
    for key in ("angle_nodes", "refinement_nodes", "refinement_replicates"):
        if type(small[key]) is not int or small[key] < 2:
            raise ValueError("quadrature controls must be integers >=2")
    if small["refinement_nodes"] <= small["angle_nodes"]:
        raise ValueError("refinement must increase node count")
    if small["prior"] != "independent_half_cauchy_center_1":
        raise ValueError("the study requires its explicitly named prior")
    if not repeated["blocks"] or any(type(n) is not int or n < 2 for n in repeated["blocks"]):
        raise ValueError("block counts must be integers >=2")
    return dict(kind="joint_sanity", raw=config, cells=small["scales"]+repeated["blocks"])


def run_preset(path, out_root, *, replicates=None, progress=print):
    config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if replicates is not None:
        config["replicates"] = replicates
    config = validate(config)["raw"]
    target = Path(out_root)/config["name"]/datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")
    records, metadata = run(config, target)
    if progress:
        progress(str(target))
    return target, dict(summary=records, **metadata)


def write_summary(records, path, probability=.95):
    lines = ["# Joint Gaussian sanity study", "",
        f"Target: {probability:.1%} intervals. Replicates per cell: {sorted(set(int(r['replicates']) for r in records))}.",
        "Binomial intervals describe Monte Carlo uncertainty; quadrature and local-Wald limitations remain.", "",
        "| Arm | Setting | Method | Quantity | Coverage | 95% Monte Carlo interval | Failures |",
        "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in records:
        lines.append(f"| {r['arm']} | {r['setting']} | {r['method']} | {r['quantity']} | "
            f"{float(r['coverage']):.2%} | {float(r['mc_low']):.2%}–{float(r['mc_high']):.2%} | {r['failures']} |")
    Path(path).write_text("\n".join(lines)+"\n", encoding="utf-8")


def quantities(masses):
    x, y = np.asarray(masses).T
    return np.array([x, y, x+y, x*y/(x+y)]).T


def profile_derivatives(masses, scatter, variances):
    """Exact one-type summed profile gradient and negative Hessian, d=1."""
    c = np.r_[masses, 1.]
    p = 1/np.asarray(variances)
    h = p*c
    gram, t = c@h, scatter@h
    value = h@t
    u, pt = h[:2], p[:2]*t[:2]
    gradient = pt/gram-value*u/gram**2
    hessian = ((p[:2, None]*scatter[:2, :2]*p[None, :2])/gram
        -(2*np.outer(pt, u)+2*np.outer(u, pt)+value*np.diag(p[:2]))/gram**2
        +4*value*np.outer(u, u)/gram**3)
    return value/(2*gram), gradient, -hessian


def fit_and_widths(samples, variances):
    scatters = [y.T@y for y in samples]

    def objective(theta):
        masses = np.exp(theta)
        parts = [profile_derivatives(masses, s, v) for s, v in zip(scatters, variances)]
        count = sum(len(y) for y in samples)
        return -sum(p[0] for p in parts)/count, -sum((p[1] for p in parts))*masses/count

    pooled = np.concatenate(samples)
    ratio = np.mean(pooled[:, :2], axis=0)/np.mean(pooled[:, 2])
    starts = [np.zeros(2), np.log(np.clip(ratio, 1e-4, 1e4))]
    fits = [optimize.minimize(objective, start, jac=True, method="L-BFGS-B",
                             bounds=[(-15, 15)]*2,
                             options={"ftol": 1e-13, "gtol": 1e-7}) for start in starts]
    best = min(fits, key=lambda f: f.fun)
    masses = np.exp(best.x)
    if not best.success or np.any(np.abs(best.x) > 14.9):
        raise ValueError("optimizer failed or reached a numerical boundary")
    sensitivity = sum(profile_derivatives(masses, s, v)[2] for s, v in zip(scatters, variances))
    meat = np.zeros((2, 2))
    for data, variance in zip(samples, variances):
        p, c = 1/np.asarray(variance), np.r_[masses, 1.]
        latent = (data@(p*c))/(c@(p*c))
        scores = latent[:, None]*p[:2]*(data[:, :2]-latent[:, None]*masses)
        meat += scores.T@scores
    np.linalg.cholesky(sensitivity)
    raw = np.linalg.inv(sensitivity)
    sandwich = raw@meat@raw
    np.linalg.cholesky(sandwich)
    total = sum(masses)
    gradients = np.array([[1., 0.], [0., 1.], [1., 1.],
                          [masses[1]**2/total**2, masses[0]**2/total**2]])
    widths = [np.sqrt(np.einsum("ai,ij,aj->a", gradients, cov, gradients))
              for cov in (raw, sandwich, np.diag(np.diag(sandwich)))]
    return quantities(masses), widths


class JointQuadrature:
    """Full positive quadrant; uniform angles = product half-Cauchy prior."""
    def __init__(self, nodes, probability):
        roots, weights = np.polynomial.legendre.leggauss(nodes)
        self.axis = np.tan((roots+1)*np.pi/4)
        self.prior_weights = np.outer(weights, weights)/4
        x, y = np.meshgrid(self.axis, self.axis, indexing="ij")
        self.masses = np.stack([x, y], axis=-1)
        self.gram = 1+x*x+y*y
        self.values = quantities(self.masses.reshape(-1, 2))
        self.orders = [np.argsort(self.values[:, j], kind="stable") for j in range(4)]
        self.ranks = [(1-probability)/2, (1+probability)/2]

    def intervals(self, data, scale):
        projection = np.einsum("...j,j->...", self.masses, data[:2])+data[2]
        log_factor = projection**2/(2*scale**2*self.gram)
        weights = self.prior_weights*np.exp(log_factor-np.max(log_factor))
        weights /= weights.sum()

        def summarize(joint):
            flat = joint.ravel()
            result = []
            for j, order in enumerate(self.orders):
                # Discrete quadrature CDF; paired grid refinement exposes its error.
                cdf = np.cumsum(flat[order])
                indices = np.minimum(np.searchsorted(cdf, self.ranks), len(order)-1)
                result.append(self.values[order[indices], j])
            return np.asarray(result)

        return summarize(weights), summarize(np.outer(weights.sum(axis=1), weights.sum(axis=0)))


def run(config, output):
    config = validate(config)["raw"]
    source = Path(__file__).read_bytes()
    executed_code_hash = hashlib.sha256(source).hexdigest()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output/"study_source.py").write_bytes(source)
    rng = np.random.default_rng(config["seed"])
    count, probability = config["replicates"], config["probability"]
    z = stats.norm.ppf((1+probability)/2)
    records, arrays = [], {}

    def record(arm, setting, method, intervals, truth, failures=0):
        covered = (intervals[:, :, 0] <= truth)&(truth <= intervals[:, :, 1])
        for j, quantity in enumerate(QUANTITIES):
            successes = int(covered[:, j].sum())
            # Exact binomial confidence interval describes Monte Carlo uncertainty.
            low = stats.beta.ppf(.025, successes, count-successes+1) if successes else 0.
            high = stats.beta.ppf(.975, successes+1, count-successes) if successes < count else 1.
            records.append(dict(arm=arm, setting=setting, method=method, quantity=quantity,
                covered=successes, replicates=count, coverage=successes/count,
                mc_low=low, mc_high=high, failures=failures,
                median_width=float(np.nanmedian(intervals[:, j, 1]-intervals[:, j, 0]))))
        arrays[f"{arm}_{setting}_{method}"] = intervals

    small = config["small_noise"]
    if small["prior"] != "independent_half_cauchy_center_1":
        raise ValueError("this reference requires its explicitly named prior")
    quad = JointQuadrature(small["angle_nodes"], probability)
    refined = JointQuadrature(small["refinement_nodes"], probability)
    refinement = []
    truth = quantities(small["masses"])
    for scale in small["scales"]:
        noise = rng.normal(size=(count, 3))
        data = small["acceleration"]*np.r_[small["masses"], 1.]+scale*noise
        arrays[f"small_{scale}_data"] = data
        results = [quad.intervals(y, scale) for y in data]
        for method, index in (("joint_law", 0), ("discard_correlation", 1)):
            record("small", scale, method, np.array([r[index] for r in results]), truth)
        changes, shifts = [], []
        for index, y in enumerate(data[:small["refinement_replicates"]]):
            finer = refined.intervals(y, scale)[0]
            coarse = results[index][0]
            changes.append(((coarse[:, 0] <= truth)&(truth <= coarse[:, 1])) !=
                           ((finer[:, 0] <= truth)&(truth <= finer[:, 1])))
            shifts.append(np.abs(finer-coarse)/(finer[:, 1]-finer[:, 0])[:, None])
        refinement.append(dict(scale=scale, paired_replicates=len(changes),
            coverage_decision_changes=np.sum(changes, axis=0).tolist(),
            max_endpoint_shift_in_interval_widths=np.max(shifts, axis=(0, 2)).tolist()))
        print(f"small-noise scale {scale} complete", flush=True)

    repeated = config["repeated"]
    truth = quantities(repeated["masses"])
    variances = np.array(repeated["variance_diagonals"])
    for blocks in repeated["blocks"]:
        intervals = np.full((3, count, 4, 2), np.nan)
        failures = 0
        all_data = []
        for rep in range(count):
            samples = [a*np.r_[repeated["masses"], 1.]+rng.normal(size=(blocks, 3))*np.sqrt(v)
                       for a, v in zip(repeated["accelerations"], variances)]
            all_data.append(samples)
            try:
                center, widths = fit_and_widths(samples, variances)
                for method, width in enumerate(widths):
                    intervals[method, rep] = np.column_stack([center-z*width, center+z*width])
            except (ValueError, np.linalg.LinAlgError, FloatingPointError):
                failures += 1
        arrays[f"repeated_{blocks}_data"] = np.asarray(all_data)
        for index, method in enumerate(("raw_local", "sandwich_local", "discard_correlation")):
            record("repeated", blocks, method, intervals[index], truth, failures)
        print(f"repeated blocks {blocks} complete; failures={failures}", flush=True)
    np.savez_compressed(output/"values.npz", **arrays)
    with (output/"summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    write_summary(records, output/"summary.md", probability)
    metadata = dict(preset=config, refinement=refinement, python=platform.python_version(),
                    numpy=np.__version__, code_sha256=executed_code_hash,
                    limitations=["Quadrature is approximate; paired refinement is not a tail certificate.",
                    "Replication intervals are local Wald diagnostics, not a global calibrated positive law.",
                    "Binomial intervals quantify simulation error, not proof of coverage."])
    (output/"run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (output/"preset.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return records, metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("preset")
    parser.add_argument("--out")
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.preset).read_text(encoding="utf-8"))
    target = args.out or str(Path("results")/config["name"]/datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S"))
    run(config, target)
    print(target)
