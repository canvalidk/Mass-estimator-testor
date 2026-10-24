"""Combine independent single-reading mass evidence under declared update rules.

One reading i supplies force and acceleration vectors with known isotropic
Gaussian SDs.  Newton II is assumed; mass is common to all readings while each
reading keeps its own latent excitation (magnitude and direction).  Write
u = log m (physical), s_i = sigma_F,i / sigma_a,i, z_i = u - log s_i and
theta_i = atan(exp z_i).  In standardized polar coordinates the original flat
reference is df dalpha dOmega = sigma_F sigma_a r dr dtheta dOmega.
Integrating a reading's own radius and direction under r dr dOmega gives its
likelihood function of mass,

    L_i(u) proportional to J1(h_i(theta_i)),   h_i = |x_i sin theta_i + y_i cos theta_i|,

where x_i, y_i are the standardized readings and J1 is the radial kernel of
``comparison_posterior``.  The single-reading flat law has log-mass density
L_1(u) sech(z_1)/2: the likelihood times a half-Cauchy (uniform-angle) prior.

Three declared ways of joining N readings are implemented.  All reproduce the
single-reading law at N = 1.

``symmetric``
    prior(u) * prod_i L_i(u).  Radius r_i, which treats standardized force and
    acceleration alike, is each reading's nuisance; the prior on mass enters
    once (default: the first reading's half-Cauchy).
``mass_matching``
    Contribution 11, eq. (18): multiply the single-reading mass densities
    p_i(m).  Equivalent to taking acceleration magnitude, with weight
    alpha dalpha, as each reading's nuisance.
``log_matching``
    Contribution 11, eq. (20): multiply the single-reading log-mass densities.

With a common s the three differ by explicit factors that grow with N:
mass_matching = symmetric * cos(theta)**(2(N-1)) and
log_matching = symmetric * (sin(2 theta)/2)**(N-1), up to normalization.
Each single-reading law carries a half-Cauchy factor; the matching rules count
it once per reading.

Summaries are continuous: the log-mass posterior is integrated with
piecewise Gauss-Legendre panels on an adaptively located finite support,
and quantiles use the polynomial antiderivatives of ``comparison_posterior``.
The ratio-of-means readout of the joint law is contribution 11, eq. (28):
the conditional law of each reading's latent pair given mass is the same
under all three references, so E[alpha_i | m] = sigma_a,i cos(theta_i)
h_i / erf(h_i / sqrt 2) applies to each.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np

from comparison_posterior import _PROBS, _antiderivative, _evaluate, _quadrature, _quantiles
from trial_archive import load_trial_archive, load_trial_manifest
from weighted_posterior import standardize


RULES = ("symmetric", "mass_matching", "log_matching")
_SQRT_HALF_PI = math.sqrt(math.pi / 2)
_Z95 = (0, 6)
_Z80 = (1, 5)


def _erf(values):
    flat = np.fromiter((math.erf(float(v)) for v in values.flat), dtype=float, count=values.size)
    return flat.reshape(values.shape)


def log_half_sech(z):
    """log(sech(z)/2), the uniform-angle density of z, stable for large |z|."""
    return -np.logaddexp(z, -z)


def trial_terms(x, y, s, u):
    """Per-reading log likelihood and E[alpha | m] / sigma_a at log masses u.

    ``x`` and ``y`` are standardized readings of shape (N, 3), ``s`` has shape
    (N,) and ``u`` any shape.  Returns arrays of shape (N,) + u.shape and the
    log of the uniform-angle factor for each reading.
    """
    x, y = np.atleast_2d(x), np.atleast_2d(y)
    s = np.asarray(s, dtype=float).reshape(-1)
    u = np.asarray(u, dtype=float)
    z = u[None, ...] - np.log(s).reshape((-1,) + (1,) * u.ndim)
    sn = np.exp(-.5 * np.logaddexp(0., -2 * z))
    cs = np.exp(-.5 * np.logaddexp(0., 2 * z))
    shape = (-1,) + (1,) * u.ndim
    p = np.sum(x * x, axis=1).reshape(shape)
    q = np.sum(y * y, axis=1).reshape(shape)
    d = np.sum(x * y, axis=1).reshape(shape)
    h2 = np.maximum(0., p * sn * sn + q * cs * cs + 2 * d * sn * cs)
    h = np.sqrt(h2)
    factor = np.ones_like(h)
    np.divide(_SQRT_HALF_PI * _erf(h / math.sqrt(2)), h, out=factor, where=h > 1e-12)
    return h2 / 2 + np.log(factor), cs * _SQRT_HALF_PI / factor, log_half_sech(z)


def _combine_terms(log_like, log_angle, u, rule, prior):
    count = log_like.shape[0]
    total = np.sum(log_like, axis=0)
    if rule == "symmetric":
        return total + (log_angle[0] if prior is None else prior(u))
    if prior is not None:
        raise ValueError("an explicit prior is defined only for the symmetric rule")
    angle = np.sum(log_angle, axis=0)
    if rule == "mass_matching":
        return total + angle + (1 - count) * u
    if rule == "log_matching":
        return total + angle
    raise ValueError(f"rule must be one of {RULES}")


def log_normal_prior(centre, width):
    """Log-density in u = log m of a normal law on log mass."""
    if not (centre > 0 and width > 0 and math.isfinite(centre) and math.isfinite(width)):
        raise ValueError("centre and width must be positive and finite")
    mu = math.log(centre)
    return lambda u: -(u - mu) ** 2 / (2 * width * width)


def combine(x, y, s, sigma_a, rule="symmetric", prior=None, *, span=32., step=.05,
            panels=48, order=96, cutoff=60., coarse=None):
    """Joint mass posterior of N readings; returns continuous summaries.

    ``sigma_a`` (N,) supplies physical acceleration SDs for the readout.
    ``prior`` is an optional callable log-density in u for the symmetric rule.
    ``coarse`` may pass precomputed (u, log_like, log_angle) on the coarse grid.
    """
    x, y = np.atleast_2d(np.asarray(x, dtype=float)), np.atleast_2d(np.asarray(y, dtype=float))
    s = np.asarray(s, dtype=float).reshape(-1)
    sigma_a = np.asarray(sigma_a, dtype=float).reshape(-1)
    count = len(x)
    if not count or x.shape != y.shape or x.shape[1] != 3 or s.shape != (count,) or sigma_a.shape != (count,):
        raise ValueError("readings must be (N,3) with one s and sigma_a per reading")
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}")
    u0 = math.log(s[0])
    if coarse is None:
        grid = u0 + np.arange(-span, span + step / 2, step)
        log_like, _, log_angle = trial_terms(x, y, s, grid)
    else:
        grid, log_like, log_angle = coarse
    lp = _combine_terms(log_like, log_angle, grid, rule, prior)
    keep = np.flatnonzero(lp > lp.max() - cutoff)
    lo = max(grid[0], grid[keep[0]] - 2 * step)
    hi = min(grid[-1], grid[keep[-1]] + 2 * step)
    edges = np.linspace(lo, hi, panels + 1)
    half = np.full(panels, (hi - lo) / (2 * panels))
    nodes, weights, vander, projection = _quadrature(order)
    u = (edges[:-1] + edges[1:])[:, None] / 2 + half[:, None] * nodes
    log_like, conditional, log_angle = trial_terms(x, y, s, u)
    lp = _combine_terms(log_like, log_angle, u, rule, prior)
    density = np.exp(lp - lp.max())
    raw_mass = np.sum(density * weights, axis=1) * half
    normalization = raw_mass.sum()
    if not math.isfinite(normalization) or normalization <= 0:
        raise FloatingPointError("posterior normalization failed")
    pdf = density / normalization
    panel_mass = raw_mass / normalization
    panel_end = np.cumsum(panel_mass)
    cdf_coeff = _antiderivative((pdf @ projection)[None]) * half[None, :, None]
    cdf_coeff[..., 0] += (panel_end - panel_mass)[None]
    quantiles = _quantiles(cdf_coeff, edges[None], panel_end[None], np.zeros(1), np.zeros(1))[0]
    w = pdf * weights * half[:, None]
    mean_log = float(np.sum(w * u))
    log_sd = math.sqrt(max(0., float(np.sum(w * u * u)) - mean_log**2))
    acceleration = np.sum(sigma_a.reshape((-1, 1, 1)) * conditional, axis=0)
    readout = float(np.sum(w * np.exp(u) * acceleration) / np.sum(w * acceleration))
    return {"rule": rule, "readings": count, "readout": readout,
            "log_quantiles": quantiles.tolist(), "probabilities": _PROBS.tolist(),
            "median": float(math.exp(quantiles[3])), "mean_log_mass": mean_log, "log_sd": log_sd,
            "support": [float(lo), float(hi)],
            "boundary_panel_mass": float(panel_mass[0] + panel_mass[-1])}


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _series_summaries(x, y, s, sigma_a, prefixes, rules, step=.05, span=32.):
    grid = math.log(s[0]) + np.arange(-span, span + step / 2, step)
    log_like, _, log_angle = trial_terms(x, y, s, grid)
    result = {}
    for count in prefixes:
        coarse = (grid, log_like[:count], log_angle[:count])
        for rule in rules:
            result[(rule, count)] = combine(x[:count], y[:count], s[:count], sigma_a[:count], rule,
                                            step=step, span=span, coarse=coarse)
    return result


def _fraction(values):
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    return {"value": mean, "mcse": float(math.sqrt(max(mean * (1 - mean), 0) / len(values)))}


def _aggregate(summaries, truth, s_ref, null):
    readout = np.array([r["readout"] for r in summaries])
    median = np.array([r["median"] for r in summaries])
    q = np.array([r["log_quantiles"] for r in summaries])
    scale = s_ref if null else truth
    y = math.log(scale)
    result = {"series": len(summaries), "reference": "noise_ratio_s" if null else "true_mass",
              "readout_ratio_quartiles": np.quantile(readout / scale, [.25, .5, .75]).tolist(),
              "median_ratio_quartiles": np.quantile(median / scale, [.25, .5, .75]).tolist(),
              "typical_95_interval_ratio": [float(np.exp(np.median(q[:, _Z95[0]]) - y)),
                                            float(np.exp(np.median(q[:, _Z95[1]]) - y))],
              "median_log_width_95": float(np.median(q[:, _Z95[1]] - q[:, _Z95[0]])),
              "log_median_spread_iqr": float(np.subtract(*np.quantile(np.log(median), [.75, .25]))),
              "max_boundary_panel_mass": float(max(r["boundary_panel_mass"] for r in summaries))}
    if not null:
        for level, (lo, hi) in (("95", _Z95), ("80", _Z80)):
            result[f"coverage_{level}"] = _fraction((q[:, lo] <= y) & (y <= q[:, hi]))
            result[f"truth_below_{level}"] = _fraction(y < q[:, lo])
            result[f"truth_above_{level}"] = _fraction(y > q[:, hi])
    return result


def analyze_run(run_dir, *, series_length=50, prefixes=(1, 2, 5, 10, 20, 50), rules=RULES, cells=None):
    """Split each saved cell into consecutive series and summarize each rule."""
    run_dir = Path(run_dir)
    if max(prefixes) > series_length or min(prefixes) < 1:
        raise ValueError("prefixes must lie between 1 and the series length")
    manifest = load_trial_manifest(run_dir)
    output = []
    for entry in manifest["entries"]:
        if cells is not None and entry["id"] not in cells:
            continue
        arrays = load_trial_archive(run_dir, entry)
        x, y, s = standardize(arrays["observation__force"], arrays["observation__acceleration"],
                              arrays["noise__supplied_force_sd"], arrays["noise__supplied_acceleration_sd"])
        sigma_a = arrays["noise__supplied_acceleration_sd"][:, 0]
        scenario = entry["scenario"]
        truth = float(scenario["true_mass"])
        null = float(scenario["true_acceleration_snr"]) == 0.
        count = len(x) // series_length
        collected = {(rule, n): [] for rule in rules for n in prefixes}
        for k in range(count):
            part = slice(k * series_length, (k + 1) * series_length)
            summaries = _series_summaries(x[part], y[part], s[part], sigma_a[part], prefixes, rules)
            for key, value in summaries.items():
                collected[key].append(value)
        output.append({"cell": entry["id"], "archive_sha256": entry["sha256"],
                       "true_mass": truth, "true_force_snr": scenario["true_force_snr"],
                       "true_acceleration_snr": scenario["true_acceleration_snr"], "null_excitation": null,
                       "series_length": series_length, "series": count,
                       "rules": {rule: {str(n): _aggregate(collected[(rule, n)], truth, float(s[0]), null)
                                        for n in prefixes} for rule in rules}})
    if not output:
        raise ValueError("no requested cells were found")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, type=Path, help="experiment attempt directory with trial archives")
    parser.add_argument("--output", required=True, type=Path, help="new JSON file")
    parser.add_argument("--cells", nargs="+", help="scenario ids to analyze (default all)")
    parser.add_argument("--series-length", type=int, default=50)
    parser.add_argument("--prefixes", nargs="+", type=int, default=[1, 2, 5, 10, 20, 50])
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; choose a new file")
    try:
        cells = analyze_run(args.run, series_length=args.series_length,
                            prefixes=tuple(args.prefixes), cells=args.cells)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    here = Path(__file__).resolve().parent
    record = {"schema_version": 1, "kind": "sequential_update_series",
              "run": str(args.run), "manifest_sha256": _sha256(args.run / "trial_manifest.json"),
              "rules": list(RULES), "prefixes": args.prefixes, "series_length": args.series_length,
              "note": ("Consecutive replicate series of one saved cell share mass and excitation "
                       "magnitude; each reading keeps its own latent excitation in the model. "
                       "Null cells report location and width relative to s, not coverage."),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__},
              "source_sha256": {name: _sha256(here / name) for name in
                                ("sequential_update.py", "weighted_posterior.py", "comparison_posterior.py", "trial_archive.py")},
              "cells": cells}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "x", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": "completed", "output": str(args.output), "cells": len(cells)}))


if __name__ == "__main__":
    main()
