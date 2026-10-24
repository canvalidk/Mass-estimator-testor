"""Exact-scatter check of ours27's R^2 -> infinity consistency claim.

Scratch only. Sources: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8,
mass_estimator_lemmas_2026-09-27.md L11/T17 and ours27 equation.
Uses a one-pair experiment with unit noise and mass 2. No latent simulation
is needed beyond the exact 2x2 noncentral-Wishart sufficient statistic.

The size mean mu_M(sqrt(H)) is approximated by sqrt(M+H). Gaussian
Poincare gives Var(norm)<=1, hence sqrt(M+H-1)<=mu<=sqrt(M+H).
The reported point has a uniform multiplicative bound of sqrt(1-1/M)
to 1/sqrt(1-1/M) relative to the exact point, before grid error.
"""
from pathlib import Path
import json
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss


def scatter_draw(rng, dimension, total_signal_squared, trials):
    # Orthogonality in reading-space puts the deterministic signal in row 1.
    n = dimension - 1
    a = rng.chisquare(n, trials)
    b = rng.normal(size=trials)
    c = rng.chisquare(n - 1, trials)
    z = rng.normal(size=(trials, 2))
    z[:, 0] += np.sqrt(total_signal_squared)
    parallel = a + z[:, 0] ** 2
    perpendicular = b * b + c + z[:, 1] ** 2
    cross = np.sqrt(a) * b + z[:, 0] * z[:, 1]
    # u0=(sin theta0,cos theta0), v0=(cos theta0,-sin theta0).
    sn, cs = 2 / np.sqrt(5), 1 / np.sqrt(5)
    p = sn * sn * parallel + cs * cs * perpendicular + 2 * sn * cs * cross
    q = cs * cs * parallel + sn * sn * perpendicular - 2 * sn * cs * cross
    d = sn * cs * (parallel - perpendicular) + (cs * cs - sn * sn) * cross
    return p, q, d


@lru_cache(None)
def quadrature(nodes):
    x, w = leggauss(nodes)
    return (x + 1) * np.pi / 4, w * np.pi / 4


def infer_many(p, q, d, dimension, nodes=2048):
    theta, weights = quadrature(nodes)
    sn, cs = np.sin(theta), np.cos(theta)
    points, widths, covers = [], [], []
    for left in range(0, len(p), 128):
        pp, qq, dd = (v[left:left + 128, None] for v in (p, q, d))
        # Remove the theta-independent trace before exponentiation.
        h_centered = .5 * (pp - qq) * (sn * sn - cs * cs) + 2 * dd * sn * cs
        probability = np.exp(.5 * (h_centered - h_centered.max(axis=1, keepdims=True))) * weights
        h = h_centered + .5 * (pp + qq)
        size_factor = np.sqrt(dimension + h)
        point = (probability * size_factor * sn).sum(axis=1) / (probability * size_factor * cs).sum(axis=1)
        cumulative = np.cumsum(probability, axis=1)
        cumulative /= cumulative[:, -1, None]
        low = np.array([np.interp(.025, cc, theta) for cc in cumulative])
        high = np.array([np.interp(.975, cc, theta) for cc in cumulative])
        points.extend(point)
        widths.extend(high - low)
        covers.extend((low < np.arctan(2)) & (high > np.arctan(2)))
    return np.asarray(points), np.asarray(widths), np.asarray(covers)


def main():
    rng = np.random.default_rng(270928)
    rows = []
    for power in (.25, .75):
        for dimension in (100, 10000, 1000000):
            signal = dimension ** power
            p, q, d = scatter_draw(rng, dimension, signal, 2000)
            points, widths, covers = infer_many(p, q, d, dimension)
            coarse = infer_many(p[:32], q[:32], d[:32], dimension, 4096)[0]
            rows.append({
                "M": dimension,
                "R_squared": signal,
                "power": power,
                "R_squared_over_sqrt_M": signal / np.sqrt(dimension),
                "mass_point_quantiles_10_50_90": np.quantile(points, [.1, .5, .9]).tolist(),
                "fraction_mass_point_within_20_percent_of_2": float(np.mean(abs(points / 2 - 1) < .2)),
                "law_95pct_coverage": float(np.mean(covers)),
                "median_law_95pct_angle_width": float(np.median(widths)),
                "quadrature_refinement_relative_error_first_32": float(np.max(abs(coarse / points[:32] - 1))),
                "size_approximation_relative_bound": float(1 / np.sqrt(1 - 1 / dimension) - 1),
            })
    result = {"source_commit": "b4bfc2040031ef010416334d0c5f4aa2ad05cbd8", "seed": 270928, "trials_per_cell": 2000, "rows": rows}
    Path(__file__).with_name("consistency_counterexample_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
