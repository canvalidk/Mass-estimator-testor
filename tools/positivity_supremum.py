"""B02's independent deterministic integrals; never imported by the estimator.

The derivation and exact-arithmetic Gauss error bound are in problems/O13.md,
"30 Oct (B02)". The bound is for discretization, not an interval certificate
for floating-point roots, weights, special functions or evaluation of the bound.

Reproduce the mathematical grid (no simulated worlds or random draws):
    python tools/positivity_supremum.py --grid results/positivity-supremum/grid.json

The hypergeometric series and Kummer transformation are DLMF 13.2.2 and
13.2.39: https://dlmf.nist.gov/13.2.E2 and https://dlmf.nist.gov/13.2.E39.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.special import hyp1f1, roots_legendre


STRENGTHS = (0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0,
             3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0)
ELLIPSES = (1.025, 1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 3.0)
ROUNDING_NOTE = ("Analytic discretization bound in exact arithmetic; floating-point "
                 "nodes, weights, special functions and evaluation of the bound "
                 "are not interval-enclosed.")


def _arguments(k, r, b):
    if type(k) is not int or k not in (1, 2, 3):
        raise ValueError("dimension must be 1, 2 or 3")
    for strength in (r, b):
        if (isinstance(strength, (bool, np.bool_)) or not np.isscalar(strength)
                or not math.isfinite(strength) or strength < 0):
            raise ValueError("channel strengths must be finite and nonnegative")
    if not math.isfinite(r * r + b * b):
        raise ValueError("squared channel strength exceeds the reference range")


def _log_majorant(k, strength_squared, ellipse):
    eta = math.pi / 8 * (ellipse - 1 / ellipse)
    c = math.cosh(eta)
    z = strength_squared / 4 * (1 + c * c)
    p = math.ceil(k / 2 + 1)
    rising, polynomial = 1.0, 1.0
    for j in range(1, p + 1):
        rising *= k / 2 + j - 1
        polynomial += math.comb(p, j) * z ** j / rising
    radial = (1.0, c, 2 * c * c)[k - 1]
    return (math.log(radial * (1 + c * c)) + math.log(polynomial)
            + strength_squared * math.sinh(eta) ** 2 / 4)


def quadrature_error_bound(k, strength_squared, nodes, ellipse):
    """Analytic Bernstein-ellipse bound, with all four angular quadrants split.

    This controls quadrature in exact arithmetic; see ROUNDING_NOTE.
    """
    constant = 4.0 if k == 1 else 4 * math.pi * k
    log_bound = (math.log(constant) + _log_majorant(k, strength_squared, ellipse)
                 - 2 * nodes * math.log(ellipse) - math.log1p(-1 / ellipse))
    if log_bound > 709:
        return math.inf
    # A representable upper value even when the analytic bound underflows.
    return math.exp(max(log_bound, -740))


def _rule(k, strength_squared, tolerance, max_nodes):
    if (isinstance(tolerance, bool) or not math.isfinite(tolerance)
            or not 0 < tolerance < 1):
        raise ValueError("tolerance must be finite and between zero and one")
    if type(max_nodes) is not int or max_nodes < 16:
        raise ValueError("max_nodes must be an integer at least 16")
    options = []
    constant = 4.0 if k == 1 else 4 * math.pi * k
    for ellipse in ELLIPSES:
        log_prefactor = (math.log(constant) + _log_majorant(k, strength_squared, ellipse)
                         - math.log1p(-1 / ellipse))
        needed = math.ceil((log_prefactor - math.log(tolerance)) / (2 * math.log(ellipse)))
        nodes = max(16, 8 * math.ceil(needed / 8))
        if nodes <= max_nodes:
            bound = quadrature_error_bound(k, strength_squared, nodes, ellipse)
            options.append((nodes, bound, ellipse))
    if not options:
        raise ValueError("requested error bound needs more nodes than max_nodes")
    return min(options)


@lru_cache(maxsize=32)
def _nodes(nodes):
    x, w = roots_legendre(nodes)
    return (x + 1) * math.pi / 4, w * math.pi / 4


def expectations(k, r, b, *, tolerance=1e-9, max_nodes=1024):
    """E Q_+ and E Q for x=r*e1+xi, y=b*e1+zeta, xi,zeta~N(0,I).

    The total norm T is integrated analytically. Only the compact Gram gap
    and eigenvector-angle integrals remain. Every quadrant is analytic,
    including at the positive-domain boundaries.
    """
    _arguments(k, r, b)
    strength_squared = float(r * r + b * b)
    nodes, bound, ellipse = _rule(k, strength_squared, tolerance, max_nodes)
    angle, weights = _nodes(nodes)
    if k == 1:
        rho, one_minus_rho, radial_weights = 1.0, 0.0, 1.0
    else:
        rho = np.sin(angle)[:, None]
        # Avoid subtracting two nearly equal numbers at the radial endpoint.
        one_minus_rho = (np.cos(angle) ** 2 / (1 + np.sin(angle)))[:, None]
        density = np.sin(angle) if k == 2 else 2 * np.sin(angle) * np.cos(angle)
        radial_weights = (weights * density)[:, None]
    delta = 2 * math.atan2(b, r)
    positive, unsigned = 0.0, 0.0
    for quadrant in range(4):
        gamma = angle + quadrant * math.pi / 2
        z = strength_squared / 4 * (1 + rho * np.cos(gamma - delta))
        loss = strength_squared / 4 * (
            one_minus_rho + rho * 2 * np.sin((gamma - delta) / 2) ** 2)
        # Kummer transformation: exp(z-lambda/2)*M(-k/2-1,k/2,-z).
        if k == 2:
            polynomial = 1 + 2 * z + z * z / 2
        else:
            polynomial = hyp1f1(-k / 2 - 1, k / 2, -z)
        kernel = np.exp(-loss) * polynomial
        if quadrant < 2:
            cost = one_minus_rho
        else:
            one_minus_h = np.sin(gamma) ** 2 / (1 + np.abs(np.cos(gamma)))
            cost = one_minus_rho + rho * one_minus_h
        mass = radial_weights * weights * kernel
        positive += float(np.sum(mass * cost))
        unsigned += float(np.sum(mass * one_minus_rho))
    scale = k / (2 * math.pi)
    return dict(positive=scale * positive, unsigned=scale * unsigned,
                quadrature_error_bound=bound, nodes=nodes, ellipse=ellipse,
                floating_point_enclosure=False)


def grid(strengths=STRENGTHS, *, tolerance=1e-9):
    """A fixed Cartesian grid; the reflected half uses the proved swap symmetry."""
    rows, maxima = [], []
    for k in (1, 2, 3):
        dimension_rows = []
        for i, r in enumerate(strengths):
            for b in strengths[:i + 1]:
                result = expectations(k, r, b, tolerance=tolerance)
                dimension_rows.append(dict(k=k, R=r, b=b, **result))
                if r != b:
                    dimension_rows.append(dict(k=k, R=b, b=r, **result))
        rows.extend(dimension_rows)
        best = max(dimension_rows, key=lambda row: row["positive"])
        boundary = next(row for row in dimension_rows
                        if row["R"] == max(strengths) and row["b"] == 0)
        maxima.append(dict(k=k, grid_maximum=best, largest_boundary=boundary,
                           maximum_bound=max(row["quadrature_error_bound"]
                                             for row in dimension_rows)))
    return dict(strengths=list(strengths), grid_points_per_dimension=len(strengths) ** 2,
                reflected_half_uses_proved_swap=True, tolerance=tolerance,
                rounding_note=ROUNDING_NOTE, numpy=np.__version__, scipy=scipy.__version__,
                maxima=maxima, rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True, help="save the deterministic grid as JSON")
    args = parser.parse_args()
    result = grid()
    args.grid.parent.mkdir(parents=True, exist_ok=True)
    args.grid.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for item in result["maxima"]:
        row = item["grid_maximum"]
        print(f"k={item['k']}: max {row['positive']:.12g} at ({row['R']}, {row['b']}); "
              f"quadrature bound {row['quadrature_error_bound']:.3g}; "
              f"nodes {row['nodes']}")
    print(ROUNDING_NOTE)


if __name__ == "__main__":
    main()
