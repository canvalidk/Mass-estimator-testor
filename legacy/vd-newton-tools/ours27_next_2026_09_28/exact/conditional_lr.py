"""Exact finite-sample conditional uncertainty benchmark for one noise pair.

Scratch derivation, not a new mass-definition law or production change.
Managed source: VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

At a tested angle theta, write A=x*sin(theta)+y*cos(theta) and
B=x*cos(theta)-y*sin(theta). Under this candidate experiment, A~N(v,I_M)
and B~N(0,I_M) independently. Conditioning on A removes the entire unknown
v, with no size prior. In particular C/sqrt(H)~N(0,1) exactly, where
H=A'A, Q=B'B, C=A'B. Its confidence-set inversion is exact but can accept
the wrong stationary branch at weak signal.

A better benchmark is the unrestricted-line likelihood ratio
L=Q-lambda_min(S)=(Q-H+sqrt((Q-H)^2+4*C^2))/2.
Conditional on H=h, rotate B into A's direction: C=sqrt(h)*Z,
Q=Z^2+W, Z~N(0,1), W~chi2_(M-1), independently. Consequently

  P(L<=ell | H=h)
   = 2*sqrt(ell)/sqrt(2*pi) * integral_0^1
       exp(-ell*t^2/2) F_chi2_(M-1)((ell+h)*(1-t^2)) dt.

This follows by rearranging L<=ell into
Z^2 <= ell*(1-W/(ell+h)). For M=1, W=0 and L~chi2_1.
For h=0, L~chi2_M. Large h gives chi2_1. At the candidate angle, the
conditional CDF is uniform, so {m>0 : conditional p(m)>alpha} has exact
1-alpha coverage for every deterministic v, finite M, correct Gaussian
noise. This is a confidence SET, not a proper probability law on masses.
It may be unbounded, disconnected, all of (0,infinity), or empty. It is
therefore NOT offered as the requested single total positive estimator.

The formula is the one-parameter conditional LR distribution familiar
from Moreira (2003), derived directly here for our known-noise experiment:
https://doi.org/10.1111/1468-0262.00438 . A current primary paper prints the
same distribution in its one-parameter corollary:
https://arxiv.org/html/2509.04144v1 . Its IV asymptotic setting is distinct;
the present finite-sample exactness follows from our independent Gaussians.

Locally, C changes with angle at rate approximately R^2 while Var(C)=
R^2+M. Thus the width scale sqrt(R^2+M)/R^2 emerges without choosing an
arbitrary candidate-size measure. This does not turn its confidence curve
into a unique joint candidate law: a normalization or density conversion
would be an additional, generally non-coverage-preserving step.
"""
from pathlib import Path
import json
import math
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss


def gamma_p(a, x):
    """Regularized lower incomplete gamma via series/continued fraction."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    small = (x > 0) & (x < a + 1)
    xs = x[small]
    if len(xs):
        delta = np.full_like(xs, 1 / a)
        total = delta.copy()
        for n in range(1, 500):
            delta *= xs / (a + n)
            total += delta
            if np.max(abs(delta / total)) < 2e-15:
                break
        else:
            raise AssertionError("Gamma series failed")
        out[small] = total * np.exp(a * np.log(xs) - xs - math.lgamma(a))
    large = x >= a + 1
    xl = x[large]
    if len(xl):
        b = xl + 1 - a
        c = np.full_like(xl, 1e300)
        d = 1 / b
        fraction = d.copy()
        for n in range(1, 500):
            an = -n * (n - a)
            b += 2
            d = an * d + b
            d = np.where(abs(d) < 1e-300, 1e-300, d)
            c = b + an / c
            c = np.where(abs(c) < 1e-300, 1e-300, c)
            d = 1 / d
            change = d * c
            fraction *= change
            if np.max(abs(change - 1)) < 2e-15:
                break
        else:
            raise AssertionError("Gamma continued fraction failed")
        out[large] = 1 - np.exp(a * np.log(xl) - xl - math.lgamma(a)) * fraction
    return np.clip(out, 0, 1)


@lru_cache(None)
def grid(nodes):
    z, w = leggauss(nodes)
    return (z + 1) / 2, w / 2


def lr_cdf(ell, h, dimension, nodes=128):
    ell, h = np.broadcast_arrays(np.asarray(ell, dtype=float), np.asarray(h, dtype=float))
    if dimension == 1:
        return gamma_p(.5, ell / 2)
    t, w = grid(nodes)
    ll, hh = ell[..., None], h[..., None]
    cdf = gamma_p((dimension - 1) / 2, (ll + hh) * (1 - t * t) / 2)
    integral = (np.exp(-ll * t * t / 2) * cdf) @ w
    return np.clip(np.sqrt(2 * ell / np.pi) * integral, 0, 1)


def statistic(h, q, c):
    d = q - h
    root = np.sqrt(d * d + 4 * c * c)
    out = (d + root) / 2
    negative = d < 0
    out[negative] = 2 * c[negative]**2 / (root[negative] - d[negative])
    return out


def main():
    rng = np.random.default_rng(2809273)
    check_x = np.geomspace(1e-8, 100, 100)
    known_gamma_half = np.array([math.erf(math.sqrt(v)) for v in check_x])
    known_gamma_two = -np.expm1(-check_x) - check_x * np.exp(-check_x)
    gamma_error = max(float(max(abs(gamma_p(.5, check_x) - known_gamma_half))),
                      float(max(abs(gamma_p(2., check_x) - known_gamma_two))))
    assert gamma_error < 2e-14
    special_errors = []
    for dimension in (3, 6, 21, 63):
        ell = np.array([.01, .2, 1., 4., 12., 30., 60., 100.])
        error = np.max(abs(lr_cdf(ell, 0., dimension, 256) - gamma_p(dimension / 2, ell / 2)))
        special_errors.append(float(error))
    assert max(special_errors) < 2e-8

    rows = []
    trials = 6000
    for dimension, energy in ((1, 0.), (6, 0.), (6, 6.), (21, 7.), (63, 21.), (63, 189.)):
        a = rng.normal(size=(trials, dimension))
        a[:, 0] += np.sqrt(energy)
        b = rng.normal(size=(trials, dimension))
        h = np.sum(a*a, axis=1)
        q = np.sum(b*b, axis=1)
        c = np.sum(a*b, axis=1)
        ell = statistic(h, q, c)
        p_values = 1 - lr_cdf(ell, h, dimension)
        refined = 1 - lr_cdf(ell[:40], h[:40], dimension, 256)
        refinement_error = float(np.max(abs(refined - p_values[:40])))
        assert refinement_error < 2e-7
        score = c / np.sqrt(h)
        coverage = float(np.mean(p_values > .05))
        assert abs(coverage - .95) < 5 * np.sqrt(.95*.05/trials)
        rows.append({
            "M": dimension, "R_squared": energy, "trials": trials,
            "conditional_LR_95pct_set_coverage": coverage,
            "exact_conditional_score_95pct_set_coverage": float(np.mean(abs(score) < 1.959963984540054)),
            "ordinary_chi1_LR_95pct_set_coverage": float(np.mean(ell < 3.841458820694124)),
            "mean_conditional_LR_p_value": float(np.mean(p_values)),
            "conditional_CDF_quadrature_refinement_error": refinement_error,
        })
    critical_values = []
    for dimension in (1, 6, 21, 63):
        for h in (0., float(dimension), 10. * dimension):
            lo, hi = 0., dimension + 15*np.sqrt(2*dimension) + 50
            for _ in range(50):
                mid = (lo + hi) / 2
                if float(lr_cdf(mid, h, dimension)) < .95:
                    lo = mid
                else:
                    hi = mid
            critical_values.append({"M": dimension, "h": h, "conditional_95pct_critical_LR": (lo + hi)/2})
    output = {"source_commit": "b4bfc2040031ef010416334d0c5f4aa2ad05cbd8", "seed": 2809273,
              "gamma_known_formula_max_error": gamma_error,
              "h_zero_equals_chi_squared_M_max_error": max(special_errors),
              "rows": rows, "critical_values": critical_values,
              "scope": "Exact confidence sets only; no proper mass law or joint candidate definition claimed."}
    Path(__file__).with_name("conditional_lr_results.json").write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
