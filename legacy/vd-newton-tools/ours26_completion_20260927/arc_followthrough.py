"""Exploration of a PROPOSED arc prior; no production model is replaced.

Same source snapshot as probe.py. Numerical examples are specified observations,
not an accuracy ranking or a study of repeated-experiment performance.
"""
from pathlib import Path
import json
import math
import numpy as np
from probe import gl, halfcauchy_grid, direct_terms, point, infer, check_near


def logshape(x, ratios, counts=None):
    ratios = np.asarray(ratios)
    counts = np.ones(len(ratios)) if counts is None else np.asarray(counts)
    return .5 * np.sqrt(np.sum(counts / np.cosh(np.asarray(x)[..., None] - np.log(ratios))**2, axis=-1))


def prior_summary(ratios, counts=None):
    x = np.linspace(np.log(min(ratios)) - 30, np.log(max(ratios)) + 30, 100001)
    shape = logshape(x, ratios, counts)
    dx = x[1] - x[0]
    cdf = np.r_[0, np.cumsum((shape[1:] + shape[:-1]) * dx / 2)]
    length = cdf[-1]
    cdf /= length
    modes = (np.diff(shape)[:-1] > 0) & (np.diff(shape)[1:] < 0)
    return {
        "ratios": list(ratios),
        "counts": [int(v) for v in (np.ones(len(ratios), dtype=int) if counts is None else counts)],
        "arc_length": float(length),
        "median": float(np.exp(np.interp(.5, cdf, x))),
        "central_80_percent": np.exp(np.interp([.1, .9], cdf, x)).tolist(),
        "log_density_modes": np.exp(x[1:-1][modes]).tolist(),
        "probability_between_half_and_twice_geometric_midpoint": float(
            np.interp(.5*np.log(min(ratios)*max(ratios))+np.log(2), x, cdf)
            - np.interp(.5*np.log(min(ratios)*max(ratios))-np.log(2), x, cdf)
        ),
    }


def arc_weights(masses, base_weights, base_center, reference_ratios, counts=None):
    # Convert probability quadrature weights for base half-Cauchy to arc dm.
    base_log_density = masses * (2/np.pi) * base_center / (base_center**2 + masses**2)
    return base_weights * logshape(np.log(masses), reference_ratios, counts) / base_log_density


def arc_point(y, sigma, reference_ratios, counts=None, order=400):
    center = np.sqrt(min(reference_ratios)*max(reference_ratios))
    masses, weights = halfcauchy_grid(center, order)
    q, a = direct_terms(y, sigma, masses)
    wa = arc_weights(masses, weights, center, reference_ratios, counts)
    return point(q, a, masses, wa)


def main():
    out = {"source_commit": "c877b1df131bc3a809ef25bc270160eb5b11b6d5"}
    out["two_ratio_log_bimodality_threshold"] = 2 + np.sqrt(3)
    out["priors"] = [prior_summary(r, n) for r, n in [
        ([1., 1.], [1, 1]),
        ([.5, 2.], [1, 1]),
        ([.25, 4.], [1, 1]),
        ([.1, 10.], [1, 1]),
        ([.1, 10.], [1, 9]),
        ([.8, 2., 5.], [1, 1, 1]),
    ]]
    checks = {}
    x = np.linspace(-7, 7, 1001)
    check_near("reciprocal_log_reflection", logshape(x, [.8, 2., 5.]), logshape(-x, [1/.8, 1/2, 1/5]), 1e-14, checks)
    check_near("uniform_replication_shape", logshape(x, [.8, 2., 5.], [7, 7, 7])/np.sqrt(7), logshape(x, [.8, 2., 5.]), 1e-14, checks)
    check_near("same_ratio_reduction", logshape(x, [2., 2., 2.])/np.sqrt(3), .5/np.cosh(x-np.log(2)), 1e-14, checks)
    sigma = np.diag([1., 1., 1., (1/5)**2, (1/2)**2, (1/.8)**2])
    examples = {}
    data = {
        "all_zero": np.zeros(6),
        "weak": np.array([1., .5, -.3, .1, .2, -.1]),
        "strong": np.array([20., 19., 21., 10., 9.8, 10.1]),
    }
    for name, y in data.items():
        arc = arc_point(y, sigma, [5., 2., .8])
        refined = arc_point(y, sigma, [5., 2., .8], order=800)
        check_near(name+"_quadrature_refinement", arc, refined, 2e-8, checks)
        examples[name] = {
            "observations_F_then_a": y.tolist(),
            "half_cauchy_c2_mass": infer(y, sigma, 2., 800),
            "arc_mass": refined,
        }
    out["ours26_three_axis_examples"] = examples

    # Uninformative additional independent push, fixed s=1 in both channels.
    # Its observed zero vector leaves likelihood and BOTH arc/half-Cauchy prior
    # exactly unchanged. Only the existing stacked-norm readout changes.
    padding = []
    for force in (1., 4.):
        base_y = np.array([force, 0.])
        base = infer(base_y, np.eye(2), 1., 600)
        cases = []
        for noise in (1., 10., 100., 1000.):
            y = np.array([force, 0., 0., 0.])
            cov = np.diag([1., noise*noise, 1., noise*noise])
            cases.append({"extra_noise_sd": noise, "mass": infer(y, cov, 1., 600)})
        theta, weights = gl(0, np.pi/2, 1000)
        ll = np.exp(-.5*force*force*np.cos(theta)**2)
        limit = np.dot(weights, np.sin(theta)*ll) / np.dot(weights, np.cos(theta)*ll)
        padding.append({"original_F": force, "original_a": 0., "original_mass": base, "added_independent_zero_push": cases, "infinite_noise_limit_mass": float(limit)})
    out["no_information_padding_existing_readout"] = padding

    # Same measurements, different physical designs. Two instruments see either
    # one common acceleration (stable pair), or two different accelerations.
    sf, sa = np.array([1., 100.]), np.ones(2)
    seff = np.sqrt(np.sum(1/sa**2) / np.sum(1/sf**2))
    out["stable_versus_independent_reference"] = {
        "force_sd": sf.tolist(), "acceleration_sd": sa.tolist(),
        "stable_pair_prior_median": float(seff),
        "independent_push_prior": prior_summary([1., 100.]),
    }
    out["checks"] = checks
    path = Path(__file__).with_name("arc_results.json")
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
