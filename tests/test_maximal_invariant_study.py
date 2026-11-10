"""Independent production anchors and numerical checks for B21's study."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.full import _h_terms, isotropic, ours30
from tools.maximal_invariant import gram_energy, profile_log_kernel, integrated_log_kernel
from tools.maximal_invariant_study import angular_cdfs, run_cell, series_statistics, validate_config


@pytest.mark.parametrize("count", [1, 5, 20])
def test_study_ours30_kernel_and_curvature_match_full_equation(count):
    rng = np.random.default_rng(6969+count)
    x, y = rng.normal(size=(2, 2, count, 3))
    x += [.6, .3, 0.]
    y += [.3, .1, 0.]
    gram, weighted, eta = series_statistics(x, y)
    angles = np.linspace(.01, math.pi/2-.01, 31)
    for b in range(2):
        readings = [isotropic(f, a, 1., 1.) for f, a in zip(x[b], y[b])]
        full = ours30(readings, 1., readout=False)
        assert eta[b] == pytest.approx(full["eta"], rel=1e-10, abs=1e-12)
        expected = full["eta"]*sum(g*_h_terms(r.pre(), 1., angles)/2
                                    for g, r in zip(full["g"], readings))
        np.testing.assert_allclose(eta[b]*gram_energy(weighted[b], angles)/2, expected,
                                   rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(gram_energy(gram[b], angles),
                                   sum(_h_terms(r.pre(), 1., angles) for r in readings),
                                   rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("kind", ["profile", "integrated", "flat", "ours30"])
def test_study_angular_cdf_matches_adaptive_integral(kind):
    rng = np.random.default_rng(6901)
    x, y = rng.normal(size=(2, 1, 5, 3))
    x += [1.5, 0., 0.]
    y += [.6, 0., 0.]
    gram, weighted, eta = series_statistics(x, y)
    angle = math.atan(1.7)
    got = angular_cdfs(gram, weighted, eta, 15, angle, 256)[kind][0]
    kernel = {"profile": lambda t: profile_log_kernel(gram_energy(gram[0], t), 15),
              "integrated": lambda t: integrated_log_kernel(gram_energy(gram[0], t), 15),
              "flat": lambda t: gram_energy(gram[0], t)/2,
              "ours30": lambda t: eta[0]*gram_energy(weighted[0], t)/2}[kind]
    peak = max(float(kernel(t)) for t in np.linspace(0., math.pi/2, 100))
    f = lambda t: math.exp(float(kernel(t))-peak)
    left = quad(f, 0., angle, epsabs=1e-10, epsrel=1e-10, limit=100)[0]
    right = quad(f, angle, math.pi/2, epsabs=1e-10, epsrel=1e-10, limit=100)[0]
    assert got == pytest.approx(left/(left+right), abs=2e-7)


def test_study_seeded_cell_replays_and_records_paired_decisions():
    config = {"seed": 69261109, "series": 4, "dimension": 3, "mass": 1., "quadrature_nodes": 32}
    job = (config, 0, .5, 5)
    summary, arrays = run_cell(job)
    other, replay = run_cell(job)
    assert summary == other
    for key in arrays:
        np.testing.assert_array_equal(arrays[key], replay[key])
    for law in ("profile", "integrated", "flat", "ours30"):
        cdf = arrays[law+"_cdf"]
        np.testing.assert_array_equal(arrays[law+"_covered"], (cdf >= .025) & (cdf <= .975))
        assert summary[law] == arrays[law+"_covered"].mean()
    assert summary["profile_minus_flat"] == pytest.approx(summary["profile"]-summary["flat"])


@pytest.mark.parametrize("field,value", [("noise_sd", 2.), ("integrated_measure", "exponential")])
def test_locked_study_rejects_misdescribed_noise_or_nuisance_measure(field, value):
    config = {"dimension": 3, "mass": 1., "prior_center": 1., "noise_sd": 1.,
              "integrated_measure": "flat_dDelta"}
    validate_config(config)
    config[field] = value
    with pytest.raises(ValueError):
        validate_config(config)
