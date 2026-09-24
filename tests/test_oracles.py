"""Declared oracles: known excitation (ceiling) and known beta (fair floor)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from met.analyst import Estimator, known_excitation_summaries, world_beta
from met.law import ReadingSet

WORLD = {"dimension": 3, "design": "new_excitation", "mass": 2.0, "acceleration_snr": 1.5,
         "direction": "random", "readings": 5, "noise": {"force_sd": 1.0, "acceleration_sd": 0.5},
         "supplied_noise": "exact"}


def test_known_excitation_is_a_truncated_normal():
    rng = np.random.default_rng(1)
    a_true = rng.normal(size=(1, 5, 3))
    f = 2.0 * a_true + 0.7 * rng.normal(size=(1, 5, 3))
    readings = ReadingSet(f, a_true + rng.normal(size=(1, 5, 3)), np.full((1, 5), 0.7), np.full((1, 5), 1.0))
    out = known_excitation_summaries(readings, a_true, ["median", "interval_95", "ratio_of_means"])
    prec = np.sum(a_true ** 2) / 0.49
    mu = np.sum(f * a_true) / 0.49 / prec
    sd = 1 / math.sqrt(prec)
    dens = lambda m: math.exp(-(m - mu) ** 2 / (2 * sd * sd))
    z = quad(dens, 0, np.inf)[0]
    mean = quad(lambda m: m * dens(m), 0, np.inf)[0] / z
    assert out["ratio_of_means"][0] == pytest.approx(mean, rel=1e-9)
    assert quad(dens, 0, out["median"][0])[0] / z == pytest.approx(0.5, abs=1e-9)
    lo, hi = out["interval_95"][0]
    assert quad(dens, 0, lo)[0] / z == pytest.approx(0.025, abs=1e-9)
    assert quad(dens, 0, hi)[0] / z == pytest.approx(0.975, abs=1e-9)


def test_world_beta():
    # omega^2 = R^2 (1 + (m/s)^2) / d = 1.5^2 (1 + (2/2)^2) / 3 = 1.5
    assert world_beta(WORLD) == pytest.approx(1.5 / 2.5)
    assert world_beta(dict(WORLD, acceleration_snr={"normal_rms": 1.5})) == pytest.approx(1.5 / 2.5)


def test_oracles_must_be_declared_by_name():
    with pytest.raises(ValueError, match="oracle must be one of"):
        Estimator({"name": "x", "oracle": "peek", "readouts": ["median"]})
    with pytest.raises(ValueError, match="known_beta"):
        Estimator({"name": "x", "oracle": "known_beta", "reduce": "none",
                   "law": {"reference": "cartesian", "nuisance": "radius"},
                   "combine": {"rule": "likelihood_product", "prior": ["flat_mass"]}, "readouts": ["median"]})
