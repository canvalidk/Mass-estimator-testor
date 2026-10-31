"""B13: the fixed-N weak-identification corner of O15, not the N-growing regime.

The minimax coverage infimum is zero; this does not imply zero limiting
coverage at each fixed finite mass. In particular exp(+/-3) is still inside
the half-Cauchy central 95% prior interval. With noisy, zero true pushes the
measured readings are random, not the all-zero measured input of T24.
"""
import math

import numpy as np

from met.analyst import posterior_summaries
from met.full import isotropic, ours30
from met.law import ReadingSet


def test_zero_push_coverage_at_fixed_mass_approaches_positive_null_floor():
    """300 paired 3D series: shrinking deterministic pushes lose planted-mass coverage.

    c=1, N=20, known Sigma_i=I_6; each true stacked pair has length `signal`.
    Each reading's direction is fixed in advance, cycling the three axes.
    Use the same Gaussian errors at every signal and both planted masses.
    This is a finite Monte Carlo illustration; the minimax result is proved
    separately in O15. It does not establish exact coverage limits from a seed.
    """
    rng = np.random.default_rng(31313)
    series, count, dimension = 300, 20, 3
    noise_force = rng.normal(size=(series, count, dimension))
    noise_acceleration = rng.normal(size=(series, count, dimension))
    directions = np.eye(dimension)[np.arange(count) % dimension]
    sd = np.ones((series, count))
    signals = np.array([30.0, 10.0, 3.0, 1.0, 0.1, 0.0])
    law = {"reference": "cartesian2", "nuisance": "radius"}
    combine = {"rule": "likelihood_product",
               "prior": [{"uniform_angle": {"center": 1.0}}],
               "misfit": "kink", "positivity": True, "calibrate": "bartlett"}
    null_intervals = []
    for mass in (math.exp(-3), math.exp(3)):
        covered = []
        for signal in signals:
            true_acceleration = signal * directions / math.sqrt(1 + mass * mass)
            true_force = mass * true_acceleration
            readings = ReadingSet(noise_force + true_force, noise_acceleration + true_acceleration,
                                  sd, sd)
            result = posterior_summaries(readings, law, combine, ["interval_95"])
            interval = result["interval_95"]
            assert not np.any(result["unresolved"])
            assert np.all(np.isfinite(interval))
            assert np.all((0 < interval[:, 0]) & (interval[:, 0] < interval[:, 1]))
            covered.append((interval[:, 0] <= mass) & (mass <= interval[:, 1]))
            if signal == 0:
                null_intervals.append(interval)

        coverage = np.mean(covered, axis=1)
        # Broad deterministic-seed checks: MC standard errors are at most .029.
        assert np.min(coverage[:3]) > 0.90
        assert coverage[3] > coverage[4] + 0.10
        assert coverage[0] > coverage[-1] + 0.30
        assert abs(coverage[-2] - coverage[-1]) < 0.05
        # The requested two fixed masses have a positive null floor, not limit zero.
        assert 0.40 < coverage[-1] < 0.70

    # At zero true push the data law is exactly the same for every planted mass.
    np.testing.assert_array_equal(null_intervals[0], null_intervals[1])
    # Anchor the fast isotropic batch route to the current full-input ours30 law.
    full = ours30([isotropic(f, a, 1.0, 1.0)
                   for f, a in zip(noise_force[0], noise_acceleration[0])],
                  1.0, readout=False)
    np.testing.assert_allclose(null_intervals[0][0], full["interval"], rtol=1e-6)
