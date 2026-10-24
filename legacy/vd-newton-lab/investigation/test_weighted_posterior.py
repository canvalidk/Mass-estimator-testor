"""Independent oracles for prior-weighted flat-law posteriors."""

import math
import unittest

import numpy as np

from comparison_posterior import infer_batch
from weighted_posterior import (LogNormalFactor, LogNormalPrior, SechTilt,
                                infer_weighted, standardize, weighted_points)


def _readings(seed, count=30, force=2., acceleration=.7):
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(count, 3)) + np.array([force, 0, 0]),
            rng.normal(size=(count, 3)) + np.array([acceleration, 0, 0]))


class WeightedPosteriorTests(unittest.TestCase):
    def test_no_factor_reproduces_existing_flat_law(self):
        force, acceleration = _readings(1, force=3., acceleration=1.)
        expected = infer_batch(force, acceleration, truth=3.)
        result = infer_weighted(force, acceleration, None, truth=3.)
        np.testing.assert_allclose(result["point"], expected["points"]["flat_joint"], rtol=1e-14)
        for key, value in expected["distributions"]["flat_joint"].items():
            np.testing.assert_allclose(result["distribution"][key], value, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(weighted_points(force, acceleration, [None])[0],
                                   expected["points"]["flat_joint"], rtol=1e-13)

    def test_sech_tilt_null_matches_readme_closed_form(self):
        # Independent oracle: at zero readings log(M) has density
        # proportional to exp(lam sech l) / cosh l (README section 4).
        zero = np.zeros((1, 3))
        l = np.linspace(-40, 40, 400001)
        for lam in (-4., 0., 4., 16., 64.):
            density = np.exp(lam * (1 / np.cosh(l) - 1)) / np.cosh(l)
            density /= np.trapezoid(density, l)
            sd = math.sqrt(np.trapezoid(l * l * density, l))
            cdf = np.cumsum(density) * (l[1] - l[0])
            upper = l[np.searchsorted(cdf, .975)]
            result = infer_weighted(zero, zero, SechTilt(lam))
            self.assertAlmostEqual(result["point"][0], 1., places=12)
            self.assertAlmostEqual(result["distribution"]["log_sd"][0], sd, places=5)
            self.assertAlmostEqual(result["distribution"]["log_upper_95"][0], upper, places=3)
            self.assertAlmostEqual(result["distribution"]["log_median"][0], 0., places=9)

    def test_lognormal_prior_is_returned_at_zero_readings(self):
        zero = np.zeros((1, 3))
        for centre, width in ((1., .5), (4., .3), (.2, 1.)):
            summary = infer_weighted(zero, zero, LogNormalPrior(centre, width))["distribution"]
            self.assertAlmostEqual(summary["mean_log_mass"][0], math.log(centre), places=9)
            self.assertAlmostEqual(summary["log_sd"][0], width, places=9)
            self.assertAlmostEqual(summary["log_upper_95"][0],
                                   math.log(centre) + 1.959963984540054 * width, places=6)

    def test_narrow_factor_pulls_point_to_its_centre(self):
        force, acceleration = _readings(3, force=2., acceleration=2.)
        point = weighted_points(force, acceleration, [LogNormalFactor(5., .01)])[0]
        np.testing.assert_allclose(point, 5., rtol=.05)

    def test_refinement_and_shared_panels(self):
        force, acceleration = _readings(2)
        factors = [SechTilt(.35), SechTilt(8., 2.), LogNormalFactor(3., .25), LogNormalPrior(.5, 2.)]
        for factor in factors:
            coarse = infer_weighted(force, acceleration, factor, order=96)
            fine = infer_weighted(force, acceleration, factor, order=192)
            np.testing.assert_allclose(coarse["point"], fine["point"], rtol=1e-9)
            np.testing.assert_allclose(coarse["distribution"]["log_median"],
                                       fine["distribution"]["log_median"], atol=1e-7)
            np.testing.assert_allclose(weighted_points(force, acceleration, [factor])[0],
                                       coarse["point"], rtol=1e-12)
        # A batch sharing one panel set agrees with separate evaluation.
        together = weighted_points(force, acceleration, factors)
        for k, factor in enumerate(factors):
            np.testing.assert_allclose(together[k], weighted_points(force, acceleration, [factor])[0], rtol=1e-9)

    def test_per_reading_centres_broadcast(self):
        force, acceleration = _readings(4, count=5)
        centres = np.array([.5, 1., 2., 4., 8.])
        joint = weighted_points(force, acceleration, [LogNormalFactor(centres, .4)])[0]
        for k, centre in enumerate(centres):
            single = weighted_points(force[k:k+1], acceleration[k:k+1], [LogNormalFactor(centre, .4)])[0, 0]
            self.assertAlmostEqual(joint[k], single, places=9)

    def test_standardize(self):
        x, y, s = standardize(np.ones((2, 3)) * 4, np.ones((2, 3)), [2., 1.], np.full((2, 3), .5))
        np.testing.assert_allclose(x[0], 2.)
        np.testing.assert_allclose(s, [4., 2.])
        with self.assertRaises(ValueError):
            standardize(np.ones((1, 3)), np.ones((1, 3)), np.array([[1., 1., 2.]]), 1.)

    def test_invalid_factors(self):
        for build in (lambda: SechTilt(float("nan")), lambda: LogNormalFactor(-1., .5),
                      lambda: LogNormalFactor(1., 0.)):
            with self.assertRaises(ValueError):
                build()


if __name__ == "__main__":
    unittest.main()
