"""Oracles for the declared multi-reading mass update rules."""

import math
import unittest

import numpy as np

from sequential_update import RULES, _combine_terms, combine, log_normal_prior, trial_terms
from weighted_posterior import infer_weighted, standardize


class SequentialUpdateTests(unittest.TestCase):
    def test_contribution_11_two_trial_readout(self):
        # Contribution 11, section 7: exact eq. (28) values for eq. (18).
        force = np.array([[2, .2, 0], [5, 0, .1]], float)
        acceleration = np.array([[1, 0, 0], [2.1, .2, 0]], float)
        sf, sa = np.array([.5, 1.]), np.array([.4, .6])
        x, y, s = standardize(force, acceleration, sf, sa)
        self.assertAlmostEqual(combine(x[:1], y[:1], s[:1], sa[:1], "mass_matching")["readout"], 1.993026469631, places=10)
        self.assertAlmostEqual(combine(x[1:], y[1:], s[1:], sa[1:], "mass_matching")["readout"], 2.373589604907, places=10)
        self.assertAlmostEqual(combine(x, y, s, sa, "mass_matching")["readout"], 2.008902161134, places=10)
        scale = np.array([[10.], [1.]])
        x, y, s = standardize(force * scale, acceleration * scale, sf * scale[:, 0], sa * scale[:, 0])
        self.assertAlmostEqual(combine(x, y, s, sa * scale[:, 0], "mass_matching")["readout"], 1.997267503687, places=10)

    def test_single_reading_is_the_flat_law_for_every_rule(self):
        rng = np.random.default_rng(5)
        force = rng.normal(size=(1, 3)) + [2, 0, 0]
        acceleration = rng.normal(size=(1, 3)) + [1, 0, 0]
        flat = infer_weighted(force, acceleration)
        for rule in RULES:
            result = combine(force, acceleration, np.ones(1), np.ones(1), rule)
            self.assertAlmostEqual(result["readout"], flat["point"][0], places=12)
            self.assertAlmostEqual(result["log_quantiles"][3], flat["distribution"]["log_median"][0], places=8)
            self.assertAlmostEqual(result["log_quantiles"][6], flat["distribution"]["log_upper_95"][0], places=8)

    def test_zero_readings(self):
        zero = np.zeros((2, 3))
        ones = np.ones(2)
        # Contribution 11, eqs. (21) and section 7: readouts s/2 and s.
        matching = combine(zero, zero, ones, ones, "mass_matching")
        self.assertAlmostEqual(matching["readout"], .5, places=10)
        cdf = lambda r: 2 / math.pi * (math.atan(r) + r / (1 + r * r))
        self.assertAlmostEqual(cdf(matching["median"]), .5, places=8)
        log_matching = combine(zero, zero, ones, ones, "log_matching")
        self.assertAlmostEqual(log_matching["readout"], 1., places=10)
        self.assertAlmostEqual(log_matching["median"], 1., places=8)
        # Readings with no information leave the symmetric posterior unchanged.
        rng = np.random.default_rng(6)
        first = (rng.normal(size=(1, 3)) + [3, 0, 0], rng.normal(size=(1, 3)) + [1, 0, 0])
        alone = combine(first[0], first[1], np.ones(1), np.ones(1), "symmetric")
        padded = combine(np.vstack([first[0], np.zeros((3, 3))]), np.vstack([first[1], np.zeros((3, 3))]),
                         np.ones(4), np.ones(4), "symmetric")
        np.testing.assert_allclose(padded["log_quantiles"], alone["log_quantiles"], atol=1e-7)
        self.assertAlmostEqual(padded["log_sd"], alone["log_sd"], places=8)

    def test_rule_identities_in_angle(self):
        rng = np.random.default_rng(7)
        count, s = 5, 2.
        x = rng.normal(size=(count, 3)) + [1.5, 0, 0]
        y = rng.normal(size=(count, 3)) + [1., 0, 0]
        u = np.linspace(-3, 4, 41)
        like, _, angle = trial_terms(x, y, np.full(count, s), u)
        # The module's own rule combination, not re-typed formulas.
        symmetric = _combine_terms(like, angle, u, "symmetric", None)
        theta = np.arctan(np.exp(u) / s)
        mass_matching = _combine_terms(like, angle, u, "mass_matching", None)
        log_matching = _combine_terms(like, angle, u, "log_matching", None)
        for rule_value, factor in ((mass_matching, 2 * np.log(np.cos(theta))),
                                   (log_matching, np.log(np.sin(2 * theta) / 2))):
            difference = rule_value - symmetric - (count - 1) * factor
            np.testing.assert_allclose(difference, difference[0], atol=1e-10)

    def test_refinement(self):
        rng = np.random.default_rng(8)
        x = rng.normal(size=(20, 3)) + [3, 0, 0]
        y = rng.normal(size=(20, 3)) + [3, 0, 0]
        for rule in RULES:
            coarse = combine(x, y, np.ones(20), np.ones(20), rule)
            fine = combine(x, y, np.ones(20), np.ones(20), rule, panels=96, order=192, step=.025)
            np.testing.assert_allclose(coarse["log_quantiles"], fine["log_quantiles"], atol=1e-6)
            self.assertAlmostEqual(coarse["readout"], fine["readout"], places=9)
            self.assertLess(coarse["boundary_panel_mass"], 1e-12)

    def test_log_normal_prior_at_zero_readings(self):
        zero = np.zeros((3, 3))
        result = combine(zero, zero, np.ones(3), np.ones(3), "symmetric", prior=log_normal_prior(4., .3))
        self.assertAlmostEqual(result["mean_log_mass"], math.log(4.), places=8)
        self.assertAlmostEqual(result["log_sd"], .3, places=8)
        with self.assertRaises(ValueError):
            combine(zero, zero, np.ones(3), np.ones(3), "mass_matching", prior=log_normal_prior(4., .3))


if __name__ == "__main__":
    unittest.main()
