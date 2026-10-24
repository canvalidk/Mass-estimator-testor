"""Independent identities and integration comparisons for the deliverable."""
import sys
from pathlib import Path
import math
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'pass to VD-docs/Mass estimator equation'))
sys.path.insert(0, str(ROOT / '03_trace_and_evaluator/mass_estimator'))
sys.path.insert(0, str(ROOT / '03_trace_and_evaluator/mass_estimator/investigation'))
from mass_estimator import estimate_mass
from estimator import estimate as reference_estimate
from practical_formulas import summarize
from comparison_posterior import infer_batch


def covariance(d, sf=1., sa=1.):
    return np.diag([sf**2] * d + [sa**2] * d)


class Checks(unittest.TestCase):
    def test_analytic_null_all_dimensions(self):
        for d in (1, 2, 3):
            r = estimate_mass(np.zeros(d), np.zeros(d), covariance(d, 2, .5))
            self.assertAlmostEqual(r.mass, 4., places=10)
            np.testing.assert_allclose(r.interval, 4*np.tan(np.pi/2*np.array([.025, .975])), rtol=1e-11)
            self.assertAlmostEqual(r.log_standard_deviation, np.pi/2, places=8)
            self.assertAlmostEqual(r.mean_force_magnitude, 2*np.sqrt(2/np.pi), places=10)

    def test_correlated_null_analytic(self):
        for c in (-.7, .6):
            cov = np.array([[4., c], [c, .25]])
            r = estimate_mass([0], [0], cov, rtol=2e-5)
            beta = np.pi/2 + np.arcsin(c)
            expected = 4*(c+np.sqrt(1-c*c)*np.tan(beta*np.array([.025, .975])-np.arcsin(c)))
            np.testing.assert_allclose(r.interval, expected, rtol=2e-5)
            self.assertAlmostEqual(r.mass, 4, places=10)
            self.assertAlmostEqual(r.log_standard_deviation, np.sqrt((np.pi**2-beta**2)/3), places=5)

    def test_known_direction_truncated_normals(self):
        def mean(x, s):
            z = x/s
            return x + s*np.exp(-z*z/2)/np.sqrt(2*np.pi)/(.5*math.erfc(-z/np.sqrt(2)))
        r = estimate_mass([-1, 7, 0], [2, -3, 0], covariance(3, 1, .5), known_direction=[1, 0, 0])
        self.assertAlmostEqual(r.mass, mean(-1, 1)/mean(2, .5), places=7)
        self.assertAlmostEqual(r.mean_force_magnitude, mean(-1, 1), places=7)

    def test_latest_isotropic_backend(self):
        for f, a in (([4, 0, 0], [1, 1, 0]), ([3, 0, 0], [-3, 0, 0]),
                     ([0, 0, 0], [8, 0, 0]), ([20, 0, 0], [20, 1, 0])):
            r = estimate_mass(f, a, np.eye(6), rtol=3e-5, max_refinements=4)
            v = infer_batch(f, a, order=192)
            law = v['distributions']['flat_joint']
            self.assertAlmostEqual(r.mass, float(v['points']['flat_joint'][0]), places=7)
            expected = np.exp([law['log_lower_95'][0], law['log_upper_95'][0]])
            np.testing.assert_allclose(r.interval, expected, rtol=3e-5)
            self.assertAlmostEqual(r.log_standard_deviation, float(law['log_sd'][0]), places=5)

    def test_general_covariance_original_reference(self):
        cov = np.array([[1., .2, .1, 0.], [.2, .7, .05, -.08],
                        [.1, .05, .4, .04], [0., -.08, .04, .6]])
        r = estimate_mass([3, 1], [.8, .5], cov, rtol=2e-5)
        old = summarize(reference_estimate([3, 1], [.8, .5], cov,
                                          direction_order=256, ratio_order=4096))
        self.assertAlmostEqual(r.mass, old.mass, places=8)
        np.testing.assert_allclose(r.interval, [old.lower, old.upper], rtol=3e-5)
        self.assertAlmostEqual(r.log_standard_deviation, old.log_standard_deviation, places=5)

    def test_original_coordinate_integral(self):
        # Independent 2D f,alpha,u integral, with no radial reduction.
        cov = np.array([[1., .1, .1, 0], [.1, .8, 0, -.05],
                        [.1, 0, .5, .1], [0, -.05, .1, .7]])
        f, a = np.array([2., .3]), np.array([.7, .8])
        x, w = np.polynomial.legendre.leggauss(72)
        mag, w = 6*(x+1), 6*w  # [0,12], negligible Gaussian tail for this case
        ff, aa = np.meshgrid(mag, mag, indexing='ij')
        ww = np.outer(w, w)
        totals = np.zeros(3)
        prec = np.linalg.inv(cov)
        for phi in np.arange(192)*2*np.pi/192:
            u = np.array([np.cos(phi), np.sin(phi)])
            residual = np.r_[f, a] - np.concatenate((ff[..., None]*u, aa[..., None]*u), axis=2)
            mass = np.exp(-.5*np.einsum('...i,ij,...j->...', residual, prec, residual))*ww
            totals += [mass.sum(), (mass*ff).sum(), (mass*aa).sum()]
        r = estimate_mass(f, a, cov, rtol=2e-5)
        np.testing.assert_allclose([r.mass, r.mean_force_magnitude, r.mean_acceleration_magnitude],
                                  [totals[1]/totals[2], totals[1]/totals[0], totals[2]/totals[0]], rtol=1e-7)

    def test_rotation_reciprocity_and_units(self):
        f, a = np.array([2., 1.]), np.array([.6, .1])
        cov = np.array([[1., .2, .1, 0.], [.2, .7, .05, -.08],
                        [.1, .05, .4, .04], [0., -.08, .04, .6]])
        base = estimate_mass(f, a, cov, rtol=2e-5)
        angle = .37
        rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        transform = np.zeros((4, 4)); transform[:2, :2] = rot; transform[2:, 2:] = rot
        r = estimate_mass(rot@f, rot@a, transform@cov@transform.T, rtol=2e-5)
        np.testing.assert_allclose([r.mass, *r.interval, r.log_standard_deviation],
                                  [base.mass, *base.interval, base.log_standard_deviation], rtol=2e-5)
        indices = [2, 3, 0, 1]
        swapped = estimate_mass(a, f, cov[np.ix_(indices, indices)], rtol=2e-5)
        np.testing.assert_allclose([swapped.mass, *swapped.interval],
                                  [base.inverse_mass, *base.inverse_mass_interval], rtol=2e-5)
        units = np.diag([1000., 1000., .01, .01])
        scaled = estimate_mass(1000*f, .01*a, units@cov@units, rtol=2e-5)
        np.testing.assert_allclose([scaled.mass, *scaled.interval],
                                  100000*np.array([base.mass, *base.interval]), rtol=2e-5)

    def test_3d_anisotropic(self):
        cov = np.diag([1., .7, 1.5, .5, .3, .8])
        cov[0, 3] = cov[3, 0] = .2
        r = estimate_mass([2, .4, -.3], [.8, .1, .2], cov)
        old = summarize(reference_estimate([2, .4, -.3], [.8, .1, .2], cov,
                                          direction_order=48, ratio_order=2048))
        np.testing.assert_allclose([r.mass, *r.interval, r.log_standard_deviation],
                                  [old.mass, old.lower, old.upper, old.log_standard_deviation], rtol=2e-4)

    def test_invalid_and_nonconvergence(self):
        bad = [dict(force=[1], acceleration=[1], covariance=[[0, 0], [0, 1]]),
               dict(force=[np.nan], acceleration=[1], covariance=np.eye(2)),
               dict(force=[1], acceleration=[1], covariance=np.eye(2), probability=1),
               dict(force=[1], acceleration=[1], covariance=np.eye(2), known_direction=[2])]
        for kw in bad:
            with self.assertRaises(ValueError): estimate_mass(**kw)
        with self.assertRaises(ArithmeticError):
            estimate_mass([20, 0, 0], [20, 1, 0], np.eye(6), ratio_order=32,
                          rtol=1e-14, max_refinements=1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
