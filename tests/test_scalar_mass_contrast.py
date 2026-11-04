"""Independent physical residual and continuous-prior T64 checks."""
import math

import numpy as np
import pytest
from scipy import integrate
from scipy.stats import norm

from met.mass_contrast_anchor import scalar_prior_projector_contrast
from met.mass_contrast_scalar import ScalarMassContrastQuadrature
from met.score_law_obstruction import ScalarScoreQuadrature, witness_constants


CONTROLS = dict(prior_center=1., prior_bounds=(.1, 10.), anchor_prefactor=.7,
                prior_nodes=96, epsabs=1e-10, epsrel=1e-9, limit=400)


def _example():
    data = [np.array([[2., 1.], [1., 2.], [-.3, .7]]),
            np.array([[1., -.2], [.4, .9], [2., 1.]])]
    covariances = [np.array([[1., .7], [.7, 1.3]]), np.diag([3., .8])]
    return data, covariances


def _physical_fields(data, covariances, mass):
    score, covariance, cost = 0., 0., 0.
    for points, sigma in zip(data, covariances):
        c, r = np.array([mass, 1.]), np.array([1., -mass])
        h = np.linalg.solve(sigma, c)
        fitted = points@h/(c@h)
        residual = points@r
        variance = r@sigma@r
        score += fitted@residual/variance
        covariance += fitted@fitted/variance
        cost += residual@residual/variance
    return score, covariance, cost


def test_grouped_scalar_fields_and_prior_variance_match_independent_physical_formulas():
    data, covariances = _example()
    q = ScalarMassContrastQuadrature([x.T@x for x in data], covariances, **CONTROLS)
    expected_u = 0.
    for points, sigma, group in zip(data, covariances, q.groups):
        coefficient = scalar_prior_projector_contrast(sigma, prior_center=1.,
            mass_bounds=(.1, 10.), epsabs=1e-12, epsrel=1e-11, limit=400)["coefficient"]
        assert group["projector_variance"] == pytest.approx(coefficient, rel=2e-12)
        expected_u += coefficient*np.einsum("ni,ij,nj->", points, np.linalg.inv(sigma), points)
    assert q.contrast_energy == pytest.approx(expected_u, rel=2e-12)
    reference_cost = _physical_fields(data, covariances, 1.)[2]
    for mass in (.15, .4, 1.2, 3., 9.):
        physical = _physical_fields(data, covariances, mass)
        score, covariance, profile = q.fields(math.atan(mass))
        assert score == pytest.approx(physical[0]*(1+mass*mass), rel=2e-13, abs=2e-14)
        assert covariance == pytest.approx(physical[1]*(1+mass*mass)**2, rel=2e-13)
        assert profile == pytest.approx(-.5*(physical[2]-reference_cost), rel=2e-13, abs=2e-14)


def test_continuous_law_cdf_matches_independent_mass_integral():
    data, covariances = _example()
    q = ScalarMassContrastQuadrature([x.T@x for x in data], covariances, **CONTROLS)
    for method in ("trace", "contrast"):
        def density(mass):
            score, covariance, cost = _physical_fields(data, covariances, mass)
            return math.exp(-.5*(score*score/covariance+q.weights[method]*cost))/(1+mass*mass)
        total = integrate.quad(density, .1, 10., epsabs=1e-12, epsrel=1e-11)[0]
        expected = integrate.quad(density, .1, .8, epsabs=1e-12, epsrel=1e-11)[0]/total
        assert q.cdf(.8, method) == pytest.approx(expected, abs=3e-11)


def test_scalar_laws_preserve_units_reciprocal_coordinates_and_exact_zero_groups():
    data, covariances = _example()
    scatters = [x.T@x for x in data]
    q = ScalarMassContrastQuadrature(scatters, covariances, **CONTROLS)
    transform = np.diag([3., .4])
    factor = 3/.4
    converted = ScalarMassContrastQuadrature([transform@s@transform.T for s in scatters],
        [transform@s@transform.T for s in covariances],
        **dict(CONTROLS, prior_center=factor, prior_bounds=(.1*factor, 10*factor)))
    reciprocal = ScalarMassContrastQuadrature([s[::-1, ::-1] for s in scatters],
        [s[::-1, ::-1] for s in covariances], **CONTROLS)
    padded = ScalarMassContrastQuadrature(scatters+[np.zeros((2, 2))],
        covariances+[np.diag([20., 1.])], **CONTROLS)
    assert converted.contrast_energy == pytest.approx(q.contrast_energy, rel=2e-13)
    assert reciprocal.contrast_energy == pytest.approx(q.contrast_energy, rel=2e-13)
    assert padded.contrast_energy == q.contrast_energy
    for method in ("trace", "contrast"):
        assert converted.cdf(.8*factor, method) == pytest.approx(q.cdf(.8, method), abs=2e-11)
        assert reciprocal.cdf(1/.8, method) == pytest.approx(1-q.cdf(.8, method), abs=2e-11)
        assert padded.cdf(.8, method) == q.cdf(.8, method)


def test_extreme_null_geometry_keeps_small_positive_mass_contrast_energy():
    ratio = 1e30
    sigma = np.diag([ratio*ratio, 1.])
    q = ScalarMassContrastQuadrature([17*sigma], [sigma], **CONTROLS)
    independent = scalar_prior_projector_contrast(sigma, prior_center=1., mass_bounds=(.1, 10.),
        epsabs=1e-12, epsrel=1e-11, limit=400)["coefficient"]
    angle_width = math.atan(10)-math.atan(.1)
    mean = .5*math.log(101/1.01)/angle_width
    second = (10-.1-angle_width)/angle_width
    expected = 34*(second-mean*mean)/(ratio*ratio)
    assert 0 < q.contrast_energy < 1e-55
    assert q.contrast_energy == pytest.approx(expected, rel=2e-12, abs=0.)
    assert q.contrast_energy == pytest.approx(34*independent, rel=2e-12, abs=0.)


def test_compact_prior_score_law_matches_existing_single_covariance_reference():
    case = witness_constants()
    scatter = 75*case["population_scatter"]
    q = ScalarMassContrastQuadrature([scatter], [case["covariance"]], **CONTROLS)
    old = ScalarScoreQuadrature(scatter, case["covariance"], prior_center=1.,
        anchor_exponent=.5, epsabs=1e-11, epsrel=1e-10, limit=400,
        anchor_weight_override=q.weights["contrast"])
    lo, hi = (old.cdf(m, "anchored_score") for m in (.1, 10.))
    assert q.cdf(.7, "contrast") == pytest.approx(
        (old.cdf(.7, "anchored_score")-lo)/(hi-lo), abs=2e-10)


def test_zero_scatter_returns_the_explicit_truncated_prior_and_rejects_invalid_inputs():
    q = ScalarMassContrastQuadrature([], [], **CONTROLS)
    expected = (math.atan(.8)-math.atan(.1))/(math.atan(10)-math.atan(.1))
    for method in ("trace", "contrast", "score"):
        assert q.cdf(.8, method) == pytest.approx(expected, abs=1e-14)
    for changes in ({"prior_bounds": (0., 10.)}, {"prior_nodes": 1}, {"anchor_prefactor": 0.}):
        with pytest.raises(ValueError):
            ScalarMassContrastQuadrature([], [], **dict(CONTROLS, **changes))


def test_weak_identifying_history_retains_the_independent_godambe_width():
    errors = []
    for count in (1e8, 1e16):
        energy = count**.6
        scatter = count*np.eye(2)+energy*np.full((2, 2), .5)
        q = ScalarMassContrastQuadrature([scatter], [np.eye(2)], **CONTROLS)
        low, high = (q.quantile(p, "contrast") for p in (.025, .975))
        # Identity-noise line angle has independent sampling SD
        # sqrt(N+Delta)/Delta. At true m=1, dm/d(angle)=2.
        sampling_sd = 2*math.sqrt(count+energy)/energy
        ratio = (high-low)/(2*norm.ppf(.975)*sampling_sd)
        errors.append(abs(ratio-1))
    assert errors[-1] < errors[0]
    assert errors[-1] < .02
