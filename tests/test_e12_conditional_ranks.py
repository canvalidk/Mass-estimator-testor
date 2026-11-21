"""Finite rank calculus and the prior completion on the scalar noise-angle arc."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import chi2, ncf, norm

from explore.e12_conditional_ranks import HALF_PI, rank_measurement


def uniform_angle(x):
    return np.asarray(x)/HALF_PI


def direct_score_curvature(rows, theta):
    u = np.array([np.sin(theta), np.cos(theta)])
    v = np.array([np.cos(theta), -np.sin(theta)])
    x, y = rows @ u, rows @ v
    size = np.linalg.norm(x)
    score = x @ y / size
    curvature = (y@y-score**2)/size-size
    return score, curvature


@pytest.mark.parametrize("data", [np.array([[2., 1.], [1., 3.], [.4, -.9]]),
                                   np.array([[2., -1.]]), np.array([[2., 0.]])])
def test_branch_increment_matches_integrated_positive_derivative(data):
    measure = rank_measurement(data, discount=False)
    def density(theta):
        z, curvature = direct_score_curvature(data, theta)
        return norm.pdf(z)*max(-curvature, 0)
    for end in [.39, .98, HALF_PI]:
        if np.linalg.matrix_rank(data) == 1:
            cusp = (measure.major_angle+HALF_PI) % np.pi
            points = [cusp] if 0 < cusp < end else []
        else:
            width = np.arctan((measure.largest/measure.smallest)**.25)
            points = [(measure.major_angle+h) % np.pi for h in (-width, width)]
            points = [p for p in points if 0 < p < end]
        reference = quad(density, 0, end, points=points)[0]
        assert measure.increasing_mass(end) == pytest.approx(reference, abs=1e-10)
    assert measure.increasing_mass(HALF_PI) <= 2*ndtr(np.sqrt(measure.largest)-
                                                       np.sqrt(measure.smallest))-1+1e-12


def test_rank_folds_and_its_true_branch_has_exact_conditional_calibration():
    rows = np.diag([3., 1.])
    ranks = [ndtr(-direct_score_curvature(rows, x)[0]) for x in [0., .4, HALF_PI]]
    assert ranks[0] == pytest.approx(ranks[-1])
    assert ranks[1] < ranks[0]-.2
    # R>J at truth. Calculate the noncentral F probability from a separate
    # integral over the central chi-square variable J using a Poisson expansion for R.
    from scipy.stats import poisson
    n, energy = 4, 3.
    k = np.arange(32)
    p = poisson.pmf(k, energy/2)
    direct = quad(lambda j: (p @ chi2.sf(j, n+2*k))*chi2.pdf(j,n-1), 0, np.inf)[0]
    assert ncf.sf((n-1)/n, n,n-1,energy) == pytest.approx(direct)
    rng = np.random.default_rng(8712)
    x = rng.normal(size=(70000, n))
    x[:, 0] += np.sqrt(energy)
    v = rng.normal(size=x.shape)
    r = (x*x).sum(axis=1)
    z = (x*v).sum(axis=1)/np.sqrt(r)
    j = (v*v).sum(axis=1)-z*z
    assert np.mean(r>j) == pytest.approx(direct, abs=.005)
    assert np.mean(ndtr(-z[r>j]) < .2) == pytest.approx(.2, abs=.006)


def test_report_is_proper_swap_equivariant_zero_neutral_and_continuous_at_tie():
    data = np.array([[1., 2.], [-4., .5], [3., 1.]])
    grid = np.linspace(0, HALF_PI, 37)
    f = rank_measurement(data).cdf(grid, uniform_angle)
    assert np.all(np.diff(f)>0)
    np.testing.assert_allclose(f[[0,-1]], [0,1])
    np.testing.assert_allclose(1-rank_measurement(data[:, ::-1]).cdf(HALF_PI-grid, uniform_angle), f)
    np.testing.assert_allclose(rank_measurement(np.vstack((-data[::-1], [0,0]))).cdf(grid, uniform_angle),f)
    for data in [[], np.eye(2), np.diag([1, 1+1e-9])]:
        np.testing.assert_allclose(rank_measurement(data).cdf(grid, uniform_angle), uniform_angle(grid),
                                   atol=1e-8, rtol=0)


def test_fixed_conflict_retains_spread_and_unmodified_rank_collapses():
    rows = np.linalg.cholesky([[5., 4.],[4., 5.]]).T
    grid = np.linspace(0, HALF_PI, 21)
    at_100 = rank_measurement(rows*100)
    at_1000 = rank_measurement(rows*1000)
    np.testing.assert_allclose(at_100.cdf(grid, uniform_angle), at_1000.cdf(grid, uniform_angle),
                               atol=1e-13)
    assert at_1000.prior_weight > 0
    close = np.pi/4+np.array([-.01,.01])
    assert np.diff(at_1000.cdf(close, uniform_angle)).item() < .3
    assert np.diff(rank_measurement(rows*1000, discount=False).cdf(close, uniform_angle)).item() > .9999


def test_axis_data_keep_prior_weight_and_identifying_local_coordinate_is_normal():
    axis = rank_measurement([[100., 0.]])
    assert axis.discount == pytest.approx(.01)
    assert axis.increasing_mass(HALF_PI) == pytest.approx(ndtr(1)-.5)
    n = 200
    rng = np.random.default_rng(173)
    theta = .6
    data = rng.normal(size=(n, 2))+40*np.array([np.sin(theta),np.cos(theta)])
    law = rank_measurement(data)
    assert law.discount == 1
    scale = (law.largest-law.smallest)/np.sqrt(law.largest)
    z = np.arange(-3,4)
    np.testing.assert_allclose(law.cdf(law.major_angle+z/scale, uniform_angle), ndtr(z), atol=2e-6)
