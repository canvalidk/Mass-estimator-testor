"""Physical reference Gaussian ratios on a candidate graph and their conflict cost."""
import numpy as np
import pytest
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad
from scipy.stats import chi2


def normal_nodes(dim):
    x, w = hermgauss(32)
    nodes = np.array(np.meshgrid(*([x*np.sqrt(2)]*dim), indexing="ij")).reshape(dim, -1).T
    weights = np.prod(np.array(np.meshgrid(*([w/np.sqrt(np.pi)]*dim), indexing="ij")), axis=0)
    return nodes, weights.ravel()


def log_divisor(y, rank, tilt):
    p = y.shape[-1]
    cost = np.sum(y[..., rank:]**2, axis=-1)
    return (-p*np.log(2*np.pi) - rank*np.log(2)
            - (p-rank)*np.log(2-tilt) - tilt*cost) / 2


def log_density(y, mean):
    return -(y.shape[-1]*np.log(2*np.pi) + np.sum((y-mean)**2, axis=-1)) / 2


@pytest.mark.parametrize("tilt", [0., .3, .8, 1.])
def test_reference_comparison_integrates_to_exact_bound_at_arbitrary_drives(tilt):
    nodes, weights = normal_nodes(3)
    a = np.array([.5, 0., 0.])
    y = nodes + a
    for b in [a, [.5, -.6, .9], [-.8, -.6, .9]]:
        b = np.asarray(b)
        mean = weights @ np.exp(log_density(y, b) - log_divisor(y, 1, tilt))
        expected = np.exp(-(a[0]-b[0])**2/4
                          -(1-tilt)*np.dot(b[1:], b[1:])/(2*(2-tilt)))
        assert mean == pytest.approx(expected, rel=2e-11)
    actual_log_score = weights @ (log_density(y, a) - log_divisor(y, 1, tilt))
    optimal_log_score = (np.log(2) - 1) / 2
    assert optimal_log_score - actual_log_score == pytest.approx(1-tilt-np.log(2-tilt))


def test_double_smoothing_recovers_matching_null_for_nonradial_competitor():
    nodes, ws = normal_nodes(2)
    x, wh = normal_nodes(1)
    center = .7
    # D = D_star / s is valid for every reference because 0 < s <= 1.
    # s depends on both graph and transverse coordinates.
    def s(u, v):
        return .7 + .2*np.cos(u + .6*v)
    half_drive = center + x[:, 0] / np.sqrt(2)
    local_means = [ws @ s(h+nodes[:, 0]/np.sqrt(2), nodes[:, 1]) for h in half_drive]
    actual = s(center+nodes[:, 0], nodes[:, 1])
    assert wh @ local_means == pytest.approx(ws @ actual, abs=1e-12)
    regret = -(ws @ np.log(actual))
    assert regret > 0.
    assert ws @ (actual-1-np.log(actual)) <= regret


def test_common_divisor_costs_same_against_each_oracle_reference():
    z, w = normal_nodes(3)
    for b in [[0., .9, -.5], [-1.2, .3, 0.]]:
        b = np.array(b)
        y = z+b
        projection = b.copy()
        projection[2:] = 0.
        oracle = log_density(y, b) - log_density(y, projection)
        reusable = log_density(y, b) - log_divisor(y, 2, 1.)
        assert w @ oracle == pytest.approx(b[2]**2/2)
        assert w @ (oracle-reusable) == pytest.approx(1-np.log(2))


def test_tempering_budget_and_candidate_order_at_finite_conflict():
    epsilon = .01
    for j in (1, 2, 17):
        cutoff = max(j, chi2.isf(2*epsilon/j, j+2))
        tail_cost = quad(lambda u: u*chi2.pdf(u, j)/2, cutoff, np.inf)[0]
        assert tail_cost <= epsilon + 1e-12
    cutoff = chi2.isf(epsilon, 4)
    rows = np.array([[2., 1.], [1., 1.]])
    eigenvalues = np.linalg.eigvalsh(rows.T @ rows)
    q = eigenvalues[0]
    assert np.diff(eigenvalues)[0] == pytest.approx(3*np.sqrt(5))
    # A log-odds difference against the equal mass candidate.
    residual = rows @ (np.array([1., -1.])/np.sqrt(2))
    difference = residual @ residual - q
    assert difference > 0.
    tempered, full = [], []
    for sigma in [.01, .001, .0001]:
        damping = min(1., cutoff*sigma**2/q)
        full.append(difference/(2*sigma**2))
        tempered.append(damping*full[-1])
    np.testing.assert_allclose(tempered, cutoff*difference/(2*q))
    assert full[-1] == pytest.approx(10000*full[0])


def test_block_optimum_is_the_product_of_single_reading_optima():
    rows = np.array([[.4, -.4], [-.7, 1.2], [.1, -.6]])
    stacked = np.r_[rows[:, 0], rows[:, 1]]
    assert log_divisor(stacked, 3, 1.) == pytest.approx(log_divisor(rows, 1, 1.).sum())
