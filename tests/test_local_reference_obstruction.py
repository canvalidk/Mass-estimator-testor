"""T66 exact finite geometry, not a simulated T65 full history."""
import math

import numpy as np
import pytest

from met.joint_mass import common_acceleration
from met.pairwise_anchor import scalar_local_projection_energy


def _density(m):
    return 2/(math.pi*(1+m*m))


def _scatter_events(scatter, covariance):
    factor = np.linalg.cholesky(scatter)
    return [common_acceleration([[row[0]]], [row[1]], covariance) for row in factor.T]


@pytest.mark.parametrize("mass", [.3, 1.2, 3.])
def test_prior_rank_derivative_matches_finite_differences_and_coordinate_change(mass):
    cov = np.array([[1., .9], [.9, 1.]])
    event = common_acceleration([[1.4]], [.7], cov)
    w = np.linalg.inv(np.linalg.cholesky(cov))
    def projection(q):
        m = math.tan(math.pi*q/2)
        v = w @ np.array([m, 1.])
        return np.outer(v, v)/(v@v) @ w @ event.data
    rank = 2/math.pi*math.atan(mass)
    derivative = (projection(rank+1e-6)-projection(rank-1e-6))/(2e-6)
    actual = scalar_local_projection_energy([event], mass, prior_density=_density(mass))
    assert actual == pytest.approx(derivative@derivative, rel=2e-9)
    # A passive force-unit change scales m and changes its prior density.
    transform = np.diag([3., 1.])
    changed = common_acceleration([[3*event.data[0]]], [event.data[1]],
                                  transform@cov@transform.T)
    assert scalar_local_projection_energy([changed], 3*mass,
        prior_density=_density(mass)/3) == pytest.approx(actual, rel=3e-14)


def test_local_state_has_summable_null_envelope_and_correct_informative_constant():
    m = 1.2
    for index in (10, 100, 1000):
        event = common_acceleration([[index]], [2.], np.diag([index**2, 1.]))
        value = scalar_local_projection_energy([event], m, prior_density=_density(m))
        assert value <= 5/(index*_density(m))**2
    cov = np.array([[1., .9], [.9, 1.]])
    mean = np.array([1.2, 1.])
    events = _scatter_events(cov+np.outer(mean, mean), cov)
    for mass in (.3, 1.2, 3.):
        q = mass*mass-1.8*mass+1
        c = 33*math.pi**2/200*(1+mass*mass)**2/q**2
        assert scalar_local_projection_energy(events, mass, prior_density=_density(mass)) == pytest.approx(c, rel=3e-14)


def test_remote_profile_reference_can_diverge_while_local_fields_stay_fixed():
    mass, count = 1.2, 400.
    cov = np.array([[1., .9], [.9, 1.]])
    mean = np.array([mass, 1.])
    events = _scatter_events(count*(cov+np.outer(mean, mean)), cov)
    left, right = mass-.5/math.sqrt(count), mass+.5/math.sqrt(count)
    def fields(data, m):
        energy, explained = 0., 0.
        for event in data:
            w = np.linalg.inv(np.linalg.cholesky(event.covariance))
            z, v = w@event.data, w@np.array([m, 1.])
            energy += z@z
            explained += (v@z)**2/(v@v)
        local = scalar_local_projection_energy(data, m, prior_density=_density(m))
        return .5*explained, 1/math.sqrt(1+local)
    baseline = [fields(events, m) for m in (left, right)]
    bounds = []
    for radius in (20., 100., 1000.):
        # e=R^-2, whitened data (R,0); this is a finite allowed input.
        remote = common_acceleration([[radius**3]], [0.], np.diag([radius**4, 1.]))
        data = events+[remote]
        values = [fields(data, m) for m in (left, right)]
        np.testing.assert_allclose([v[1] for v in values], [v[1] for v in baseline], rtol=1e-3)
        lleft, wleft = values[0]
        lright, wright = values[1]
        assert wright > wleft
        # Sup L is at least its boundary limit; for a positive weight
        # difference this gives an UPPER bound on right/left log odds.
        boundary = .5*sum((e.data[0]**2/e.covariance[0, 0]) for e in [remote])
        bound = -(wright-wleft)*boundary+wright*lright-wleft*lleft
        bounds.append(bound)
    assert bounds[2] < bounds[1] < bounds[0]
    assert bounds[-1] < -100.
