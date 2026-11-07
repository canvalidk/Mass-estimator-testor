"""Forward eta32 keeps the production residual and stationary-line semantics."""
import math

import numpy as np
import pytest

from met.forward import _isotropic_eta32, forward, from_invariants
from met.full import isotropic
from met.ours31 import _per_reading
from met.ours32 import bartlett_eta


@pytest.mark.parametrize("c", [0.1, 1.0, 10.0])
def test_high_push_forward_eta_matches_general_engine(c):
    truths = from_invariants(3, 1e4, math.atan(2.0), readings=20)
    fast = forward(truths, draws=3, seed=90210, nu_declared=10.0, c=c, engine="isotropic")
    general = forward(truths, draws=3, seed=90210, nu_declared=10.0, c=c, engine="general")
    np.testing.assert_allclose(fast["eta32"], general["eta32"], rtol=2e-7, atol=2e-8)
    assert not any(key.startswith("_") for key in fast)


def test_forward_eta_matches_exact_symmetric_pair():
    theta, strength, delta, nu = math.atan(2.0), 1000.0, 3.0, 10.0
    along = np.array([math.sin(theta), math.cos(theta)])
    normal = np.array([math.cos(theta), -math.sin(theta)])
    observed = np.array([strength * along + sign * delta * normal for sign in (-1, 1)])
    x, y = observed[None, :, :1], observed[None, :, 1:]
    gram = ((x*x).sum(-1), (y*y).sum(-1), (x*y).sum(-1))
    actual = _isotropic_eta32(gram, np.zeros((1, 2)), 1, nu, observations=(x, y))[0]
    expected = (nu - delta**2) / (nu * delta**2) - (nu + delta**2) / (nu * strength**2)
    assert actual == pytest.approx(expected, rel=1e-10)


def test_large_antiparallel_forward_factor_stays_finite():
    x, y = np.array([[[2e8]]]), np.array([[[-1e8]]])
    gram = ((x*x).sum(-1), (y*y).sum(-1), (x*y).sum(-1))
    with np.errstate(all="raise"):
        actual = _isotropic_eta32(gram, np.array([[1e16]]), 1, 0.5, observations=(x, y))[0]
    reading = isotropic(x[0, 0], y[0, 0], 1.0, 1.0)
    reference = bartlett_eta(_per_reading([reading], 1.0, [0.5]), 1.0)[0]
    assert math.isfinite(actual) and actual == pytest.approx(reference)
