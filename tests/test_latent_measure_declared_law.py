"""G19's latent measure, inside the case study where how well sigma is known is handed over as q
(Claude, 17 Nov 2026; problems/O18.md; R-2026-11-17-3: q is not an input of the definition).

In noise units a reading is y = t*u + sqrt(v)*eps: t the push along the candidate's
direction u, eps standard normal in k+1 dimensions, v the noise size with the declared
law q = IG(nu/2, nu/2). B17's family puts a weight h(v) = v**b on the push's measure.

What is checked here:
1. with the weight, the reading's law given the true push is the (h*q) mixture, not the
   declared one; b = 0 is the declared Student law, whose flat-push marginal has the
   exponent (nu + k)/2 that ours33 uses;
2. the pressure test: a proper push whose spread scales with the noise size (the conjugate
   normal-inverse-gamma prior's shape) leaves the noise size's overall law at q exactly and
   still gives b = -1/2. So "the noise-size marginal is q" does not force b = 0; what does
   is that the reading's law given the true pair is the declared error law.
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats

NU, K_DIM = 10.0, 2


def _q(v, nu=NU):
    return stats.invgamma.pdf(v, a=nu / 2, scale=nu / 2)


def _quad(f):
    # the integrands are smooth and vanish faster than any power at 0; split where q peaks
    parts = [integrate.quad(f, lo, hi, epsabs=0, epsrel=1e-11, limit=400)
             for lo, hi in ((0, 0.5), (0.5, 2.0), (2.0, 50.0), (50.0, np.inf))]
    return sum(p[0] for p in parts)


def _kernel(cost, b, k=K_DIM, nu=NU):
    """The push integrated against h(v) dt = v**b dt, then v against the declared q."""
    return _quad(lambda v: v ** b * (2 * math.pi * v) ** (-k / 2) * math.exp(-cost / (2 * v)) * _q(v, nu))


@pytest.mark.parametrize("b", [0.0, -0.5, -1.0, 0.5])
def test_the_weight_sets_the_exponent(b):
    for cost in (0.3, 4.0, 60.0):
        got = math.log(_kernel(cost, b) / _kernel(0.0, b))
        assert got == pytest.approx(-(NU + K_DIM - 2 * b) / 2 * math.log1p(cost / NU), rel=1e-8)


def test_with_no_weight_the_factor_is_the_declared_student_laws_marginal():
    # the declared error law is the multivariate Student law with nu degrees of freedom;
    # its k-dimensional marginal is the Student law with the same nu
    for cost in (0.3, 4.0, 60.0):
        r = np.zeros(K_DIM)
        r[0] = math.sqrt(cost)
        student = stats.multivariate_t.pdf(r, loc=np.zeros(K_DIM), shape=np.eye(K_DIM), df=NU)
        assert _kernel(cost, 0.0) == pytest.approx(student, rel=1e-8)


def test_a_weight_replaces_the_declared_noise_size_law_given_the_true_push():
    # h = 1/v: given any true push, the noise size follows h*q = IG(nu/2 + 1, nu/2), not q
    for v in (0.4, 1.0, 2.5):
        assert _q(v) / v == pytest.approx(stats.invgamma.pdf(v, a=NU / 2 + 1, scale=NU / 2), rel=1e-12)
    assert _quad(lambda v: _q(v) / v) == pytest.approx(1.0, rel=1e-9)                  # h*q is normalised: the declared mean precision is 1
    assert _quad(lambda v: _q(v) / v / v) == pytest.approx(1 + 2 / NU, rel=1e-9)       # mean precision; declared: 1


def test_a_push_scaled_by_the_noise_size_keeps_q_overall_and_still_changes_the_factor():
    # t | v uniform on [-T sqrt(v), T sqrt(v)], v ~ q: a proper joint law whose noise-size
    # marginal is q by construction. The whole reading (the residual and the coordinate along u):
    T, y_along = 1e6, 0.7

    def reading_density(cost):
        def f(v):
            s = math.sqrt(v)
            push = (stats.norm.cdf(T - y_along / s) - stats.norm.cdf(-T - y_along / s)) / (2 * T * s)
            return _q(v) * (2 * math.pi * v) ** (-K_DIM / 2) * math.exp(-cost / (2 * v)) * push
        return _quad(f)

    for cost in (0.3, 4.0, 60.0):
        got = math.log(reading_density(cost) / reading_density(0.0))
        assert got == pytest.approx(-(NU + K_DIM + 1) / 2 * math.log1p(cost / NU), rel=1e-7)     # b = -1/2
        assert abs(got + (NU + K_DIM) / 2 * math.log1p(cost / NU)) > 0.01                        # not b = 0

    # and there the noise size's law given the true push is q(v) v**-0.5 renormalised, not q:
    norm = _quad(lambda v: _q(v) * v ** -0.5)
    assert _quad(lambda v: _q(v) * v ** -1.5) / norm == pytest.approx(1 + 1 / NU, rel=1e-9)     # mean precision
