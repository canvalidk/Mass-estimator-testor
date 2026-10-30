"""B02: exact controls, independent integrals and checks of the supremum proof.

The proof is in O13; finite grids and seeded Monte Carlo are checks, not proofs
of a supremum. No candidate weight or production equation is changed here.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import hyp0f1
from scipy.stats import norm

from tools.positivity_supremum import expectations, quadrature_error_bound


NULL = {
    1: (0.5 - 1 / math.pi, 0.0),
    2: (1.5 - math.pi / 4, 2 - math.pi / 2),
    3: (2 - 2 / math.pi, 1.0),
}


@pytest.mark.parametrize("k", (1, 2, 3))
def test_exact_null_controls(k):
    result = expectations(k, 0.0, 0.0, tolerance=1e-11)
    positive, unsigned = NULL[k]
    allowance = result["quadrature_error_bound"] + 2e-13
    assert result["positive"] == pytest.approx(positive, abs=allowance)
    assert result["unsigned"] == pytest.approx(unsigned, abs=allowance)
    assert not result["floating_point_enclosure"]


@pytest.mark.parametrize("r,b", ((0., 0.), (0., 2.), (1., .5), (4., .25), (3., 2.)))
def test_one_dimension_against_independent_normal_survival_integral(r, b):
    # E min(X²,Y²) 1{XY<0} by the tail integral; no Gram-angle formula.
    weak = min(r, b)
    limit = weak + 12
    integrand = lambda u: 2 * u * (
        norm.cdf(r - u) * norm.cdf(-b - u)
        + norm.cdf(-r - u) * norm.cdf(b - u))
    reference, integration_error = quad(integrand, 0, limit, epsabs=1e-12, epsrel=1e-12)
    # Tail <= E[(weak+Z)^2 1{|weak+Z|>limit}], an explicit upper bound.
    tail = ((1 + weak * weak) * (norm.sf(limit - weak) + norm.sf(limit + weak))
            + (limit + weak) * norm.pdf(limit - weak)
            + (limit - weak) * norm.pdf(limit + weak))
    result = expectations(1, r, b, tolerance=1e-11)
    assert result["positive"] == pytest.approx(
        reference, abs=result["quadrature_error_bound"] + integration_error + tail + 2e-12)
    assert result["unsigned"] == 0


@pytest.mark.parametrize("k", (1, 2, 3))
@pytest.mark.parametrize("r,b", ((0., 2.), (.25, 4.), (1., 3.), (16., .5)))
def test_swap_symmetry_is_checked_by_separate_quadratures(k, r, b):
    left, right = expectations(k, r, b), expectations(k, b, r)
    allowance = left["quadrature_error_bound"] + right["quadrature_error_bound"] + 2e-11
    for quantity in ("positive", "unsigned"):
        assert left[quantity] == pytest.approx(right[quantity], abs=allowance)


@pytest.mark.parametrize("k", (1, 2, 3))
def test_fixed_total_strength_is_worst_at_the_boundary(k):
    length = 4.0
    boundary = expectations(k, length, 0.0, tolerance=1e-11)
    values = []
    for alpha in np.linspace(0, math.pi / 4, 9):
        result = expectations(k, length * math.cos(alpha), length * math.sin(alpha),
                              tolerance=1e-11)
        values.append(result["positive"])
        # Eigenvalues forget the orientation of the two-channel mean.
        assert result["unsigned"] == pytest.approx(boundary["unsigned"], abs=2e-10)
    assert np.all(np.diff(values) < 1e-10)
    assert values[-1] < values[0]


@pytest.mark.parametrize("k", (1, 2, 3))
def test_boundary_is_increasing_and_approaches_both_suprema(k):
    strengths = (0., .5, 1., 2., 4., 8., 16., 32.)
    values = [expectations(k, r, 0.) for r in strengths]
    assert np.all(np.diff([row["positive"] for row in values]) > -2e-10)
    assert np.all(np.diff([row["unsigned"] for row in values]) > -2e-10)
    last = values[-1]
    rounding = 2e-10
    assert last["positive"] <= k - .5 + last["quadrature_error_bound"] + rounding
    assert last["positive"] > k - .5 - .002
    if k == 1:
        assert all(row["unsigned"] == 0 for row in values)
    else:
        assert last["unsigned"] < k - 1
        assert last["unsigned"] > k - 1 - .003


@pytest.mark.parametrize("k", (1, 2, 3))
def test_circular_arc_mass_decreases_away_from_kernel_centre(k):
    # The elementary layer-cake step of the proof, on a nontrivial kernel.
    kernel = lambda s: hyp0f1(k / 2, (0.2 + 2.8 * math.cos(s) ** 2))
    for radius in (.05, .2, .6):
        arc_masses = [quad(kernel, d - radius, d + radius, epsabs=1e-11)[0]
                      for d in np.linspace(math.pi / 4, math.pi / 2, 13)]
        assert np.all(np.diff(arc_masses) < 2e-10)


@pytest.mark.parametrize("k", (1, 2, 3))
def test_strict_boundary_gap_contains_an_open_event(k):
    x, y = np.zeros(k), np.zeros(k)
    x[0], y[0] = 1.0, -2.0
    a, b, c = x @ x, y @ y, x @ y
    ray_cost = b - max(c, 0) ** 2 / a
    positive_cost = min(a, b)
    assert c < 0 and a < b
    assert positive_cost < ray_cost
    # These inequalities persist in an open neighbourhood, of positive
    # probability under every finite Gaussian mean.


@pytest.mark.parametrize("k", (1, 2, 3))
def test_seeded_matrix_draws_crosscheck_the_integral(k):
    # Monte Carlo is only an independent finite-cell cross-check.
    rng = np.random.default_rng(20261030 + k)
    count, r, b = 120000, 3.0, .5
    x, y = rng.normal(size=(2, count, k))
    x[:, 0] += r
    y[:, 0] += b
    singular = np.linalg.svd(np.stack((x, y), axis=-1), compute_uv=False)
    unsigned = singular[:, -1] ** 2 if k > 1 else np.zeros(count)
    a, bb, c = (x * x).sum(1), (y * y).sum(1), (x * y).sum(1)
    positive = np.where(c > 0, unsigned, np.minimum(a, bb))
    result = expectations(k, r, b, tolerance=1e-11)
    for sample, quantity in ((positive, "positive"), (unsigned, "unsigned")):
        standard_error = sample.std(ddof=1) / math.sqrt(count)
        assert abs(sample.mean() - result[quantity]) < 7 * standard_error + 2e-10


@pytest.mark.parametrize("k", (1, 2, 3))
def test_discretization_certificate_and_refinement(k):
    coarse = expectations(k, 8., 3., tolerance=1e-7)
    fine = expectations(k, 8., 3., tolerance=1e-12)
    assert coarse["quadrature_error_bound"] <= 1e-7
    assert fine["quadrature_error_bound"] <= 1e-12
    assert fine["nodes"] >= coarse["nodes"]
    assert abs(coarse["positive"] - fine["positive"]) <= (
        coarse["quadrature_error_bound"] + fine["quadrature_error_bound"] + 2e-10)
    error = quadrature_error_bound(k, 73., 40, 1.4)
    refined = quadrature_error_bound(k, 73., 48, 1.4)
    assert refined == pytest.approx(error * 1.4 ** -16)


def test_reference_refuses_an_unmet_error_budget():
    with pytest.raises(ValueError, match="more nodes"):
        expectations(3, 32., 32., tolerance=1e-12, max_nodes=16)
