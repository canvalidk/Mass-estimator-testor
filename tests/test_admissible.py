"""B12: analytic cut-precision classifications and independent integral checks.

Finite numerical sign checks deliberately retain an unknown analytic verdict.
All per-reading calculations concern strong signal and whole-line Q.
"""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import gamma

from met.admissible import (
    check_eta, check_form, laplace, linear_gaussian, numerical_check, weight,
)
from met.full import misfit_weight


@pytest.mark.parametrize("nu", [0, 1, 2])
@pytest.mark.parametrize("alpha", [1.5, 3.0, 10.0])
@pytest.mark.parametrize("normalization", ["covariance", "precision"])
def test_inverse_gamma_control_matches_independent_precision_integrals(nu, alpha, normalization):
    rate = alpha - 1 if normalization == "covariance" else alpha
    tilted = lambda u: gamma.pdf(u, alpha + nu / 2, scale=1 / rate)
    result = check_form("invgamma", nu, alpha=alpha, normalization=normalization)
    assert result["cm"] == result["proper"] == "pass"
    mean_v, mean_u = rate / (alpha - 1), alpha / rate
    assert result["mean_v"] == pytest.approx(mean_v)
    assert result["mean_precision"] == pytest.approx(mean_u)
    assert result[normalization] == "pass"
    assert result["precision" if normalization == "covariance" else "covariance"] == "fail"
    for q in (0.0, 0.7, 5.0):
        denominator = quad(lambda u: math.exp(-q * u / 2) * tilted(u),
                           0, np.inf, epsabs=1e-11)[0]
        numerator = quad(lambda u: u * math.exp(-q * u / 2) * tilted(u),
                         0, np.inf, epsabs=1e-11)[0]
        assert laplace("invgamma", q, nu, alpha=alpha, normalization=normalization) == pytest.approx(
            denominator, rel=2e-9, abs=1e-11)
        assert weight("invgamma", q, nu, alpha=alpha, normalization=normalization) == pytest.approx(
            numerator / denominator, rel=2e-9)


@pytest.mark.parametrize("nu", [0, 1, 2])
@pytest.mark.parametrize("normalization", ["covariance", "precision"])
def test_mellin_moment_identity_including_zero_degrees_of_freedom(nu, normalization):
    alpha = 3.0
    rate = alpha - 1 if normalization == "covariance" else alpha
    transform = lambda q: laplace("invgamma", q, nu, alpha=alpha, normalization=normalization)
    if nu == 0:
        mean_v = quad(transform, 0, np.inf, epsabs=1e-10)[0] / (2 * transform(0))
        h = 1e-5
        derivative = (-3 * transform(0) + 4 * transform(h) - transform(2 * h)) / (2 * h)
        mean_u = -2 * derivative / transform(0)
    else:
        r = nu / 2
        denominator = quad(lambda q: q ** (r - 1) * transform(q),
                           0, np.inf, epsabs=1e-10)[0]
        mean_v = quad(lambda q: q ** r * transform(q), 0, np.inf, epsabs=1e-10)[0]
        mean_v /= 2 * r * denominator
        mean_u = quad(lambda q: q ** (r - 1) * transform(q) * weight(
            "invgamma", q, nu, alpha=alpha, normalization=normalization),
            0, np.inf, epsabs=1e-10)[0] / denominator
    assert mean_v == pytest.approx(rate / (alpha - 1), rel=2e-9)
    assert mean_u == pytest.approx(alpha / rate, rel=2e-9)


@pytest.mark.parametrize("nu", [0, 1, 2])
def test_fixed_gaussian_passes_both_normalizations(nu):
    result = check_form("gaussian", nu)
    assert [result[key] for key in ("cm", "proper", "covariance", "precision")] == ["pass"] * 4
    assert result["mean_v"] == result["mean_precision"] == 1
    for q in (0, 0.5, 20):
        assert weight("gaussian", q, nu) == 1
        assert laplace("gaussian", q, nu) == pytest.approx(math.exp(-q / 2))


@pytest.mark.parametrize("form,nu,threshold", [
    ("ours29", 1, 1), ("ours29", 2, 2),
    ("ours30", 0, 0.5), ("ours30", 1, 1.5), ("ours30", 2, 2.5),
    ("floor_mle", 1, 1), ("floor_mle", 2, 2),
    ("floor_map", 0, 6), ("floor_map", 1, 7), ("floor_map", 2, 8),
])
def test_capped_forms_fail_bernstein_at_second_derivative_jump(form, nu, threshold):
    result = check_form(form, nu, exponent=3)
    assert result["cm"] == result["proper"] == "fail"
    assert result["mean_v"] is None and result["mean_precision"] is None
    assert weight(form, threshold / 2, nu, exponent=3) == 1
    assert weight(form, 2 * threshold, nu, exponent=3) == pytest.approx(0.5)
    transform = lambda q: laplace(form, q, nu, exponent=3)
    middle = math.exp(-threshold / 2)
    assert transform(threshold) == pytest.approx(middle)
    assert transform(2 * threshold) == pytest.approx(middle * 2 ** (-threshold / 2))
    def jump(h):
        left = (transform(threshold) - 2 * transform(threshold - h)
                + transform(threshold - 2 * h)) / h ** 2
        right = (transform(threshold + 2 * h) - 2 * transform(threshold + h)
                 + transform(threshold)) / h ** 2
        return right - left
    h = threshold * 1e-4
    # Extrapolate away the first-order error of the one-sided differences.
    assert 2 * jump(h / 2) - jump(h) == pytest.approx(middle / (2 * threshold), rel=1e-4)


@pytest.mark.parametrize("form", ["ours29", "smooth"])
def test_zero_nu_implementation_guard_is_gaussian(form):
    result = check_form(form, 0)
    assert [result[key] for key in ("cm", "proper", "covariance", "precision")] == ["pass"] * 4
    for q in (0, 1, 10):
        assert weight(form, q, 0) == misfit_weight(q, 0, "smooth" if form == "smooth" else "kink") == 1
        assert laplace(form, q, 0) == pytest.approx(math.exp(-q / 2))


@pytest.mark.parametrize("nu", [1, 2])
def test_smooth_transform_is_gamma_but_inverse_tilt_is_not_a_probability(nu):
    r = nu / 2
    result = check_form("smooth", nu)
    assert result["cm"] == "pass" and result["proper"] == "fail"
    assert result["mean_v"] is None and result["mean_precision"] is None
    for q in (0, 0.5, 7):
        direct = quad(lambda u: math.exp(-q * u / 2) * gamma.pdf(u, r, scale=1 / r),
                      0, np.inf, epsabs=1e-10)[0]
        assert laplace("smooth", q, nu) == pytest.approx(direct, rel=2e-9)
    # Inverting the u^r tilt leaves u^-1 exp(-r*u). Infinitely many disjoint
    # shells near zero each have mass approaching log(2), so Z is infinite.
    for epsilon in (1e-2, 1e-4, 1e-6):
        shell = quad(lambda t: math.exp(-r * epsilon * t) / t, 1, 2)[0]
        assert shell >= math.exp(-2 * r * epsilon) * math.log(2)
        assert shell > 0.67


@pytest.mark.parametrize("nu", [0, 1, 2])
@pytest.mark.parametrize("exponent", [1.5, 3.0])
def test_exact_floor_cut_matches_precision_integrals_and_reports_both_moments(nu, exponent):
    result = check_form("floor_exact", nu, exponent=exponent)
    assert result["cm"] == result["proper"] == "pass"
    assert result["covariance"] == result["precision"] == "fail"
    expected_v = (exponent - 1) / (exponent - 2) if exponent > 2 else math.inf
    assert result["mean_v"] == pytest.approx(expected_v)
    assert result["mean_precision"] == pytest.approx((exponent - 1) / exponent)
    shape = exponent + nu / 2 - 1
    for q in (0, 0.7, 12):
        denominator = quad(lambda u: u ** (shape - 1) * math.exp(-q * u / 2),
                           0, 1, epsabs=1e-10)[0]
        numerator = quad(lambda u: u ** shape * math.exp(-q * u / 2),
                         0, 1, epsabs=1e-10)[0]
        assert laplace("floor_exact", q, nu, exponent=exponent) == pytest.approx(
            shape * denominator, rel=2e-9)
        assert weight("floor_exact", q, nu, exponent=exponent) == pytest.approx(
            numerator / denominator, rel=2e-9)


@pytest.mark.parametrize("nu", [1, 2])
def test_improper_floor_prior_can_have_a_finite_cut_without_a_probability_law(nu):
    result = check_form("floor_exact", nu, exponent=1)
    assert result["cm"] == "pass" and result["proper"] == "fail"
    assert result["mean_v"] is None and result["mean_precision"] is None


def test_zero_nu_improper_floor_has_no_finite_cut_transform():
    result = check_form("floor_exact", 0, exponent=1)
    assert result["cm"] == result["proper"] == "fail"
    for function in (weight, laplace):
        with pytest.raises(ValueError):
            function("floor_exact", 1, 0, exponent=1)


def test_zero_nu_floor_mle_has_no_unique_plugin():
    assert check_form("floor_mle", 0)["cm"] == "unknown"
    for function in (weight, laplace):
        with pytest.raises(ValueError):
            function("floor_mle", 0, 0)


@pytest.mark.parametrize("function", [weight, laplace])
@pytest.mark.parametrize("form,q,nu,kwargs", [
    ("unknown", 1, 1, {}), ("gaussian", -1, 1, {}), ("gaussian", 1, -1, {}),
    ("invgamma", 1, 1, {"alpha": 0}),
    ("invgamma", 1, 1, {"alpha": 1, "normalization": "covariance"}),
    ("invgamma", 1, 1, {"normalization": "unknown"}),
])
def test_analytic_forms_reject_invalid_inputs(function, form, q, nu, kwargs):
    with pytest.raises(ValueError):
        function(form, q, nu, **kwargs)


def test_finite_grid_success_does_not_prove_complete_monotonicity_or_properness():
    result = numerical_check(lambda q: 1.0, 1)
    assert result["cm"] == result["proper"] == "unknown"
    assert result["sampled_cm"] == "pass"
    assert result["mean_v"] is None and result["mean_precision"] is None


def test_finite_grid_check_can_detect_a_sign_violation():
    result = numerical_check(lambda q: 1 + q, 1)
    assert result["cm"] == "unknown"
    assert result["sampled_cm"] == "fail"


@pytest.mark.parametrize("n", [2, 5])
def test_equal_precision_eta_is_the_capped_reciprocal_residual(n):
    x = np.arange(n, dtype=float) * 2
    c = np.full(n, 3.0)
    result = linear_gaussian(x, c)
    theta = float(x.mean())
    residual = x - theta
    s = float(3 * np.sum(residual ** 2))
    assert result["theta"] == pytest.approx(theta)
    assert result["H"] == pytest.approx(3 * n)
    assert result["S"] == pytest.approx(s)
    assert result["J"] == pytest.approx(3 * s)
    assert result["nu"] == n - 1
    assert result["eta"] == pytest.approx(min(1, n / s))


def test_heterogeneous_precisions_make_eta_depend_on_more_than_shared_scale_residual():
    c = [1, 2, 3]
    first = linear_gaussian([2 * math.sqrt(2), -math.sqrt(2), 0], c)
    second = linear_gaussian([3, 0, -1], c)
    for result in (first, second):
        assert result["theta"] == pytest.approx(0, abs=1e-15)
        assert result["H"] == pytest.approx(6)
        assert result["S"] == pytest.approx(12)
    assert first["J"] == pytest.approx(16)
    assert second["J"] == pytest.approx(18)
    assert first["eta"] == pytest.approx(3 / 8)
    assert second["eta"] == pytest.approx(1 / 3)


@pytest.mark.parametrize("n", [1, 2, 5])
def test_series_eta_cut_classification_including_one_reading(n):
    result = check_eta(n)
    expected = "pass" if n == 1 else "fail"
    assert [result[key] for key in ("cm", "proper", "covariance", "precision")] == [expected] * 4


@pytest.mark.parametrize("x,c", [([], []), ([1], [0]), ([1, 2], [1]), ([1], [-1])])
def test_linear_limit_rejects_invalid_readings(x, c):
    with pytest.raises(ValueError):
        linear_gaussian(x, c)


def test_recorded_preset_reproduces_full_grid_and_cli_distinguishes_undefined_from_infinite(capsys):
    from pathlib import Path
    from met.admissible import main, run_preset

    path = Path(__file__).resolve().parents[1] / "presets" / "studies" / "admissible_tempering.yaml"
    result = run_preset(path)
    assert result == run_preset(path)
    assert result["seed"] == 20261030 and result["sampling"] is False
    assert len(result["rows"]) == 57
    assert [sum(row["k"] == k for row in result["rows"]) for k in (1, 2, 3)] == [19] * 3
    assert [row["proper"] for row in result["eta"]] == ["pass", "fail", "fail"]
    assert main([str(path)]) == 0
    output = capsys.readouterr().out
    assert "seed=20261030; deterministic (no sampling)" in output
    assert "undefined" in output and "infinity" in output
    assert sum(line.startswith(("| 1 |", "| 2 |", "| 3 |")) for line in output.splitlines()) == 57
