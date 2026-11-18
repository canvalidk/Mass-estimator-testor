"""Finite witnesses for explore/e12-report.md. No grid truncation of mass endpoints."""

import math

import numpy as np
import pytest
from scipy import integrate, optimize, special

from explore.e12_scalar_witness import pooled_witness
from met.conditional_score import gaussian_score_accepts, gaussian_score_test
from met.full import Reading


P = 0.95
Z = special.ndtri((1.0 + P) / 2.0)


def _rows(points):
    return [Reading([f], [a], np.eye(2)) for f, a in points]


def _statistic(gram, theta):
    u = np.array([math.sin(theta), math.cos(theta)])
    v = np.array([math.cos(theta), -math.sin(theta)])
    variance = u @ gram @ u
    if variance == 0:
        return 0.0
    return (v @ gram @ u) ** 2 / variance


def _ends(gram):
    """Strict scores at the axes; all examples below are off the threshold."""
    cross_squared = gram[0, 1] ** 2
    return cross_squared < Z**2 * gram[1, 1], cross_squared < Z**2 * gram[0, 0]


def _central_angles(density, q):
    total = integrate.quad(density, 0, math.pi / 2)[0]
    angles = []
    for rank in ((1 - q) / 2, (1 + q) / 2):
        def residual(x):
            return integrate.quad(density, 0, x)[0] - rank * total
        angles.append(optimize.brentq(residual, 0, math.pi / 2))
    return np.array(angles)


def test_all_zero_and_single_perpendicular_candidate_accept():
    for observations in ([], _rows([(0, 0)]), _rows([(1, -1)])):
        test = gaussian_score_test(observations, 1.0)
        assert test["degenerate"] and test["statistic"] == 0.0
    assert gaussian_score_accepts(_rows([(0, 0)]), [0.01, 1.0, 100.0],
                                 probability=P).all()


def test_connected_hull_can_be_whole_despite_a_hole_in_the_score_set():
    readings = _rows([(4, 0), (0, 1)])
    gram = np.diag([16.0, 1.0])
    assert _ends(gram) == (True, True)
    test = gaussian_score_test(readings, 1.0)
    assert test["statistic"] == pytest.approx(225 / 34)
    assert test["statistic"] > Z**2


@pytest.mark.parametrize("theta", [0.3, math.pi / 4, 1.1])
def test_strict_overcoverage_neighborhood_at_an_interior_truth(theta):
    # An interior point of the open event where the score rejects the true
    # candidate yet accepts both boundaries, including after a small extra row.
    for points in ([(0.95 * Z, -0.9 * Z)], [(0.95 * Z, -0.9 * Z), (0.01, 0.02)]):
        rows = np.array(points)
        assert _ends(rows.T @ rows) == (True, True)
        result = gaussian_score_test(_rows(points), math.tan(theta))
        assert not result["degenerate"]
        assert result["statistic"] > Z**2


def test_isotropic_blocks_lose_the_state_for_an_updated_interval_hull():
    extra = (3, -1)
    ends = []
    for scale in (1, 3):
        root = math.sqrt(scale)
        points = [(root, 0), (0, root)]
        gram = np.eye(2) * scale
        assert _ends(gram) == (True, True)
        for theta in (0.25, 0.8, 1.4):
            assert _statistic(gram, theta) == pytest.approx(0, abs=1e-28)
        updated = gram + np.outer(extra, extra)
        ends.append(_ends(updated))
        accepted = gaussian_score_accepts(_rows(points + [extra]), [1e-6, 1e6],
                                          probability=P)
        assert tuple(accepted) == ends[-1]
    assert ends == [(False, True), (True, True)]


def test_two_finite_reports_on_the_same_compatible_reading_are_distinct():
    radius = 4.0
    centre = math.pi / 4
    score_arc = centre + np.array([-1, 1]) * math.asin(Z / radius)
    reading = _rows([(radius / math.sqrt(2), radius / math.sqrt(2))])
    for angle in score_arc:
        actual = gaussian_score_test(reading, math.tan(angle))["statistic"]
        assert actual == pytest.approx(Z**2)
    # A single compatible reading at face value, c=s, with the angular prior.
    def density(angle):
        return math.exp(-0.5 * radius**2 * math.sin(angle - centre)**2)
    hulls = []
    for k in (0, 1000):
        q = P + (1 - P) * k / (1 + k + radius**2)
        inner = _central_angles(density, q)
        hulls.append([min(inner[0], score_arc[0]), max(inner[1], score_arc[1])])
    assert 0 < hulls[1][0] < hulls[0][0] < hulls[0][1] < hulls[1][1] < math.pi / 2


@pytest.mark.parametrize("share", [0.03, 0.5, 1.0])
def test_local_score_and_adaptive_raw_interval_endpoints_coalesce(share):
    # Reachable scatters with the same interior principal angle, including a
    # rank-one case and a small eigengap relative to the leading eigenvalue.
    centre = math.pi / 4
    for information in (100.0, 6400.0):
        leading = information / share**2
        trailing = leading * (1.0 - share)
        gram = np.array([[leading + trailing, leading - trailing],
                         [leading - trailing, leading + trailing]]) / 2.0
        score_lo = optimize.brentq(lambda x: _statistic(gram, x) - Z**2, 0, centre)
        q = P + (1 - P) * 3 / (1 + 3 + information)

        def density(angle):
            return math.exp(-0.5 * information * math.sin(angle - centre)**2)

        raw_lo, _ = _central_angles(density, q)
        scaled = math.sqrt(information) * (centre - np.array([score_lo, raw_lo]))
        tolerance = 0.10 if information == 100 else 0.002
        np.testing.assert_allclose(scaled, Z, atol=tolerance, rtol=0)


def test_conflicting_fixed_readings_shrink_the_score_arc_but_retain_raw_spread():
    physical = np.array([[0.2, 1], [1, 1], [1.4, 1]])
    gram = physical.T @ physical
    np.testing.assert_allclose(gram, [[3, 2.6], [2.6, 3]])
    centre = math.pi / 4
    # The pooled raw angular density is independent of sigma in the cap's tail.
    def density(angle):
        u = np.array([math.sin(angle), math.cos(angle)])
        return math.exp((3 / (2 * 0.4)) * (u @ gram @ u - 5.6))
    interval = _central_angles(density, P)
    score_widths = []
    for sd in (0.05, 0.005):
        whitened = gram / sd**2
        assert _ends(whitened) == (False, False)
        left = optimize.brentq(lambda x: _statistic(whitened, x) - Z**2, 0, centre)
        right = 2 * centre - left
        assert interval[0] < left < centre < right < interval[1]
        readings = [Reading([f], [a], np.eye(2) * sd**2) for f, a in physical]
        assert gaussian_score_test(readings, math.tan(left))["statistic"] == pytest.approx(Z**2)
        score_widths.append(right - left)
    assert score_widths[1] < 0.11 * score_widths[0]


@pytest.mark.parametrize("n,whole", [(1, 0.90250), (5, 0.91402), (20, 0.92996)])
def test_exact_finite_null_whole_domain_and_active_weight_probabilities(n, whole):
    # Quantile-transform the independent angular Beta variable to uniform for
    # stable quadrature at its singular lower endpoint.
    def pass_given_angle(u):
        t = 1.0 if n == 1 else special.betaincinv(0.5, (n - 1) / 2, u)
        return special.gammainc(n / 2, Z**2 / (2 * t))
    first = integrate.quad(pass_given_angle, 0, 1, epsabs=1e-10)[0]
    both = integrate.quad(lambda u: pass_given_angle(u)**2, 0, 1, epsabs=1e-10)[0]
    neither = integrate.quad(lambda u: (1 - pass_given_angle(u))**2, 0, 1,
                             epsabs=1e-10)[0]
    assert first == pytest.approx(P, abs=1e-8)
    assert both == pytest.approx(whole, abs=5e-6)
    assert P**2 - 1e-12 <= both <= P
    assert neither == pytest.approx(1 - 2 * P + both, abs=1e-8)
    assert 0 < neither <= 1 - P


def test_limit_whole_score_set_and_whole_hull_probabilities_differ():
    # E2: score process A*sin(2theta)+B*cos(2theta); its supremum in absolute
    # value is sqrt(A²+B²). Both endpoints have absolute value |B|.
    radial = -math.expm1(-Z**2 / 2)
    axial = special.erf(Z / math.sqrt(2))
    assert radial == pytest.approx(0.8535, abs=5e-5)
    assert axial == pytest.approx(0.95)
    assert radial < axial


def _paired_readout(row, rank, scale):
    # The prior is uniform in noise angle (physical half-Cauchy scale c=s).
    fit = pooled_witness([row], alpha=1 / 8)

    def density(angle):
        return math.exp(-0.5 * fit.precision * math.sin(angle - fit.angle)**2)

    angles = _central_angles(density, 1.0 - 2.0 * rank)
    endpoints = scale * np.tan(angles)
    return math.sqrt(endpoints[0] * endpoints[1])


@pytest.mark.parametrize("rank", [0.1, 0.25, 0.5])
def test_paired_point_returns_zero_center_and_respects_unit_and_channel_changes(rank):
    scale = 7.0
    assert _paired_readout([0, 0], rank, scale) == pytest.approx(scale)
    row = math.sqrt(8) * np.array([math.sin(0.63), math.cos(0.63)])
    point = _paired_readout(row, rank, scale)
    assert _paired_readout(row, rank, scale * 10) == pytest.approx(10 * point)
    assert _paired_readout(row[::-1], rank, 1 / scale) == pytest.approx(1 / point)


def test_paired_quartiles_and_median_differ_on_one_reachable_raw_law():
    centre = 0.63
    row = math.sqrt(8) * np.array([math.sin(centre), math.cos(centre)])
    assert _paired_readout(row, 0.25, 1) == pytest.approx(0.7470473930)
    assert _paired_readout(row, 0.5, 1) == pytest.approx(0.7639974014)

    # Certify the ordering using coarse quantile brackets and the derivative
    # bound |f'| <= 4 on the unnormalized angular density.
    n = 16384
    axis = math.pi / 2

    def midpoint(x):
        angles = (np.arange(n) + 0.5) * (x / n)
        return x * np.mean(np.exp(-4 * np.sin(angles - centre)**2))

    total = midpoint(axis)
    for rank, mass, sign in [(0.25, 0.45, 1), (0.5, 0.760, -1), (0.75, 1.254, 1)]:
        x = math.atan(mass)
        difference = midpoint(x) - rank * total
        error_bound = (x**2 + rank * axis**2) / n
        assert sign * difference > error_bound
    assert 0.45 * 1.254 < 0.760**2


def test_paired_quartiles_and_median_agree_at_first_order_on_aligned_arrays():
    direction = np.array([math.sin(0.63), math.cos(0.63)])
    scaled_differences = []
    for radius in (8, 16, 32):
        row = radius * direction
        quartiles = _paired_readout(row, 0.25, 1)
        median = _paired_readout(row, 0.5, 1)
        scaled_differences.append(radius * abs(math.log(quartiles / median)))
    assert scaled_differences[2] < 0.6 * scaled_differences[1]
    assert scaled_differences[1] < 0.6 * scaled_differences[0]
