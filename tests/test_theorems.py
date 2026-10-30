"""Theorems of definition/ that had no test, checked against met/full.py where they say
something numeric about it. One section per theorem; each theorem file names its test.

Not here, with the reason in their files: T21 (a derivation about which readouts are
admissible; its representation step is O5) and T27 (own-world calibration needs a declared
world of signals, which ours29 does not have, L10).
"""
import math

import numpy as np
import pytest
from scipy import integrate, stats
from scipy.special import exp1, gammainc

from met.full import (Reading, isotropic, ours29, bartlett_eta, reading_state, _h_terms,
                      _log_law, _log_prior, _series_terms)
from tests.test_full import _series_full, _full_cov

SF, SA = 0.3, 0.1
S = SF / SA


def _phys_c(reading, m):
    """C_m = (m N; N), the candidate map in physical units."""
    n = reading.subspace
    return np.vstack([m * n, n])


def _H_full(reading, c, m):
    """H at mass m through met.full."""
    return float(_h_terms(reading.pre(), c, np.array([math.atan(m / c)]))[0])


# --- T2: candidates are codirectional pairs F* = m a* --------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_t2_explained_energy_is_the_best_codirectional_fit(seed):
    """H(m) is the whitened energy explained by the best pair (m N alpha, N alpha) over alpha:
    the candidates are exactly the codirectional pairs with ratio m."""
    rng = np.random.default_rng(seed)
    reading = _series_full(seed, 3, 1)[0]
    y = np.r_[reading.force, reading.acceleration]
    w = np.linalg.cholesky(np.linalg.inv(reading.cov)).T          # whitening, w^T w = Sigma^-1
    for m in rng.uniform(0.05, 20.0, size=5):
        cm = _phys_c(reading, m)
        alpha, *_ = np.linalg.lstsq(w @ cm, w @ y, rcond=None)
        residual = np.sum((w @ (y - cm @ alpha)) ** 2)
        explained = y @ np.linalg.solve(reading.cov, y) - residual
        for c in (0.5, 2.0, 7.0):                                   # any parametrisation
            assert _H_full(reading, c, m) == pytest.approx(explained, rel=1e-9, abs=1e-12)
        f_star, a_star = (cm @ alpha)[:3], (cm @ alpha)[3:]
        np.testing.assert_allclose(f_star, m * a_star, rtol=1e-12)


# --- T3: separation, and the readout in use meets it ---------------------------------------

def _ratio_of_means(x, y):
    return x.mean() / y.mean()


def _median_ratio(x, y):
    return np.median(x / y)


def test_t3_ratio_of_means_separates_and_composes():
    """R(X,Y) = R(X,1) R(1,Y) with lumping (L5) and common force (L14), for the ratio of means
    on any joint cloud; the median of X/Y fails lumping (the property has teeth)."""
    rng = np.random.default_rng(7)
    x1, x2 = rng.lognormal(0, 1, 50_000), rng.gamma(2.0, 1.5, 50_000)
    y1, y2 = rng.lognormal(0.3, 0.6, 50_000), rng.gamma(3.0, 0.5, 50_000) + 0.1 * x1
    one = np.ones_like(x1)
    R = _ratio_of_means
    assert R(x1 + x2, y1) == pytest.approx(R(x1, y1) + R(x2, y1), rel=1e-12)
    assert 1 / R(x1, y1 + y2) == pytest.approx(1 / R(x1, y1) + 1 / R(x1, y2), rel=1e-12)
    assert R(x1, y1) == pytest.approx(R(x1, one) * R(one, y1), rel=1e-12)
    lumped, parts = _median_ratio(x1 + x2, y1), _median_ratio(x1, y1) + _median_ratio(x2, y1)
    assert abs(lumped / parts - 1) > 0.01


# --- T4: the readout is E[X]/E[Y] = E[alpha m]/E[alpha] on the joint law ----------------------

def test_t4_readout_is_the_ratio_of_means_of_the_joint_cloud():
    """Draw m from the law and each reading's candidate alpha | m ~ N(G^-1 b, G^-1); the ratio
    of the mean stacked |F*| to the mean stacked |a*| is ours29's readout."""
    readings = _series_full(11, 3, 3, mass=2.0, signal=0.8)
    r = ours29(readings, 1.0)
    rng = np.random.default_rng(12)
    m_nodes = np.exp(r["lm"])
    idx = rng.choice(len(m_nodes), size=200_000, p=r["law"])
    nodes, counts = np.unique(idx, return_counts=True)
    xs, ys = [], []
    for node, count in zip(nodes, counts):
        m = m_nodes[node]
        stack = []
        for rd in readings:
            p = rd.pre()
            g = m * m * p["k_ff"] + m * p["k_s"] + p["k_aa"]
            mu = np.linalg.solve(g, m * p["u"] + p["v"])
            chol = np.linalg.cholesky(np.linalg.inv(g))
            stack.append(mu[None] + rng.normal(size=(count, len(mu))) @ chol.T)
        length = np.linalg.norm(np.hstack(stack), axis=1)
        xs.append(m * length)
        ys.append(length)
    x, y = np.concatenate(xs), np.concatenate(ys)
    est = x.mean() / y.mean()
    # delta-method standard error of the ratio
    se = est * math.sqrt(np.var(x / x.mean() - y / y.mean()) / len(x))
    assert abs(est - r["readout"]) < 4 * se
    assert se / est < 0.01


# --- T6, T7: the Gaussian family in noise units, its flat endpoint, prior once --------------

def _log_marginal(y, cov, cm):
    """log N(y; 0, cov) evaluated directly (dense linear algebra, no shortcut)."""
    sign, logdet = np.linalg.slogdet(2 * math.pi * cov)
    assert sign > 0
    return -0.5 * logdet - 0.5 * y @ np.linalg.solve(cov, y)


@pytest.mark.parametrize("omega", [0.3, 1.0, 5.0, 100.0])
def test_t6_t7_noise_unit_gaussian_push_gives_omega_tempered_h(omega):
    """alpha ~ N(0, omega^2 G_m^-1) (the push Gaussian in noise units): the marginal of y is
    omega^2/(1+omega^2) * H(m)/2 + an m-free constant, for every omega (T7: the determinant is
    m-free), and the flat weight that ours29 uses is its omega -> infinity endpoint (T6)."""
    reading = _series_full(21, 3, 1)[0]
    y = np.r_[reading.force, reading.acceleration]
    si = np.linalg.inv(reading.cov)
    ms = [0.1, 0.7, 2.0, 9.0]
    resid = []
    for m in ms:
        cm = _phys_c(reading, m)
        gm = cm.T @ si @ cm
        lm = _log_marginal(y, reading.cov + omega ** 2 * cm @ np.linalg.solve(gm, cm.T), cm)
        resid.append(lm - 0.5 * omega ** 2 / (1 + omega ** 2) * _H_full(reading, 1.0, m))
    assert np.ptp(resid) < 1e-9


def test_t7_a_push_weight_in_physical_units_puts_m_into_the_reading():
    """The same Gaussian in physical units, alpha ~ N(0, omega^2 I), leaves an m-dependent
    determinant: something attached to the reading depends on m, which T7 forbids."""
    reading = _series_full(21, 3, 1)[0]
    y = np.r_[reading.force, reading.acceleration]
    resid = []
    for m in [0.1, 0.7, 2.0, 9.0]:
        cm = _phys_c(reading, m)
        lm = _log_marginal(y, reading.cov + 25.0 * cm @ cm.T, cm)
        resid.append(lm - 0.5 * 25.0 / 26.0 * _H_full(reading, 1.0, m))
    assert np.ptp(resid) > 1.0


# --- T15: equal energies give s --------------------------------------------------------------

def test_t15_equal_energies_give_s():
    """P_g = Q_g (noise units, g-weighted) and c = s: median = readout = s. One 3D reading;
    1D readings with equal sums (eta < 1 here); 3D readings each with P = Q and one g < 1."""
    rng = np.random.default_rng(3)
    x, y = rng.normal(size=3), rng.normal(size=3)
    y *= np.linalg.norm(x) / np.linalg.norm(y)
    series = {"one 3D": [isotropic(SF * x, SA * y, SF, SA)]}
    xs, ys = 2 * rng.normal(size=5), 2 * rng.normal(size=5)
    ys *= math.sqrt((xs ** 2).sum() / (ys ** 2).sum())
    series["1D sums"] = [isotropic([SF * a], [SA * b], SF, SA) for a, b in zip(xs, ys)]
    each = []
    for _ in range(4):
        x, y = rng.normal(size=3), rng.normal(size=3)
        each.append(isotropic(SF * x, SA * y * np.linalg.norm(x) / np.linalg.norm(y), SF, SA))
    series["3D each"] = each
    results = {k: ours29(v, S) for k, v in series.items()}
    assert results["1D sums"]["eta"] < 0.9
    assert min(results["3D each"]["g"]) < 0.9
    for r in results.values():
        assert r["median"] == pytest.approx(S, rel=1e-9)
        assert r["readout"] == pytest.approx(S, rel=1e-9)


def test_t15_the_sums_are_g_weighted():
    """Equal unweighted sums are not enough when a reading's misfit weight is below 1."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        x, y = rng.normal(size=(4, 3)), rng.normal(size=(4, 3))
        y *= np.linalg.norm(x) / np.linalg.norm(y)
        readings = [isotropic(SF * x[i], SA * y[i], SF, SA) for i in range(4)]
        g = [min(1.0, 2.0 / reading_state(rd, S)["Q"]) for rd in readings]
        p, q = (x ** 2).sum(1), (y ** 2).sum(1)
        if min(g) < 0.7 and abs(np.dot(g, q - p)) > 0.5:
            break
    r = ours29(readings, S)
    assert abs(r["median"] / S - 1) > 1e-3


# --- T16: a tilt is an unseen reading ------------------------------------------------------

@pytest.mark.parametrize("tilt_centre, lam, law_centre", [(1.7, 0.8, 1.7), (1.7, 0.8, 0.4),
                                                        (0.3, -1.5, 2.0)])
def test_t16_a_tilt_is_one_reading_with_equal_energies(tilt_centre, lam, law_centre):
    """exp(lam sin 2 theta_c) times the law (raw, eta = 1) is the law with one more 1D reading of
    noise ratio c, P = Q, D = 2 lam; whatever the starting law's centre."""
    readings = _series_full(5, 3, 3)
    x = math.sqrt(2 * abs(lam))
    extra = isotropic([tilt_centre * 0.5 * x], [math.copysign(0.5 * x, lam)],
                      tilt_centre * 0.5, 0.5)
    g0 = ours29(readings, law_centre, calibrate=False, readout=False)["g"]
    g1 = ours29(readings + [extra], law_centre, calibrate=False, readout=False)["g"]
    assert g1 == g0 + [1.0]
    lm = np.linspace(-6, 6, 241)
    diff = _log_law(readings + [extra], law_centre, g1, 1.0, lm) - _log_law(readings, law_centre, g0, 1.0, lm)
    theta = np.arctan(np.exp(lm) / tilt_centre)
    assert np.ptp(diff - lam * np.sin(2 * theta)) < 1e-12


# --- T17: convergence ------------------------------------------------------------------------

def test_t17_readout_converges_to_the_mass_as_the_signal_grows():
    """Fixed noise draw, full covariances, signal scaled up: the readout tends to m."""
    rng = np.random.default_rng(31)
    covs = [_full_cov(rng, 3) for _ in range(4)]
    pushes = [rng.normal(size=3) for _ in range(4)]
    noise = [np.linalg.cholesky(cv) @ rng.normal(size=6) for cv in covs]
    m = 2.5
    errors = []
    for scale in [1.0, 10.0, 100.0, 1000.0]:
        rd = []
        for cv, push, e in zip(covs, pushes, noise):
            y = scale * np.r_[m * push, push] + e
            rd.append(Reading(y[:3], y[3:], cv))
        errors.append(abs(ours29(rd, 1.0)["readout"] - m))
    assert all(b < a for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 2e-3


# --- T18: the flat law's width error, and eta as its measured inverse -----------------------

@pytest.mark.parametrize("r", [1.0, 2.0])
def test_t18_width_ratio_is_one_plus_n_over_sum_r2(r):
    """One noise pair, N readings of signal r (noise units): the true variance of the best line
    over the flat law's variance is 1 + 1/r^2 (T18 with one pair). Closed form of equation.md."""
    rng = np.random.default_rng(int(10 * r))
    n, reps, theta0 = 400, 4000, 0.6
    x = r * math.sin(theta0) + rng.normal(size=(reps, n))
    y = r * math.cos(theta0) + rng.normal(size=(reps, n))
    total = ((y ** 2 - x ** 2) + 2j * x * y).sum(axis=1)
    ratio = np.var(np.angle(total) / 2) * np.mean(np.abs(total))
    assert ratio == pytest.approx(1 + 1 / r ** 2, rel=0.08)


def _closed_form_eta(x, y):
    """eta for 1D readings of one noise pair in noise units (equation.md, closed form)."""
    c = (y ** 2 - x ** 2) + 2j * x * y
    total = c.sum()
    d1 = -0.5 * np.imag(np.conj(c) * np.exp(1j * np.angle(total)))
    return abs(total) / np.sum(d1 ** 2)


@pytest.mark.parametrize("r", [1.0, 2.0])
def test_t18_eta_is_the_inverse_on_long_series(r):
    """eta is the measured inverse of T18's ratio: averaged over long series it is
    r^2 / (1 + r^2); met.full's eta is the closed form on each series. One series of 400 is
    not enough on its own: J's heavy tail gives eta a spread of about 0.09 there (O14)."""
    theta0, n = 0.6, 400
    etas = []
    for seed in range(200):
        rng = np.random.default_rng(seed)
        x = r * math.sin(theta0) + rng.normal(size=n)
        y = r * math.cos(theta0) + rng.normal(size=n)
        etas.append(_closed_form_eta(x, y))
        if seed == 0:
            readings = [isotropic([SF * a], [SA * b], SF, SA) for a, b in zip(x, y)]
            assert bartlett_eta(readings, S, [1.0] * n)[0] == pytest.approx(etas[0], rel=1e-9)
    assert np.mean(etas) == pytest.approx(r * r / (1 + r * r), abs=0.04)


# --- T19: no mean ----------------------------------------------------------------------------

def test_t19_the_law_has_no_mean():
    """The law's density in log m falls as 1/m at large m (density in m as 1/m^2), so the
    truncated mean grows by the same amount per unit of log m without end."""
    readings = _series_full(41, 3, 3)
    r = ours29(readings, 1.0, grid=(60.0, 6001), readout=False)
    lm, lp = r["lm"], np.log(r["law"])
    top = lm > lm[-1] - 10
    slope = np.polyfit(lm[top], lp[top], 1)[0]
    assert slope == pytest.approx(-1.0, abs=1e-3)
    m_p = np.exp(lm) * r["law"]
    def truncated(upper):
        return m_p[lm <= upper].sum()
    hi = lm[-1]
    step1 = truncated(hi) - truncated(hi - 10)
    step2 = truncated(hi - 10) - truncated(hi - 20)
    assert step1 > 0 and step1 == pytest.approx(step2, rel=1e-3)


# --- T20: pi is declared; c enters only through it ------------------------------------------

def test_t20_the_centre_enters_only_through_the_starting_law():
    """ours29 has no default centre, and changing c changes the law by the ratio of the two
    starting laws alone: g and eta do not depend on c (T20, the declaration)."""
    with pytest.raises(TypeError):
        ours29(_series_full(1, 3, 2))                           # c is required
    readings = _series_full(62, 3, 4, signal=0.7)
    ra, rb = ours29(readings, 0.6, readout=False), ours29(readings, 4.0, readout=False)
    np.testing.assert_allclose(ra["g"], rb["g"], rtol=1e-9)
    assert ra["eta"] == pytest.approx(rb["eta"], rel=1e-7)
    assert ra["eta"] < 1.0
    lm = np.linspace(-5, 5, 201)
    la = _log_law(readings, 0.6, ra["g"], ra["eta"], lm) - _log_prior(lm, 0.6)
    lb = _log_law(readings, 4.0, rb["g"], rb["eta"], lm) - _log_prior(lm, 4.0)
    assert np.ptp(la - lb) < 1e-8


# --- T23: the scope of the noise-model consistency and zero-reading routes -----------------

def _zero_reading_factor(pdf, phi):
    """Evidence of an exact zero 1D reading at the line phi (c = s), flat weight in noise units:
    int f(u1 b) f(u2 b) db with (u1, u2) the unit candidate direction in noise units."""
    u1, u2 = math.sin(phi), math.cos(phi)
    return integrate.quad(lambda b: pdf(u1 * b) * pdf(u2 * b), -np.inf, np.inf,
                          epsabs=1e-14, epsrel=1e-12)[0]


def test_t23_only_gaussian_noise_leaves_a_zero_reading_alone():
    """With independent errors on the two channels, an exact zero reading's factor is m-free for
    Gaussian noise and tilts the law for this Student t5 product (max/min 1.077).
    These examples do not characterize every independent-channel density; B05's
    non-Gaussian counterexample below has the same constant-factor property."""
    phis = np.linspace(0.0, math.pi / 2, 9)
    gauss = np.array([_zero_reading_factor(stats.norm.pdf, p) for p in phis])
    t5 = np.array([_zero_reading_factor(lambda z: stats.t.pdf(z, 5), p) for p in phis])
    assert np.ptp(gauss) / gauss.mean() < 1e-10
    assert t5.max() / t5.min() == pytest.approx(1.0775, abs=5e-4)
    # A t law elliptical in the whitened plane is also m-free. Maxwell's
    # characterization requires joint spherical symmetry as well as independence;
    # the zero-line integral alone is a weaker condition.
    def elliptic_factor(phi):
        u = np.array([math.sin(phi), math.cos(phi)])
        dens = lambda b: (1 + np.sum((u * b) ** 2) / 5) ** -3.5
        return integrate.quad(dens, -np.inf, np.inf, epsabs=1e-14, epsrel=1e-12)[0]
    ell = np.array([elliptic_factor(p) for p in phis])
    assert np.ptp(ell) / ell.mean() < 1e-10
    # and met.full's zero reading is m-free (H = 0 everywhere)
    zero = isotropic([0.0], [0.0], SF, SA)
    assert np.all(_h_terms(zero.pre(), S, phis) == 0.0)


def _t23_ours30_weight(f, a):
    """Exact 1D positive-mass geometric weight, with unit declared covariance."""
    q_plus = 0.0 if f * a >= 0 else min(f * f, a * a)
    return 1.0 if q_plus <= 0.5 else 0.5 / q_plus


def test_t23_gaussian_ours30_fails_literal_total_signal_consistency():
    """An analytic population-score bound, not a large-N simulation.

    The world and the model are both Gaussian with covariance I. A fixed
    weak signal with true mass 2 has a nonzero limiting score at that mass.
    The positive-mass misfit weights cause this; a positive Bartlett power
    cannot move the population optimum back to the true mass. O4/T23 give
    the law-concentration argument as N and hence total true signal grow.
    """
    kappa = 0.5
    threshold = math.sqrt(kappa)
    b_integral = 2 / math.pi * integrate.quad(
        lambda v: (v - kappa / v) * math.exp(-v * v),
        threshold, np.inf, epsabs=1e-13, epsrel=1e-12)[0]
    b = (math.exp(-kappa) - kappa * exp1(kappa)) / math.pi
    assert b_integral == pytest.approx(b, abs=2e-14)
    assert b == pytest.approx(0.10397397063916228, abs=2e-14)

    theta0 = math.atan(2.0)
    r = 0.01
    # Under the centered Gaussian, E[g X Y]=b. For a Gaussian mean shift
    # of length r, Cauchy-Schwarz bounds expectation changes by
    # sqrt(E_0 f²) sqrt(exp(r²)-1). The directional score obeys E_0 f²<=2.
    null_score = b * math.cos(2 * theta0)
    score_upper_bound = null_score + math.sqrt(2 * math.expm1(r * r))
    cross_lower_bound = b - math.sqrt(math.expm1(r * r))
    assert null_score == pytest.approx(-3 * b / 5)
    assert score_upper_bound < -0.048
    assert cross_lower_bound > 0.09
    # C/N>=2 E[gXY]>0; J/N is positive and bounded by E||Y||^4/4.
    fourth_moment_bound = (r ** 4 + 8 * r * r + 8) / 4
    eta_lower_bound = min(1.0, 2 * cross_lower_bound / fourth_moment_bound)
    assert eta_lower_bound > 0.09
    # The bound uses exactly the production g, including its positive-domain cut.
    from met.full import ours30
    for f, a in ((1, 1), (1, -1), (0.2, -0.3), (3, -2)):
        measured = ours30([isotropic([f], [a], 1, 1)], 1, readout=False)
        assert measured["g"][0] == pytest.approx(_t23_ours30_weight(f, a), abs=1e-14)


def _t23_independent_non_gaussian_pdf(x):
    """B05's standardized b=5 Mellin counterexample, with variance exactly one."""
    scale = math.sqrt(14 / 5)
    z = scale * float(x)
    if abs(z) < 1e-3:
        h = stats.norm.pdf(0) * (0.8 - z * z / 7 + z ** 6 / 264)
    else:
        h = 13.5 * gammainc(2.5, z * z / 2) / abs(z) ** 5 - stats.norm.pdf(z)
    return scale * 0.8 * h


def test_t23_independent_non_gaussian_errors_can_leave_zero_readings_alone():
    """Independent channels + a constant zero-line factor do not force Gaussianity."""
    pdf = _t23_independent_non_gaussian_pdf
    assert 2 * integrate.quad(pdf, 0, np.inf, epsabs=1e-12)[0] == pytest.approx(1, abs=2e-11)
    assert 2 * integrate.quad(lambda x: x * x * pdf(x), 0, np.inf,
                              epsabs=1e-12)[0] == pytest.approx(1, abs=2e-11)
    assert pdf(0) > 1.05 * stats.norm.pdf(0)
    expected = math.sqrt(14 / 5) * 0.64 / math.sqrt(2 * math.pi)
    for phi in (0, 0.2, 0.7, math.pi / 4, 1.3, math.pi / 2):
        assert _zero_reading_factor(pdf, phi) == pytest.approx(expected, abs=3e-10)


def test_t23_elliptical_non_gaussian_scale_mixture_has_correct_covariance_and_zero_factor():
    """A bounded Gaussian-scale mixture is spherical, non-Gaussian, and Z-neutral."""
    variances = np.array([0.5, 1.5])
    probabilities = np.array([0.5, 0.5])
    assert probabilities @ variances == 1  # Each channel's declared variance.
    assert 3 * probabilities @ variances ** 2 == 3.75  # Not the Gaussian fourth moment 3.
    expected = sum(p / math.sqrt(2 * math.pi * v) for p, v in zip(probabilities, variances))
    for phi in (0, 0.4, 1.1, math.pi / 2):
        direction = np.array([math.sin(phi), math.cos(phi)])
        def density_on_line(t):
            radius2 = float(np.dot(t * direction, t * direction))
            return sum(p * math.exp(-radius2 / (2 * v)) / (2 * math.pi * v)
                       for p, v in zip(probabilities, variances))
        actual = integrate.quad(density_on_line, -np.inf, np.inf, epsabs=1e-12)[0]
        assert actual == pytest.approx(expected, abs=2e-12)


@pytest.mark.parametrize("mixture", [False, True], ids=["Gaussian", "bounded-elliptical-mixture"])
def test_t23_fixed_n_strong_signal_consistency_does_not_select_gaussian(mixture):
    """Frozen finite noise, with active g and eta, for two declared-covariance models.

    This is the fixed-N amplitude limit, distinct from L11's unrestricted
    total-signal growth. Each alternative factor is raised to the unchanged
    production g and eta; its latent integral is the residual-space mixture.
    """
    from met.full import ours30
    errors = np.array([[0.0, 8.0], [40.0, -6.0], [-30.0, 4.0]])
    deviations = []
    for strength in (100.0, 1000.0, 10000.0):
        readings = [isotropic(2 * np.array([strength, 0.0]) + error,
                              [strength, 0.0], 1, 1) for error in errors]
        oracle = ours30(readings, 1.0, readout=False)
        assert max(oracle["g"]) < 1
        assert 0 < oracle["eta"] < 0.1
        midpoint = (oracle["lm"][0] + oracle["lm"][-1]) / 2
        halfspan = oracle["lm"][-1] - oracle["lm"][0]
        lm = np.linspace(midpoint - halfspan, midpoint + halfspan, 6001)
        mass = np.exp(lm)
        log_law = -np.logaddexp(lm, -lm)  # Half-Cauchy prior in log mass.
        for reading, g in zip(readings, oracle["g"]):
            residual = reading.force[None] - mass[:, None] * reading.acceleration[None]
            k = np.sum(residual * residual, axis=1) / (1 + mass * mass)
            if mixture:
                # k_i=2: integrating the 2D flat latent vector leaves v^-1.
                components = [-math.log(2 * v) - k / (2 * v) for v in (0.5, 1.5)]
                log_factor = np.logaddexp(*components)
            else:
                log_factor = -k / 2
            log_law += oracle["eta"] * g * log_factor
        density = np.exp(log_law - log_law.max())
        cdf = integrate.cumulative_trapezoid(density, lm, initial=0)
        median = float(np.interp(cdf[-1] / 2, cdf, mass))
        if not mixture:
            assert median == pytest.approx(oracle["median"], rel=2e-6)
        deviations.append(abs(median - 2))
    assert deviations[2] < deviations[1] < deviations[0]
    assert deviations[-1] < 0.001


# --- T26: per-reading rules cannot fix the width; pure noise adds score, not curvature ------

def test_t26_a_pure_noise_reading_adds_score_noise_and_no_curvature():
    """Pure-noise 1D readings at a fixed line: mean l'' ~ 0 while mean l'^2 ~ 1, so any
    positive per-reading weight makes the law too narrow; eta, from the series, falls instead."""
    rng = np.random.default_rng(61)
    n = 4000
    x, y = rng.normal(size=n), rng.normal(size=n)
    readings = [isotropic([SF * a], [SA * b], SF, SA) for a, b in zip(x, y)]
    phi = np.array([0.4])
    _, per = _series_terms(readings, S, [1.0] * n, phi, 2)
    d1 = np.array([p[1][0] for p in per])
    d2 = np.array([p[2][0] for p in per])
    assert abs(d2.mean()) < 4 * d2.std() / math.sqrt(n)
    assert np.mean(d1 ** 2) == pytest.approx(1.0, rel=0.06)
    short = readings[:400]
    r = ours29(short, S)
    assert r["eta"] < 0.2
    def theta_width(interval):
        return math.atan(interval[1] / S) - math.atan(interval[0] / S)
    # quoted: about the starting law's own 95% (0.95 pi / 2); raw: five times narrower
    assert theta_width(r["interval"]) == pytest.approx(0.95 * math.pi / 2, abs=0.05)
    assert theta_width(r["raw_interval"]) < 0.25 * theta_width(r["interval"])


# --- T28: mass and push size are orthogonal; the reference law is the half-Cauchy -----------

def test_t28_fisher_orthogonal_and_uniform_in_theta():
    """One noise pair, 1D: mean of (F, a) = (m alpha, alpha). In (theta, r), theta = atan(m/s),
    r the signal length in noise units, the Fisher matrix is diag(r^2, 1): orthogonal, and
    Jeffreys for theta at fixed r is uniform for every r. Uniform theta on (0, pi/2) is the
    half-Cauchy at s, which is met.full's starting law at c = s."""
    def mean(theta, r):
        return np.array([SF * r * math.sin(theta), SA * r * math.cos(theta)])
    w = np.diag([1 / SF, 1 / SA])
    for theta, r in [(0.3, 0.5), (0.9, 2.0), (1.4, 7.0)]:
        h = 1e-6
        jt = (mean(theta + h, r) - mean(theta - h, r)) / (2 * h)
        jr = (mean(theta, r + h) - mean(theta, r - h)) / (2 * h)
        j = np.column_stack([jt, jr])
        fisher = j.T @ w @ w @ j
        np.testing.assert_allclose(fisher, np.diag([r * r, 1.0]), atol=1e-6)
    lm = np.linspace(-8, 8, 161)
    m = np.exp(lm)
    uniform_theta_in_log_m = (2 / math.pi) * (S / (S * S + m * m)) * m   # d theta / d log m
    np.testing.assert_allclose(np.exp(_log_prior(lm, S)) * 2 / math.pi, uniform_theta_in_log_m,
                               rtol=1e-12)


# --- T34: the gain is the geometric mean of the pairings (27 Oct) -----------------------

def _pairings(force, acc):
    force, acc = np.asarray(force, float), np.asarray(acc, float)
    dot = force @ acc
    return dot / (acc @ acc), (force @ force) / dot, math.sqrt(force @ force) / math.sqrt(acc @ acc)


T34_PAIRS = [([2.0], [1.0]), ([3.0, 0.0], [1.5, 0.0]), ([2.0, 0.5], [1.0, 0.1]),
             ([1.0, 1.0, 0.2], [0.3, 0.9, 0.0]), ([5.0, -1.0, 2.0], [1.0, 0.2, 0.3])]


@pytest.mark.parametrize("force, acc", T34_PAIRS)
def test_t34_gain_is_the_geometric_mean_of_the_pairings(force, acc):
    from met.full import exact_mass
    m1, m2, gain = _pairings(force, acc)
    assert m1 * m2 == pytest.approx(gain ** 2, rel=1e-13)
    assert m1 <= gain * (1 + 1e-14) and gain <= m2 * (1 + 1e-14)
    cos = (np.asarray(force) @ np.asarray(acc)) / (np.linalg.norm(force) * np.linalg.norm(acc))
    assert m1 / gain == pytest.approx(cos, rel=1e-13) and gain / m2 == pytest.approx(cos, rel=1e-13)
    if exact_mass(force, acc)[0] == "one":           # codirectional: T1's one mass is all three
        assert exact_mass(force, acc)[1] == pytest.approx(gain, rel=1e-12)
        assert m1 == pytest.approx(m2, rel=1e-12)


def test_t34_swap():
    rng = np.random.default_rng(34)
    for _ in range(200):
        force, acc = rng.normal(size=3), rng.normal(size=3)
        if force @ acc <= 0:
            acc = -acc
        m1, m2, gain = _pairings(force, acc)
        s1, s2, sgain = _pairings(acc, force)
        assert s1 == pytest.approx(1 / m2, rel=1e-12)
        assert s2 == pytest.approx(1 / m1, rel=1e-12)
        assert sgain == pytest.approx(1 / gain, rel=1e-12)
