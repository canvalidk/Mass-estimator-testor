"""The cartesian reference (24 September): flat in the pair as a vector.

Each closed form is checked against an integral done a different way:
the radial kernel and mean radius by direct radial quadrature; the whole
law against brute force in the 21 September coordinates (f, alpha) with
weight r^(d-2); and in 1D against brute force in the record's own variables
(true acceleration a*, mass m). The record's section 3 formulas are checked
term by term.
"""

import math
import zlib

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import i0e, roots_legendre

from met.analyst import Estimator, joint_log_density, posterior_summaries, tail_slopes
from met.law import ReadingSet, cartesian_mean_radius, mean_radius, radial_terms

CART = {"reference": "cartesian", "nuisance": "radius"}
FLAT = {"reference": "flat", "nuisance": "radius"}
SINGLE = {"rule": "single", "prior": []}
READOUTS = ["ratio_of_means", "median", "log_sd", "interval_95"]


def _series(force, acceleration, sf, sa):
    force, acceleration = np.atleast_2d(np.asarray(force, float)), np.atleast_2d(np.asarray(acceleration, float))
    n = force.shape[0]
    return ReadingSet(force[None], acceleration[None], np.broadcast_to(sf, (1, n)), np.broadcast_to(sa, (1, n)))


def _angular_scaled(z, d):
    """< exp(z u) >_u * exp(-z), stable for large z."""
    if d == 1:
        return 0.5 * (1 + math.exp(-2 * z))
    if d == 2:
        return i0e(z)
    return (1 - math.exp(-2 * z)) / (2 * z) if z > 1e-9 else 1.0


# ---------------------------------------------------------------- radial closed forms

@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("h", [0.0, 0.3, 1.0, 3.0, 8.0, 20.0])
def test_kernel_and_mean_radius_match_direct_integration(d, h):
    """K = integral r^(d-1) exp(-r^2/2) <exp(r h.u)> dr is exp(h^2/2) times a constant;
    E[r | theta] is the noncentral chi mean."""
    def integrand(r, k):
        return r**k * math.exp(-(r - h) ** 2 / 2) * _angular_scaled(r * h, d)

    k0 = quad(integrand, 0, np.inf, args=(d - 1,), limit=400, epsabs=0, epsrel=1e-13)[0]
    k1 = quad(integrand, 0, np.inf, args=(d,), limit=400, epsabs=0, epsrel=1e-13)[0]
    constant = {1: math.sqrt(math.pi / 2), 2: 1.0, 3: math.sqrt(math.pi / 2)}[d]
    log_k, radius = radial_terms(h, d, "cartesian")
    assert math.log(k0 / constant) + h * h / 2 == pytest.approx(float(log_k), abs=1e-11)
    assert k1 / k0 == pytest.approx(float(radius), rel=1e-11)


def test_two_dimensions_is_the_21_september_law():
    """In 2D the two references are the same measure: every curve, every nuisance."""
    rng = np.random.default_rng(20)
    readings = ReadingSet(rng.normal(size=(3, 4, 2)) * 2, rng.normal(size=(3, 4, 2)),
                          rng.uniform(0.5, 2, size=(3, 4)), rng.uniform(0.5, 2, size=(3, 4)))
    u = np.tile(np.linspace(-6, 6, 121), (3, 1))
    for nuisance in ("radius", "acceleration", "force"):
        flat = readings.reading_terms(u, nuisance, reference="flat")
        cart = readings.reading_terms(u, nuisance, reference="cartesian")
        for a, b in zip(flat, cart):
            assert np.allclose(a, b, rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("d", [1, 3])
def test_mean_radius_is_continuous_at_zero(d):
    hs = np.array([0.0, 1e-9, 1e-7, 1e-6, 1.0001e-6, 1e-5])
    at_zero = {1: math.sqrt(2 / math.pi), 3: 2 * math.sqrt(2 / math.pi)}[d]
    assert np.allclose(cartesian_mean_radius(hs, d), at_zero, rtol=1e-9)


# ---------------------------------------------------------------- the record's section 3

@pytest.mark.parametrize("d", [1, 2, 3])
def test_single_reading_law_is_the_records_closed_form(d):
    """p(m) ~ 1/(sF^2 + m^2 sa^2) exp[-|F - m a|^2 / (2 (sF^2 + m^2 sa^2))], and
    E[alpha | m] = sigma_e chi_d(|a_hat| / sigma_e) with the record's a_hat, sigma_e."""
    rng = np.random.default_rng(30 + d)
    for _ in range(3):
        f, a = rng.normal(size=d) * 2, rng.normal(size=d)
        sf, sa = rng.uniform(0.4, 2.0, size=2)
        u = np.linspace(-5, 5, 201)[None]
        m = np.exp(u[0])
        lp, acc = joint_log_density(_series(f, a, sf, sa), u, CART, SINGLE)
        var = sf**2 + m**2 * sa**2
        resid = np.sum((f[None] - m[:, None] * a[None]) ** 2, axis=1)
        expected = -np.log(var) - resid / (2 * var) + u[0]          # density in u = log m
        diff = lp[0] - expected
        assert np.ptp(diff) < 1e-11
        a_hat = (m[:, None] * sa**2 * f[None] + sf**2 * a[None]) / var[:, None]
        sigma_e = sf * sa / np.sqrt(var)
        chi = cartesian_mean_radius(np.linalg.norm(a_hat, axis=1) / sigma_e, d)
        assert np.allclose(acc[0], sigma_e * chi, rtol=1e-12)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_symmetric_rule_adds_residuals(d):
    """Many readings under the symmetric rule: exp[-sum |F_i - m a_i|^2 / (2 (sF^2 + m^2 sa^2))]
    times the half-Cauchy prior, in every dimension."""
    rng = np.random.default_rng(40 + d)
    n, sf, sa = 6, 1.3, 0.6
    f, a = rng.normal(size=(n, d)) * 2 + 1, rng.normal(size=(n, d)) + 0.5
    combine = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
    u = np.linspace(-5, 5, 201)[None]
    m = np.exp(u[0])
    lp, _ = joint_log_density(_series(f, a, sf, sa), u, CART, combine)
    var = sf**2 + m**2 * sa**2
    resid = np.sum((f[None] - m[:, None, None] * a[None]) ** 2, axis=(1, 2))
    expected = -np.log(var) - resid / (2 * var) + u[0]
    assert np.ptp(lp[0] - expected) < 1e-10


def test_a_d_dimensional_reading_is_d_one_dimensional_readings():
    """Record section 3: under the symmetric rule a 3D reading is exactly three 1D readings."""
    rng = np.random.default_rng(50)
    f, a = rng.normal(size=(4, 3)) * 2 + 1, rng.normal(size=(4, 3)) + 0.5
    combine = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": 1.0}}]}
    u = np.linspace(-5, 5, 201)[None]
    lp3, _ = joint_log_density(_series(f, a, 1.0, 1.0), u, CART, combine)
    lp1, _ = joint_log_density(_series(f.reshape(-1, 1), a.reshape(-1, 1), 1.0, 1.0), u, CART, combine)
    assert np.ptp(lp3 - lp1) < 1e-10


# ---------------------------------------------------------------- brute force

def _brute_force_flat_coordinates(f_obs, a_obs, sf, sa, weight_power, nodes=400, radius=60.0):
    """E[f]/E[alpha] and the law of theta, integrating the 21 September measure
    df d(alpha) dOmega times r^weight_power numerically in (r, theta), with the
    direction averaged in closed form. No radial closed form is used."""
    d = len(f_obs)
    x, w = roots_legendre(nodes)
    theta = (x + 1) * math.pi / 4
    wt = w * math.pi / 4
    r = (x + 1) * radius / 2
    wr = w * radius / 2
    R, T = np.meshgrid(r, theta, indexing="ij")
    f, alpha = sf * R * np.sin(T), sa * R * np.cos(T)
    vec = f[..., None] * np.asarray(f_obs) / sf**2 + alpha[..., None] * np.asarray(a_obs) / sa**2
    z = np.linalg.norm(vec, axis=-1)
    if d == 1:
        log_avg = z + np.log(0.5 * (1 + np.exp(-2 * z)))
    elif d == 2:
        log_avg = z + np.log(i0e(z))
    else:
        safe = np.where(z > 1e-12, z, 1.0)
        log_avg = np.where(z > 1e-12, z + np.log(-np.expm1(-2 * safe) / (2 * safe)), 0.0)
    log_int = -0.5 * (R**2) + log_avg + np.log(np.maximum(R, 1e-300)) * (1 + weight_power)
    log_int = log_int - log_int.max()
    dens = np.exp(log_int) * wr[:, None] * wt[None, :]
    return float(np.sum(dens * f) / np.sum(dens * alpha)), theta, np.sum(dens, axis=0) / wt


@pytest.mark.parametrize("d", [1, 3])
def test_cartesian_is_the_flat_law_reweighted_by_r_to_the_d_minus_2(d):
    """The isolation identity: m_hat(cartesian) = E_flat[f r^(d-2)] / E_flat[alpha r^(d-2)]."""
    rng = np.random.default_rng(zlib.crc32(f"iso{d}".encode()))
    cases = [(rng.normal(size=d) * k, rng.normal(size=d) * k, *rng.uniform(0.5, 2.0, size=2))
             for k in (0.4, 1.0, 3.0)]
    for f, a, sf, sa in cases:
        brute, _, _ = _brute_force_flat_coordinates(f, a, sf, sa, weight_power=d - 2)
        mine = posterior_summaries(_series(f, a, sf, sa), CART, SINGLE, ["ratio_of_means"])
        assert mine["ratio_of_means"][0] == pytest.approx(brute, rel=1e-9)
        # and the same brute force with no weight is the 21 September estimator
        brute_flat, _, _ = _brute_force_flat_coordinates(f, a, sf, sa, weight_power=0)
        flat = posterior_summaries(_series(f, a, sf, sa), FLAT, SINGLE, ["ratio_of_means"])
        assert flat["ratio_of_means"][0] == pytest.approx(brute_flat, rel=1e-9)


def test_one_dimension_against_the_records_variables():
    """1D brute force in (a*, m) with measure da* dm (1 + m^2/s^2)^(-1/2): flat in the
    true acceleration, completed to uniform theta. Inner integral over a* adaptive
    (split at the kink of |a*|), outer over theta = atan(m/s) by Gauss-Legendre."""
    f_obs, a_obs, sf, sa = 1.7, 0.4, 1.2, 0.5
    s = sf / sa
    x, w = roots_legendre(200)
    num = den = 0.0
    for th, wt in zip((x + 1) * math.pi / 4, w * math.pi / 4):
        m = s * math.tan(th)
        dm = s / math.cos(th) ** 2 * wt

        def g(a_star, m=m):
            return math.exp(-(f_obs - m * a_star) ** 2 / (2 * sf**2) - (a_obs - a_star) ** 2 / (2 * sa**2)
                            - 0.5 * math.log1p(m * m / s**2)) * abs(a_star)

        # where to integrate: +-40 conditional SDs around the conditional peak (placement only)
        var = sf**2 + m * m * sa**2
        centre, width = (m * sa**2 * f_obs + sf**2 * a_obs) / var, 40 * sf * sa / math.sqrt(var)
        edges = sorted({centre - width, centre + width} | ({0.0} if abs(centre) < width else set()))
        inner = sum(quad(g, lo, hi, epsabs=0, epsrel=1e-13, limit=200)[0] for lo, hi in zip(edges, edges[1:]))
        num += m * inner * dm
        den += inner * dm
    mine = posterior_summaries(_series([f_obs], [a_obs], sf, sa), CART, SINGLE, ["ratio_of_means"])
    assert mine["ratio_of_means"][0] == pytest.approx(num / den, rel=1e-11)


def test_plain_flat_in_mass_completion_by_brute_force():
    """3D: d^3 a* dm (flat in the vector AND flat in m) is likelihood_product with
    nuisance acceleration and a flat_mass prior. Brute force in the flat coordinates
    with weight alpha^(d-2) = alpha."""
    f, a, sf, sa = np.array([2.0, 0.5, -0.3]), np.array([0.8, 0.1, 0.2]), 1.0, 0.7
    law = {"reference": "cartesian", "nuisance": "acceleration"}
    combine = {"rule": "likelihood_product", "prior": ["flat_mass"]}
    mine = posterior_summaries(_series(f, a, sf, sa), law, combine, ["ratio_of_means"])
    x, w = roots_legendre(400)
    theta, wt = (x + 1) * math.pi / 4, w * math.pi / 4
    r, wr = (x + 1) * 30.0, w * 30.0
    R, T = np.meshgrid(r, theta, indexing="ij")
    fm, al = sf * R * np.sin(T), sa * R * np.cos(T)
    z = np.linalg.norm(fm[..., None] * f / sf**2 + al[..., None] * a / sa**2, axis=-1)
    log_int = -0.5 * R**2 + z + np.log(-np.expm1(-2 * z) / (2 * z)) + np.log(R) + np.log(al)
    dens = np.exp(log_int - log_int.max()) * wr[:, None] * wt[None, :]
    assert mine["ratio_of_means"][0] == pytest.approx(float(np.sum(dens * fm) / np.sum(dens * al)), rel=1e-9)


# ---------------------------------------------------------------- symmetries and nulls

@pytest.mark.parametrize("d", [1, 2, 3])
def test_null_law_is_the_same_half_cauchy(d):
    sf, sa = 2.0, 0.5
    out = posterior_summaries(_series(np.zeros(d), np.zeros(d), sf, sa), CART, SINGLE, READOUTS)
    assert out["ratio_of_means"][0] == pytest.approx(sf / sa, rel=1e-12)
    assert out["log_sd"][0] == pytest.approx(math.pi / 2, rel=1e-12)


@pytest.mark.parametrize("d", [1, 3])
def test_reciprocal_symmetry_and_units(d):
    rng = np.random.default_rng(60 + d)
    for _ in range(4):
        f, a = rng.normal(size=d) * 3, rng.normal(size=d)
        sf, sa = rng.uniform(0.3, 2, size=2)
        r = posterior_summaries(_series(f, a, sf, sa), CART, SINGLE, READOUTS)
        w = posterior_summaries(_series(a, f, sa, sf), CART, SINGLE, READOUTS)
        assert r["ratio_of_means"][0] * w["ratio_of_means"][0] == pytest.approx(1, rel=1e-9)
        assert r["median"][0] * w["median"][0] == pytest.approx(1, rel=1e-9)
        k = 7.3
        scaled = posterior_summaries(_series(f * k, a, sf * k, sa), CART, SINGLE, READOUTS)
        assert scaled["ratio_of_means"][0] == pytest.approx(k * r["ratio_of_means"][0], rel=1e-9)


# ---------------------------------------------------------------- properness now depends on d

@pytest.mark.parametrize("d,nuisance,rule,prior,n", [
    (d, nu, rule, prior, n)
    for d in (1, 2, 3)
    for nu in ("radius", "acceleration", "force")
    for rule in ("single", "likelihood_product", "posterior_mass", "posterior_log")
    for prior in ([], ["flat_mass"], ["flat_inverse_mass"], [{"uniform_angle": {"center": "first_reading"}}],
                  [{"angle_power": {"sin": 0.5, "cos": -1.0, "center": 1.3}}])
    for n in (1, 3)
    if not (rule == "single" and n != 1)
])
def test_static_tail_slopes_match_the_computed_density(d, nuisance, rule, prior, n):
    combine = {"rule": {"posterior_mass": "posterior_product", "posterior_log": "posterior_product"}.get(rule, rule),
               "prior": prior}
    if rule == "posterior_mass":
        combine["coordinate"] = "mass"
    if rule == "posterior_log":
        combine["coordinate"] = "log_mass"
    rng = np.random.default_rng(zlib.crc32(f"{d}{nuisance}{rule}{n}{prior}".encode()))
    readings = ReadingSet(rng.normal(size=(1, n, d)) * 2 + 1, rng.normal(size=(1, n, d)) + 0.5,
                          rng.uniform(0.5, 2, size=(1, n)), rng.uniform(0.5, 2, size=(1, n)))
    law = {"reference": "cartesian", "nuisance": nuisance}
    u = np.array([[-80.0, -70.0, 70.0, 80.0]])
    lp, _ = joint_log_density(readings, u, law, combine)
    expected_plus, expected_minus, _ = tail_slopes(law, combine, n, d)
    assert (lp[0, 3] - lp[0, 2]) / 10 == pytest.approx(expected_plus, abs=1e-6)
    assert (lp[0, 1] - lp[0, 0]) / 10 == pytest.approx(expected_minus, abs=1e-6)


def test_flat_in_mass_completion_is_refused_in_1d_only():
    """Record section 2: alpha^(d-2) df d(alpha) dOmega leaves the 1D single-reading law improper."""
    est = Estimator({"name": "plain", "reduce": "none",
                     "law": {"reference": "cartesian", "nuisance": "acceleration"},
                     "combine": {"rule": "likelihood_product", "prior": ["flat_mass"]},
                     "readouts": ["ratio_of_means"]})
    world = {"readings": 1, "design": "same_pair"}
    with pytest.raises(ValueError, match="cannot be normalised"):
        est.check_world(dict(world, dimension=1))
    est.check_world(dict(world, dimension=2))
    est.check_world(dict(world, dimension=3))


def test_unknown_reference_is_refused():
    with pytest.raises(ValueError, match="reference must be one of"):
        Estimator({"name": "x", "reduce": "none", "law": {"reference": "tube", "nuisance": "radius"},
                   "combine": SINGLE, "readouts": ["median"]})


def test_mean_radius_2d_is_shared():
    h = np.linspace(0, 10, 11)
    assert np.array_equal(cartesian_mean_radius(h, 2), mean_radius(h, 2))


def test_angle_power_contains_uniform_angle_and_converts_between_completions():
    """{sin: 1, cos: 1} is uniform_angle; and the flat-in-m completion of the premise
    (alpha^(d-2) df d(alpha) dOmega) is the cartesian law times cos(theta)^(d-2)."""
    f, a, sf, sa = np.array([2.0, 0.5, -0.3]), np.array([0.8, 0.1, 0.2]), 1.0, 0.7
    readings = _series(f, a, sf, sa)
    u = np.linspace(-5, 5, 101)[None]
    one = joint_log_density(readings, u, FLAT, {"rule": "likelihood_product",
                                                "prior": [{"uniform_angle": {"center": 2.0}}]})[0]
    two = joint_log_density(readings, u, FLAT, {"rule": "likelihood_product",
                                                "prior": [{"angle_power": {"sin": 1, "cos": 1, "center": 2.0}}]})[0]
    assert np.ptp(one - two) < 1e-12
    plain = posterior_summaries(readings, {"reference": "cartesian", "nuisance": "acceleration"},
                                {"rule": "likelihood_product", "prior": ["flat_mass"]}, READOUTS)
    tilted = posterior_summaries(readings, CART, {"rule": "single", "prior": [
        {"angle_power": {"sin": 0, "cos": 1, "center": "first_reading"}}]}, READOUTS)
    for key in READOUTS:
        assert np.allclose(plain[key], tilted[key], rtol=1e-9), key
