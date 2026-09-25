"""The two cartesian references (24 September): flat in the pair as a vector.

cartesian2 (7081a44) is checked first, cartesian1 (665a321) after it, each by the
checks its own line wrote; then the two are checked against each other.

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
from scipy.integrate import dblquad, quad
from scipy.special import i0e, roots_legendre

from met.analyst import Estimator, estimate, joint_log_density, posterior_summaries, tail_slopes
from met.law import ReadingSet, cartesian1_mean_radius, cartesian2_mean_radius, mean_radius, radial_terms

CART = {"reference": "cartesian2", "nuisance": "radius"}
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
    log_k, radius = radial_terms(h, d, "cartesian2")
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
        cart = readings.reading_terms(u, nuisance, reference="cartesian2")
        for a, b in zip(flat, cart):
            assert np.allclose(a, b, rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("d", [1, 3])
def test_mean_radius_is_continuous_at_zero(d):
    hs = np.array([0.0, 1e-9, 1e-7, 1e-6, 1.0001e-6, 1e-5])
    at_zero = {1: math.sqrt(2 / math.pi), 3: 2 * math.sqrt(2 / math.pi)}[d]
    assert np.allclose(cartesian2_mean_radius(hs, d), at_zero, rtol=1e-9)


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
        chi = cartesian2_mean_radius(np.linalg.norm(a_hat, axis=1) / sigma_e, d)
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
    law = {"reference": "cartesian2", "nuisance": "acceleration"}
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
    law = {"reference": "cartesian2", "nuisance": nuisance}
    u = np.array([[-80.0, -70.0, 70.0, 80.0]])
    lp, _ = joint_log_density(readings, u, law, combine)
    expected_plus, expected_minus, _ = tail_slopes(law, combine, n, d)
    assert (lp[0, 3] - lp[0, 2]) / 10 == pytest.approx(expected_plus, abs=1e-6)
    assert (lp[0, 1] - lp[0, 0]) / 10 == pytest.approx(expected_minus, abs=1e-6)


def test_flat_in_mass_completion_is_refused_in_1d_only():
    """Record section 2: alpha^(d-2) df d(alpha) dOmega leaves the 1D single-reading law improper."""
    est = Estimator({"name": "plain", "reduce": "none",
                     "law": {"reference": "cartesian2", "nuisance": "acceleration"},
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
    assert np.array_equal(cartesian2_mean_radius(h, 2), mean_radius(h, 2))


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
    plain = posterior_summaries(readings, {"reference": "cartesian2", "nuisance": "acceleration"},
                                {"rule": "likelihood_product", "prior": ["flat_mass"]}, READOUTS)
    tilted = posterior_summaries(readings, CART, {"rule": "single", "prior": [
        {"angle_power": {"sin": 0, "cos": 1, "center": "first_reading"}}]}, READOUTS)
    for key in READOUTS:
        assert np.allclose(plain[key], tilted[key], rtol=1e-9), key


def test_the_new_mass_law_is_the_21_september_2d_law_in_every_dimension():
    """The cartesian kernel is exp(h^2/2) for every d, the flat one only in 2D. So a
    reading's cartesian mass law (not its E[alpha | m]) equals the 21 September law of a
    2D reading with the same |x|, |y| and x.y."""
    rng = np.random.default_rng(70)
    f3, a3 = rng.normal(size=3) * 2, rng.normal(size=3)
    x, y = np.linalg.norm(f3), np.linalg.norm(a3)
    cos_xy = f3 @ a3 / (x * y)
    f2, a2 = np.array([x, 0.0]), np.array([y * cos_xy, y * math.sqrt(1 - cos_xy**2)])
    u = np.linspace(-6, 6, 241)[None]
    lp3, _ = joint_log_density(_series(f3, a3, 1.0, 1.0), u, CART, SINGLE)
    lp2, _ = joint_log_density(_series(f2, a2, 1.0, 1.0), u, FLAT, SINGLE)
    assert np.ptp(lp3 - lp2) < 1e-10
    # and an aligned reading gives the same law in 1D and 3D
    aligned3 = joint_log_density(_series([x, 0, 0], [y, 0, 0], 1.0, 1.0), u, CART, SINGLE)[0]
    aligned1 = joint_log_density(_series([x], [y], 1.0, 1.0), u, CART, SINGLE)[0]
    assert np.ptp(aligned3 - aligned1) < 1e-10


@pytest.mark.parametrize("d", [1, 3])
def test_instrument_free_completion_by_brute_force(d):
    """(f alpha)^((d-2)/2) df d(alpha) dOmega, the sigma-free member of the forced class,
    is the cartesian law times angle_power {sin: (d-2)/2, cos: (d-2)/2} at the instrument ratio."""
    rng = np.random.default_rng(80 + d)
    f, a = rng.normal(size=d) * 1.5, rng.normal(size=d)
    sf, sa = 1.3, 0.6
    k = (d - 2) / 2
    combine = {"rule": "single", "prior": [{"angle_power": {"sin": k, "cos": k, "center": "first_reading"}}]}
    mine = posterior_summaries(_series(f, a, sf, sa), CART, combine, ["ratio_of_means"],
                               {"span": 200.0} if d == 1 else None)
    # theta = (pi/2) sin^2(pi t/2) clusters nodes at both ends, where (f alpha)^k is singular in 1D
    x, w = roots_legendre(300)
    t = (x + 1) / 2
    theta = (math.pi / 2) * np.sin(math.pi * t / 2) ** 2
    wt = (w / 2) * (math.pi**2 / 2) * np.sin(math.pi * t / 2) * np.cos(math.pi * t / 2)
    r, wr = (x + 1) * 30.0, w * 30.0
    R, T = np.meshgrid(r, theta, indexing="ij")
    fm, al = sf * R * np.sin(T), sa * R * np.cos(T)
    z = np.linalg.norm(fm[..., None] * f / sf**2 + al[..., None] * a / sa**2, axis=-1)
    log_avg = z + (np.log(0.5 * (1 + np.exp(-2 * z))) if d == 1 else np.log(-np.expm1(-2 * z) / (2 * z)))
    log_int = -0.5 * R**2 + log_avg + np.log(R) + k * np.log(fm * al)
    dens = np.exp(log_int - log_int.max()) * wr[:, None] * wt[None, :]
    assert mine["ratio_of_means"][0] == pytest.approx(float(np.sum(dens * fm) / np.sum(dens * al)), rel=1e-9)


def test_axis_splitting_holds_for_cartesian_and_fails_for_flat_and_the_sigma_free_completion():
    """Axis splitting: one 3D reading, and its three axes treated as 1D readings under the
    symmetric rule (prior counted once), give the same law of mass. Only the cartesian law
    with the uniform-angle completion passes."""
    rng = np.random.default_rng(90)
    f, a = rng.normal(size=3) * 2 + 1, rng.normal(size=3) + 0.5
    u = np.linspace(-5, 5, 201)[None]
    prior = [{"uniform_angle": {"center": 1.0}}]

    def gap(law, prior3, prior1):
        one = joint_log_density(_series(f, a, 1.0, 1.0), u, law, {"rule": "likelihood_product", "prior": prior3})[0]
        three = joint_log_density(_series(f[:, None], a[:, None], 1.0, 1.0), u, law,
                                  {"rule": "likelihood_product", "prior": prior1})[0]
        return np.ptp(one - three)

    assert gap(CART, prior, prior) < 1e-10
    assert gap(FLAT, prior, prior) > 0.1
    free3 = prior + [{"angle_power": {"sin": 0.5, "cos": 0.5, "center": 1.0}}]
    free1 = prior + [{"angle_power": {"sin": -0.5, "cos": -0.5, "center": 1.0}}]
    assert gap(CART, free3, free1) > 0.1


# ---------------------------------------------------------------- cartesian1: checks from its own line (665a321)
# Written independently of the checks above: the mean radius by 2D quadrature about
# w's axis, the whole-estimator 2D identity, the single-reading law by per-component
# brute force, and the von Mises state being exact under the cartesian1 reference.

CART1 = {"reference": "cartesian1", "nuisance": "radius"}


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("h", [0.0, 0.4, 1.5, 4.0])
def test_mean_radius_is_the_noncentral_chi_mean(d, h):
    """E|v| for v ~ N(w, I_d), |w| = h, by direct integration in polar form about w's axis."""
    if d == 1:
        ref = quad(lambda v: abs(v) * math.exp(-(v - h) ** 2 / 2), -np.inf, np.inf)[0] / math.sqrt(2 * math.pi)
    elif d == 2:
        ref = dblquad(lambda phi, r: r * r * math.exp(-(r * r - 2 * r * h * math.cos(phi) + h * h) / 2),
                      0, 12 + h, 0, 2 * math.pi)[0] / (2 * math.pi)
    else:
        ref = dblquad(lambda c, r: r ** 3 * math.exp(-(r * r - 2 * r * h * c + h * h) / 2),
                      0, 12 + h, -1, 1)[0] * 2 * math.pi / (2 * math.pi) ** 1.5
    assert float(cartesian1_mean_radius(h, d)) == pytest.approx(ref, rel=1e-7)


def test_cartesian_equals_flat_in_two_dimensions():
    for seed in range(3):
        rng = np.random.default_rng(seed)
        f, a = rng.normal(size=2) * 3, rng.normal(size=2)
        law_c = estimate(f, a, 0.8, 1.1, READOUTS, law=CART1)
        law_f = estimate(f, a, 0.8, 1.1, READOUTS, law=FLAT)
        for key in READOUTS:
            assert np.allclose(law_c[key], law_f[key], rtol=1e-10), key


@pytest.mark.parametrize("d", [1, 3])
def test_cartesian_differs_from_flat_outside_two_dimensions(d):
    f, a = np.full(d, 1.0), np.full(d, 0.6)
    assert estimate(f, a, 1, 1, law=CART1)["ratio_of_means"] != pytest.approx(
        estimate(f, a, 1, 1, law=FLAT)["ratio_of_means"], rel=1e-6)


@pytest.mark.parametrize("d", [1, 3])
def test_single_reading_law_by_brute_force(d):
    """Integrate exp(-|x - v sin t|^2/2 - |y - v cos t|^2/2) over v in R^d directly, on a theta grid."""
    rng = np.random.default_rng(40 + d)
    x, y = rng.normal(size=d) * 2 + 1, rng.normal(size=d) + 0.5
    theta = np.linspace(0.05, math.pi / 2 - 0.05, 7)
    brute = []
    for t in theta:
        w = x * math.sin(t) + y * math.cos(t)
        # the v-integral factorises over components; do each 1D integral numerically
        val = 1.0
        for k in range(d):
            val *= quad(lambda v: math.exp(-(x[k] - v * math.sin(t)) ** 2 / 2 - (y[k] - v * math.cos(t)) ** 2 / 2),
                        -np.inf, np.inf)[0]
        brute.append(math.log(val))
    readings = ReadingSet(x[None, None], y[None, None], [[1.0]], [[1.0]])
    u = np.log(np.tan(theta))[None]
    like, _, _ = readings.reading_terms(u, "radius", reference="cartesian1")
    diff = np.array(brute) - like[0, 0]
    assert np.ptp(diff) < 1e-8


@pytest.mark.parametrize("d", [1, 3])
def test_vonmises_state_is_exact_under_cartesian_in_every_dimension(d):
    rng = np.random.default_rng(50 + d)
    n, b = 7, 3
    v = rng.normal(size=(b, n, d))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    r = ReadingSet(2.0 * v + rng.normal(size=(b, n, d)), v + 0.6 * rng.normal(size=(b, n, d)),
                   np.full((b, n), 1.0), np.full((b, n), 0.6))
    old = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
    a = posterior_summaries(r, CART1, {"rule": "sequential", "old": old, "carry": "curve_only"}, READOUTS)
    b_ = posterior_summaries(r, CART1, {"rule": "sequential", "old": old, "carry": "vonmises_state"}, READOUTS)
    for key in READOUTS:
        assert np.allclose(a[key], b_[key], rtol=1e-9), key


# ---------------------------------------------------------------- cartesian1 against cartesian2
# The two implementations were written independently on 24 September, each
# believing it was the cartesian reference. Where the maths says they agree they
# must agree to rounding; where it says they differ, the difference must be
# exactly the predicted factor, and nothing else.

def _random_readings(seed, b, n, d):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(b, n, d))
    return ReadingSet(1.7 * v + rng.normal(size=(b, n, d)), v + 0.8 * rng.normal(size=(b, n, d)),
                      rng.uniform(0.5, 2.0, size=(b, n)), rng.uniform(0.5, 2.0, size=(b, n)))


@pytest.mark.parametrize("d", [1, 2, 3])
def test_the_two_mean_radii_agree(d):
    h = np.concatenate([[0.0, 1e-9, 1e-7, 1e-6, 2e-6], np.linspace(0.0, 40.0, 4001)])
    assert np.allclose(cartesian1_mean_radius(h, d), cartesian2_mean_radius(h, d), rtol=1e-14, atol=0)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_radius_nuisance_terms_agree(d):
    r = _random_readings(60 + d, 4, 6, d)
    u = np.log(r.s[:, :1]) + np.linspace(-8, 8, 801)[None, :]
    one = r.reading_terms(u, "radius", reference="cartesian1")
    two = r.reading_terms(u, "radius", reference="cartesian2")
    for a, b in zip(one, two):
        assert np.allclose(a, b, rtol=1e-14, atol=1e-14)


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("nuisance", ["radius", "acceleration", "force"])
def test_single_reading_law_agrees_under_every_split(d, nuisance):
    r = _random_readings(70 + d, 4, 3, d)
    u = np.log(r.s[:, :1]) + np.linspace(-8, 8, 801)[None, :]
    l1, f1, c1 = r.reading_terms(u, nuisance, reference="cartesian1")
    l2, f2, c2 = r.reading_terms(u, nuisance, reference="cartesian2")
    assert np.allclose(l1 + f1, l2 + f2, rtol=1e-12, atol=1e-12)
    assert np.allclose(c1, c2, rtol=1e-14, atol=0)


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("nuisance", ["acceleration", "force"])
def test_the_splits_differ_by_exactly_the_predicted_factor(d, nuisance):
    """cartesian1's likelihood is cartesian2's times cos^(2-d) (acceleration) or sin^(2-d) (force)."""
    r = _random_readings(80 + d, 4, 5, d)
    u = np.log(r.s[:, :1]) + np.linspace(-8, 8, 801)[None, :]
    l1 = r.reading_terms(u, nuisance, reference="cartesian1")[0]
    l2 = r.reading_terms(u, nuisance, reference="cartesian2")[0]
    z = u[:, None, :] - np.log(r.s)[:, :, None]
    log_trig = -0.5 * np.logaddexp(0.0, 2.0 * z) if nuisance == "acceleration" else -0.5 * np.logaddexp(0.0, -2.0 * z)
    assert np.allclose(l1 - l2, (2 - d) * log_trig, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("d", [1, 3])
def test_combined_laws_agree_for_radius_and_part_for_the_other_splits(d):
    """With several readings: identical under the radius nuisance (and in 2D always);
    different under the acceleration and force nuisances outside 2D."""
    r = _random_readings(90 + d, 3, 6, d)
    prior = [{"uniform_angle": {"center": "first_reading"}}]
    combine = {"rule": "likelihood_product", "prior": prior}
    law = lambda ref, nuisance: {"reference": ref, "nuisance": nuisance}
    one = posterior_summaries(r, law("cartesian1", "radius"), combine, READOUTS)
    two = posterior_summaries(r, law("cartesian2", "radius"), combine, READOUTS)
    for key in READOUTS:
        assert np.allclose(one[key], two[key], rtol=1e-12), key
    flat_mass = {"rule": "likelihood_product", "prior": ["flat_mass"]}
    one = posterior_summaries(r, law("cartesian1", "acceleration"), flat_mass, ["median"])
    two = posterior_summaries(r, law("cartesian2", "acceleration"), flat_mass, ["median"])
    assert np.all(np.isfinite(one["median"])) and np.all(np.isfinite(two["median"]))
    assert not np.allclose(one["median"], two["median"], rtol=1e-6)
    r2 = _random_readings(99, 3, 6, 2)
    one = posterior_summaries(r2, law("cartesian1", "acceleration"), flat_mass, READOUTS)
    two = posterior_summaries(r2, law("cartesian2", "acceleration"), flat_mass, READOUTS)
    for key in READOUTS:
        assert np.all(np.isfinite(one[key])), key
        assert np.allclose(one[key], two[key], rtol=1e-12), key


def test_the_old_name_is_refused_with_the_reason():
    base = {"name": "e", "reduce": "none", "readouts": ["median"],
            "combine": {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}}
    with pytest.raises(ValueError, match="cartesian1 .*cartesian2"):
        Estimator(dict(base, law={"reference": "cartesian", "nuisance": "radius"}))
