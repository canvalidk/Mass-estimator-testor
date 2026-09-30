"""ours29 (30 September): the cartesian likelihood product with each reading's log-likelihood
weighted by its own misfit (misfit: kink | smooth) and the product raised to the Bartlett
power eta = min(1, curv/J) (calibrate: bartlett).

Checks: the closed form against the numerical path and against a brute-force grid on the
whole projective line; Z (exact zero readings change nothing); one reading gives eta = 1;
rotation, swap and units; the misfit weight's two ends (g = 1 in 1D, the impossible pair's
bounded tilt); the relation to ours25's sandwich; validation.
"""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from met.analyst import Estimator, bartlett_weight, estimate, joint_log_density, misfit_weights, sandwich_weight
from met.law import ReadingSet

CART = {"reference": "cartesian2", "nuisance": "radius"}
PRIOR = [{"uniform_angle": {"center": "first_reading"}}]
OURS29 = {"rule": "likelihood_product", "prior": PRIOR, "misfit": "kink", "calibrate": "bartlett"}
READOUTS = ["ratio_of_means_rss", "median", "interval_95"]


def _series(rng, n, d, m=1.0, snr=0.0, sf=1.0, sa=1.0):
    u = rng.normal(size=(n, d))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    a_star = snr * sa * u
    force = m * a_star + sf * rng.normal(size=(n, d))
    acceleration = a_star + sa * rng.normal(size=(n, d))
    return force, acceleration


def _set(force, acceleration, sf=1.0, sa=1.0):
    n = force.shape[0]
    return ReadingSet(force[None], acceleration[None], np.full((1, n), sf), np.full((1, n), sa))


def _brute_eta(force, acceleration, sf, sa, g):
    """eta by a fine grid on the whole projective line, independent of the closed form."""
    phi = np.linspace(-np.pi / 2, np.pi / 2, 200001, endpoint=False)
    s = sf / sa
    cp, sp = np.cos(phi), np.sin(phi)
    x, y = force / sf, acceleration / sa
    # h_i^2 = |x sin + y cos|^2 in the angle theta with m = s tan(theta); l_i = g_i h_i^2 / 2
    hx = np.einsum("nd,k->nkd", x, sp) + np.einsum("nd,k->nkd", y, cp)
    ell = 0.5 * g[:, None] * np.sum(hx ** 2, axis=2)
    total = ell.sum(axis=0)
    j = int(np.argmax(total))
    h = phi[1] - phi[0]
    psi = (ell[:, (j + 1) % len(phi)] - ell[:, j - 1]) / (2 * h)
    curv = -np.sum(ell[:, (j + 1) % len(phi)] - 2 * ell[:, j] + ell[:, j - 1]) / h ** 2
    return min(1.0, curv / np.sum(psi ** 2))


@pytest.mark.parametrize("d", [1, 2, 3])
def test_closed_form_matches_numeric_path_and_brute_force(d):
    rng = np.random.default_rng(10 + d)
    force, acceleration = _series(rng, 60, d, m=2.0, snr=0.3, sf=1.3, sa=0.7)
    closed = _set(force, acceleration, 1.3, 0.7)
    eta = bartlett_weight(closed, OURS29)[0]
    g = misfit_weights(closed, OURS29)[0]
    # per-axis SDs that differ by 1e-9 force the numerical path
    tiny = 1.3 * (1 + 1e-9 * np.arange(d))
    numeric = ReadingSet(force[None], acceleration[None], np.broadcast_to(tiny, (1, 60, d)), np.full((1, 60, d), 0.7))
    assert bartlett_weight(numeric, dict(OURS29))[0] == pytest.approx(eta, abs=1e-6)
    assert np.allclose(misfit_weights(numeric, dict(OURS29))[0], g, atol=1e-6)
    assert _brute_eta(force, acceleration, 1.3, 0.7, g) == pytest.approx(eta, rel=1e-4)
    assert 0 < eta <= 1


def test_misfit_weight_is_one_in_1d_and_equals_the_gram_eigenvalue_rule():
    rng = np.random.default_rng(3)
    f1, a1 = _series(rng, 40, 1)
    assert np.all(misfit_weights(_set(f1, a1), OURS29) == 1.0)
    f2, a2 = _series(rng, 40, 2)
    q = np.array([np.linalg.eigvalsh(np.array([[fi @ fi, fi @ ai], [fi @ ai, ai @ ai]]))[0] for fi, ai in zip(f2, a2)])
    g = misfit_weights(_set(f2, a2), OURS29)[0]
    assert np.allclose(g, np.where(q > 1, 1 / q, 1.0))
    smooth = misfit_weights(_set(f2, a2), dict(OURS29, misfit="smooth"))[0]
    assert np.allclose(smooth, 1 / (1 + q))


@pytest.mark.parametrize("d", [1, 3])
def test_exact_zero_readings_change_nothing(d):
    """Z: an exact zero reading has zero score and zero curvature, and no misfit."""
    rng = np.random.default_rng(20 + d)
    force, acceleration = _series(rng, 25, d, m=1.5, snr=0.5)
    base = estimate(force, acceleration, 1.0, 1.0, READOUTS[1:], law=CART, combine=OURS29)
    zf = np.vstack([force, np.zeros((40, d))])
    za = np.vstack([acceleration, np.zeros((40, d))])
    more = estimate(zf, za, 1.0, 1.0, READOUTS[1:], law=CART, combine=OURS29)
    assert bartlett_weight(_set(zf, za), OURS29)[0] == pytest.approx(bartlett_weight(_set(force, acceleration), OURS29)[0], abs=1e-12)
    assert more["median"] == pytest.approx(base["median"], rel=1e-9)
    assert np.allclose(more["interval_95"], base["interval_95"], rtol=1e-9)


def test_sandwich_is_not_z_invariant_but_bartlett_is():
    """ours25's Bessel-corrected J counts readings; appending zeros moves its weight."""
    rng = np.random.default_rng(18)
    force, acceleration = _series(rng, 12, 3, m=1.0, snr=0.0)            # pure noise, peak inside m > 0
    zf, za = np.vstack([force, np.zeros((12, 3))]), np.vstack([acceleration, np.zeros((12, 3))])
    sand = {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"}
    bart = {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "bartlett"}
    w0, w1 = sandwich_weight(_set(force, acceleration), CART, sand)[0], sandwich_weight(_set(zf, za), CART, sand)[0]
    b0, b1 = bartlett_weight(_set(force, acceleration), bart)[0], bartlett_weight(_set(zf, za), bart)[0]
    assert b0 < 1 and 0.01 < w0 < 1
    assert b0 == pytest.approx(b1, abs=1e-12)
    assert w1 / w0 == pytest.approx((12 / 11) / (24 / 23), rel=1e-3)   # the Bessel factor N/(N-1) moved


def test_one_reading_gives_eta_one():
    rng = np.random.default_rng(5)
    for d in (1, 2, 3):
        f, a = _series(rng, 1, d, snr=0.2)
        assert bartlett_weight(_set(f, a), OURS29)[0] == 1.0


def test_rotation_swap_and_units():
    rng = np.random.default_rng(6)
    force, acceleration = _series(rng, 30, 3, m=2.0, snr=0.4, sf=1.0, sa=0.5)
    base = estimate(force, acceleration, 1.0, 0.5, READOUTS, law=CART, combine=OURS29)
    eta = bartlett_weight(_set(force, acceleration, 1.0, 0.5), OURS29)[0]
    rots = Rotation.random(30, random_state=7).as_matrix()               # each reading its own frame
    rf, ra = np.einsum("nij,nj->ni", rots, force), np.einsum("nij,nj->ni", rots, acceleration)
    rot = estimate(rf, ra, 1.0, 0.5, READOUTS, law=CART, combine=OURS29)
    assert bartlett_weight(_set(rf, ra, 1.0, 0.5), OURS29)[0] == pytest.approx(eta, abs=1e-12)
    assert rot["median"] == pytest.approx(base["median"], rel=1e-9)
    assert rot["ratio_of_means_rss"] == pytest.approx(base["ratio_of_means_rss"], rel=1e-9)
    swap = estimate(acceleration, force, 0.5, 1.0, READOUTS, law=CART, combine=OURS29)
    assert bartlett_weight(_set(acceleration, force, 0.5, 1.0), OURS29)[0] == pytest.approx(eta, abs=1e-12)
    assert swap["median"] * base["median"] == pytest.approx(1.0, rel=1e-6)
    units = estimate(3.0 * force, acceleration, 3.0, 0.5, READOUTS, law=CART, combine=OURS29)
    assert units["median"] / base["median"] == pytest.approx(3.0, rel=1e-6)
    assert units["ratio_of_means_rss"] / base["ratio_of_means_rss"] == pytest.approx(3.0, rel=1e-6)


def test_impossible_pair_tilt_stays_bounded():
    """F = (2, 0), a = (0, 1): Q = 1/sigma^2 (the smaller Gram eigenvalue), |C| = 3/sigma^2, so the
    untempered tilt over m > 0 is |C|/2 = 1.5/sigma^2 and the kink form's is exactly
    nu a/(1 - a) = 1.5 nats (a = |C|/T = 0.6) once Q >= nu (readings-are-evidence record, T30)."""
    u = np.linspace(-14, 14, 40001)[None]
    plain = {"rule": "likelihood_product", "prior": PRIOR}
    kink = dict(plain, misfit="kink")
    for sigma in (1.0, 0.5, 0.1, 0.005):
        rs = _set(np.array([[2.0, 0.0]]), np.array([[0.0, 1.0]]), sigma, sigma)
        prior = joint_log_density(rs, u, CART, plain)[0][0] - \
            np.sum(rs.reading_terms(u, "radius", reference="cartesian2")[0], axis=1)[0]
        tempered = np.ptp(joint_log_density(rs, u, CART, kink)[0][0] - prior)
        untempered = np.ptp(joint_log_density(rs, u, CART, plain)[0][0] - prior)
        assert tempered == pytest.approx(1.5, abs=2e-3)
        assert untempered == pytest.approx(1.5 / sigma ** 2, rel=2e-3)


def test_bartlett_and_sandwich_agree_up_to_bessel_where_the_peak_is_interior():
    rng = np.random.default_rng(8)
    force, acceleration = _series(rng, 50, 3, m=1.0, snr=0.6)
    rs = _set(force, acceleration)
    sand = sandwich_weight(rs, CART, {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "sandwich"})[0]
    bart = bartlett_weight(rs, {"rule": "likelihood_product", "prior": PRIOR, "calibrate": "bartlett"})[0]
    assert sand < 1 and bart < 1
    assert sand == pytest.approx(bart * 49 / 50, rel=1e-3)


def test_misfit_leaves_the_readout_weight_alone():
    """g and eta change how much a reading counts, not the size of its candidate."""
    rng = np.random.default_rng(9)
    force, acceleration = _series(rng, 20, 2, snr=0.5)
    rs = _set(force, acceleration)
    u = np.linspace(-3, 3, 61)[None]
    _, acc0 = joint_log_density(rs, u, CART, {"rule": "likelihood_product", "prior": PRIOR})
    _, acc1 = joint_log_density(rs, u, CART, OURS29)
    assert np.array_equal(acc0, acc1)


def test_ours29_is_refused_where_it_cannot_apply():
    base = {"name": "x", "reduce": "none", "readouts": ["median"]}
    world = {"readings": 5, "design": "new_excitation", "dimension": 3}
    for law in ({"reference": "flat", "nuisance": "radius"}, {"reference": "cartesian2", "nuisance": "acceleration"}):
        prior = PRIOR if law["nuisance"] == "radius" else ["flat_mass"]
        est = Estimator(dict(base, law=law, combine={"rule": "likelihood_product", "prior": prior, "calibrate": "bartlett"}))
        with pytest.raises(ValueError, match="ours29"):
            est.check_world(world)
    with pytest.raises(ValueError, match="misfit must be one of"):
        Estimator(dict(base, law=CART, combine={"rule": "likelihood_product", "prior": PRIOR, "misfit": "magic"}))
    with pytest.raises(ValueError):
        Estimator(dict(base, law=CART, combine={"rule": "posterior_product", "coordinate": "mass", "prior": [],
                                                "misfit": "kink"}))
    Estimator(dict(base, law=CART, combine=OURS29)).check_world(dict(world, readings=1))
