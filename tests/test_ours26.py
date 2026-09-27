"""ours26: ours25b with a different noise SD on each axis.

The law is the per-axis Beale sum (axis splitting), the readout keeps its definition
E|F*_stack| / E|a*_stack| with a weighted noncentral chi mean as its weight, and the
prior's centre is declared, because per-axis noise does not force it.

Checks, each against a computation done a different way: the weighted mean against
the closed forms, polar-coordinate integration and nested quadrature; the whole readout
against an independent quadrature of one 2D reading; the law against the Beale formula.
Then what should hold: isotropic data give ours25b exactly, a vector reading and its
axes agree, units and the reciprocal symmetry, and the refusals.
"""

import math

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from met.analyst import Estimator, estimate, known_excitation_summaries, posterior_summaries
from met.law import ReadingSet, StackedState, stacked_mean_radius, weighted_stacked_mean
from met.world import generate, per_axis_noise, series_rngs, validate_cell

CART2 = {"reference": "cartesian2", "nuisance": "radius"}
CART1 = {"reference": "cartesian1", "nuisance": "radius"}
READOUTS = ["ratio_of_means_rss", "median", "interval_95"]


def _combine(center):
    return {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": center}}]}


GEO = _combine("first_reading_geometric")
SF, SA = np.array([0.9, 0.4, 1.5]), np.array([0.5, 0.5, 0.2])      # three different ratios


def _series(seed, n, d=3):
    rng = np.random.default_rng(seed)
    return rng.normal(1.0, 1.0, (n, d)), rng.normal(0.4, 0.6, (n, d))


# ---------------------------------------------------------------- the weighted mean

@pytest.mark.parametrize("k", [1, 2, 3, 7, 40])
@pytest.mark.parametrize("h", [0.0, 0.3, 4.0, 60.0])
@pytest.mark.parametrize("c", [1e-3, 1.0, 50.0])
def test_one_group_is_the_noncentral_chi_mean(k, h, c):
    got = weighted_stacked_mean(np.array([c * c]), np.array([h * h]), np.array([k]))
    assert got == pytest.approx(c * stacked_mean_radius(np.array([h]), k)[0], rel=1e-9)


def _polar_mean(c1, c2, w1, w2, nr=200, nphi=256):
    """E sqrt(c1^2 v1^2 + c2^2 v2^2), v ~ N(w, I_2), in polar coordinates about the origin
    (Gauss-Legendre in the radius, the periodic trapezoid in the angle)."""
    top = math.hypot(w1, w2) + 14.0
    x, wx = leggauss(nr)
    rho, wr = (x + 1) * top / 2, wx * top / 2
    phi = np.arange(nphi) * 2 * math.pi / nphi
    r, f = np.meshgrid(rho, phi, indexing="ij")
    dens = np.exp(-0.5 * ((r * np.cos(f) - w1) ** 2 + (r * np.sin(f) - w2) ** 2)) / (2 * math.pi)
    g = r * np.sqrt(c1 ** 2 * np.cos(f) ** 2 + c2 ** 2 * np.sin(f) ** 2)
    return float(np.sum(wr[:, None] * (2 * math.pi / nphi) * g * dens * r))


@pytest.mark.parametrize("c1,c2,w1,w2", [(1, 3, 0.7, 2.0), (0.2, 1, 3, 0), (2, 0.5, 5.4, 1.4)])
def test_two_groups_by_polar_integration(c1, c2, w1, w2):
    got = weighted_stacked_mean(np.array([c1 * c1, c2 * c2]), np.array([w1 * w1, w2 * w2]), np.array([1, 1]))
    assert got == pytest.approx(_polar_mean(c1, c2, w1, w2), rel=1e-9)


def test_very_unequal_weights_by_nested_quadrature():
    def normal(v, mean):
        return math.exp(-0.5 * (v - mean) ** 2) / math.sqrt(2 * math.pi)

    def inner(v1):
        return quad(lambda v2: math.sqrt(v1 * v1 + 1e-4 * v2 * v2) * normal(v2, 10.0), -5, 25,
                    epsabs=0, epsrel=1e-13, limit=200)[0]
    ref = quad(lambda v1: inner(v1) * normal(v1, 0.0), -14, 14, points=[0], epsabs=0, epsrel=1e-13, limit=400)[0]
    got = weighted_stacked_mean(np.array([1.0, 1e-4]), np.array([0.0, 100.0]), np.array([1, 1]))
    assert got == pytest.approx(ref, rel=1e-9)


# ---------------------------------------------------------------- the law and the readout

def test_law_is_the_per_axis_beale_sum():
    f, a = _series(1, 4)
    r = ReadingSet(f[None], a[None], np.broadcast_to(SF, (1, 4, 3)), np.broadcast_to(SA, (1, 4, 3)))
    u = np.log(1.6) + np.linspace(-4, 4, 301)[None]
    like, ref, acc = r.reading_terms(u, "radius", reference="cartesian2")
    m = np.exp(u[0])
    beale = -np.sum((f[:, :, None] - m * a[:, :, None]) ** 2
                    / (2 * (SF[None, :, None] ** 2 + m ** 2 * SA[None, :, None] ** 2)), axis=(0, 1))
    assert np.ptp(like[0].sum(axis=0) - beale) < 1e-11
    assert ref is None and acc is None


def test_readout_against_independent_quadrature():
    """One 2D reading, two ratios: the posterior by Gauss-Legendre in log m, the weight by
    polar integration at each node."""
    f, a = np.array([[6.0, 2.0]]), np.array([[3.0, 1.2]])
    sf, sa = np.array([0.5, 1.0]), np.array([0.3, 0.2])
    got = estimate(f, a, sf, sa, ["ratio_of_means_rss"], CART2, GEO)["ratio_of_means_rss"]
    x, y, s = f[0] / sf, a[0] / sa, sf / sa
    c0 = np.mean(np.log(s))

    def log_post(u):
        th = np.arctan(np.exp(u)[:, None] / s)
        return np.sum(-0.5 * (x * np.cos(th) - y * np.sin(th)) ** 2, axis=1) - np.logaddexp(u - c0, c0 - u)
    wide = np.linspace(-20, 20, 40001)
    lw = log_post(wide)
    keep = wide[lw > lw.max() - 45]
    t, wt = leggauss(300)
    u = keep[0] + (t + 1) * (keep[-1] - keep[0]) / 2
    dens = wt * np.exp(log_post(u) - lw.max())
    th = np.arctan(np.exp(u)[:, None] / s)
    w, c = x * np.sin(th) + y * np.cos(th), sa * np.cos(th)
    acc = np.array([_polar_mean(c[i, 0], c[i, 1], w[i, 0], w[i, 1]) for i in range(len(u))])
    assert got == pytest.approx(np.sum(dens * np.exp(u) * acc) / np.sum(dens * acc), rel=1e-8)


def test_spline_carries_the_weight_to_the_grid():
    f, a = _series(2, 4)
    r = ReadingSet(f[None], a[None], np.broadcast_to(SF, (1, 4, 3)), np.broadcast_to(SA, (1, 4, 3)))
    u = np.log(1.6) + np.linspace(-3, 3, 2001)[None]
    assert np.allclose(r.stacked_cond_alpha(u, reference="cartesian2"), r._weighted_cond_alpha_at(u), rtol=2e-8)


# ---------------------------------------------------------------- what should hold

@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("n", [1, 5])
@pytest.mark.parametrize("law", [CART2, CART1])
def test_isotropic_data_give_ours25b(d, n, law):
    """With one SD on every axis, ours26 (per-axis SDs, geometric centre) is ours25b."""
    f, a = _series(10 * d + n, n, d)
    ours25b = estimate(f, a, 0.9, 0.5, READOUTS, law, _combine("first_reading"))
    ours26 = estimate(f, a, [0.9] * d, [0.5] * d, READOUTS, law, GEO)
    for key in READOUTS:
        assert np.allclose(ours26[key], ours25b[key], rtol=1e-12), key


def test_unequal_sds_between_readings_use_the_weighted_mean():
    """Per-reading SDs that differ within a series (ours25b refused these): the readout's
    weight is the weighted mean, which here is checked by splitting the stacked vector
    into its two noise groups by hand."""
    f, a = _series(3, 2)
    r = ReadingSet(f[None], a[None], np.array([[1.0, 2.0]]), np.array([[0.5, 0.5]]))
    u = np.log(2.0) + np.linspace(-2, 2, 7)[None]
    got = r.stacked_cond_alpha(u, reference="cartesian2")
    s = np.array([2.0, 4.0])
    th = np.arctan(np.exp(u[0])[:, None] / s)
    x, y = f / np.array([1.0, 2.0])[:, None], a / 0.5
    w2 = np.sum((x[None] * np.sin(th)[:, :, None] + y[None] * np.cos(th)[:, :, None]) ** 2, axis=2)
    ref = weighted_stacked_mean((0.5 * np.cos(th)) ** 2, w2, np.array([3, 3]))
    assert np.allclose(got[0], ref, rtol=1e-12)


@pytest.mark.parametrize("center", ["first_reading_geometric", 1.7])
def test_vector_and_axes_agree(center):
    """Axis splitting: a 3D series with per-axis noise and the same series as 3N one-dimensional
    readings give one ours26 answer (the first-reading centre means the same reading after splitting)."""
    rng = np.random.default_rng(5)
    f, a = rng.normal(size=(2, 4, 3)), rng.normal(0.3, 1.0, size=(2, 4, 3))
    readings = ReadingSet(f, a, np.broadcast_to(SF, (2, 4, 3)), np.broadcast_to(SA, (2, 4, 3)))
    spec = {"name": "e", "law": CART2, "combine": _combine(center), "readouts": READOUTS}
    vector = Estimator(dict(spec, reduce="none")).evaluate(readings)
    axes = Estimator(dict(spec, reduce="split_axes")).evaluate(readings)
    assert np.allclose(axes["ratio_of_means_rss"], vector["ratio_of_means_rss"], rtol=1e-9)
    assert np.allclose(axes["median"], vector["median"], rtol=1e-8)


def test_series_enters_only_through_the_sums_of_each_noise_pair():
    """Mixing the readings by a rotation, separately on each axis (each axis is one noise
    pair here), keeps every pair's sum x.x, y.y, x.y and so the answer; mixing all axes
    together would not."""
    rng = np.random.default_rng(12)
    f, a = _series(11, 3)
    rot = [np.linalg.qr(rng.normal(size=(3, 3)))[0] for _ in range(3)]
    f2 = np.stack([rot[k] @ f[:, k] for k in range(3)], axis=1)
    a2 = np.stack([rot[k] @ a[:, k] for k in range(3)], axis=1)
    r1 = estimate(f, a, SF, SA, READOUTS, CART2, _combine(2.0))
    r2 = estimate(f2, a2, SF, SA, READOUTS, CART2, _combine(2.0))
    for key in ("ratio_of_means_rss", "median"):
        assert r2[key] == pytest.approx(r1[key], rel=1e-10)


def test_zero_readings_median_is_the_centre_but_the_readout_is_not():
    out = estimate(np.zeros((1, 3)), np.zeros((1, 3)), [1.0, 1.0, 1.0], [0.2, 0.5, 1.25],
                   ["ratio_of_means_rss", "median"], CART2, GEO)
    assert out["median"] == pytest.approx(2.0, rel=1e-10)
    assert abs(out["ratio_of_means_rss"] / 2.0 - 1) > 0.05


def test_units_and_reciprocal():
    f, a = _series(6, 4)
    base = estimate(f, a, SF, SA, READOUTS, CART2, GEO)
    scaled = estimate(7.0 * f, a, 7.0 * SF, SA, READOUTS, CART2, GEO)
    swapped = estimate(a, f, SA, SF, READOUTS, CART2, GEO)
    for key in ("ratio_of_means_rss", "median"):
        assert scaled[key] == pytest.approx(7.0 * base[key], rel=1e-9)
        assert swapped[key] * base[key] == pytest.approx(1.0, rel=1e-9)


def test_pooling_is_per_axis_inverse_variance():
    rng = np.random.default_rng(7)
    f, a = rng.normal(size=(1, 3, 2)), rng.normal(size=(1, 3, 2))
    sf = np.array([[[1.0, 2.0], [0.5, 2.0], [1.0, 4.0]]])
    r = ReadingSet(f, a, sf, 1.0).pooled()
    w = 1 / sf ** 2
    assert np.allclose(r.force[0, 0], np.sum(w * f, axis=1)[0] / np.sum(w, axis=1)[0])
    assert np.allclose(r.force_sd_axes[0, 0], 1 / np.sqrt(np.sum(w, axis=1)[0]))


def test_known_excitation_with_per_axis_noise():
    rng = np.random.default_rng(8)
    true_a = rng.normal(size=(2, 3, 3))
    f = 2.0 * true_a + SF * rng.normal(size=(2, 3, 3))
    r = ReadingSet(f, true_a, np.broadcast_to(SF, (2, 3, 3)), np.broadcast_to(SA, (2, 3, 3)))
    got = known_excitation_summaries(r, true_a, ["median"])["median"]
    w = 1 / SF ** 2
    mu = np.sum(w * f * true_a, axis=(1, 2)) / np.sum(w * true_a ** 2, axis=(1, 2))
    assert np.allclose(got, mu, rtol=1e-6)     # far from zero, the truncation does not move the median


# ---------------------------------------------------------------- the world

def _world(**changes):
    world = {"dimension": 3, "design": "new_excitation", "mass": 2.0, "acceleration_snr": 3.0,
             "direction": "random", "readings": 4,
             "noise": {"force_sd": {"per_axis": [1.0, 1.0, 1.0]}, "acceleration_sd": {"per_axis": [0.2, 0.5, 1.25]}},
             "supplied_noise": "exact"}
    world.update(changes)
    return world


def test_world_draws_per_axis_noise():
    world = _world(acceleration_snr=0.0, readings=50)
    validate_cell(world)
    readings, truth = generate(world, series_rngs(3, world, 40))
    assert np.allclose(readings.acceleration_sd_axes, [0.2, 0.5, 1.25])
    spread = np.std(readings.acceleration.reshape(-1, 3), axis=0)
    assert np.allclose(spread, [0.2, 0.5, 1.25], rtol=0.05)
    assert per_axis_noise(world) == (True, True)


def test_world_snr_is_against_the_rms_axis_sd():
    world = _world(direction="fixed", design="same_pair", readings=1)
    _, truth = generate(world, series_rngs(3, world, 2))
    rms = math.sqrt((0.2 ** 2 + 0.5 ** 2 + 1.25 ** 2) / 3)
    assert np.allclose(np.linalg.norm(truth["true_acceleration"], axis=2), 3.0 * rms)


def test_world_refuses_a_wrong_axis_count():
    with pytest.raises(ValueError, match="3 entries"):
        validate_cell(_world(noise={"force_sd": 1.0, "acceleration_sd": {"per_axis": [0.2, 0.5]}}))


def test_isotropic_world_draws_are_unchanged():
    """A per-axis setting with equal SDs draws the same numbers as the plain setting."""
    plain = _world(noise={"force_sd": 1.0, "acceleration_sd": 0.5})
    axes = _world(noise={"force_sd": {"per_axis": [1.0] * 3}, "acceleration_sd": {"per_axis": [0.5] * 3}})
    r1, _ = generate(plain, series_rngs(3, plain, 5))
    r2, _ = generate(axes, series_rngs(3, plain, 5))
    assert np.array_equal(r1.force, r2.force) and np.array_equal(r1.acceleration, r2.acceleration)
    assert r2.isotropic and np.array_equal(r1.force_sd, r2.force_sd)


# ---------------------------------------------------------------- refusals

def _spec(**changes):
    spec = {"name": "e", "reduce": "none", "law": CART2, "combine": GEO, "readouts": ["ratio_of_means_rss", "median"]}
    spec.update(changes)
    return spec


def test_world_level_refusals():
    world = _world()
    Estimator(_spec()).check_world(world)                              # ours26 itself
    Estimator(_spec(combine=_combine(2.0))).check_world(world)         # a numeric centre
    with pytest.raises(ValueError, match="not forced"):
        Estimator(_spec(combine=_combine("first_reading"))).check_world(world)
    with pytest.raises(ValueError, match="built for reference cartesian1 or cartesian2"):
        Estimator(_spec(law={"reference": "flat", "nuisance": "radius"}, readouts=["median"])).check_world(world)
    with pytest.raises(ValueError, match="built for reference cartesian1 or cartesian2"):
        Estimator(_spec(law={"reference": "cartesian2", "nuisance": "acceleration"},
                        combine={"rule": "likelihood_product", "prior": ["flat_log_mass"]},
                        readouts=["median"])).check_world(world)
    with pytest.raises(ValueError, match="ratio_of_means_rss"):
        Estimator(_spec(readouts=["ratio_of_means"])).check_world(world)
    # split first, and every one-dimensional law is available again (the centre still is not forced)
    Estimator(_spec(reduce="split_axes", law={"reference": "flat", "nuisance": "radius"},
                    readouts=["ratio_of_means", "median"])).check_world(world)
    with pytest.raises(ValueError, match="not forced"):
        Estimator(_spec(reduce="split_axes", combine=_combine("first_reading"))).check_world(world)


def test_one_ratio_on_every_axis_keeps_the_first_reading_centre():
    """sigma_F,k / sigma_a,k the same on every axis, the SDs themselves not: the angle is common,
    so the centre is forced again, and ours26 with first_reading runs."""
    world = _world(noise={"force_sd": {"per_axis": [0.4, 1.0, 2.5]}, "acceleration_sd": {"per_axis": [0.2, 0.5, 1.25]}})
    assert per_axis_noise(world) == (True, False)
    est = Estimator(_spec(combine=_combine("first_reading")))
    est.check_world(world)
    readings, _ = generate(world, series_rngs(3, world, 3))
    out = est.evaluate(readings)
    geo = Estimator(_spec()).evaluate(readings)
    assert np.allclose(out["ratio_of_means_rss"], geo["ratio_of_means_rss"], rtol=1e-12)


def test_data_level_refusals():
    f, a = _series(9, 3)
    r = ReadingSet(f[None], a[None], np.broadcast_to(SF, (1, 3, 3)), np.broadcast_to(SA, (1, 3, 3)))
    u = np.zeros((1, 5))
    with pytest.raises(ValueError, match="per-axis noise"):
        r.reading_terms(u, "radius", reference="flat")
    with pytest.raises(ValueError, match="per-axis noise"):
        r.force_sd
    with pytest.raises(ValueError, match="not forced"):
        posterior_summaries(r, CART2, _combine("first_reading"), ["median"])
    with pytest.raises(ValueError, match="use ratio_of_means_rss|ours26's readout"):
        posterior_summaries(r, CART2, GEO, ["ratio_of_means"])
    with pytest.raises(ValueError, match="not built for per-axis noise"):
        posterior_summaries(r, CART2, {"rule": "posterior_product", "coordinate": "log_mass",
                                       "prior": []}, ["median"])
    with pytest.raises(ValueError, match="stacked state"):
        StackedState.from_readings(r)
