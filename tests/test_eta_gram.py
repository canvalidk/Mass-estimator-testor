"""B22: opt-in Gram eta, preserved defaults, premises and fixed-signal limits.

Population checks concern fixed positive signal, not unrestricted L11 or coverage.
"""
import math

import numpy as np
import pytest

import met.full as full


GRID = (10., 1201)


def readings(seed=2211):
    rng = np.random.default_rng(seed)
    basis, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    answer = []
    for _ in range(3):
        z = rng.normal(size=(6, 6))
        covariance = z @ z.T + np.eye(6)
        data = rng.normal(size=6)*3
        answer.append(full.Reading(data[:3], data[3:], covariance, basis[:, :2]))
    return answer


def estimate(rs, positive=True, **kwargs):
    return full.ours29(rs, 1.7, positivity=positive, readout=False, grid=GRID, **kwargs)


@pytest.mark.parametrize("positive", [False, True])
def test_default_is_exactly_explicit_bartlett_and_gram_calls_the_existing_ratio(positive):
    rs = readings()
    implicit = estimate(rs, positive)
    explicit = estimate(rs, positive, eta_form="bartlett")
    for key in implicit:
        np.testing.assert_equal(implicit[key], explicit[key])
    gram = estimate(rs, positive, eta_form="gram")
    eta, q, h, phi = full.gram_eta(rs, 1.7, gram["g"])
    assert gram["eta_form"] == "gram"
    assert (gram["eta"], gram["gram_Q"], gram["gram_H"], gram["phi_hat"]) == (eta, q, h, phi)
    assert np.isnan(gram["curv"]) and np.isnan(gram["J"])
    np.testing.assert_array_equal(gram["g"], implicit["g"])
    assert gram["phi_hat"] == pytest.approx(implicit["phi_hat"], abs=1e-12)


def test_ours30_forwards_option_and_raw_bypasses_both_calibrations(monkeypatch):
    rs = [full.isotropic([2.], [1.], 1., 1.), full.isotropic([1.], [3.], 1., 1.)]
    expected = full.ours29(rs, 1., positivity=True, eta_form="gram", readout=False)
    actual = full.ours30(rs, 1., eta_form="gram", readout=False)
    for key in expected:
        np.testing.assert_equal(actual[key], expected[key])
    def forbidden(*args):
        raise AssertionError("calibration must be bypassed")
    monkeypatch.setattr(full, "bartlett_eta", forbidden)
    # The Gram route must not secretly compute a Bartlett correction.
    full.ours30(rs, 1., eta_form="gram", readout=False)
    monkeypatch.setattr(full, "gram_eta", forbidden)
    raw = [full.ours30(rs, 1., eta_form=form, calibrate=False, readout=False)
           for form in ("bartlett", "gram")]
    assert all(result["eta"] == 1. and result["eta_form"] == "none" for result in raw)
    np.testing.assert_array_equal(raw[0]["law"], raw[1]["law"])
    for calibrate in (True, False):
        with pytest.raises(ValueError, match="unknown eta form"):
            full.ours30(rs, 1., calibrate=calibrate, eta_form="invalid")


@pytest.mark.parametrize("positive", [False, True])
def test_gram_z_preserves_full_law_and_point_exactly(positive):
    rs = readings()
    zero = full.Reading(np.zeros(3), np.zeros(3), np.eye(6))
    original = estimate(rs, positive, eta_form="gram")
    padded = estimate(rs+[zero], positive, eta_form="gram")
    for key in ("law", "lm", "point", "interval", "eta", "gram_Q", "gram_H"):
        np.testing.assert_equal(original[key], padded[key])
    null = estimate([zero], positive, eta_form="gram")
    empty = estimate([], positive, eta_form="gram")
    np.testing.assert_array_equal(null["law"], empty["law"])
    assert null["point"] == pytest.approx(1.7, rel=1e-7)


@pytest.mark.parametrize("positive", [False, True])
@pytest.mark.parametrize("operation", ["units", "swap", "frame"])
def test_gram_l3_swap_and_frame_transform_the_full_declared_inputs(positive, operation):
    rs, c = readings(), 1.7
    if operation == "units":
        scale = 3.7/.4
        transform = np.diag([3.7]*3+[.4]*3)
        other_c, rotate = c*scale, np.eye(3)
    elif operation == "swap":
        scale = 1.
        transform = np.block([[np.zeros((3,3)), np.eye(3)], [np.eye(3), np.zeros((3,3))]])
        other_c, rotate = 1/c, np.eye(3)
    else:
        scale = 1.
        rotate, _ = np.linalg.qr(np.random.default_rng(2212).normal(size=(3,3)))
        transform = np.kron(np.eye(2), rotate)
        other_c = c
    transformed = []
    for r in rs:
        y = transform @ np.r_[r.force, r.acceleration]
        transformed.append(full.Reading(y[:3], y[3:], transform @ r.cov @ transform.T,
                                        rotate @ r.subspace))
    before = estimate(rs, positive, eta_form="gram")
    after = full.ours29(transformed, other_c, positivity=positive, eta_form="gram",
                        readout=False, grid=GRID)
    assert after["eta"] == pytest.approx(before["eta"], abs=2e-12)
    np.testing.assert_allclose(after["g"], before["g"], rtol=2e-11, atol=2e-12)
    if operation == "swap":
        np.testing.assert_allclose(after["law"], before["law"][::-1], rtol=1e-10, atol=1e-13)
        assert after["point"] == pytest.approx(1/before["point"], rel=1e-7)
        np.testing.assert_allclose(after["interval"], 1/np.array(before["interval"])[::-1], rtol=1e-7)
    else:
        np.testing.assert_allclose(after["law"], before["law"], rtol=1e-10, atol=1e-13)
        assert after["point"] == pytest.approx(scale*before["point"], rel=1e-10)


def evidence_span(f, a, sd, positive):
    rs = [full.isotropic(f, a, sd, sd)]
    result = full.ours29(rs, 1., positivity=positive, eta_form="gram", readout=False, grid=GRID)
    phi = np.linspace(0., math.pi/2, 1001)
    h = full._h_terms(rs[0].pre(), 1., phi)
    span = result["eta"]*result["g"][0]*np.ptp(h)/2
    return span, math.log(result["interval"][1]/result["interval"][0])


@pytest.mark.parametrize("f,a", [
    ([2.,0.], [math.cos(math.pi/6),math.sin(math.pi/6)]),
    ([math.sqrt(3),1.,0.], [1.,0.,0.]),
    ([-2.], [1.]), ([-2.,0.,0.], [1.,0.,0.]),
])
def test_ours30_gram_t30_keeps_a_finite_nonvanishing_incompatible_pair_contrast(f,a):
    values = np.array([evidence_span(f,a,sd,True) for sd in (.03,.01,.003)])
    np.testing.assert_allclose(values, np.broadcast_to(values[0], values.shape), rtol=1e-8, atol=1e-8)
    assert values[0,0] > .1 and values[0,1] > .05


def test_sign_blind_gram_base_retains_t30_antiparallel_failure():
    first = evidence_span([-2.], [1.], .1, False)[0]
    second = evidence_span([-2.], [1.], .01, False)[0]
    assert second/first == pytest.approx(100., rel=1e-10)


def test_full_type_can_return_negative_gram_eta_without_an_unrequested_floor():
    # Three different correlated instruments: the common-instrument eigenvalue
    # interpretation no longer bounds Q/H by one. This is not a search artifact.
    rs = []
    for theta in (.2,.75,1.3):
        u = np.array([math.sin(theta),math.cos(theta)])
        n = np.array([math.cos(theta),-math.sin(theta)])
        rs.append(full.Reading([u[0]],[u[1]],np.outer(u,u)+.01*np.outer(n,n)))
    result = full.ours30(rs,1.,eta_form="gram",readout=False)
    eta,q,h,phi = full.gram_eta(rs,1.,result["g"])
    np.testing.assert_allclose(result["g"],1.,atol=0.)
    assert result["eta"] == eta == pytest.approx(-.8521796125593,abs=2e-11)
    assert h == pytest.approx(1.0518271664203,abs=2e-11)
    assert q == pytest.approx(3-h,abs=2e-12)
    dense = sum(full._h_terms(r.pre(),1.,np.linspace(-math.pi/2,math.pi/2,20001)) for r in rs)
    assert dense.max() == pytest.approx(h,abs=5e-7)
    lm = np.log(np.tan([phi,.01]))
    evidence = full._log_law(rs,1.,result["g"],eta,lm)-full._log_prior(lm,1.)
    assert evidence[0] < evidence[1]  # negative eta reverses the raw-evidence ordering
    assert np.all(np.isfinite(result["law"])) and np.all(result["law"] >= 0)
    assert np.sum(result["law"]) == pytest.approx(1.)
    assert 0 < result["interval"][0] < result["point"] < result["interval"][1] < np.inf


@pytest.mark.parametrize("k,signal", [(1,1.),(1,2.),(3,1.),(3,2.)])
def test_fixed_signal_population_optima_and_positive_eta_limits(k,signal):
    # Same Gaussian cells as test_weight_consistency, truth away from prior centre.
    # Positive scalar eta changes precision, never the quadratic optimum. The
    # positive empirical limits check the fixed-signal argument, not full L11.
    rng = np.random.default_rng(2210)
    mean = signal*np.array([2.,1.])/math.sqrt(5)
    x, y = rng.normal(size=(400_000,k)), rng.normal(size=(400_000,k))
    x[:,0] += mean[0]; y[:,0] += mean[1]
    a, b, cross = (x*x).sum(1), (y*y).sum(1), (x*y).sum(1)
    q = (a+b)/2-np.hypot((a-b)/2,cross)
    for positive in (False, True):
        misfit = np.where(cross>0,q,np.minimum(a,b)) if positive else q
        count = k-.5 if positive else k-1
        g = np.ones_like(q) if count == 0 else np.minimum(1., count/np.maximum(misfit,1e-300))
        gram = np.array([[(g*a).mean(),(g*cross).mean()],[(g*cross).mean(),(g*b).mean()]])
        eigenvalues, vectors = np.linalg.eigh(gram)
        u = vectors[:,1]
        mass = u[0]/u[1]
        gram_eta = 1-eigenvalues[0]/eigenvalues[1]
        scores = g*(u[0]*u[1]*(a-b)+(u[1]**2-u[0]**2)*cross)
        bartlett = min(1.,np.diff(eigenvalues)[0]/np.mean(scores*scores))
        assert 0 < gram_eta <= 1 and 0 < bartlett <= 1
        if positive:
            assert mass < (1.87 if signal == 1 else 1.99)
        else:
            assert mass == pytest.approx(2.,abs=.05)
