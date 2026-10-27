"""Where the cloud readout's Z failure lives (O2, 27 Oct, Claude).

Per candidate X = m Y, so T10's readout is the law tilted by A(m) = E[Y | m], and only A's
m-dependence matters. Two facts: the series join of T9 (L12) is second order, and null readings
drive A toward a shape set by the instrument alone, (1 + m^2/s^2)^(-1/2).
"""
import numpy as np
import pytest

from met.full import isotropic, ours29

SF, SA = 0.2, 0.1
S = SF / SA
F0 = [[4.0, 0.3, 0.0], [6.1, -0.2, 0.1], [1.0, 0.8, 0.0], [3.0, 0.1, -0.2]]
A0 = [[2.0, 0.0, 0.1], [2.0, 0.1, 0.0], [0.6, -0.4, 0.0], [1.0, 0.0, 0.0]]


def _series(scale, nulls=0):
    sig = [isotropic([scale * v for v in f], [scale * v for v in x], SF, SA) for f, x in zip(F0, A0)]
    return sig + [isotropic([0.0] * 3, [0.0] * 3, SF, SA) for _ in range(nulls)]


def _join_readout(result, scale, p, draws=800, seed=0):
    """T10 with the series size joined as an l_p norm, by Monte Carlo on a_i* | m."""
    m = np.exp(result["lm"]); law = result["law"]
    v = 1 / (1 / SA**2 + m**2 / SF**2)
    z = np.random.default_rng(seed).standard_normal((len(F0), draws, 3))
    sizes = []
    for i, (f, x) in enumerate(zip(F0, A0)):
        mu = (scale * np.asarray(x) / SA**2 + m[:, None] * scale * np.asarray(f) / SF**2) * v[:, None]
        sizes.append(np.linalg.norm(mu[:, None, :] + np.sqrt(v)[:, None, None] * z[i][None], axis=2))
    sizes = np.stack(sizes, axis=1)
    y = sizes.max(axis=1) if p == np.inf else (sizes ** p).sum(axis=1) ** (1 / p)
    weight = y.mean(axis=1)
    return float((m * weight * law).sum() / (weight * law).sum())


def test_series_join_is_second_order():
    """At weak signal (scale 0.05), l_1, l_2 and l_inf give readouts within 0.2%; l_2 is T9."""
    result = ours29(_series(0.05), S)
    readouts = {p: _join_readout(result, 0.05, p) for p in (1, 2, np.inf)}
    assert readouts[2] == pytest.approx(result["readout"], rel=2e-3)
    assert max(readouts.values()) / min(readouts.values()) - 1 < 2e-3


@pytest.mark.parametrize("scale, limit", [(0.05, 2.3930), (0.02, 2.0826)])
def test_null_readings_drive_the_weight_to_the_instrument_shape(scale, limit):
    """200 null readings (1000 under MET_SLOW, to 3e-4) leave the law unchanged (T24) and replace the readout's weight by
    (1 + m^2/s^2)^(-1/2), a function of the declared noise alone."""
    base = ours29(_series(scale), S)
    m = np.exp(base["lm"]); law = base["law"]
    w = (1 + m**2 / S**2) ** -0.5
    tilted = float((m * w * law).sum() / (w * law).sum())
    assert tilted == pytest.approx(limit, abs=2e-4)
    import os
    nulls, tol = (1000, 3e-4) if os.environ.get("MET_SLOW") else (200, 1e-3)
    padded = ours29(_series(scale, nulls=nulls), S)
    assert padded["readout"] == pytest.approx(tilted, abs=tol)
    assert abs(padded["readout"] - S) < abs(base["readout"] - S)      # pulled toward s
