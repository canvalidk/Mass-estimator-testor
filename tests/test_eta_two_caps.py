"""Eta's two jobs as an exact product, and a cap on each (Claude, 17 Nov 2026; problems/O14.md, "B22 read").

On one instrument the Bartlett eta is the Gram eta times expected-over-observed squared scores:
gap/J = (gap/H_max) * (H_max/J). The first factor is the weak-signal share; the second is the discount for
readings that disagree, and it exceeds 1 whenever they agree better than the declared noise predicts.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from met.full import Reading, ours29

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import eta_two_caps  # noqa: E402
import ours30_eta_gram_coverage as b22  # noqa: E402


def _scalar_series(seed, n, signal):
    rng = np.random.default_rng(seed)
    u = np.array([math.sin(math.pi / 4), math.cos(math.pi / 4)])
    rows = signal * u + rng.normal(size=(n, 2))
    return [Reading(np.array([y[0]]), np.array([y[1]]), np.eye(2)) for y in rows]


@pytest.mark.parametrize("seed,n,signal", [(22, 6, 1.2), (23, 12, 0.4), (24, 4, 3.0)])
def test_the_bartlett_eta_is_the_gram_eta_times_expected_over_observed_squared_scores(seed, n, signal):
    rd = _scalar_series(seed, n, signal)
    b, g = ours29(rd, 1.0, misfit=False), ours29(rd, 1.0, misfit=False, eta_form="gram")
    gap = g["gram_H"] - g["gram_Q"]
    assert b["curv"] == pytest.approx(gap, rel=1e-9)                       # the fitted curvature is the eigengap
    eta_gram = gap / g["gram_H"]
    assert g["eta"] == pytest.approx(eta_gram, rel=1e-9)
    assert b["curv"] / b["J"] == pytest.approx(eta_gram * g["gram_H"] / b["J"], rel=1e-9)
    assert b["eta"] == pytest.approx(min(1.0, eta_gram * g["gram_H"] / b["J"]), rel=1e-9)


def test_the_joint_cap_lets_good_agreement_narrow_the_width_in_most_null_series():
    # 20 scalar readings of pure noise (B22's cell): in most series the readings agree better at the fitted
    # line than the declared noise predicts, the second factor exceeds 1, and the Bartlett eta sits above
    # the Gram eta. Capping each factor at 1 gives the smaller of the two etas.
    res = eta_two_caps.cell(1, 20, 0.0, series=300)["ours29"]
    assert res["bartlett_above_gram"] > 0.75
    assert res["bartlett"] == pytest.approx(0.843, abs=0.005) and res["gram"] == pytest.approx(0.943, abs=0.005)
    assert res["two_caps"] >= max(res["bartlett"], res["gram"])


def test_two_caps_keeps_the_discount_the_gram_eta_drops():
    # Can's far pair in B22's summaries: readings of one mass at 5 and 6 (sigma_F = 0.1, c = s)
    force = np.array([[50.0], [60.0]])            # F / sigma_F, sigma_F = 0.1
    acceleration = np.array([[500.0], [500.0]])   # a / sigma_a, a = 1 and sigma_a = 0.002
    _, _, etas = b22.summaries(force, acceleration, positivity=False)
    assert etas["gram"] > 0.99 and etas["bartlett"] < 0.2
    assert min(etas["gram"], etas["bartlett"]) == etas["bartlett"]
