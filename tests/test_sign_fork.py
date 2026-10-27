"""The sign of F·a: what the misfit weight sees and what T30 requires (O13, 27 Oct sort).

T30 requires bounded confidence on a fixed pair that no positive mass fits, as the noise
shrinks. T32's Q is taken on the whole projective line, so it measures the k - 1
codirectionality conditions and not L1's m > 0. These tests pin both facts down.
"""
import math

import numpy as np
import pytest

import met.full as full

LM = np.linspace(-12.0, 12.0, 4001)


def _evidence_span(force, acceleration, sd, c=1.0):
    """Span in nats over m > 0 of the tempered evidence eta * sum g_i H_i / 2 (the law
    without the starting law), with Q, g and eta."""
    readings = [full.isotropic(force, acceleration, sd, sd)]
    state = full.reading_state(readings[0], c)
    g = [full.misfit_weight(state["Q"], state["nu"])]
    eta = full.bartlett_eta(readings, c, g)[0]
    evidence = full._log_law(readings, c, g, eta, LM) - full._log_prior(LM, c)
    return float(evidence.max() - evidence.min()), state["Q"], g[0]


@pytest.mark.parametrize("sd", [0.1, 0.03, 0.01])
def test_misfit_weight_cannot_tell_anti_parallel_from_aligned(sd):
    """(2,0) against (1,0) and against (-1,0): same Q, same g, same evidence span."""
    aligned = _evidence_span([2.0, 0.0], [1.0, 0.0], sd)
    anti = _evidence_span([2.0, 0.0], [-1.0, 0.0], sd)
    assert anti[0] == pytest.approx(aligned[0], rel=1e-9)
    assert anti[1] == pytest.approx(aligned[1], abs=1e-9)
    assert anti[2] == aligned[2] == 1.0


def test_t30_bounded_tilt_holds_for_a_noncollinear_pair():
    """T30 as intended: the 30-degree pair's evidence freezes at 9.0825 nats."""
    a = [math.cos(math.pi / 6), math.sin(math.pi / 6)]
    spans = [_evidence_span([2.0, 0.0], a, sd)[0] for sd in (0.1, 0.03, 0.01)]
    assert max(spans) - min(spans) < 1e-6
    assert spans[0] == pytest.approx(9.0825, abs=1e-3)


@pytest.mark.xfail(strict=True, reason="T30 requires bounded confidence on a fixed pair no positive "
                   "mass fits. The anti-parallel pair (1,0)/(-1,0) has Q = 0, g = 1, and its evidence "
                   "span grows as SD^-2: 49.9994 nats at SD 0.1, 4999.94 at SD 0.01. O13, 27 Oct sort: "
                   "the wrong premise is L8's count, which omits L1's sign condition.")
def test_t30_bounded_tilt_on_an_anti_parallel_pair():
    wide = _evidence_span([1.0, 0.0], [-1.0, 0.0], 0.1)[0]
    narrow = _evidence_span([1.0, 0.0], [-1.0, 0.0], 0.01)[0]
    assert narrow <= 2 * wide
