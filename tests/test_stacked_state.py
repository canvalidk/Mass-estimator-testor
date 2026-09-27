"""ours25b's exact carried state (C, T, M) and the combination equation.

The state is three running sums per series; combining sets of readings is adding
their states. Checked against the rule applied to all readings (the whole law and the
ratio_of_means_rss readout), against its own algebra (order and grouping), and against
the per-reading likelihood it replaces.
"""

import numpy as np
import pytest

from met.analyst import Estimator, combined_stacked_state, posterior_summaries
from met.law import ReadingSet, StackedState

CART2 = {"reference": "cartesian2", "nuisance": "radius"}
CART1 = {"reference": "cartesian1", "nuisance": "radius"}
OLD = {"rule": "likelihood_product", "prior": [{"uniform_angle": {"center": "first_reading"}}]}
SEQ = {"rule": "sequential", "carry": "stacked_state", "old": OLD}
READOUTS = ["ratio_of_means_rss", "median", "log_sd", "interval_95"]


def _readings(seed, n, d, b=3, sf=1.0, sa=0.6, scale=1.0):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(b, n, d))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    return ReadingSet(2.0 * scale * v + sf * rng.normal(size=(b, n, d)), scale * v + sa * rng.normal(size=(b, n, d)),
                      np.full((b, n), sf), np.full((b, n), sa))


@pytest.mark.parametrize("d", [1, 2, 3])
@pytest.mark.parametrize("n", [2, 6, 21])
@pytest.mark.parametrize("law", [CART2, CART1])
def test_carry_equals_the_rule_on_all_readings(d, n, law):
    r = _readings(10 * d + n, n, d)
    full = posterior_summaries(r, law, OLD, READOUTS)
    carried = posterior_summaries(r, law, SEQ, READOUTS)
    for key in READOUTS:
        assert np.allclose(carried[key], full[key], rtol=1e-12), key


@pytest.mark.parametrize("d", [1, 3])
def test_state_is_order_and_grouping_free(d):
    r = _readings(5, 7, d)
    whole = StackedState.from_readings(r)
    parts = [StackedState.from_readings(r.take_readings(slice(i, i + 1))) for i in range(7)]
    left = parts[0]
    for p in parts[1:]:
        left = left.add(p)
    right = parts[-1]
    for p in reversed(parts[:-1]):
        right = p.add(right)
    grouped = StackedState.from_readings(r.take_readings(slice(0, 3))).add(
        StackedState.from_readings(r.take_readings(slice(3, 7))))
    for st in (left, right, grouped):
        assert np.allclose(st.C, whole.C, rtol=1e-13) and np.allclose(st.T, whole.T, rtol=1e-13)
        assert st.M == whole.M == 7 * d


@pytest.mark.parametrize("d", [1, 2, 3])
def test_state_reproduces_the_per_reading_likelihood_and_weight(d):
    r = _readings(6, 5, d)
    u = np.log(r.s[:, :1]) + np.linspace(-3, 3, 41)[None, :]
    like, _, _ = r.reading_terms(u, "radius", reference="cartesian2")
    st = StackedState.from_readings(r)
    total = np.sum(like, axis=1)
    diff = st.log_likelihood(u) - total
    assert np.ptp(diff, axis=1) == pytest.approx(0.0, abs=1e-10)          # equal up to a constant per series
    assert np.allclose(st.cond_alpha(u), r.stacked_cond_alpha(u, reference="cartesian2"), rtol=1e-12)


def test_combined_state_is_old_plus_new():
    r = _readings(7, 6, 3)
    assert np.allclose(combined_stacked_state(r).C, StackedState.from_readings(r).C, rtol=1e-13)


def test_states_with_different_instruments_do_not_add():
    a = StackedState.from_readings(_readings(8, 2, 3, sf=1.0, sa=0.6))
    b = StackedState.from_readings(_readings(9, 2, 3, sf=1.0, sa=0.3))
    with pytest.raises(ValueError, match="one noise ratio"):
        a.add(b)


def _spec(**changes):
    spec = {"name": "e", "reduce": "none", "law": CART2, "combine": SEQ, "readouts": ["ratio_of_means_rss", "median"]}
    spec.update(changes)
    return spec


def test_refusals():
    with pytest.raises(ValueError, match="cartesian1 or cartesian2"):
        Estimator(_spec(law={"reference": "flat", "nuisance": "radius"}))
    with pytest.raises(ValueError, match="cartesian1 or cartesian2"):
        Estimator(_spec(combine=dict(SEQ, old={"rule": "posterior_product", "coordinate": "mass", "prior": []})))
    with pytest.raises(ValueError, match="does not carry eq. 28"):
        Estimator(_spec(readouts=["ratio_of_means", "median"]))
    Estimator(_spec())   # the valid form is accepted
