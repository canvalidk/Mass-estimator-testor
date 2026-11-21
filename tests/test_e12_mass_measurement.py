"""Finite consequences of persistent mass readouts and the exchange-phase bound."""
from fractions import Fraction

import numpy as np
import pytest


IDENTITY = np.eye(2)
EXCHANGE = np.array([[0., 1.], [1., 0.]])
PLUS = np.array([1., 1.]) / np.sqrt(2)
MINUS = np.array([1., -1.]) / np.sqrt(2)


def reciprocal_standard(pair, load):
    x, y = pair
    return 1 / (x + load) + 1 / y


def permitted_transitions(tags):
    return np.array([[a == b for b in tags] for a in tags])


def test_classical_protected_values_allow_only_transitions_inside_joint_level_sets():
    pairs = [(Fraction(x), Fraction(y)) for x in (1, 2, 3) for y in (1, 2, 3)]
    symmetric = [(x+y, x*y/(x+y)) for x, y in pairs]
    moves = permitted_transitions(symmetric)
    swap = np.array([[p == q[::-1] for q in pairs] for p in pairs])
    np.testing.assert_array_equal(moves, np.eye(len(pairs), dtype=bool) | swap)
    # A nontrivial stochastic readout can randomly exchange the names while retaining both tags.
    kernel = .3*np.eye(len(pairs)) + .7*swap
    np.testing.assert_allclose(kernel.sum(axis=1), 1)
    assert not np.any(kernel[~moves])
    i, j = pairs.index((1, 3)), pairs.index((3, 1))
    assert kernel[i, j] == pytest.approx(.7)
    assert reciprocal_standard(pairs[i], 1) != reciprocal_standard(pairs[j], 1)

    labelled = [(reciprocal_standard(pair, 1), reciprocal_standard(pair, 2)) for pair in pairs]
    np.testing.assert_array_equal(permitted_transitions(labelled), np.eye(len(pairs), dtype=bool))


def outcome_law(outcomes, density):
    return np.array([sum(np.trace(k @ density @ k.conj().T).real for k in kraus)
                     for kraus in outcomes])


def after_readout(outcomes, density):
    return sum(k @ density @ k.conj().T for kraus in outcomes for k in kraus)


def tag_errors(outcomes):
    state0, state1 = np.diag([1., 0.]), np.diag([0., 1.])
    return after_readout(outcomes, state0)[1, 1].real, after_readout(outcomes, state1)[0, 0].real


def contrast(outcomes):
    law_plus = outcome_law(outcomes, np.outer(PLUS, PLUS))
    law_minus = outcome_law(outcomes, np.outer(MINUS, MINUS))
    np.testing.assert_allclose([law_plus.sum(), law_minus.sum()], 1., atol=1e-12)
    return np.abs(law_plus - law_minus).sum()/2


def tag_bound(e0, e1):
    return np.sqrt((1-e0)*e1) + np.sqrt(e0*(1-e1))


def test_exchange_phase_is_detectable_while_protecting_only_the_symmetric_composites():
    outcomes = [[(IDENTITY + EXCHANGE)/2], [(IDENTITY - EXCHANGE)/2]]
    np.testing.assert_allclose(outcome_law(outcomes, np.outer(PLUS, PLUS)), [1, 0], atol=1e-15)
    np.testing.assert_allclose(outcome_law(outcomes, np.outer(MINUS, MINUS)), [0, 1], atol=1e-15)
    np.testing.assert_allclose(outcome_law(outcomes, IDENTITY/2), [.5, .5])
    pairs = [(Fraction(1), Fraction(3)), (Fraction(3), Fraction(1))]
    totals = [x+y for x, y in pairs]
    reduced = [x*y/(x+y) for x, y in pairs]
    loaded = [reciprocal_standard(pair, 1) for pair in pairs]
    assert totals == [4, 4] and reduced == [Fraction(3, 4)]*2
    assert loaded == [Fraction(5, 6), Fraction(5, 4)]
    state0 = np.diag([1., 0.])
    np.testing.assert_allclose(after_readout(outcomes, state0), IDENTITY/2)
    np.testing.assert_allclose(tag_errors(outcomes), [.5, .5])


@pytest.mark.parametrize("v", [0., .05, .6, .9, 1.])
def test_equal_error_square_root_instrument_attains_visibility_bound(v):
    # Square roots of (I ± v X)/2 via the two exchange eigenspaces.
    left = np.sqrt((1+v)/2)
    right = np.sqrt((1-v)/2)
    a, b = (left+right)/2, (left-right)/2
    outcomes = [[a*IDENTITY + b*EXCHANGE], [a*IDENTITY - b*EXCHANGE]]
    errors = tag_errors(outcomes)
    np.testing.assert_allclose(errors, (1-np.sqrt(1-v**2))/2, atol=1e-15)
    assert contrast(outcomes) == pytest.approx(v)
    assert contrast(outcomes) == pytest.approx(tag_bound(*errors))


@pytest.mark.parametrize("e0,e1", [(0., .3), (.03, .41), (.49, .07), (.5, .5)])
def test_unequal_tag_error_bound_is_sharp(e0, e1):
    outcomes = [[np.array([[np.sqrt(1-e0), s*np.sqrt(e1)],
                           [np.sqrt(e0), s*np.sqrt(1-e1)]])/np.sqrt(2)] for s in [1, -1]]
    complete = sum(k.T @ k for kraus in outcomes for k in kraus)
    np.testing.assert_allclose(complete, IDENTITY, atol=1e-15)
    np.testing.assert_allclose(tag_errors(outcomes), [e0, e1], atol=1e-15)
    assert contrast(outcomes) == pytest.approx(tag_bound(e0, e1))


@pytest.mark.parametrize("seed", [13, 21, 71])
def test_complex_unresolved_kraus_outcomes_obey_tag_bound(seed):
    rng = np.random.default_rng(seed)
    random = rng.normal(size=(12, 2)) + 1j*rng.normal(size=(12, 2))
    isometry, _ = np.linalg.qr(random)
    operators = np.sqrt(.4)*isometry.reshape(6, 2, 2)
    outcomes = [[np.sqrt(.6)*IDENTITY, *operators[:2]], list(operators[2:5]), [operators[5]]]
    complete = sum(k.conj().T @ k for kraus in outcomes for k in kraus)
    np.testing.assert_allclose(complete, IDENTITY, atol=1e-14)
    e0, e1 = tag_errors(outcomes)
    assert 0 < e0 < .5 and 0 < e1 < .5
    assert 0.001 < contrast(outcomes) <= tag_bound(e0, e1) + 1e-14
    epsilon = max(e0, e1)
    assert tag_bound(e0, e1) <= 2*np.sqrt(epsilon*(1-epsilon)) + 1e-14
