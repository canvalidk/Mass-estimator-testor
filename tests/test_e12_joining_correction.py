"""Finite consequences of associativity for reciprocal disagreement factors."""
import math

import numpy as np
import pytest

from explore.e12_tangent_pairings import tangent_pairings
from met.full import isotropic


def state(weight, displacement):
    return np.array([weight, weight * displacement])


def discrepancy(a, b):
    ha, sa = a
    hb, sb = b
    return (ha * hb / (ha + hb)) * (sa / ha - sb / hb) ** 2


def tree_increments(a, b, c):
    return (discrepancy(a, b), discrepancy(a + b, c),
            discrepancy(b, c), discrepancy(a, b + c))


@pytest.mark.parametrize("p,q", [(.2, 3.), (1., 1.), (3., .4), (11., .02)])
def test_three_singletons_realize_an_arbitrary_split(p, q):
    t = math.sqrt((p+q)/q) - 1.
    change = math.sqrt(p * (1.+t)/t)
    a, b, c = state(1., change), state(t, 0.), state(1., 0.)
    assert tree_increments(a, b, c) == pytest.approx((p, q, 0., p+q))
    # Ordinary Gaussian readings realize these tangent statistics at mass 1.
    u = np.array([1., 1.]) / math.sqrt(2)
    v = np.array([1., -1.]) / math.sqrt(2)
    row_states = []
    for weight, weighted_displacement in (a, b, c):
        row = math.sqrt(weight) * (u + weighted_displacement / weight * v)
        block = tangent_pairings([isotropic([row[0]], [row[1]], 1., 1.)], 1.)
        row_states.append([block.weight, block.score])
    assert np.array(tree_increments(*map(np.array, row_states))) == pytest.approx([p, q, 0., p+q])


def test_tree_constant_has_a_zero_anchor_inside_positive_mass():
    chosen = math.pi / 4
    u, v = np.array([1, 1])/math.sqrt(2), np.array([1, -1])/math.sqrt(2)
    rows = np.array([u + 4*v, math.sqrt(3)*(u + 5*v), math.sqrt(2)*(u + 7*v)])
    root = math.atan(-rows[1, 1]/rows[1, 0])
    assert 0. < root < chosen
    errors = []
    for step in [1e-3, 1e-5, 1e-7]:
        th = root + step
        parallel = rows @ [math.sin(th), math.cos(th)]
        transverse = rows @ [math.cos(th), -math.sin(th)]
        states = np.column_stack((parallel**2, parallel * transverse))
        left1, left2, right1, right2 = tree_increments(*states)
        errors.append(max(abs(left1 - right1), abs(left2 - right2)))
        if step == 1e-7:
            assert left1 == pytest.approx(rows[1] @ rows[1], abs=.001)
    assert errors[2] < .001 and errors[2] < errors[1] < errors[0]


def test_nonlinearity_of_new_discrepancy_alone_breaks_grouping_but_stateful_potential_works():
    a, b, c = state(1., 3.), state(2., -1.), state(1., -1.)
    p, q, r, s = tree_increments(a, b, c)
    f = math.log1p
    assert abs((f(p)+f(q)) - (f(r)+f(s))) > .2
    # Give the three starting blocks some measured internal disagreement, too.
    u, v, w = .2, 1.4, 2.
    def joined(x, y, increment):
        return f(x+y+increment) - f(x) - f(y)
    left = joined(u,v,p) + joined(u+v+p,w,q)
    right = joined(v,w,r) + joined(u,v+w+r,s)
    assert left == pytest.approx(right)
    assert left == pytest.approx(f(u+v+w+p+q) - f(u) - f(v) - f(w))


def test_balanced_conflict_is_identically_uninformative_for_all_J_potentials():
    for scale in [1., 2., 100.]:
        rows = [isotropic([scale], [0.], 1., 1.), isotropic([0.], [scale], 1., 1.)]
        for mass in [.09, .5, 1., 2.3, 9.]:
            block = tangent_pairings(rows, mass)
            assert block.between == pytest.approx(scale**2)
            assert block.common == pytest.approx(0., abs=1e-22)
            assert block.within == pytest.approx(0., abs=1e-22)


def test_the_gaussian_pressure_cannot_be_canceled_by_a_function_of_J():
    gap_ratio = 1.1
    determinant = 225.
    j0 = 9.25
    levels = j0 * gap_ratio ** np.arange(3)
    assert 9 < levels[0] < levels[2] < 225/17
    for f in [lambda x: -.5*x, np.sqrt, lambda x: np.sin(x) + .001*x*x]:
        t = 7.
        def contrast(scale, a, b):
            return scale*determinant/2*(1/a-1/b) + f(scale*a) - f(scale*b)
        difference = contrast(gap_ratio*t, *levels[:2]) - contrast(t, *levels[1:])
        assert difference == pytest.approx(t*determinant/(2*j0)*(1-1/gap_ratio)*(gap_ratio-1/gap_ratio))


@pytest.mark.parametrize("c,r", [(-1., 4.), (0., 2.), (1., 1.2), (100., 1.01)])
def test_conflicting_acute_pair_gives_unique_concentrating_fit(c, r):
    angles = np.linspace(0, math.pi / 2, 401)
    largest, smallest = (r+1)**2, (r-1)**2
    h = smallest + (largest-smallest)*np.cos(angles-math.pi/4)**2
    log_score = .5*h + c * largest * smallest/h
    middle = len(angles)//2
    assert np.all(np.diff(log_score[:middle+1]) > 0)
    assert np.all(np.diff(log_score[middle:]) < 0)
    for scale in (1., 100.):
        h_scaled = scale*h
        j_scaled = scale**2 * largest*smallest/h_scaled
        np.testing.assert_allclose(scale*log_score, .5*h_scaled + c*j_scaled)


def test_arbitrary_singleton_four_list_cancellation_has_nonconstant_contrast():
    x, y, z = np.array([[1.,1.],[2.,1.],[1.,2.]])
    increments = []
    for sin, cos in ((0.,1.), (1.,0.)):
        def observed(row):
            u = row @ [sin, cos]
            return np.array([u*u, u * (row @ [cos, -sin])])
        X, Y, Z = map(observed, (x,y,z))
        increments.append([discrepancy(X+Y,X), discrepancy(X+Z,X)])
    np.testing.assert_allclose(increments, [[1/6, 2/15], [2/15, 1/6]])


def test_four_observations_force_the_unequal_precision_coefficient():
    h, k, x, y = 2., 5., -1.25, .4
    a, b, c, d = state(h, x), state(h, y), state(k, x), state(k, y)

    def tree(cost, a, b, c, d):
        return cost(a, b) + cost(c, d) + cost(a+b, c+d)

    assert tree(discrepancy, a, b, c, d) == pytest.approx(tree(discrepancy, a, d, b, c))

    def precision_dependent_cost(u, v):
        ratio = (u[0]-v[0])/(u[0]+v[0])
        return (1. + ratio**2) * discrepancy(u, v)

    homogeneous = tree(precision_dependent_cost, a, b, c, d)
    crossed = tree(precision_dependent_cost, a, d, b, c)
    assert abs(homogeneous-crossed) > 1.


def test_continuous_flat_singletons_collapse_on_a_nonzero_axis():
    angles = (np.arange(10001) + .5) * math.pi / (2 * 10001)
    sin, cos = np.sin(angles), np.cos(angles)

    def law(force, acceleration):
        log_weights = .5 * (force*sin + acceleration*cos)**2
        weights = np.exp(log_weights - log_weights.max())
        return weights / weights.sum()

    # At each fixed precision the positive-singleton law fixes its axis limit.
    for force in (2., 8., 32.):
        axis = law(force, 0.)
        distances = [np.abs(law(force, a) - axis).sum() for a in (.1, .001, 1e-6)]
        assert distances[0] > distances[1] > distances[2]
        assert distances[2] < 2e-6

    # Along the axis shrinking the errors moves the law to the excluded endpoint.
    outside = angles < math.pi/2 - .2
    masses = [law(force, 0.)[outside].sum() for force in (2., 8., 32.)]
    assert masses[0] > .5 and masses[0] > masses[1] > masses[2]
    assert masses[2] < 1e-8
