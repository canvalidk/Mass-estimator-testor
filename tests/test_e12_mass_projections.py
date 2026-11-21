"""Finite identities in the mass projection and adjustable-standard proof."""
from collections import Counter
from fractions import Fraction
from itertools import product
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j1


def reciprocal_composite(pair, standard):
    x, y = pair
    return 1 / (x + standard) + 1 / y


def test_known_labelled_standard_distinguishes_marginally_identical_ordinary_lumps():
    pairs = [(Fraction(1), Fraction(2)), (Fraction(2), Fraction(3)), (Fraction(3), Fraction(1))]
    transposed = [(y, x) for x, y in pairs]
    for observable in (lambda p: p[0], lambda p: p[1], sum, lambda p: reciprocal_composite(p, 0)):
        assert Counter(map(observable, pairs)) == Counter(map(observable, transposed))
    assert Counter(reciprocal_composite(p, 1) for p in pairs)[1] == 1
    assert Counter(reciprocal_composite(p, 1) for p in transposed)[1] == 0


def masses_from_paired_outputs(first, second, b, c):
    x = math.sqrt((b - c)**2 / 4 + (c - b) / (first - second)) - (b + c) / 2
    return x, 1 / (first - 1 / (x + b))


def test_two_paired_outputs_identify_the_pair_but_two_marginals_allow_a_checkerboard():
    b, c = 0.0, 1.0
    marginals_a = [(2.0, 1.1), (2.1, 1.2)]
    marginals_b = [(2.0, 1.2), (2.1, 1.1)]
    joint_mass_laws = []
    for outputs in (marginals_a, marginals_b):
        atoms = [masses_from_paired_outputs(*record, b, c) for record in outputs]
        assert np.min(atoms) > 0
        for j, standard in enumerate((b, c)):
            observed = [reciprocal_composite(p, standard) for p in atoms]
            np.testing.assert_allclose(observed, [record[j] for record in outputs])
        joint_mass_laws.append(atoms)
    assert sorted(x[0] for x in marginals_a) == sorted(x[0] for x in marginals_b)
    assert sorted(x[1] for x in marginals_a) == sorted(x[1] for x in marginals_b)
    assert set(joint_mass_laws[0]).isdisjoint(joint_mass_laws[1])


def test_any_fixed_finite_linear_design_has_a_signed_null_measure():
    # Each signed two-point difference vanishes in one measurement;
    # convolution makes their product vanish in every chosen direction.
    directions = ((1, 0), (0, 1), (1, 1), (2, 5))
    signed = Counter()
    for bits in product((0, 1), repeat=len(directions)):
        x = 10 - sum(bit * b for bit, (_, b) in zip(bits, directions))
        y = 10 + sum(bit * a for bit, (a, _) in zip(bits, directions))
        signed[x, y] += (-1)**sum(bits)

    def push(a, b):
        values = Counter()
        for (x, y), weight in signed.items():
            values[a * x + b * y] += weight
        return {k: v for k, v in values.items() if v}

    assert all(min(atom) > 0 for atom in signed)
    assert all(not push(a, b) for a, b in directions)
    assert push(1, 3)


@pytest.mark.parametrize("standard,t", [(0.15, .4), (1., 2.3), (8., 5.)])
def test_sweep_factors_as_a_laplace_transform_of_the_joint_signature(standard, t):
    atoms = np.array([[0.8, 3.], [4., 1.], [1.2, 0.5]])
    weights = np.array([.2, .35, .45])
    x, z = atoms[:, 0], 1 / atoms[:, 1]
    observed = weights @ (np.exp(-t*z) * -np.expm1(-t / (standard + x)))

    def integrand(u):
        # The entire series in the proof is sqrt(t/u) J_1(2 sqrt(tu)).
        kernel = t if u == 0 else math.sqrt(t/u) * j1(2 * math.sqrt(t*u))
        joint_signature = weights @ np.exp(-u*x - t*z)
        return math.exp(-standard*u) * kernel * joint_signature

    reconstructed, error = quad(integrand, 0, np.inf, epsabs=1e-11)
    assert error < 1e-8
    assert reconstructed == pytest.approx(observed, abs=2e-10)


def test_reusing_an_unknown_mass_and_using_an_independent_copy_give_different_laws():
    # Same law for each copy, independent of the second body of known unit mass.
    masses = (1, 3)
    reused = Counter(2*x + 1 for x in masses)
    copies = Counter(x + independent_x + 1 for x, independent_x in product(masses, repeat=2))
    assert reused[5] == 0
    assert copies[5] == 2


def test_rapid_joint_variations_can_disappear_from_the_measured_sweep():
    from scipy.integrate import trapezoid

    grid = np.linspace(1, 2, 481)
    outputs = np.linspace(1, 3, 481)
    def perturbation(x, z, frequency):
        # A compact C^3 bump suffices for this finite check.
        sx = np.sin(math.pi * (x - 1))
        sz = np.sin(math.pi * (z - 1))
        cz = np.cos(math.pi * (z - 1))
        phase = frequency * (x + z)
        derivative = sz**4 * np.cos(phase) + 4*math.pi*sz**3*cz * np.sin(phase)/frequency
        return np.where((z >= 1) & (z <= 2), sx**4 * derivative, 0)

    frequency = 96
    signed = perturbation(grid[None, :], grid[:, None], frequency)
    joint_norm = trapezoid(trapezoid(np.abs(signed), grid, axis=1), grid)
    assert joint_norm > 0.08
    for standard in (0, 0.2, 1, 12, 100):
        z = outputs[:, None] - 1/(grid[None, :] + standard)
        observed = trapezoid(perturbation(grid[None, :], z, frequency), grid, axis=1)
        assert trapezoid(np.abs(observed), outputs) < 0.001 * joint_norm
