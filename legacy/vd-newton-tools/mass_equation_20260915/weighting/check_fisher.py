"""Supplement direct Fisher-volume algebra with coordinate and quadrature checks."""

import math
import numpy as np


def check_determinants():
    rng = np.random.default_rng(9152026)
    worst = 0.0
    for d in (2, 3, 4):
        for _ in range(100):
            q, _ = np.linalg.qr(rng.normal(size=(d, d)))
            u, e = q[:, 0], q[:, 1:]
            f, alpha = np.exp(rng.uniform(-1, 1, size=2))
            m, v = f / alpha, alpha * u
            b = rng.normal(size=(2*d, 2*d))
            g = b.T @ b + np.eye(2*d)

            j_pair = np.zeros((2*d, d+1))
            j_pair[:d, 0] = u
            j_pair[d:, 1] = u
            j_pair[:d, 2:] = f * e
            j_pair[d:, 2:] = alpha * e
            j_mass = np.zeros((2*d, d+1))
            j_mass[:d, 0] = v
            j_mass[:d, 1:] = m * np.eye(d)
            j_mass[d:, 1:] = np.eye(d)
            old = np.sqrt(np.linalg.det(j_pair.T @ g @ j_pair))
            new = np.sqrt(np.linalg.det(j_mass.T @ g @ j_mass))
            expected = old * alpha**(2-d)
            error = abs(new/expected - 1)
            worst = max(worst, error)
            assert error < 5e-12

            sf, sa = np.exp(rng.uniform(-1, 1, size=2))
            g_iso = np.diag([sf**-2]*d + [sa**-2]*d)
            observed = np.sqrt(np.linalg.det(j_pair.T @ g_iso @ j_pair))
            expected = (f*f/sf**2 + alpha*alpha/sa**2)**((d-1)/2)/(sf*sa)
            error = abs(observed/expected - 1)
            worst = max(worst, error)
            assert error < 5e-12
    return worst


def check_sequential_example():
    nodes, weights = np.polynomial.legendre.leggauss(128)
    radius = (nodes+1)*7
    wr = weights*7
    theta = (nodes+1)*math.pi/4
    wt = weights*math.pi/4
    # Standardized polar coordinates for combined precisions A=10, B=2.
    f = radius[:, None]*np.sin(theta)/math.sqrt(10)
    alpha = radius[:, None]*np.cos(theta)/math.sqrt(2)
    base = np.exp(-radius[:, None]**2/2)*radius[:, None]*wr[:, None]*wt
    results = {}
    for label, a, b in (('fixed_first_prior', 1., 1.), ('recomputed_prior', 10., 2.)):
        weighted = base*(a*f*f+b*alpha*alpha)
        estimate = (weighted*f).sum()/(weighted*alpha).sum()
        expected = math.sqrt(2/10)*(2*a/10+b/2)/(a/10+2*b/2)
        assert abs(estimate/expected-1) < 5e-13
        results[label] = float(estimate)
    assert abs(results['fixed_first_prior']/results['recomputed_prior']-7/11) < 5e-13
    return results


if __name__ == '__main__':
    print('Worst relative determinant discrepancy:', check_determinants())
    print('Sequential design example:', check_sequential_example())
    print('All assertions passed.')
