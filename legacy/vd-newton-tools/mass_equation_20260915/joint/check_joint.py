"""Checks for the common-mass Fisher determinant and null pooling identities."""
import math
import numpy as np


def density(m, vectors, sf, sa):
    c = m*m/sf**2 + 1/sa**2
    d = sf**2 + m*m*sa**2
    schur = np.sum(np.sum(vectors*vectors, axis=1)/d)
    return math.sqrt(schur)*np.prod(c**(vectors.shape[1]/2))


def check_fisher():
    rng = np.random.default_rng(9152026)
    worst = 0.
    for n in (1, 2, 4):
        for d in (1, 2, 3):
            for _ in range(50):
                m = math.exp(rng.uniform(-1, 1))
                vectors = rng.normal(size=(n, d))
                sf, sa = np.exp(rng.uniform(-1, 1, size=(2, n)))
                fisher = np.zeros((1+n*d, 1+n*d))
                for i, v in enumerate(vectors):
                    j = np.zeros((2*d, 1+n*d))
                    j[:d, 0] = v
                    j[:d, 1+i*d:1+(i+1)*d] = m*np.eye(d)
                    j[d:, 1+i*d:1+(i+1)*d] = np.eye(d)
                    g = np.diag([sf[i]**-2]*d+[sa[i]**-2]*d)
                    fisher += j.T@g@j
                direct = math.exp(np.linalg.slogdet(fisher)[1]/2)
                analytic = density(m, vectors, sf, sa)
                error = abs(direct/analytic-1)
                worst = max(worst, error)
                assert error < 2e-12
                reciprocal = density(1/m, m*vectors, sa, sf)*m**(n*d-2)
                assert abs(reciprocal/analytic-1) < 2e-12
    return worst


def trigamma_integer_or_half(x):
    # N/2 only. Recurrence starts at the exactly known integer or half-integer value.
    base = .5 if x % 1 else 1.
    result = math.pi**2/2 if base == .5 else math.pi**2/6
    while base < x:
        result -= 1/base**2
        base += 1
    return result


def check_log_pooling():
    nodes, weights = np.polynomial.legendre.leggauss(384)
    ell, w = (nodes+1)*20, weights*20
    results = {}
    for n in (1, 2, 4, 16):
        constant = math.gamma((n+1)/2)/(math.sqrt(math.pi)*math.gamma(n/2))
        pdf = constant/np.cosh(ell)**n
        normalization = 2*np.sum(w*pdf)
        variance = 2*np.sum(w*pdf*ell**2)
        expected = .5*trigamma_integer_or_half(n/2)
        assert abs(normalization-1) < 2e-11
        assert abs(variance-expected) < 2e-11
        results[n] = float(variance)
    return results


def check_noncentral_chi_integral():
    nodes, weights = np.polynomial.legendre.leggauss(192)
    beta, w = (nodes+1)*math.pi/4, weights*math.pi/4
    sine2 = np.sin(beta)**2
    cosine_log = np.log(np.cos(beta))
    worst = 0.
    for k in (1, 3, 12, 128):
        for t in (0., .5, 3., 12., 64.):
            integral = math.sqrt(2/math.pi)*np.sum(
                w*(-np.expm1(k*cosine_log-t*t*sine2/2))/sine2)
            expected = None
            if t == 0:
                expected = math.sqrt(2)*math.exp(math.lgamma((k+1)/2)-math.lgamma(k/2))
            elif k == 1:
                expected = math.sqrt(2/math.pi)*math.exp(-t*t/2)+t*math.erf(t/math.sqrt(2))
            elif k == 3:
                expected = (math.sqrt(2/math.pi)*math.exp(-t*t/2)
                            +(t+1/t)*math.erf(t/math.sqrt(2)))
            if expected is not None:
                error = abs(integral/expected-1)
                worst = max(worst, error)
                assert error < 3e-11
    return worst


if __name__ == '__main__':
    print('Worst relative determinant discrepancy:', check_fisher())
    print('Variances from multiplying null log densities:', check_log_pooling())
    print('Joint Jeffreys null log variance remains:', math.pi**2/4)
    print('Worst noncentral chi integral relative discrepancy:', check_noncentral_chi_integral())
    print('All assertions passed.')
