"""Independent mathematical checks; does not modify the production estimator.

Uses only NumPy and the standard library. Run from any directory.
"""
from pathlib import Path
import json
import math
import numpy as np


def trigamma(x):
    answer = 0.0
    while x < 30:
        answer += 1 / x**2
        x += 1
    return answer + 1/x + 1/(2*x*x) + 1/(6*x**3) - 1/(30*x**5) + 1/(42*x**7) - 1/(30*x**9)


def integrate_grid(values, step):
    return step/3 * (values[0] + values[-1] + 4*values[1:-1:2].sum() + 2*values[2:-1:2].sum())


def power_log_density(ell, k):
    nu = k + 1
    log_cosh = np.logaddexp(ell, -ell) - math.log(2)
    log_beta = 2*math.lgamma(nu/2) - math.lgamma(nu)
    return np.exp(-k*math.log(2) - log_beta - nu*log_cosh)


def main():
    checks = []
    ell = np.linspace(-100, 100, 100001)
    step = ell[1] - ell[0]
    for k in (-.5, 0., 1., 2., 9., 99.):
        density = power_log_density(ell, k)
        normalization = integrate_grid(density, step)
        mean = integrate_grid(ell*density, step)
        variance = integrate_grid(ell**2*density, step)
        exact_variance = .5 * trigamma((k+1)/2)
        assert abs(normalization - 1) < 1e-10
        assert abs(mean) < 1e-10
        assert abs(variance/exact_variance - 1) < 1e-9
        cdf = np.r_[0., np.cumsum((density[1:] + density[:-1]) * step / 2)]
        lower, upper = np.exp(np.interp([.025, .975], cdf, ell))
        checks.append(dict(k=k, point=1., log_variance=exact_variance,
                           central_95=[float(lower), float(upper)],
                           normalization_error=float(normalization-1),
                           variance_error=float(variance-exact_variance)))
    # These two CDFs come from direct integrals in mass, independent of the
    # gamma/trigamma formulas and the integration on log mass above.
    assert np.allclose(checks[1]['central_95'],
                       [math.tan(math.pi*.025/2), math.tan(math.pi*.975/2)], rtol=2e-6)
    assert np.allclose(checks[2]['central_95'],
                       [math.sqrt(.025/.975), math.sqrt(.975/.025)], rtol=2e-6)

    # Direct f/alpha integration for finite, nonzero readings under a known
    # orientation and independent Gaussian channels. This checks the universal
    # high-mass asymptotic and its normalization without using the polar code.
    z, weights = np.polynomial.legendre.leggauss(256)
    f, weights = 10*(z+1), 10*weights  # [0,20], beyond these moderate Gaussian data.
    force, accel, sf, sa = 3., 2., 1., .5
    zf = sf*math.sqrt(2*math.pi)*.5*math.erfc(-force/(math.sqrt(2)*sf))
    za = sa*math.sqrt(2*math.pi)*.5*math.erfc(-accel/(math.sqrt(2)*sa))
    normalization = zf*za
    boundary_integral = float(weights @ (f*np.exp(-.5*((f-force)/sf)**2 - .5*(accel/sa)**2)))
    tail_constant = boundary_integral/normalization
    masses = [1e2, 1e3, 1e4, 1e5, 1e6]
    scaled_densities = []
    for mass in masses:
        scaled_density = float(weights @ (f*np.exp(-.5*((f-force)/sf)**2 - .5*((f/mass-accel)/sa)**2))) / normalization
        scaled_densities.append(scaled_density)
    assert abs(scaled_densities[-1]/tail_constant-1) < 3e-5

    # Null acceleration-size-biased law: q_R(r)=(1+r^2)^(-3/2).
    theta = np.linspace(0, math.pi/2, 100001)
    dtheta = theta[1]-theta[0]
    assert abs(integrate_grid(np.cos(theta), dtheta)-1) < 1e-12
    assert abs(integrate_grid(np.sin(theta), dtheta)-1) < 1e-12
    # Truncated second moment is asinh(T)-T/sqrt(1+T^2), asymptotic log(2T)-1.
    truncation = [10., 100., 1e3, 1e6]
    second = [math.asinh(t)-t/math.sqrt(1+t*t) for t in truncation]
    assert abs(second[-1]-(math.log(2*truncation[-1])-1)) < 1e-10
    q_interval = [p/math.sqrt(1-p*p) for p in (.025,.975)]
    result = dict(power_null=checks,
                  finite_reading_tail=dict(readings=[force,accel], sigmas=[sf,sa],
                                           constant=tail_constant,masses=masses,
                                           mass_squared_times_density=scaled_densities),
                  size_biased_null=dict(mean=1.,median=1/math.sqrt(3),central_95=q_interval,
                                        second_moment_truncations=truncation,
                                        truncated_second_moments=second),
                  checked="normalization, log moments, two exact CDFs, non-null tail, size-biased mean and divergent second moment")
    out = Path(__file__).with_name('verification.json')
    out.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
