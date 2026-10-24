"""Continuous posterior summaries for the isotropic three-dimensional study.

The two channel standard deviations are one.  ``flat_joint`` uses
df d(alpha) dOmega; ``tube_joint`` uses (f**2 + alpha**2) times that measure.
Radial and direction integrals are analytic.  Only log mass is integrated
numerically, using piecewise Gauss--Legendre quadrature and polynomial
antiderivatives.  Quantiles and CRPS describe a continuous distribution,
not a distribution with atoms at quadrature nodes.

This is an investigation implementation, not a general covariance estimator.
In particular neither posterior's arithmetic mean of mass is reported: the
positive mass tail makes that mean divergent under these measures.
"""

from functools import lru_cache
import math

import numpy as np


_LIMIT = 32.0
_PROBS = np.array([.025, .1, .25, .5, .75, .9, .975])


@lru_cache(None)
def _quadrature(order):
    # ``order`` is a resolution budget; there are about 20 smooth panels.
    degree = max(8, int(math.ceil(order / 6)))
    x, w = np.polynomial.legendre.leggauss(degree)
    v = np.polynomial.legendre.legvander(x, degree)
    projection = v[:, :degree] * w[:, None]
    projection *= (2 * np.arange(degree) + 1) / 2
    return x, w, v, projection


def _radial(h2, peak):
    """J1 and J3 with exp(-peak/2) common scaling, plus first moments."""
    h = np.sqrt(np.maximum(h2, 0))
    exponential = np.exp((h2 - peak) / 2)
    erf = np.fromiter((math.erf(float(v) / math.sqrt(2)) for v in h.flat),
                      dtype=float, count=h.size).reshape(h.shape)
    factor = np.ones_like(h)
    np.divide(math.sqrt(math.pi / 2) * erf, h, out=factor, where=h > 1e-12)
    j1 = exponential * factor
    return (j1, (1 + h2) * j1 + np.exp(-peak / 2),
            exponential, (h2 + 3) * exponential, factor)


def _antiderivative(coeff):
    """Legendre coefficients of the integral from -1, last axis polynomial."""
    result = np.zeros(coeff.shape[:-1] + (coeff.shape[-1] + 1,))
    result[..., 0] = coeff[..., 0]
    result[..., 1] = coeff[..., 0]
    for k in range(1, coeff.shape[-1]):
        term = coeff[..., k] / (2 * k + 1)
        result[..., k + 1] += term
        result[..., k - 1] -= term
    return result


def _evaluate(coeff, x):
    return np.polynomial.legendre.legval(x, np.moveaxis(coeff, -1, 0), tensor=False)


def _quantiles(cdf_coeff, edges, panel_end, left_mass, right_mass):
    n = len(edges)
    probability = np.broadcast_to(_PROBS, (n, len(_PROBS)))
    index = np.sum(panel_end[:, :, None] < probability[:, None, :], axis=1)
    index = np.minimum(index, panel_end.shape[1] - 1)
    rows = np.arange(n)[:, None]
    selected = cdf_coeff[rows, index]
    lo = np.full(probability.shape, -1.)
    hi = np.full(probability.shape, 1.)
    # 34 bisections give <4e-9 absolute log error even in the largest panel.
    for _ in range(34):
        middle = (lo + hi) / 2
        below = _evaluate(selected, middle) < probability
        lo = np.where(below, middle, lo)
        hi = np.where(below, hi, middle)
    lower = edges[rows, index]
    upper = edges[rows, index + 1]
    value = (lower + upper) / 2 + (upper - lower) * (lo + hi) / 4
    # Analytic endpoint tails also cover quantiles outside the numeric range.
    lm, rm = left_mass[:, None], right_mass[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        value = np.where(probability < lm,
                         -_LIMIT + np.log(probability / lm), value)
        value = np.where(probability > 1 - rm,
                         _LIMIT - np.log((1 - probability) / rm), value)
    return value


def _summary(density, edges, z, half, weights, vander, projection,
             left_coefficient, right_coefficient, y, truth_log_kernel):
    raw_panel_mass = np.sum(density * weights, axis=2) * half
    raw_left = left_coefficient * math.exp(-_LIMIT)
    raw_right = right_coefficient * math.exp(-_LIMIT)
    normalization = raw_left + raw_panel_mass.sum(axis=1) + raw_right
    if np.any(~np.isfinite(normalization)) or np.any(normalization <= 0):
        raise FloatingPointError('Posterior normalization failed')
    pdf = density / normalization[:, None, None]
    left_mass, right_mass = raw_left / normalization, raw_right / normalization
    panel_mass = raw_panel_mass / normalization[:, None]
    panel_end = left_mass[:, None] + np.cumsum(panel_mass, axis=1)
    panel_start = panel_end - panel_mass
    pdf_coeff = pdf @ projection
    cdf_coeff = _antiderivative(pdf_coeff) * half[:, :, None]
    cdf_coeff[..., 0] += panel_start
    cdf = np.clip(cdf_coeff @ vander.T, 0, 1)
    integral_cdf_coeff = _antiderivative(cdf_coeff) * half[:, :, None]
    panel_cdf_integral = _evaluate(integral_cdf_coeff, np.ones(half.shape))

    log_first = np.sum(np.sum(pdf * z * weights, axis=2) * half, axis=1)
    log_first += (_LIMIT + 1) * (right_mass - left_mass)
    log_second = np.sum(np.sum(pdf * z*z * weights, axis=2) * half, axis=1)
    log_second += (_LIMIT**2 + 2*_LIMIT + 2) * (left_mass + right_mass)
    quantiles = _quantiles(cdf_coeff, edges, panel_end, left_mass, right_mass)

    rows = np.arange(len(edges))
    index = np.clip(np.sum(edges <= y, axis=1) - 1, 0, half.shape[1] - 1)
    position = (y - (edges[rows, index] + edges[rows, index + 1])/2) / half[rows, index]
    position = np.clip(position, -1, 1)
    cdf_truth = np.clip(_evaluate(cdf_coeff[rows, index], position), 0, 1)
    previous_integral = np.cumsum(panel_cdf_integral, axis=1) - panel_cdf_integral
    integral_to_truth = (left_mass + previous_integral[rows, index]
                         + _evaluate(integral_cdf_coeff[rows, index], position))
    if y < -_LIMIT:
        cdf_truth = left_mass * math.exp(y + _LIMIT)
        integral_to_truth = cdf_truth.copy()
    elif y > _LIMIT:
        tail_factor = math.exp(_LIMIT - y)
        cdf_truth = 1 - right_mass * tail_factor
        integral_to_truth = (left_mass + panel_cdf_integral.sum(axis=1)
                             + y - _LIMIT - right_mass * (1 - tail_factor))

    # .5 E|Z-Z'| = E[Z(2F(Z)-1)], with analytic endpoint tail contributions.
    pair_half = np.sum(np.sum(z * (2*cdf - 1) * pdf * weights, axis=2) * half, axis=1)
    pair_half += ((_LIMIT + 1) * (left_mass + right_mass)
                  - (_LIMIT + .5) * (left_mass**2 + right_mass**2))
    crps = log_first - y + 2*integral_to_truth - pair_half
    density_score = np.log(normalization) - truth_log_kernel
    result = {'cdf_truth': cdf_truth, 'log_crps': crps,
              'log_density_score': density_score,
              'mean_log_mass': log_first,
              'log_sd': np.sqrt(np.maximum(0, log_second - log_first**2)),
              'log_median': quantiles[:, 3],
              'quadrature_tail_probability': left_mass + right_mass}
    for coverage, lower, upper in ((50, 2, 4), (80, 1, 5), (95, 0, 6)):
        result[f'log_lower_{coverage}'] = quantiles[:, lower]
        result[f'log_upper_{coverage}'] = quantiles[:, upper]

    roots = []
    for power in (.5, -.5):
        moment = np.sum(np.sum(pdf * np.exp(power*z) * weights, axis=2) * half, axis=1)
        moment += (left_coefficient * math.exp(-(1+power)*_LIMIT)/(1+power)
                   + right_coefficient * math.exp(-(1-power)*_LIMIT)/(1-power)) / normalization
        roots.append(moment)
    return result, roots[0] / roots[1]


def infer_batch(force, acceleration, truth=1., order=96):
    """Return batched point estimates and continuous log-mass diagnostics.

    ``force`` and ``acceleration`` have shape (n,3), or shape (3,) for one
    reading.  ``truth`` is one strictly positive finite mass, used only for
    CDF and scoring; it never selects estimator values or quadrature panels.
    ``order`` >= 24 controls within-panel resolution (96 means degree 16).
    Increase to 192 for convergence checks.  Fixed endpoint tails outside
    [-32,32] use their leading exponential, whose omitted corrections are
    negligible in the intended observed signal range (norms up to about 25).
    """
    force = np.atleast_2d(np.asarray(force, dtype=float))
    acceleration = np.atleast_2d(np.asarray(acceleration, dtype=float))
    if force.shape != acceleration.shape or force.ndim != 2 or force.shape[1] != 3:
        raise ValueError('force and acceleration must have identical (n,3) shapes')
    if not np.all(np.isfinite(force)) or not np.all(np.isfinite(acceleration)):
        raise ValueError('readings must be finite')
    if not np.isfinite(truth) or truth <= 0 or order < 24:
        raise ValueError('truth must be finite and positive; order must be >= 24')
    if not len(force):
        raise ValueError('at least one reading is required')
    p, q = np.sum(force**2, axis=1), np.sum(acceleration**2, axis=1)
    dot = np.sum(force*acceleration, axis=1)
    peak = np.where(dot >= 0, (p+q+np.hypot(p-q, 2*dot))/2, np.maximum(p, q))
    centre = .5*np.log((p+1)/(q+1))
    scale = np.sqrt(1/(p+1) + 1/(q+1))
    fixed = np.array([-_LIMIT, -20., -12., -8., -6., -4., -2., 0.,
                      2., 4., 6., 8., 12., 20., _LIMIT])
    # Narrow panels straddle the resolved interior mode.
    adaptive = centre[:, None] + scale[:, None] * np.array([-2., 0., 2.])
    adaptive = np.clip(adaptive, -_LIMIT + 1e-8, _LIMIT - 1e-8)
    edges = np.sort(np.concatenate((np.broadcast_to(fixed, (len(p), len(fixed))), adaptive), axis=1), axis=1)
    # Exact coincidences create zero panels, which are harmless for integrals;
    # the truth-panel lookup selects the last edge and skips such panels.
    half = np.diff(edges, axis=1) / 2
    midpoint = (edges[:, 1:] + edges[:, :-1]) / 2
    x, weights, vander, projection = _quadrature(order)
    z = midpoint[:, :, None] + half[:, :, None]*x
    theta = np.arctan(np.exp(z))
    sn, cs = np.sin(theta), np.cos(theta)
    jacobian = 1/(2*np.cosh(z))
    h2 = np.maximum(0, p[:, None, None]*sn**2 + q[:, None, None]*cs**2
                    + 2*dot[:, None, None]*sn*cs)
    radial = _radial(h2, peak[:, None, None])
    endpoints = _radial(np.stack((q, p), axis=1), peak[:, None])
    y = math.log(truth)
    # Stable sine/cosine for extreme finite truth, without squaring truth.
    ytheta = math.atan(truth)
    ys, yc = math.sin(ytheta), math.cos(ytheta)
    yh2 = np.maximum(0, p*ys**2 + q*yc**2 + 2*dot*ys*yc)
    yradial = _radial(yh2, peak)
    log_jacobian = -abs(y) - math.log1p(math.exp(-2*abs(y)))
    log_j1 = (yh2-peak)/2 + np.log(yradial[4])
    truth_log_kernel = (log_j1 + log_jacobian,
                        np.logaddexp(np.log1p(yh2)+log_j1, -peak/2) + log_jacobian)

    points, distributions = {}, {}
    for k, name in enumerate(('flat_joint', 'tube_joint')):
        first = radial[k+2]
        numerator = np.sum(np.sum(first*sn*jacobian*weights, axis=2)*half, axis=1)
        denominator = np.sum(np.sum(first*cs*jacobian*weights, axis=2)*half, axis=1)
        numerator += endpoints[k+2][:, 1]*math.exp(-_LIMIT)
        denominator += endpoints[k+2][:, 0]*math.exp(-_LIMIT)
        points[name] = numerator / denominator
        summary, root_point = _summary(radial[k]*jacobian, edges, z, half, weights,
                                      vander, projection, endpoints[k][:, 0], endpoints[k][:, 1],
                                      y, truth_log_kernel[k])
        # Exchange symmetry is an exact median oracle.  For equally strong
        # adverse readings, CDF is numerically indistinguishable from .5 over
        # a wide gap between modes, so direct inversion loses that information.
        summary['log_median'] = np.where(p == q, 0., summary['log_median'])
        distributions[name] = summary
        if k == 0:
            points['flat_geometric'] = np.exp(summary['mean_log_mass'])
            points['flat_median'] = np.exp(summary['log_median'])
            points['flat_reciprocal_root'] = root_point
    return {'points': points, 'distributions': distributions}
