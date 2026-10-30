"""B04: fitted-score bias and an unadopted effective-reading Bartlett power.

The exact linear identity assumes fixed nonnegative curvatures. The ours30
curvatures below are signed, data-dependent derivatives; the identity is a
comparison target, not an assumption about their repeated-world distribution.
Only isotropic readings sharing one noise ratio per series are supported.
"""
import numpy as np

from .analyst import _cdf, _quantiles, misfit_weights


# Keep this object alive: ReadingSet caches analyst results by its identity.
_OURS30 = {"misfit": "kink", "positivity": True}


def linear_moments(c):
    """Exact covariance of fitted Gaussian scores, allowing zero-weight rows.

    For independent x_i with variance 1/c_i on the positive-c_i rows, the
    fitted scores have covariance diag(c)-c*c.T/C. A zero-c_i row is an
    exact-null contribution and needs no fictitious infinite-variance draw.
    """
    c = np.asarray(c, dtype=float)
    if c.ndim != 1 or c.size == 0 or not np.isfinite(c).all() or np.any(c < 0):
        raise ValueError("c must be a nonempty finite nonnegative vector")
    total = float(c.sum())
    if not total > 0 or not np.isfinite(total):
        raise ValueError("the total curvature must be finite and positive")
    covariance = np.diag(c) - np.outer(c, c) / total
    square_sum = float(np.dot(c, c))
    return dict(C=total, n_eff=total ** 2 / square_sum,
                expected_J=total - square_sum / total,
                variance_J=2 * float(np.sum(covariance ** 2)),
                score_covariance=covariance)


def _validate_readings(readings):
    if readings.series < 1 or readings.readings < 1 or not readings.isotropic:
        raise ValueError("eta_effective needs nonempty isotropic reading series")
    # Declaring sf=c*sa can differ by a few ulps after the division sf/sa.
    # This tolerance covers that arithmetic representation of one ratio only.
    if not np.allclose(readings.s, readings.s[:, :1],
                       rtol=16 * np.finfo(float).eps, atol=0.0):
        raise ValueError("each series must share one declared noise ratio c=s")
    if not all(np.isfinite(value).all() for value in
               (readings.P, readings.Q, readings.D, readings.s)):
        raise ValueError("standardized readings and noise ratios must be finite")


def statistics(readings):
    """Ours30 per-reading derivatives at the exact whole-line best angle.

    Return B-vectors plus (B,N) scores, curvatures and fixed misfit weights.
    n_eff=C²/sum(c_i²), even when some c_i are negative; it is not necessarily
    a count >=1. If C=sum(c_i²)=0, n_eff and linear_target are NaN. If C=0
    but the individual curvatures do not vanish, n_eff=0 and the target is
    -infinity. J/C is NaN at J=C=0 and infinity at J>0=C. These degeneracies
    are exposed rather than silently included in finite repeated-world means.
    Both powers use full.py's C<=0 or J<=1e-300 fallback to 1. The candidate
    keeps eta when n_eff<=1 or n_eff is undefined.
    """
    _validate_readings(readings)
    g = misfit_weights(readings, _OURS30).copy()
    direction = (readings.Q - readings.P) + 2j * readings.D
    total = np.sum(g * direction, axis=1)
    phi = 0.5 * np.angle(total)
    phase = np.exp(2j * phi)
    per = g * np.conj(direction) * phase[:, None]
    curvatures = np.real(per)
    scores = -0.5 * np.imag(per)
    C = np.abs(total)
    J = np.sum(scores ** 2, axis=1)
    sum_c2 = np.sum(curvatures ** 2, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        n_eff = np.divide(C ** 2, sum_c2, out=np.full_like(C, np.nan), where=sum_c2 > 0)
        linear_target = 1 - sum_c2 / C ** 2
        ratio = J / C
        eta = np.where((C > 0) & (J > 1e-300), np.minimum(1.0, C / J), 1.0)
        eta_prime = np.where((C > 0) & (J > 1e-300) & (n_eff > 1),
                             np.minimum(1.0, C / J * (1 - 1 / n_eff)), eta)
    return dict(C=C, J=J, sum_c2=sum_c2, n_eff=n_eff, linear_target=linear_target,
                ratio=ratio, eta=eta, eta_prime=eta_prime, phi=phi, g=g,
                curvatures=curvatures, scores=scores,
                negative_curvature_fraction=np.mean(curvatures < 0, axis=1),
                n_eff_le_one=(n_eff <= 1), degenerate=(C <= 0),
                weighted_p=np.sum(g * readings.P, axis=1),
                weighted_q=np.sum(g * readings.Q, axis=1),
                weighted_d=np.sum(g * readings.D, axis=1))


def _log_density(x, a, b, power):
    """Log density in x=log(m/c), with an m-independent constant omitted.

    sum_i g_i H_i/2 = const + (a*cos(2phi)+b*sin(2phi))/4,
    cos(2phi)=-tanh(x), sin(2phi)=sech(x). The prior in x is sech(x).
    """
    return (-np.logaddexp(x, -x)
            + power[:, None] / 4 * (-a[:, None] * np.tanh(x) + b[:, None] / np.cosh(x)))


def _intervals_for_power(a, b, power, center, grid_points):
    # Match full.py's support search, but evaluate only the exact three-sum law.
    coarse = np.broadcast_to(np.linspace(-25.0, 25.0, 4001), (len(a), 4001))
    lp = _log_density(coarse, a, b, power)
    keep = lp > lp.max(axis=1, keepdims=True) - 40.0
    step = coarse[0, 1] - coarse[0, 0]
    left = np.argmax(keep, axis=1)
    right = keep.shape[1] - 1 - np.argmax(keep[:, ::-1], axis=1)
    rows = np.arange(len(a))
    lo, hi = coarse[rows, left] - 2 * step, coarse[rows, right] + 2 * step
    fraction = np.linspace(0.0, 1.0, grid_points)
    x = lo[:, None] + (hi - lo)[:, None] * fraction
    density_log = _log_density(x, a, b, power)
    density = np.exp(density_log - density_log.max(axis=1, keepdims=True))
    h = (hi - lo) / (grid_points - 1)
    cdf = _cdf(density, h)
    area = h * (density.sum(axis=1) - 0.5 * (density[:, 0] + density[:, -1]))
    quantiles = _quantiles(x, density / area[:, None], cdf, [0.025, 0.975])
    # The same plateau convention as full.py, for a quantile between modes.
    from .full import _plateau_centre
    for i in range(len(a)):
        for j, probability in enumerate((0.025, 0.975)):
            quantiles[i, j] = _plateau_centre(x[i], cdf[i], probability, quantiles[i, j])
    return center[:, None] * np.exp(quantiles)


def intervals(readings, stats, *, grid_points=2001):
    """Paired equal-tail 95% intervals for eta and the study-only eta_prime.

    Both powers temper exactly the same fixed-g summed likelihood. The
    half-Cauchy centre is the series' declared common ratio c=s. Numerical
    support and inversion match full.py; N affects three sums, not grid cost.
    """
    _validate_readings(readings)
    if isinstance(grid_points, bool) or not isinstance(grid_points, (int, np.integer)) or grid_points < 201:
        raise ValueError("grid_points must be an integer >=201")
    keys = ("weighted_p", "weighted_q", "weighted_d", "eta", "eta_prime")
    if any(np.shape(stats[key]) != (readings.series,) for key in keys):
        raise ValueError("statistics must contain one value per reading series")
    a, b = stats["weighted_q"] - stats["weighted_p"], 2 * stats["weighted_d"]
    center = readings.s[:, 0]
    original = _intervals_for_power(a, b, stats["eta"], center, grid_points)
    candidate = _intervals_for_power(a, b, stats["eta_prime"], center, grid_points)
    return dict(interval_lo=original[:, 0], interval_hi=original[:, 1],
                prime_lo=candidate[:, 0], prime_hi=candidate[:, 1])
