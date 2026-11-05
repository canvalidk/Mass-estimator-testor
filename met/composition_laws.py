"""B18: full-tail mass laws and independent propagation, for studies only.

The historical ours29 composition calculation does not import this module.
Masses, CDFs and survival probabilities remain in logarithmic coordinates;
normal-score recalibration never exponentiates an extreme probability rank.
The grid uses linear density within each log-mass cell, supplemented by the
half-Cauchy exponential tail limits. Midpoint integration bounds below are
conditional on these numerical marginal CDFs, not a certification of their
discretization or of the underlying physical/statistical model.
"""
import math

import numpy as np
from scipy import special


LOG_TWO = math.log(2.)
LOG_HALF = -LOG_TWO
LOG_025 = math.log(.025)
LOG_975 = math.log(.975)


def _log1mexp(logp):
    """log(1-exp(logp)), stable on the whole closed log-probability range."""
    logp = np.asarray(logp, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return np.where(logp < LOG_HALF, np.log1p(-np.exp(logp)), np.log(-np.expm1(logp)))


def _logdiff(loga, logb):
    """log(exp(loga)-exp(logb)); equality, including both -inf, gives -inf."""
    loga, logb = np.broadcast_arrays(loga, logb)
    with np.errstate(invalid="ignore"):
        difference = np.where(loga == logb, 0., logb-loga)
    return loga+_log1mexp(np.minimum(difference, 0.))


def _normal_rank(logp, logsf):
    logp, logsf = np.broadcast_arrays(np.asarray(logp, float), np.asarray(logsf, float))
    lower = logp <= logsf
    value = np.empty(logp.shape)
    value[lower] = special.ndtri_exp(logp[lower])
    value[~lower] = -special.ndtri_exp(logsf[~lower])
    return value


def _as_output(value):
    return float(value) if np.ndim(value) == 0 else value


class LogMassLaw:
    """A proper numerical law on all real log masses.

    ``log_density`` is a density with respect to log mass. Between supplied
    nodes the ordinary density is linear. Below/above the grid it has log
    density ``left_log_constant+x`` / ``right_log_constant-x``. If constants
    are omitted these tails are matched to the endpoint densities. ``warp``
    selects the ours33 normal-score quoted CDF; a non-warped law already
    contains any likelihood power and has identical quoted/raw methods.
    """

    def __init__(self, lm, log_density, *, eta=1., warp=False,
                 left_log_constant=None, right_log_constant=None):
        self.lm = np.asarray(lm, dtype=float)
        self.log_density = np.asarray(log_density, dtype=float)
        self.eta, self.warp = float(eta), bool(warp)
        if (self.lm.ndim != 1 or len(self.lm) < 2 or self.log_density.shape != self.lm.shape
                or np.any(~np.isfinite(self.lm)) or np.any(~np.isfinite(self.log_density))
                or np.any(np.diff(self.lm) <= 0)):
            raise ValueError("finite densities and increasing log-mass nodes are required")
        if not math.isfinite(self.eta) or not 0 < self.eta <= 1:
            raise ValueError("eta must be finite in (0,1]")
        self.left_log_constant = float(self.log_density[0]-self.lm[0] if left_log_constant is None
                                       else left_log_constant)
        self.right_log_constant = float(self.log_density[-1]+self.lm[-1] if right_log_constant is None
                                        else right_log_constant)
        if not all(math.isfinite(x) for x in (self.left_log_constant, self.right_log_constant)):
            raise ValueError("tail constants must be finite")
        self._log_cells = (np.log(np.diff(self.lm))
                           + np.logaddexp(self.log_density[:-1], self.log_density[1:])-LOG_TWO)
        left = self.left_log_constant+self.lm[0]
        right = self.right_log_constant-self.lm[-1]
        self._log_z = float(special.logsumexp(np.r_[left, self._log_cells, right]))
        self._lcdf = np.r_[left, np.logaddexp(left, np.logaddexp.accumulate(self._log_cells))]-self._log_z
        self._lsf = np.r_[np.logaddexp(right, np.logaddexp.accumulate(self._log_cells[::-1]))[::-1], right]-self._log_z

    def _raw_logcdf(self, logmass):
        value = np.asarray(logmass, dtype=float)
        if np.any(np.isnan(value)):
            raise ValueError("log mass must not be NaN")
        flat = value.ravel()
        lower, upper = np.empty_like(flat), np.empty_like(flat)
        left, right = flat <= self.lm[0], flat >= self.lm[-1]
        middle = ~(left | right)
        lower[left] = self.left_log_constant+flat[left]-self._log_z
        upper[left] = _log1mexp(lower[left])
        upper[right] = self.right_log_constant-flat[right]-self._log_z
        lower[right] = _log1mexp(upper[right])
        if np.any(middle):
            x = flat[middle]
            index = np.searchsorted(self.lm, x, side="right")-1
            h = self.lm[index+1]-self.lm[index]
            fraction = (x-self.lm[index])/h
            with np.errstate(divide="ignore"):
                partial = (np.log(h)+np.log(fraction)
                           + np.logaddexp(np.log1p(-fraction/2)+self.log_density[index],
                                          np.log(fraction/2)+self.log_density[index+1])-self._log_z)
                remaining = (np.log(h)+np.log1p(-fraction)
                             + np.logaddexp(np.log((1+fraction)/2)+self.log_density[index+1],
                                            np.log((1-fraction)/2)+self.log_density[index])-self._log_z)
            lower[middle] = np.logaddexp(self._lcdf[index], partial)
            upper[middle] = np.logaddexp(self._lsf[index+1], remaining)
        # Retain the directly accumulated smaller tail and obtain its
        # complement stably. Prefix/suffix rounding cannot make F+SF drift.
        small_lower = lower <= upper
        upper[small_lower] = _log1mexp(lower[small_lower])
        lower[~small_lower] = _log1mexp(upper[~small_lower])
        return lower.reshape(value.shape), upper.reshape(value.shape)

    def logcdf(self, logmass, *, quoted=True):
        """Return (log CDF, log survival), accepting scalar or array log masses."""
        lower, upper = self._raw_logcdf(logmass)
        if quoted and self.warp and self.eta < 1:
            z = math.sqrt(self.eta)*_normal_rank(lower, upper)
            lower, upper = special.log_ndtr(z), special.log_ndtr(-z)
        return _as_output(lower), _as_output(upper)

    def _inverse_tail(self, logrank, *, upper):
        nodes = -self.lm[::-1] if upper else self.lm
        density = self.log_density[::-1] if upper else self.log_density
        accumulated = self._lsf[::-1] if upper else self._lcdf
        cell = self._log_cells[::-1] if upper else self._log_cells
        constant = self.right_log_constant if upper else self.left_log_constant
        ranks = np.asarray(logrank, dtype=float)
        flat = ranks.ravel()
        output = np.empty_like(flat)
        indices = np.searchsorted(accumulated, flat, side="right")-1
        tail = indices < 0
        output[tail] = flat[tail]+self._log_z-constant
        inside = ~tail
        if np.any(inside):
            j = np.minimum(indices[inside], len(nodes)-2)
            desired = _logdiff(flat[inside], accumulated[j])
            fraction_mass = np.exp(desired-(cell[j]-self._log_z))
            # This clips only cell-inversion rounding, never a probability
            # rank or an extreme-law quantile.
            fraction_mass = np.clip(fraction_mass, 0., 1.)
            scale = np.maximum(density[j], density[j+1])
            d0, d1 = np.exp(density[j]-scale), np.exp(density[j+1]-scale)
            radicand = np.maximum(0., d0*d0+(d1-d0)*fraction_mass*(d0+d1))
            denominator = d0+np.sqrt(radicand)
            fraction = np.divide(fraction_mass*(d0+d1), denominator,
                                 out=np.zeros_like(denominator), where=denominator > 0)
            output[inside] = nodes[j]+fraction*(nodes[j+1]-nodes[j])
        result = output.reshape(ranks.shape)
        return -result if upper else result

    def logquantile(self, logp, logsf, *, quoted=True):
        """Invert a paired log probability/survival rank without tail clipping."""
        logp, logsf = np.broadcast_arrays(np.asarray(logp, float), np.asarray(logsf, float))
        if (np.any(np.isnan(logp)) or np.any(np.isnan(logsf))
                or np.any(logp > 0) or np.any(logsf > 0)
                or np.any((logp == -np.inf) & (logsf == -np.inf))):
            raise ValueError("valid paired log probability/survival ranks are required")
        if quoted and self.warp and self.eta < 1:
            z = _normal_rank(logp, logsf)/math.sqrt(self.eta)
            upper = z > 0
            logsmall = special.log_ndtr(-np.abs(z))
        else:
            upper = logsf < logp
            logsmall = np.minimum(logp, logsf)
        # Each inverse uses the smaller tail, whose rank is <= 1/2.
        result = np.empty(logsmall.shape)
        result[upper] = self._inverse_tail(logsmall[upper], upper=True)
        result[~upper] = self._inverse_tail(logsmall[~upper], upper=False)
        return _as_output(result)

    def log_interval(self, *, quoted=True):
        values = self.logquantile(np.array([LOG_025, LOG_975]),
                                  np.array([LOG_975, LOG_025]), quoted=quoted)
        return float(values[0]), float(values[1])

    def median_log(self, *, quoted=True):
        return float(self.logquantile(LOG_HALF, LOG_HALF, quoted=quoted))


def _isotropic_residual(reading, c, phi):
    """Exact isotropic shortcut, or None for the full QR reference path."""
    d = reading.dimension
    vf, va = float(reading.cov[0, 0]), float(reading.cov[d, d])
    diagonal = np.r_[np.full(d, vf), np.full(d, va)]
    if not np.array_equal(reading.cov, np.diag(diagonal)):
        return None
    force = reading.subspace.T @ reading.force/math.sqrt(vf)
    acceleration = reading.subspace.T @ reading.acceleration/math.sqrt(va)
    uf, ua = c*np.sin(phi)/math.sqrt(vf), np.cos(phi)/math.sqrt(va)
    length = np.hypot(uf, ua)
    residual = (ua/length)[:, None]*force-(uf/length)[:, None]*acceleration
    return np.sum(residual*residual, axis=1)


def _fast_isotropic_law(readings, c, equation, nu_eff):
    """Same continuous factors/search, in their common declared noise angle.

    This path applies only to exactly block-isotropic covariance with c equal
    to every force/acceleration noise ratio. The observations supply all
    sufficient statistics; no investigator truth or simulation parameter is
    passed to the fit. Other declarations use the full probe reference.
    """
    from met.full import _log_prior
    from met.ours32 import bartlett_eta

    data = []
    for reading in readings:
        d = reading.dimension
        vf, va = float(reading.cov[0, 0]), float(reading.cov[d, d])
        diagonal = np.r_[np.full(d, vf), np.full(d, va)]
        if (not np.array_equal(reading.cov, np.diag(diagonal))
                or c != math.sqrt(vf/va)):
            return None
        basis, _ = np.linalg.qr(reading.subspace, mode="reduced")
        f = basis.T @ reading.force/math.sqrt(vf)
        a = basis.T @ reading.acceleration/math.sqrt(va)
        if np.any(~np.isfinite(f)) or np.any(~np.isfinite(a)):
            raise ValueError("reading observations must be finite")
        ff, aa, fa = float(f@f), float(a@a), float(f@a)
        energy, gap = ff+aa, math.hypot(ff-aa, 2*fa)
        phi = .5*math.atan2(2*fa, aa-ff)
        q = max(0., .5*(energy-gap))
        qplus = q if fa >= 0 else min(ff, aa)
        data.append((reading, f, a, ff, aa, fa, energy, phi, q, qplus))
    if equation in ("ours32", "ours33"):
        if nu_eff is None:
            raise ValueError("ours32 and ours33 require nu_eff")
        degrees = ([float(nu_eff)]*len(readings) if np.ndim(nu_eff) == 0 else
                   [float(value) for value in nu_eff])
        if len(degrees) != len(readings) or not all(math.isfinite(v) and v > 0 for v in degrees):
            raise ValueError("nu_eff must contain one positive finite value per reading")
    else:
        degrees = [None]*len(readings)
    weights = []
    if equation in ("ours29", "ours30"):
        for reading, _, _, _, _, _, _, _, q, qplus in data:
            numerator = reading.k-(.5 if equation == "ours30" else 1.)
            misfit = qplus if equation == "ours30" else q
            weights.append(1. if numerator <= 0 or misfit <= 0 else min(1., numerator/misfit))

    def terms(unused_states, unused_c, phi, order):
        sine, cosine = np.sin(phi)[:, None], np.cos(phi)[:, None]
        answer = []
        for item, nu in zip(data, degrees):
            reading, f, a = item[:3]
            residual = f*cosine-a*sine
            cost = np.sum(residual*residual, axis=1)
            coefficient = .5*(nu+reading.k-(equation == "ours32"))
            value = -coefficient*np.log1p(cost/nu)
            if order == 0:
                answer.append((value,))
            else:
                along = f*sine+a*cosine
                h1 = 2*np.sum(residual*along, axis=1)
                h2 = 2*(cost-np.sum(along*along, axis=1))
                denominator = nu+cost
                answer.append((value, coefficient*h1/denominator,
                               coefficient*(h2/denominator+(h1/denominator)**2)))
        return answer

    if not data:
        eta = 1.
    elif equation in ("ours32", "ours33"):
        states = [(item[0], nu, dict(phi=item[7], E=item[6]), None)
                  for item, nu in zip(data, degrees)]
        eta = bartlett_eta(states, c, terms=terms)[0]
    else:
        ag, bg, cg = (math.fsum(w*item[index] for w, item in zip(weights, data))
                      for index in (3, 4, 5))
        phi = .5*math.atan2(2*cg, bg-ag)
        sine, cosine = math.sin(phi), math.cos(phi)
        first, second = [], []
        for item, weight in zip(data, weights):
            _, f, a = item[:3]
            residual, along = f*cosine-a*sine, f*sine+a*cosine
            first.append(weight*float(residual@along))
            second.append(weight*float(residual@residual-along@along))
        curvature, j = -math.fsum(second), math.fsum(value*value for value in first)
        eta = 1. if j <= 1e-300 or curvature <= 0 else min(1., curvature/j)

    def factors(phi):
        if equation in ("ours32", "ours33"):
            return ((eta if equation == "ours32" else 1.)
                    * sum((part[0] for part in terms(None, c, phi, 0)), np.zeros_like(phi)))
        sine, cosine = np.sin(phi)[:, None], np.cos(phi)[:, None]
        total = np.zeros_like(phi)
        for item, weight in zip(data, weights):
            _, f, a = item[:3]
            residual = f*cosine-a*sine
            total -= .5*eta*weight*np.sum(residual*residual, axis=1)
        return total

    center = math.log(c)
    coarse = center+np.linspace(-25., 25., 4001)
    phi = np.arctan(np.exp(coarse-center))
    value = _log_prior(coarse, c)+factors(phi)
    keep = np.flatnonzero(value > value.max()-40.)
    step = coarse[1]-coarse[0]
    native = np.linspace(coarse[keep[0]]-2*step, coarse[keep[-1]]+2*step, 2001)
    nodes = np.unique(np.r_[center+np.linspace(-35., 35., 7001), native])
    nodes = nodes[np.r_[True, np.diff(nodes) > 1e-11]]
    shifted = nodes-center
    phi = np.arctan2(np.exp(np.minimum(shifted, 0.)), np.exp(np.minimum(-shifted, 0.)))
    density = _log_prior(nodes, c)+factors(phi)
    endpoints = factors(np.array([0., math.pi/2]))
    return LogMassLaw(nodes, density, eta=eta, warp=(equation == "ours33"),
                      left_log_constant=float(endpoints[0]-center),
                      right_log_constant=float(endpoints[1]+center))


def analyse(readings, c, equation, nu_eff=None):
    """Build the chosen probe's law, with the half-Cauchy prior at c.

    ours29/30/32 already quote a likelihood-powered law; their raw/quoted
    methods coincide here. ours33 retains its unpowered raw law and applies
    its normal-score warp only through the public quoted methods.
    """
    from met.full import _log_prior, ours29, ours30
    from met.ours31 import _per_reading
    from met.ours32 import _residual_energy, bartlett_eta, ours32
    from met.ours33 import _l_terms, ours33

    readings, c = list(readings), float(c)
    if not math.isfinite(c) or c <= 0:
        raise ValueError("c must be finite and positive")
    if equation not in ("ours29", "ours30", "ours32", "ours33"):
        raise ValueError("unsupported mass equation")
    fast = _fast_isotropic_law(readings, c, equation, nu_eff)
    if fast is not None:
        return fast
    if equation in ("ours29", "ours30"):
        result = (ours29 if equation == "ours29" else ours30)(readings, c, readout=False)
        eta, weights = result["eta"], result["g"]
        degrees = [None]*len(readings)
    else:
        if nu_eff is None:
            raise ValueError("ours32 and ours33 require nu_eff")
        result = (ours32 if equation == "ours32" else ours33)(
            readings, c, nu_eff, calibrate=(equation == "ours32"))
        degrees = result["nu_eff"]
        if equation == "ours33":
            eta = bartlett_eta(_per_reading(readings, c, degrees), c, terms=_l_terms)[0]
        else:
            eta = result["eta"]
        weights = None
    center = math.log(c)
    regular = center+np.linspace(-35., 35., 7001)
    nodes = np.unique(np.r_[regular, result["lm"]])
    nodes = nodes[np.r_[True, np.diff(nodes) > 1e-11]]

    def factors(phi):
        value = np.zeros_like(phi, dtype=float)
        for index, (reading, nu) in enumerate(zip(readings, degrees)):
            residual = _isotropic_residual(reading, c, phi)
            if residual is None:
                residual = _residual_energy(reading, c, phi)
            if equation in ("ours29", "ours30"):
                value -= .5*eta*weights[index]*residual
            else:
                coefficient = .5*(nu+reading.k-(equation == "ours32"))
                value -= (eta if equation == "ours32" else 1.)*coefficient*np.log1p(residual/nu)
        return value

    shifted = nodes-center
    phi = np.arctan2(np.exp(np.minimum(shifted, 0.)), np.exp(np.minimum(-shifted, 0.)))
    density = _log_prior(nodes, c)+factors(phi)
    endpoints = factors(np.array([0., math.pi/2]))
    return LogMassLaw(nodes, density, eta=eta, warp=(equation == "ours33"),
                      left_log_constant=float(endpoints[0]-center),
                      right_log_constant=float(endpoints[1]+center))


class _Reciprocal:
    def __init__(self, law):
        self.law = law

    def logcdf(self, logmass, *, quoted=True):
        lower, upper = self.law.logcdf(-np.asarray(logmass), quoted=quoted)
        return upper, lower

    def logquantile(self, logp, logsf, *, quoted=True):
        return -np.asarray(self.law.logquantile(logsf, logp, quoted=quoted))


def composed_logcdf(p1, p2, logvalue, composition, *, quoted=True, nodes=20000):
    """Independent sum/parallel CDF, SF and log absolute midpoint error bound.

    The bound applies to the one-dimensional monotone midpoint integration
    with these supplied numerical marginal CDFs. It does not cover their
    discretization error. No simulated draws, truth-dependent tail clipping,
    dependence model, or post-composition eta is introduced here.
    """
    if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or nodes < 1:
        raise ValueError("nodes must be a positive integer")
    value = float(logvalue)
    if math.isnan(value):
        raise ValueError("log value must not be NaN")
    if composition in ("L14", "parallel", "common_force"):
        lower, upper, error = composed_logcdf(_Reciprocal(p1), _Reciprocal(p2), -value,
                                            "L5", quoted=quoted, nodes=nodes)
        return upper, lower, error
    if composition not in ("L5", "sum", "lumping"):
        raise ValueError("composition must be L5/sum or L14/parallel")
    if value == -math.inf:
        return -math.inf, 0., -math.inf
    if value == math.inf:
        return 0., -math.inf, -math.inf
    logu = np.log((np.arange(nodes)+.5)/nodes)
    logn = math.log(nodes)
    median1 = float(p1.logquantile(LOG_HALF, LOG_HALF, quoted=quoted))
    median2 = float(p2.logquantile(LOG_HALF, LOG_HALF, quoted=quoted))
    median_sum = float(np.logaddexp(median1, median2))
    if value < median_sum:
        f1, _ = p1.logcdf(value, quoted=quoted)
        f2, _ = p2.logcdf(value, quoted=quoted)
        ranks = float(f1)+logu
        quantiles = p1.logquantile(ranks, _log1mexp(ranks), quoted=quoted)
        residual = _logdiff(value, quantiles)
        conditional, _ = p2.logcdf(residual, quoted=quoted)
        lower = float(f1)+float(special.logsumexp(conditional))-logn
        bound = float(f1)+float(f2)-LOG_TWO-logn
        return lower, float(_log1mexp(lower)), bound

    half = value-LOG_TWO
    f1h, s1h = (float(x) for x in p1.logcdf(half, quoted=quoted))
    f2h, s2h = (float(x) for x in p2.logcdf(half, quoted=quoted))
    _, s1t = p1.logcdf(value, quoted=quoted)
    _, s2t = p2.logcdf(value, quoted=quoted)

    def piece(first, second, logf):
        ranks = logf+logu
        quantiles = first.logquantile(ranks, _log1mexp(ranks), quoted=quoted)
        residual = _logdiff(value, quantiles)
        _, conditional = second.logcdf(residual, quoted=quoted)
        return logf+float(special.logsumexp(conditional))-logn

    upper = float(special.logsumexp([s1h+s2h, piece(p1, p2, f1h), piece(p2, p1, f2h)]))
    bound = float(np.logaddexp(f1h+float(_logdiff(s2h, s2t)),
                              f2h+float(_logdiff(s1h, s1t)))-LOG_TWO-logn)
    return float(_log1mexp(upper)), upper, bound
