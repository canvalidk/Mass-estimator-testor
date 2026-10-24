"""Full-joint isotropic inference with independent variance calibration.

The physical reference measure is df d(alpha) dOmega.  Each channel has
one unknown isotropic variance, with scale prior p(sigma**2) proportional
to 1/sigma**2 and independent calibration nu*s**2/sigma**2 ~ chi-square(nu).
The supplied SD is for the observed mean, not an individual repeat. Only
that mean and the external calibration statistic enter this likelihood.

Calibration-only precisions are independent Gamma(nu/2, nu*s**2/2).
Writing the two scaled precisions as t*u and t*(1-u) integrates their common
Gamma scale exactly. A shifted Beta quadrature integrates the remaining
ratio, INCLUDING measurement evidence. This is not an equal mixture of
normalized conditional posteriors. Radial/direction integrals are analytic.
The mass arithmetic mean still diverges; E[f]/E[alpha] exists.

The incomplete-beta identities used below are NIST DLMF 8.17.1, 8.17.2,
8.17.4 and 8.17.22: https://dlmf.nist.gov/8.17 .
"""

from functools import lru_cache
import math

import numpy as np

from comparison_posterior import _LIMIT, _quadrature, _summary


@lru_cache(None)
def _beta_nodes(a, b, order):
    """Gauss-Jacobi nodes, normalized for a Beta(a,b) expectation."""
    alpha, beta = b - 1., a - 1.
    k = np.arange(order, dtype=float)
    diagonal = ((beta - alpha) * (beta + alpha)
                / ((2*k + alpha + beta) * (2*k + alpha + beta + 2)))
    j = k[1:]
    off = (2 / (2*j + alpha + beta)
           * np.sqrt(j*(j+alpha)*(j+beta)*(j+alpha+beta)
                     / ((2*j+alpha+beta-1)*(2*j+alpha+beta+1))))
    nodes, vectors = np.linalg.eigh(np.diag(diagonal) + np.diag(off, 1)
                                   + np.diag(off, -1))
    weights = vectors[0]**2
    weights /= weights.sum()
    return (nodes + 1)/2, weights


def _beta_fraction(a, b, x):
    """Modified-Lentz incomplete-beta continued fraction, vectorized."""
    tiny = 1e-300
    def nonzero(value):
        return np.where(np.abs(value) < tiny,
                        np.where(value < 0, -tiny, tiny), value)
    qab, qap, qam = a+b, a+1, a-1
    c = np.ones_like(x)
    d = 1 / nonzero(1 - qab*x/qap)
    value = d.copy()
    for m in range(1, 257):
        m2 = 2*m
        coefficient = m*(b-m)*x / ((qam+m2)*(a+m2))
        d = 1/nonzero(1+coefficient*d)
        c = nonzero(1+coefficient/c)
        value *= d*c
        coefficient = -(a+m)*(qab+m)*x / ((a+m2)*(qap+m2))
        d = 1/nonzero(1+coefficient*d)
        c = nonzero(1+coefficient/c)
        change = d*c
        value *= change
        if np.all(np.abs(change-1) < 8*np.finfo(float).eps):
            return value
    raise FloatingPointError('Incomplete-beta continued fraction did not converge')


def _log_beta_regularized(a, b, x):
    """Log regularized incomplete beta for strictly interior arguments."""
    log_front = (math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)
                 + a*np.log(x) + b*np.log1p(-x))
    # With b large and x tiny, evaluating the complement at 1-x loses
    # significant digits in its continued fraction. The direct fraction
    # remains well conditioned over this moderate b*x range.
    direct = (x < (a+1)/(a+b+2)) | ((x < .01) & (b*x < 20))
    result = np.empty_like(x)
    if np.any(direct):
        result[direct] = (log_front[direct]
                          + np.log(_beta_fraction(a,b,x[direct])) - math.log(a))
    if np.any(~direct):
        complement = np.exp(log_front[~direct]) * _beta_fraction(b,a,1-x[~direct])/b
        result[~direct] = np.log1p(-np.clip(complement, 0., 1.))
    return result


def _log_radial(h2, residual, power):
    """Log int_0^1 [1+residual/2+h2*v**2/2]**(-power) dv."""
    h2, residual = np.broadcast_arrays(np.asarray(h2,float), np.asarray(residual,float))
    log_c = np.log1p(residual/2)
    ratio = h2/(2+residual)
    result = -power*log_c
    small = power*ratio < 1e-5
    if np.any(small):
        r = ratio[small]
        # Expansion integrated termwise in v: integral v**(2j)=1/(2j+1).
        correction = (-power*r/3 + power*(power+1)*r*r/10
                      - power*(power+1)*(power+2)*r*r*r/42)
        result[small] += np.log1p(correction)
    if np.any(~small):
        r = ratio[~small]
        if power > 1000:
            # Large beta shapes make the complementary continued fraction
            # sensitive to subtraction at 1-x. Scale v by sqrt(power*r):
            # the resulting integrand has a resolved, near-Gaussian width.
            # For power>1000, omitting t>12 has absolute integral error
            # below 1e-59, bounded by power/[24*(power-1)] times
            # (1+144/power)**(1-power). GL64 resolves the retained integral.
            width = np.sqrt(power*r)
            upper = np.minimum(width,12.)
            nodes,weights = _radial_nodes()
            integral = np.zeros_like(width)
            for node,weight in zip(nodes,weights):
                integral += weight*np.exp(-power*np.log1p((upper*node)**2/power))
            result[~small] += np.log(integral*upper)-np.log(width)
        else:
            x = r/(1+r)
            if np.any(x >= 1):
                raise FloatingPointError('Calibration radial argument exceeds floating-point range')
            result[~small] += (-math.log(2)-.5*np.log(r)
                              + math.lgamma(.5)+math.lgamma(power-.5)-math.lgamma(power)
                              + _log_beta_regularized(.5,power-.5,x))
    return result


@lru_cache(None)
def _radial_nodes():
    nodes,weights = np.polynomial.legendre.leggauss(64)
    return (nodes+1)/2,weights/2


def _log_jacobian(z):
    return -np.abs(z)-np.log1p(np.exp(-2*np.abs(z)))


def _degrees(value):
    array = np.asarray(value, dtype=float)
    if array.shape == ():
        array = np.repeat(array, 2)
    if array.shape != (2,) or np.any(~np.isfinite(array)) or np.any(array <= 0):
        raise ValueError('degrees must be positive finite scalar or channel pair')
    return tuple(float(v) for v in array)


def infer_batch(force, acceleration, truth=1., *, force_sd, acceleration_sd,
                degrees, order=96, calibration_order=24):
    """Return physical mass points and continuous posterior diagnostics.

    Readings have shape (n,3), SDs scalar or (n,), degrees scalar or (force,
    acceleration) pair. Calibration SDs are mean-observation scales. Truth
    only evaluates scores; it does not choose points, panels, or densities.
    Increase both order and calibration_order for numerical refinement.
    """
    force = np.atleast_2d(np.asarray(force,dtype=float))
    acceleration = np.atleast_2d(np.asarray(acceleration,dtype=float))
    if (force.shape != acceleration.shape or force.ndim != 2
            or force.shape[1] != 3 or len(force) == 0
            or np.any(~np.isfinite(force)) or np.any(~np.isfinite(acceleration))):
        raise ValueError('Finite force and acceleration arrays must share nonempty (n,3) shape')
    size = len(force)
    standard_deviations = []
    for value in (force_sd, acceleration_sd):
        array = np.asarray(value,dtype=float)
        if array.shape not in ((),(size,)) or np.any(~np.isfinite(array)) or np.any(array <= 0):
            raise ValueError('Calibration SDs must be positive finite scalars or (n,) arrays')
        standard_deviations.append(np.broadcast_to(array,(size,)))
    physical_sf, physical_sa = standard_deviations
    mass_shift = np.log(physical_sf)-np.log(physical_sa)
    # Work in the fixed observed-calibration units. This is a scalar change
    # of units in each channel; its physical df d(alpha) Jacobian is constant
    # and cancels from the conditional posterior. A relative log-mass grid
    # avoids making the numerical tails depend on the chosen physical units.
    force = force/physical_sf[:,None]
    acceleration = acceleration/physical_sa[:,None]
    sf,sa = np.ones(size),np.ones(size)
    truth = np.asarray(truth,dtype=float)
    if truth.shape not in ((),(size,)) or np.any(~np.isfinite(truth)) or np.any(truth <= 0):
        raise ValueError('truth must be positive finite scalar or (n,) array')
    if (not isinstance(order,(int,np.integer)) or isinstance(order,bool) or order < 24
            or not isinstance(calibration_order,(int,np.integer))
            or isinstance(calibration_order,bool) or calibration_order < 8):
        raise ValueError('order must be integer >=24 and calibration_order integer >=8')
    nuf,nua = _degrees(degrees)
    kf,ka = nuf/2,nua/2
    bf,ba = kf*sf**2,ka*sa**2
    if np.any(~np.isfinite(bf)) or np.any(~np.isfinite(ba)) or np.any(bf<=0) or np.any(ba<=0):
        raise ValueError('Calibration precision rates must be positive and finite')
    p,q = np.sum(force**2,axis=1),np.sum(acceleration**2,axis=1)
    dot = np.sum(force*acceleration,axis=1)
    sp,sq = p/sf**2,q/sa**2
    if np.any(~np.isfinite(sp)) or np.any(~np.isfinite(sq)):
        raise ValueError('Standardized squared reading norms must be finite')
    centre = .5*(np.log1p(sp)-np.log1p(sq))+np.log(sf)-np.log(sa)
    scale = np.sqrt(1/(sp+1)+1/(sq+1))
    fixed = np.array([-_LIMIT,-20.,-12.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,12.,20.,_LIMIT])
    adaptive = np.clip(centre[:,None]+scale[:,None]*np.array([-8.,-4.,-2.,0.,2.,4.,8.]),
                       -_LIMIT+1e-8,_LIMIT-1e-8)
    edges = np.sort(np.concatenate((np.broadcast_to(fixed,(size,len(fixed))),adaptive),axis=1),axis=1)
    half = np.diff(edges,axis=1)/2
    midpoint = (edges[:,1:]+edges[:,:-1])/2
    x,weights,vander,projection = _quadrature(order)
    z = midpoint[:,:,None]+half[:,:,None]*x
    y = np.broadcast_to(np.log(truth),(size,))-mass_shift
    log_density = np.full_like(z,-np.inf)
    log_first_f = np.full_like(z,-np.inf)
    log_first_a = np.full_like(z,-np.inf)
    log_left = np.full(size,-np.inf)
    log_right = np.full(size,-np.inf)
    log_f_tail = np.full(size,-np.inf)
    log_a_tail = np.full(size,-np.inf)
    log_truth = np.full(size,-np.inf)
    nodes,beta_weights = _beta_nodes(kf+1,ka+1,calibration_order)
    power = kf+ka+2
    first_power = kf+ka+1.5
    # One calibration-ratio node at a time keeps memory independent of
    # calibration_order, with no samples x mass x scale-ratio tensor.
    for u,weight in zip(nodes,beta_weights):
        if weight == 0:
            continue
        log_w = math.log(weight)
        pp,qq = u*p/bf,(1-u)*q/ba
        dd = math.sqrt(u*(1-u))*dot/np.sqrt(bf*ba)
        log_r = .5*(math.log1p(-u)-math.log(u)+np.log(bf)-np.log(ba))
        standardized_z = z-log_r[:,None,None]
        theta = np.arctan(np.exp(standardized_z))
        sn,cs = np.sin(theta),np.cos(theta)
        h2 = np.maximum(0,pp[:,None,None]*sn*sn+qq[:,None,None]*cs*cs+2*dd[:,None,None]*sn*cs)
        residual = np.maximum(0,pp[:,None,None]*cs*cs+qq[:,None,None]*sn*sn-2*dd[:,None,None]*sn*cs)
        log_jac = _log_jacobian(standardized_z)
        log_density = np.logaddexp(log_density,log_w+_log_radial(h2,residual,power)+log_jac)
        first = log_w-first_power*np.log1p(residual/2)+log_jac
        log_first_f = np.logaddexp(log_first_f,first+.5*(np.log(bf)-math.log(u))[:,None,None]+np.log(sn))
        log_first_a = np.logaddexp(log_first_a,first+.5*(np.log(ba)-math.log1p(-u))[:,None,None]+np.log(cs))
        log_left = np.logaddexp(log_left,log_w+_log_radial(qq,pp,power)-log_r)
        log_right = np.logaddexp(log_right,log_w+_log_radial(pp,qq,power)+log_r)
        log_f_tail = np.logaddexp(log_f_tail,log_w+.5*(np.log(bf)-math.log(u))-first_power*np.log1p(qq/2)+log_r)
        log_a_tail = np.logaddexp(log_a_tail,log_w+.5*(np.log(ba)-math.log1p(-u))-first_power*np.log1p(pp/2)-log_r)
        yt = np.arctan(np.exp(y-log_r))
        ys,yc = np.sin(yt),np.cos(yt)
        yh = np.maximum(0,pp*ys*ys+qq*yc*yc+2*dd*ys*yc)
        yr = np.maximum(0,pp*yc*yc+qq*ys*ys-2*dd*ys*yc)
        log_truth = np.logaddexp(log_truth,log_w+_log_radial(yh,yr,power)+_log_jacobian(y-log_r))
    shift = np.max(log_density,axis=(1,2))
    density = np.exp(log_density-shift[:,None,None])
    summary,_ = _summary(density,edges,z,half,weights,vander,projection,
                         np.exp(log_left-shift),np.exp(log_right-shift),y,log_truth-shift)
    first_shift = np.maximum(np.max(log_first_f,axis=(1,2)),np.max(log_first_a,axis=(1,2)))
    numerator = np.sum(np.sum(np.exp(log_first_f-first_shift[:,None,None])*weights,axis=2)*half,axis=1)
    denominator = np.sum(np.sum(np.exp(log_first_a-first_shift[:,None,None])*weights,axis=2)*half,axis=1)
    numerator += np.exp(log_f_tail-first_shift-_LIMIT)
    denominator += np.exp(log_a_tail-first_shift-_LIMIT)
    point = np.exp(np.log(numerator)-np.log(denominator)+mass_shift)
    if np.any(~np.isfinite(point)) or np.any(point <= 0):
        raise FloatingPointError('Calibrated physical first-moment ratio failed')
    # Exact channel exchange symmetry, including arbitrary relative units.
    equal = (sp == sq) & (nuf == nua)
    summary['log_median'] = np.where(equal,0.,summary['log_median'])
    for name in summary:
        if name in ('mean_log_mass','log_median') or name.startswith(('log_lower_','log_upper_')):
            summary[name] += mass_shift
    return {'points': {'calibrated_joint': point,
                       'calibrated_geometric': np.exp(summary['mean_log_mass']),
                       'calibrated_median': np.exp(summary['log_median'])},
            'distributions': {'calibrated_joint': summary}}
