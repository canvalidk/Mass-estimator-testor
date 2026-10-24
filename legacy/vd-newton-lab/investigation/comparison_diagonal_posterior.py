"""The existing flat/tube laws with independent unequal coordinate errors.

Only scalar channel rescaling is used; separate coordinate whitening would
change the common-direction model and its reference measure. In RMS-scaled
channels write (f,a)=r(sin(theta),cos(theta)), x=r*u. Then
df da dOmega = dtheta d^3x/||x||. At fixed theta the likelihood in x is a
diagonal Gaussian. Its inverse-radius and radius moments reduce the angular
integral to one smooth Laplace integral, with no direction grid.

This is the same physical flat df da dOmega law as estimator.py. The tube
sensitivity law uses (f/sigma_F_RMS)^2+(a/sigma_a_RMS)^2 as its extra weight.
"""
from functools import lru_cache
import math

import numpy as np

from comparison_posterior import _LIMIT, _quadrature, _summary, infer_batch as isotropic_infer


@lru_cache(None)
def _laplace_nodes(order):
    x, w = np.polynomial.legendre.leggauss(order)
    phi = (x + 1) * (math.pi / 4)
    sine2, cosine = np.sin(phi)**2, np.cos(phi)
    return sine2, cosine**2, w * (math.pi / 4) * cosine


def gaussian_radius_moments(mean, variance, order=64):
    """E[1/R], E[R], E[R²] for diagonal Gaussian vectors in dimension 3.

    For L(t)=E exp(-t R²), E[1/R]=int t^-1/2 L(t) dt / sqrt(pi).
    E[R] uses the same integral times the tilted Gaussian E[R²], avoiding
    subtracting nearly equal numbers. With scale²=E[R²]/3, use
    t=tan(phi)²/(2 scale²), phi in [0,pi/2]. The angular mapping resolves
    the small-variance endpoint more accurately than an unpanelled compact
    linear mapping when coordinate SDs differ substantially.
    """
    mean, variance = np.broadcast_arrays(np.asarray(mean, float), np.asarray(variance, float))
    if (mean.ndim == 0 or mean.shape[-1] != 3 or np.any(variance <= 0)
            or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(variance))):
        raise ValueError("Three finite means and positive diagonal variances are required")
    if not isinstance(order, int) or order < 8:
        raise ValueError("Radial quadrature order must be an integer >=8")
    radius2 = np.sum(mean * mean + variance, axis=-1)
    if not np.all(np.isfinite(radius2)):
        raise ValueError("Gaussian second radius moment must be finite")
    scale2 = radius2 / 3
    c, mean2 = variance / scale2[..., None], mean * mean
    scaled_mean2 = mean2 / scale2[..., None]
    inverse, radius = np.zeros_like(radius2), np.zeros_like(radius2)
    sine2, cosine2, weights = _laplace_nodes(order)
    # Loop over quadrature nodes rather than materializing a fourth tensor
    # dimension for all trials, mass nodes, radial nodes and coordinates.
    for u2, v, weight in zip(sine2, cosine2, weights):
        denominator = v + c * u2
        log_kernel = -.5 * (np.sum(np.log(denominator), axis=-1)
                            + u2 * np.sum(scaled_mean2 / denominator, axis=-1))
        weighted = weight * np.exp(log_kernel)
        tilted_radius2 = np.sum(variance * v / denominator + mean2 * v*v / denominator**2, axis=-1)
        inverse += weighted
        radius += weighted * tilted_radius2
    factor = math.sqrt(2 / math.pi) / np.sqrt(scale2)
    return inverse * factor, radius * factor, radius2


def _kernel(force, acceleration, vf, va, theta, radial_order):
    sn, cs = np.sin(theta)[..., None], np.cos(theta)[..., None]
    if np.ndim(vf) == 2:
        vf = vf[:, None, :]
    if np.ndim(va) == 2:
        va = va[:, None, :]
    precision = sn * sn / vf + cs * cs / va
    variance = 1 / precision
    h = sn * force[:, None, :] / vf + cs * acceleration[:, None, :] / va
    mean = variance * h
    log_gaussian = .5 * np.sum(h * mean + np.log(variance), axis=-1)
    inverse, radius, radius2 = gaussian_radius_moments(mean, variance, radial_order)
    return log_gaussian, inverse, radius, radius2


def _rescale(result, mass_scale):
    for name in result['points']:
        result['points'][name] *= mass_scale
    shift = np.log(mass_scale)
    for distribution in result['distributions'].values():
        for key in distribution:
            if key in ('mean_log_mass', 'log_median') or key.startswith(('log_lower_', 'log_upper_')):
                distribution[key] += shift
    return result


def _infer_chunk(force, acceleration, vf, va, truth, order):
    # The profile location only places integration panels; no simulation truth
    # or true direction enters either the panel selection or the density.
    from comparison_diagonal_baselines import covariance_profile
    profile = covariance_profile(force, acceleration, np.sqrt(vf), np.sqrt(va))
    p, q = np.sum(force**2 / vf, axis=1), np.sum(acceleration**2 / va, axis=1)
    ordinary_centre = .5 * np.log((np.sum(force**2, axis=1) + 1) / (np.sum(acceleration**2, axis=1) + 1))
    centre = ordinary_centre.copy()
    good = np.isfinite(profile) & (profile > 0)
    centre[good] = np.log(profile[good])
    scale = np.sqrt(1 / (p + 1) + 1 / (q + 1))
    fixed = np.r_[-_LIMIT, -20., -12., np.arange(-8., 9.), 12., 20., _LIMIT]
    adaptive = centre[:, None] + scale[:, None] * np.array([-8., -4., -2., 0., 2., 4., 8.])
    adaptive = np.clip(adaptive, -_LIMIT + 1e-8, _LIMIT - 1e-8)
    edges = np.sort(np.concatenate((np.broadcast_to(fixed, (len(force), len(fixed))), adaptive), axis=1), axis=1)
    half = np.diff(edges, axis=1) / 2
    midpoint = (edges[:, 1:] + edges[:, :-1]) / 2
    x, weights, vander, projection = _quadrature(order)
    z = midpoint[:, :, None] + half[:, :, None] * x
    theta = np.arctan(np.exp(z))
    radial_order = max(32, order // 2)
    flat_shape = (len(force), -1)
    kernel = _kernel(force, acceleration, vf, va, theta.reshape(flat_shape), radial_order)
    kernel = [value.reshape(z.shape) for value in kernel]
    endpoints = _kernel(force, acceleration, vf, va,
                        np.broadcast_to([0., math.pi/2], (len(force), 2)), radial_order)
    peak = np.maximum(np.max(kernel[0], axis=(1, 2)), np.max(endpoints[0], axis=1))
    gaussian = np.exp(kernel[0] - peak[:, None, None])
    gaussian_end = np.exp(endpoints[0] - peak[:, None])
    sn, cs = np.sin(theta), np.cos(theta)
    jacobian = 1 / (2 * np.cosh(z))
    y = np.log(truth)
    at_truth = _kernel(force, acceleration, vf, va,
                       np.broadcast_to(np.arctan(truth), (len(force),))[:, None], radial_order)
    log_jacobian = -abs(y) - np.log1p(np.exp(-2 * abs(y)))
    points, distributions = {}, {}
    for index, name in enumerate(('flat_joint', 'tube_joint')):
        moment = gaussian if index == 0 else gaussian * kernel[3]
        moment_end = gaussian_end if index == 0 else gaussian_end * endpoints[3]
        numerator = np.sum(np.sum(moment * sn * jacobian * weights, axis=2) * half, axis=1)
        denominator = np.sum(np.sum(moment * cs * jacobian * weights, axis=2) * half, axis=1)
        numerator += moment_end[:, 1] * math.exp(-_LIMIT)
        denominator += moment_end[:, 0] * math.exp(-_LIMIT)
        points[name] = numerator / denominator
        density = gaussian * kernel[index + 1] * jacobian
        endpoint_density = gaussian_end * endpoints[index + 1]
        truth_log_kernel = (at_truth[0][:, 0] - peak + np.log(at_truth[index + 1][:, 0]) + log_jacobian)
        summary, root_point = _summary(density, edges, z, half, weights, vander, projection,
                                      endpoint_density[:, 0], endpoint_density[:, 1], y, truth_log_kernel)
        # In adverse symmetric data the CDF can round to 1/2 throughout a
        # deep gap between modes. Use the exact exchange-symmetry median only
        # where the covariance and squared coordinates warrant that identity.
        row_vf, row_va = np.broadcast_to(vf, force.shape), np.broadcast_to(va, acceleration.shape)
        same_covariance = np.all(row_vf == row_va, axis=1)
        same_coordinates = np.all(force*force == acceleration*acceleration, axis=1)
        isotropic = np.all(row_vf == row_vf[:, :1], axis=1)
        same_norm = np.sum(force*force, axis=1) == np.sum(acceleration*acceleration, axis=1)
        symmetric = same_covariance & (same_coordinates | (isotropic & same_norm))
        summary['log_median'] = np.where(symmetric, 0., summary['log_median'])
        distributions[name] = summary
        if index == 0:
            points['flat_geometric'] = np.exp(summary['mean_log_mass'])
            points['flat_median'] = np.exp(summary['log_median'])
            points['flat_reciprocal_root'] = root_point
    return {'points': points, 'distributions': distributions}


def infer_batch(force, acceleration, truth=1., *, force_sd, acceleration_sd, order=96):
    """Physical mass estimates and log-mass probability summaries, batch API.

    Each channel's supplied instrument SDs have shape (3,) or (N,3). The
    latter permits independent calibration estimates in each trial; the
    inference treats the supplied values as fixed, without integrating their
    calibration uncertainty. The latent direction is integrated out.
    Truth is used only in distribution scoring. Internal
    chunks cap the auxiliary quadrature memory independently of runner batches.
    """
    force, acceleration = np.atleast_2d(np.asarray(force, float)), np.atleast_2d(np.asarray(acceleration, float))
    sf, sa = np.asarray(force_sd, float), np.asarray(acceleration_sd, float)
    if force.shape != acceleration.shape or force.ndim != 2 or force.shape[1] != 3 or not len(force):
        raise ValueError('force and acceleration must have matching nonempty (N,3) shapes')
    if (sf.shape not in ((3,), force.shape) or sa.shape not in ((3,), force.shape)
            or np.any(sf <= 0) or np.any(sa <= 0)):
        raise ValueError('Supply positive SDs of shape (3,) or (N,3) per channel')
    if not all(np.isfinite(value).all() for value in (force, acceleration, sf, sa)):
        raise ValueError('Readings and SDs must be finite')
    if not math.isfinite(truth) or truth <= 0 or not isinstance(order, int) or order < 24:
        raise ValueError('Positive finite truth and integer order >=24 are required')
    if sf.ndim == 2 or sa.ndim == 2:
        sf, sa = np.broadcast_to(sf, force.shape), np.broadcast_to(sa, force.shape)
        rms_f = np.array([math.hypot(*row)/math.sqrt(3) for row in sf])
        rms_a = np.array([math.hypot(*row)/math.sqrt(3) for row in sa])
        mass_scale = rms_f / rms_a
        if not np.all(np.isfinite(mass_scale)) or np.any(mass_scale <= 0):
            raise ValueError('Noise scales must have a representable positive ratio')
        force, acceleration = force/rms_f[:, None], acceleration/rms_a[:, None]
        vf, va = (sf/rms_f[:, None])**2, (sa/rms_a[:, None])**2
        if np.all(sf == sf[:, :1]) and np.all(sa == sa[:, :1]):
            return _rescale(isotropic_infer(force, acceleration, truth/mass_scale, order=order), mass_scale)
        chunks = [_infer_chunk(force[i:i+16], acceleration[i:i+16], vf[i:i+16], va[i:i+16],
                               truth/mass_scale[i:i+16], order)
                  for i in range(0, len(force), 16)]
        return _rescale(_join_chunks(chunks), mass_scale)
    rms_f, rms_a = math.hypot(*sf)/math.sqrt(3), math.hypot(*sa)/math.sqrt(3)
    mass_scale = rms_f / rms_a
    force, acceleration = force/rms_f, acceleration/rms_a
    vf, va = (sf/rms_f)**2, (sa/rms_a)**2
    if np.all(sf == sf[0]) and np.all(sa == sa[0]):
        return _rescale(isotropic_infer(force, acceleration, truth/mass_scale, order=order), mass_scale)
    chunks = [_infer_chunk(force[i:i+16], acceleration[i:i+16], vf, va, truth/mass_scale, order)
              for i in range(0, len(force), 16)]
    return _rescale(_join_chunks(chunks), mass_scale)


def _join_chunks(chunks):
    return {'points': {name: np.concatenate([chunk['points'][name] for chunk in chunks])
                         for name in chunks[0]['points']},
              'distributions': {name: {key: np.concatenate([chunk['distributions'][name][key] for chunk in chunks])
                                        for key in distribution}
                                for name, distribution in chunks[0]['distributions'].items()}}
