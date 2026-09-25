"""The per-reading law: evidence about mass from one force/acceleration reading.

Model for one reading, in d = 1, 2 or 3 dimensions:

    observed force        F = f u + e_F,   e_F ~ N(0, sigma_F^2 I_d)
    observed acceleration a = alpha u + e_a, e_a ~ N(0, sigma_a^2 I_d)

with latent positive magnitudes f, alpha and a common unit direction u.
Newton II with m > 0 is assumed: m = f / alpha.

The flat reference measure is  df d(alpha) dOmega(u)  (flat in the two
physical magnitudes, uniform in direction). In standardized polar
coordinates  f / sigma_F = r sin(theta),  alpha / sigma_a = r cos(theta):

    df d(alpha) = sigma_F sigma_a r dr d(theta),    m = s tan(theta),
    s = sigma_F / sigma_a.

Integrating the direction and the radius r at fixed theta gives the radial
kernel K_d(h), with h = |x sin(theta) + y cos(theta)| and x = F / sigma_F,
y = a / sigma_a:

    K_d(h) = integral_0^inf r exp(-r^2/2) < exp(r h . u) >_u dr

    d = 1:  1 + sqrt(pi/2) h exp(h^2/2) erf(h/sqrt2)
    d = 2:  exp(h^2/2)
    d = 3:  sqrt(pi/2) exp(h^2/2) erf(h/sqrt2) / h

and the conditional mean radius  E[r | theta] = M_d(h) / K_d(h),  with
M_d = integral r^2 exp(-r^2/2) < exp(r h . u) > dr:

    d = 1:  (1 + h^2) sqrt(pi/2) / (exp(-h^2/2) + sqrt(pi/2) h erf(h/sqrt2))
    d = 2:  sqrt(pi/2) [(1 + h^2/2) I0e(h^2/4) + (h^2/2) I1e(h^2/4)]   (Rice mean)
    d = 3:  h / erf(h/sqrt2)

The same single-reading law can be split into "likelihood of mass" times
"reference factor on mass" in more than one way, according to which latent
quantity is treated as the reading's nuisance. The product is the same for
one reading; the split matters only when readings are combined.

    nuisance     nuisance measure        log_like(u)             log_ref(u), u = log m
    radius       r dr dOmega             log K                   log(sech(z)/2), z = u - log s
    acceleration alpha d(alpha) dOmega   log K + 2 log cos       u      (flat in m)
    force        f df dOmega             log K + 2 log sin       -u     (flat in 1/m)

(all up to additive constants that do not depend on u).

The cartesian reference (24 September, excitation_weighting_premise_2026-09-24)
------------------------------------------------------------------------------
Weight the true pair flat as a VECTOR: d^d v d(theta), with v = r u the pair in
noise units (r^2 = (f/sigma_F)^2 + (alpha/sigma_a)^2) and theta uniform. Against
the flat reference,

    d^d v d(theta) = r^(d-1) dr dOmega d(theta) = r^(d-2) df d(alpha) dOmega / (sigma_F sigma_a),

so the cartesian law is the flat law reweighted by r^(d-2): identical in 2D,
x r in 3D, x 1/r in 1D. The factor splits as r^(d-2) = alpha^(d-2) x
(sec(theta) / sigma_a)^(d-2): the excitation part alpha^(d-2) that the premise
forces (flat in the true acceleration vector at fixed mass), and a factor in
m alone that keeps the single-reading law uniform in theta. The radial
integral is now a full Gaussian integral over R^d:

    K_d(h) = integral_0^inf r^(d-1) exp(-r^2/2) < exp(r h . u) >_u dr  ~  exp(h^2/2)   (every d)

and E[r | theta] is the mean length of N(w, I_d) with |w| = h (noncentral chi mean):

    d = 1:  sqrt(2/pi) exp(-h^2/2) + h erf(h/sqrt2)                 (folded normal)
    d = 2:  the Rice mean (as for the flat reference)
    d = 3:  sqrt(2/pi) exp(-h^2/2) + (h + 1/h) erf(h/sqrt2)

The nuisance splits keep their names; the nuisance measures become the
d-dimensional volumes (radius: d^d v; acceleration: d^d a* = alpha^(d-1) d(alpha)
dOmega; force: d^d F*), so the powers 2 become d:

    nuisance     log_like(u)             log_ref(u)
    radius       h^2/2                   log(sech(z)/2)             (uniform angle)
    acceleration h^2/2 + d log cos       u + (2 - d) log cos
    force        h^2/2 + d log sin       -u + (2 - d) log sin

In every dimension, the radius-nuisance likelihood is exactly
-|F - m a|^2 / (2 (sigma_F^2 + m^2 sigma_a^2)) plus a constant in m.
"""

import math

import numpy as np
from scipy.special import erf, ive

SQRT_HALF_PI = math.sqrt(math.pi / 2)
SQRT2 = math.sqrt(2.0)
SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)
DIMENSIONS = (1, 2, 3)
NUISANCES = ("radius", "acceleration", "force")
REFERENCES = ("flat", "cartesian")
_SMALL_H = 1e-6


def log_radial_kernel(h, d):
    """log K_d(h) for h >= 0 (any array shape)."""
    h = np.abs(np.asarray(h, dtype=float))
    half_h2 = 0.5 * h * h
    if d == 2:
        return half_h2
    if d == 1:
        return half_h2 + np.log(np.exp(-half_h2) + SQRT_HALF_PI * h * erf(h / SQRT2))
    if d == 3:
        safe = np.where(h > _SMALL_H, h, 1.0)
        ratio = np.where(h > _SMALL_H, SQRT_HALF_PI * erf(safe / SQRT2) / safe, 1.0 - h * h / 6.0)
        return half_h2 + np.log(ratio)
    raise ValueError(f"dimension must be one of {DIMENSIONS}")


def mean_radius(h, d):
    """E[r | theta] = M_d(h) / K_d(h) for h >= 0."""
    h = np.abs(np.asarray(h, dtype=float))
    h2 = h * h
    if d == 1:
        return (1.0 + h2) * SQRT_HALF_PI / (np.exp(-0.5 * h2) + SQRT_HALF_PI * h * erf(h / SQRT2))
    if d == 2:
        x = 0.25 * h2
        return SQRT_HALF_PI * ((1.0 + 0.5 * h2) * ive(0, x) + 0.5 * h2 * ive(1, x))
    if d == 3:
        safe = np.where(h > _SMALL_H, h, 1.0)
        return np.where(h > _SMALL_H, safe / erf(safe / SQRT2), SQRT_HALF_PI * (1.0 + h2 / 6.0))
    raise ValueError(f"dimension must be one of {DIMENSIONS}")


def cartesian_mean_radius(h, d):
    """E[r | theta] under the cartesian reference: the mean of |N(w, I_d)| with |w| = h."""
    h = np.abs(np.asarray(h, dtype=float))
    h2 = h * h
    if d == 1:
        return SQRT_2_OVER_PI * np.exp(-0.5 * h2) + h * erf(h / SQRT2)
    if d == 2:
        return mean_radius(h, 2)
    if d == 3:
        safe = np.where(h > _SMALL_H, h, 1.0)
        big = SQRT_2_OVER_PI * np.exp(-0.5 * h2) + (safe + 1.0 / safe) * erf(safe / SQRT2)
        return np.where(h > _SMALL_H, big, 2.0 * SQRT_2_OVER_PI * (1.0 + h2 / 6.0))
    raise ValueError(f"dimension must be one of {DIMENSIONS}")


def radial_terms(h, d, reference):
    """(log K, E[r | theta]) for the named reference measure."""
    if reference == "flat":
        return log_radial_kernel(h, d), mean_radius(h, d)
    if reference == "cartesian":
        if d not in DIMENSIONS:
            raise ValueError(f"dimension must be one of {DIMENSIONS}")
        h = np.abs(np.asarray(h, dtype=float))
        return 0.5 * h * h, cartesian_mean_radius(h, d)
    raise ValueError(f"reference must be one of {REFERENCES}")


def nuisance_power(d, reference):
    """The power of the latent magnitude in the acceleration/force nuisance measures:
    2 for the flat reference (alpha d(alpha)), d for the cartesian one (d^d a*)."""
    return 2 if reference == "flat" else d


def standardize(force, acceleration, force_sd, acceleration_sd):
    """Standardized readings and their summary statistics.

    force, acceleration: arrays (..., d). force_sd, acceleration_sd: arrays (...),
    per-coordinate standard deviations the analyst was given.
    Returns P = |x|^2, Q = |y|^2, D = x.y, s = sigma_F/sigma_a, all shape (...).
    """
    force = np.asarray(force, dtype=float)
    acceleration = np.asarray(acceleration, dtype=float)
    sf = np.asarray(force_sd, dtype=float)
    sa = np.asarray(acceleration_sd, dtype=float)
    x = force / sf[..., None]
    y = acceleration / sa[..., None]
    return (np.sum(x * x, axis=-1), np.sum(y * y, axis=-1), np.sum(x * y, axis=-1), sf / sa)


class ReadingSet:
    """A batch of series of readings, as the analyst sees them.

    force, acceleration: (B, N, d). force_sd, acceleration_sd: (B, N).
    """

    def __init__(self, force, acceleration, force_sd, acceleration_sd):
        self.force = np.asarray(force, dtype=float)
        self.acceleration = np.asarray(acceleration, dtype=float)
        if self.force.ndim != 3 or self.force.shape != self.acceleration.shape:
            raise ValueError("readings must have shape (series, readings, dimension)")
        shape = self.force.shape[:2]
        self.force_sd = np.broadcast_to(np.asarray(force_sd, dtype=float), shape).copy()
        self.acceleration_sd = np.broadcast_to(np.asarray(acceleration_sd, dtype=float), shape).copy()
        if not (np.all(self.force_sd > 0) and np.all(self.acceleration_sd > 0)):
            raise ValueError("supplied standard deviations must be positive")
        self.dimension = self.force.shape[2]
        if self.dimension not in DIMENSIONS:
            raise ValueError(f"dimension must be one of {DIMENSIONS}")
        # Per-series quantities an analyst stage computed and later stages reuse
        # (sliced along with the series by subset()).
        self.series_data = {}
        self.P, self.Q, self.D, self.s = standardize(self.force, self.acceleration,
                                                     self.force_sd, self.acceleration_sd)

    def subset(self, index):
        """The series selected by `index` (a slice or index array)."""
        out = ReadingSet(self.force[index], self.acceleration[index],
                         self.force_sd[index], self.acceleration_sd[index])
        out.series_data = {k: v[index] for k, v in self.series_data.items()}
        return out

    def take_readings(self, index):
        """The same series, keeping only the readings selected by `index` (a slice)."""
        out = ReadingSet(self.force[:, index], self.acceleration[:, index],
                         self.force_sd[:, index], self.acceleration_sd[:, index])
        out.series_data = dict(self.series_data)
        return out

    @property
    def series(self):
        return self.force.shape[0]

    @property
    def readings(self):
        return self.force.shape[1]

    def reading_terms(self, u, nuisance, *, reference):
        """Per-reading curves on log-mass grids u of shape (B, G).

        reference: 'flat' (21 Sept, df d(alpha) dOmega) or 'cartesian' (24 Sept,
        flat in the vector). Returns (log_like, log_ref, cond_alpha), each (B, N, G):
          log_like   log-likelihood of mass with this reading's nuisance integrated,
          log_ref    the mass factor the single-reading reference implies (density in u),
          cond_alpha E[alpha | m] for this reading.
        """
        if nuisance not in NUISANCES:
            raise ValueError(f"nuisance must be one of {NUISANCES}")
        if reference not in REFERENCES:
            raise ValueError(f"reference must be one of {REFERENCES}")
        u = np.asarray(u, dtype=float)[:, None, :]
        z = u - np.log(self.s)[:, :, None]
        log_sin = -0.5 * np.logaddexp(0.0, -2.0 * z)
        log_cos = -0.5 * np.logaddexp(0.0, 2.0 * z)
        sn, cs = np.exp(log_sin), np.exp(log_cos)
        P, Q, D = (v[:, :, None] for v in (self.P, self.Q, self.D))
        h = np.sqrt(np.maximum(0.0, P * sn * sn + Q * cs * cs + 2.0 * D * sn * cs))
        log_k, radius = radial_terms(h, self.dimension, reference)
        cond_alpha = self.acceleration_sd[:, :, None] * cs * radius
        k = nuisance_power(self.dimension, reference)
        if nuisance == "radius":
            log_like, log_ref = log_k, -np.logaddexp(z, -z)
        elif nuisance == "acceleration":
            log_like, log_ref = log_k + k * log_cos, u + (2 - k) * log_cos
        else:
            log_like, log_ref = log_k + k * log_sin, -u + (2 - k) * log_sin
        log_ref = np.broadcast_to(log_ref, log_k.shape)
        return log_like, log_ref, cond_alpha

    def pooled(self):
        """Inverse-variance mean of each series' readings, as one reading per series.

        Exact sufficient reduction when every reading of a series measures the
        same latent vector pair with independent Gaussian errors (contribution 11,
        eq. 5). The channels are independent, so each is pooled separately.
        """
        wf = 1.0 / self.force_sd**2
        wa = 1.0 / self.acceleration_sd**2
        force = np.sum(wf[:, :, None] * self.force, axis=1) / np.sum(wf, axis=1)[:, None]
        acceleration = np.sum(wa[:, :, None] * self.acceleration, axis=1) / np.sum(wa, axis=1)[:, None]
        return ReadingSet(force[:, None, :], acceleration[:, None, :],
                          (1.0 / np.sqrt(np.sum(wf, axis=1)))[:, None],
                          (1.0 / np.sqrt(np.sum(wa, axis=1)))[:, None])
