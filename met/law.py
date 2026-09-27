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

Two cartesian implementations: cartesian1 and cartesian2
---------------------------------------------------------
Two sessions implemented the 24 September premise (flat in the true pair as a
vector) independently on the same day. Both are kept, under their own names,
so the investigation can check them against each other (tests/test_cartesian.py,
"cartesian1 against cartesian2"). They share the measure d^d v d(theta), the
radial kernel exp(h^2/2) and E[r | theta]; so they give the same single-reading
law, and the same combined law under the radius nuisance. They differ only in
how the acceleration and force nuisances split the single-reading law when
readings are combined: cartesian1 keeps the flat reference's split (powers 2),
cartesian2 integrates over the d-dimensional nuisance volume (powers d). Per
reading, cartesian1's likelihood is cartesian2's times cos(theta)^(2-d)
(acceleration) or sin(theta)^(2-d) (force), and its reference factor is
cartesian2's divided by the same.

cartesian1 (Line A, 665a321, 24 September 18:13)
------------------------------------------------
The cartesian1 reference replaces the flat-magnitude measure by one flat in
each component of the common latent vector v = r u (standardised), and
uniform in theta:  d^d v d(theta) = r^(d-2) x (flat measure). A d-dimensional
reading is then exactly d independent one-dimensional projection readings
sharing the mass, and the radial integral is a plain Gaussian integral:

    K(h) = exp(h^2/2)  in every dimension,   E[r | theta] = mean of a noncentral chi_d.

In 2D the two references coincide. The nuisance splits keep the same mass
factors (uniform angle, flat in m, flat in 1/m) with K replaced.

cartesian2 (Line B, 7081a44, 24 September 23:27; excitation_weighting_premise_2026-09-24)
----------------------------------------------------------------------------------------
Weight the true pair flat as a VECTOR: d^d v d(theta), with v = r u the pair in
noise units (r^2 = (f/sigma_F)^2 + (alpha/sigma_a)^2) and theta uniform. Against
the flat reference,

    d^d v d(theta) = r^(d-1) dr dOmega d(theta) = r^(d-2) df d(alpha) dOmega / (sigma_F sigma_a),

so the cartesian2 law is the flat law reweighted by r^(d-2): identical in 2D,
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

Per-axis noise (ours26, 26 September)
-------------------------------------
A reading may have a different noise SD on each axis: e_F,k ~ N(0, sigma_F,k^2),
e_a,k ~ N(0, sigma_a,k^2). Axis splitting (a vector reading is its axes read
separately) then fixes the law: each axis is a one-dimensional cartesian reading
with its own ratio s_k = sigma_F,k / sigma_a,k and angle tan(theta_k) = m / s_k, and
the radius-nuisance likelihood is the sum over axes,

    log L(m) = sum_i sum_k -(F_ik - m a_ik)^2 / (2 (sigma_F,ik^2 + m^2 sigma_a,ik^2))  (+ const),

which is the isotropic law term by term when the SDs agree. The stacked-length
readout keeps its definition, E|F*_stack| / E|a*_stack|; its weight is now a weighted
noncentral chi mean (`weighted_stacked_mean`), computed by quadrature:

    E|a*_stack| = E sqrt( sum_c sigma_a,c^2 cos^2(theta_c) v_c^2 ),   v_c ~ N(w_c, 1),

over all N d components c, with w_c = x_c sin(theta_c) + y_c cos(theta_c). What is no
longer forced is the centre of the uniform-angle prior: with several ratios s_k there
is no single instrument ratio, so 'first_reading' is refused and the centre is declared.
"""

import math

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import erf, gammaln, ive

SQRT_HALF_PI = math.sqrt(math.pi / 2)
SQRT2 = math.sqrt(2.0)
SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)
DIMENSIONS = (1, 2, 3)
NUISANCES = ("radius", "acceleration", "force")
REFERENCES = ("flat", "cartesian1", "cartesian2")
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


def cartesian1_mean_radius(h, d):
    """E|v| for v ~ N(w, I_d) with |w| = h: the mean of a noncentral chi variable.

    Under the cartesian1 reference the latent vector v is Gaussian around w at
    fixed theta, so this is E[r | theta].
    """
    h = np.abs(np.asarray(h, dtype=float))
    if d == 1:
        return math.sqrt(2 / math.pi) * np.exp(-0.5 * h * h) + h * erf(h / SQRT2)
    if d == 2:
        return mean_radius(h, 2)
    if d == 3:
        safe = np.where(h > _SMALL_H, h, 1.0)
        return np.where(h > _SMALL_H,
                        math.sqrt(2 / math.pi) * np.exp(-0.5 * h * h) + (safe + 1.0 / safe) * erf(safe / SQRT2),
                        2 * math.sqrt(2 / math.pi) * (1.0 + h * h / 6.0))
    raise ValueError(f"dimension must be one of {DIMENSIONS}")


def cartesian2_mean_radius(h, d):
    """E[r | theta] under the cartesian2 reference: the mean of |N(w, I_d)| with |w| = h."""
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


_LOG_SQRT2 = 0.5 * math.log(2.0)
_STACKED_TABLES = {}   # k -> (h_max, spline in asinh(h))


def _noncentral_chi_mean_series(k, h):
    """E|v| for v ~ N(w, I_k), |w| = h, by the Poisson mixture of central chi means:
    |v|^2 is chi^2 with k + 2J degrees of freedom, J ~ Poisson(h^2/2). Exact to rounding
    (about 1e-13 relative, checked against arbitrary precision for k up to 1200)."""
    out = np.empty(len(h))
    for i, hh in enumerate(h):
        lam = 0.5 * hh * hh
        if lam == 0.0:
            out[i] = math.exp(_LOG_SQRT2 + gammaln((k + 1) / 2) - gammaln(k / 2))
            continue
        sd = math.sqrt(lam)
        j = np.arange(max(0, int(lam - 12 * sd - 30)), int(lam + 12 * sd + 30) + 1, dtype=float)
        n = k + 2 * j
        log_terms = j * math.log(lam) - lam - gammaln(j + 1) + _LOG_SQRT2 + gammaln((n + 1) / 2) - gammaln(n / 2)
        out[i] = np.sum(np.exp(log_terms))
    return out


def stacked_mean_radius(h, k):
    """E|v| for v ~ N(w, I_k) with |w| = h: the noncentral chi mean in k dimensions.

    k = 1, 2, 3 use the closed forms; larger k (the stacked vector of N readings
    in d dimensions, k = N d) use a cubic spline in asinh(h) over the exact series,
    built once per k and extended when a larger h is asked for (about 5e-9 relative).
    """
    h = np.abs(np.asarray(h, dtype=float))
    if k in DIMENSIONS:
        return cartesian2_mean_radius(h, k)
    if k < 1 or int(k) != k:
        raise ValueError("k must be a positive integer")
    k = int(k)
    need = float(np.max(h)) if h.size else 0.0
    table = _STACKED_TABLES.get(k)
    if table is None or need > table[0]:
        h_max = max(64.0, 2.0 * need, 2.0 * table[0] if table else 0.0)
        a = np.linspace(0.0, math.asinh(h_max), int(600 * math.asinh(h_max)) + 2)
        table = (h_max, CubicSpline(a, _noncentral_chi_mean_series(k, np.sinh(a))))
        _STACKED_TABLES[k] = table
    return table[1](np.arcsinh(h))


_INV_2SQRTPI = 1.0 / (2.0 * math.sqrt(math.pi))
WEIGHTED_STEP = 0.4     # quadrature step in log t
WEIGHTED_REACH = 36.0   # nodes run this far either side of -log E[X], in log t
WEIGHTED_SPACING = 0.02  # ours26's stacked weight: spline nodes at most this far apart in log m
WEIGHTED_MIN_POINTS = 41


def weighted_stacked_mean(c2, w2, dof):
    """E sqrt(X), X = sum_g c2_g Y_g, Y_g a noncentral chi-square with dof_g degrees of
    freedom and noncentrality w2_g, independent (the last axis runs over g).

    Components sharing one noise pair (sigma_F, sigma_a) share their weight c2 and are
    grouped: their squared pushes add to one noncentral chi-square. The mean is taken
    through the Laplace transform of X, which is closed:

        sqrt(X) = (1 / (2 sqrt(pi))) int_0^inf (1 - exp(-t X)) t^(-3/2) dt,
        E exp(-t X) = prod_g (1 + 2 t c2_g)^(-dof_g/2) exp(-t c2_g w2_g / (1 + 2 t c2_g)),

    integrated by the trapezoid rule in log t (the integrand is analytic there and
    decays like exp(-|log t|/2)), centred on t = 1 / E[X], with the two tails added in
    closed form to first order. About 1e-10 relative, checked against the closed forms
    (one group) and against direct double integration (two groups).
    """
    c2 = np.asarray(c2, dtype=float)
    w2 = np.asarray(w2, dtype=float)
    dof = np.broadcast_to(np.asarray(dof, dtype=float), c2.shape)
    mean_x = np.sum(c2 * (dof + w2), axis=-1)
    tau0 = -np.log(np.where(mean_x > 0, mean_x, 1.0))
    last = int(WEIGHTED_REACH / WEIGHTED_STEP)
    total = np.zeros(mean_x.shape)
    for j in range(-last, last + 1):
        tau = tau0 + j * WEIGHTED_STEP
        t = np.exp(tau)[..., None]
        a = 2.0 * t * c2
        log_l = np.sum(-0.5 * dof * np.log1p(a) - t * c2 * w2 / (1.0 + a), axis=-1)
        term = -np.expm1(log_l) * np.exp(-0.5 * tau)
        total += 0.5 * term if abs(j) == last else term
    total *= WEIGHTED_STEP
    lo, hi = tau0 - last * WEIGHTED_STEP, tau0 + last * WEIGHTED_STEP
    total += 2.0 * mean_x * np.exp(0.5 * lo) + 2.0 * np.exp(-0.5 * hi)
    return np.where(mean_x > 0, _INV_2SQRTPI * total, 0.0)


def radial_terms(h, d, reference):
    """(log K, E[r | theta]) for the named reference measure."""
    if reference == "flat":
        return log_radial_kernel(h, d), mean_radius(h, d)
    if reference == "cartesian1":
        return 0.5 * h * h, cartesian1_mean_radius(h, d)
    if reference == "cartesian2":
        if d not in DIMENSIONS:
            raise ValueError(f"dimension must be one of {DIMENSIONS}")
        h = np.abs(np.asarray(h, dtype=float))
        return 0.5 * h * h, cartesian2_mean_radius(h, d)
    raise ValueError(reference_error(reference))


def reference_error(reference):
    """The refusal for an unknown reference; the retired name says why it was split."""
    if reference == "cartesian":
        return ("reference 'cartesian' was split on 2026-09-25 into two independent implementations, "
                "cartesian1 (665a321) and cartesian2 (7081a44); name the one you mean")
    return f"reference must be one of {REFERENCES}"


def nuisance_power(d, reference):
    """The power of the latent magnitude in the acceleration/force nuisance measures:
    2 for flat and cartesian1 (alpha d(alpha)), d for cartesian2 (d^d a*)."""
    return d if reference == "cartesian2" else 2


class StackedState:
    """The exact carried state of ours25b: three running sums per series.

    With one noise ratio s = sigma_F / sigma_a (and one sigma_a) for the whole series,
    and x = F / sigma_F, y = a / sigma_a, each reading contributes

        c_i = ((Q_i - P_i) + 2i D_i) / 4,   T_i = P_i + Q_i,   d components,

    (P = |x|^2, Q = |y|^2, D = x.y), and the series is carried by their sums C, T, M.
    Everything ours25b computes depends on the data only through

        H(theta)^2 = T/2 + 2 Re[C exp(-2i theta)],     m = s tan(theta):

        law      log L(theta) = H^2 / 2   (+ const; the prior is added once, elsewhere)
        readout  E[|a*| | theta] = sigma_a cos(theta) mu_M(H),  mu_M = stacked_mean_radius.

    Combining two sets of readings is adding their states (`add`), in any order and
    grouping. This is the combination equation, exact in every dimension, under either
    cartesian reference with the radius nuisance and the prior counted once.
    Arrays are per series, shape (B,).
    """

    def __init__(self, C, T, M, s, sigma_a):
        self.C = np.asarray(C, dtype=complex)
        self.T = np.asarray(T, dtype=float)
        self.M = int(M)
        self.s = np.asarray(s, dtype=float)
        self.sigma_a = np.asarray(sigma_a, dtype=float)

    @classmethod
    def from_readings(cls, readings):
        """The state of every reading in a ReadingSet (one s and one sigma_a per series)."""
        if not readings.isotropic:
            raise ValueError("the stacked state needs one force SD and one acceleration SD for every "
                             "component of a series; with per-axis noise (ours26) the state would be one "
                             "(C, T, M) per noise pair, which is not built")
        for name, sd in (("force", readings.force_sd), ("acceleration", readings.acceleration_sd)):
            if not np.allclose(sd, sd[:, :1], rtol=1e-12, atol=0.0):
                raise ValueError(f"the stacked state needs one {name} SD for every reading of a series")
        C = np.sum(((readings.Q - readings.P) + 2j * readings.D) / 4.0, axis=1)
        T = np.sum(readings.P + readings.Q, axis=1)
        return cls(C, T, readings.readings * readings.dimension, readings.s[:, 0], readings.acceleration_sd[:, 0])

    def add(self, other):
        """The state of the union of two sets of readings of the same series."""
        if not (np.allclose(self.s, other.s, rtol=1e-12) and np.allclose(self.sigma_a, other.sigma_a, rtol=1e-12)):
            raise ValueError("states can be added only with one noise ratio and one sigma_a")
        return StackedState(self.C + other.C, self.T + other.T, self.M + other.M, self.s, self.sigma_a)

    def _angle(self, u):
        z = np.asarray(u, dtype=float) - np.log(self.s)[:, None]
        cos2 = -np.tanh(z)                                  # cos 2theta, tan theta = exp(z)
        sin2 = 1.0 / np.cosh(np.clip(z, -700, 700))         # sin 2theta
        cos = np.exp(-0.5 * np.logaddexp(0.0, 2.0 * z))
        return cos2, sin2, cos

    def h2(self, u):
        """H(theta)^2 on log-mass grids u of shape (B, G)."""
        cos2, sin2, _ = self._angle(u)
        re = self.C.real[:, None] * cos2 + self.C.imag[:, None] * sin2
        return np.maximum(0.0, self.T[:, None] / 2.0 + 2.0 * re)

    def log_likelihood(self, u):
        """sum_i log K(h_i) = H^2 / 2 (the radius-nuisance likelihood of all readings)."""
        return 0.5 * self.h2(u)

    def cond_alpha(self, u):
        """E[|a*| | m] for the stacked true accelerations: sigma_a cos(theta) mu_M(H)."""
        _, _, cos = self._angle(u)
        return self.sigma_a[:, None] * cos * stacked_mean_radius(np.sqrt(self.h2(u)), self.M)


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


PER_AXIS_ERROR = ("this reading set has per-axis noise (different SDs on different axes of a reading); "
                  "only ours26 reads it directly: reference cartesian1 or cartesian2, nuisance radius, "
                  "rule likelihood_product (or reduce split_axes first)")


def _axis_sds(sd, shape):
    """Supplied SDs as (B, N, d): a (B, N)-broadcastable array is per reading (every axis
    the same); a three-dimensional array is per axis."""
    sd = np.asarray(sd, dtype=float)
    if sd.ndim == 3:
        return np.broadcast_to(sd, shape).copy()
    return np.broadcast_to(np.broadcast_to(sd, shape[:2])[..., None], shape).copy()


class ReadingSet:
    """A batch of series of readings, as the analyst sees them.

    force, acceleration: (B, N, d). force_sd, acceleration_sd: (B, N) per reading, or
    (B, N, d) per axis (ours26). force_sd_axes, acceleration_sd_axes are always (B, N, d).
    A set is `isotropic` when every reading has one SD on all its axes; only then are the
    per-reading force_sd, acceleration_sd, P, Q, D and s defined.
    """

    def __init__(self, force, acceleration, force_sd, acceleration_sd):
        self.force = np.asarray(force, dtype=float)
        self.acceleration = np.asarray(acceleration, dtype=float)
        if self.force.ndim != 3 or self.force.shape != self.acceleration.shape:
            raise ValueError("readings must have shape (series, readings, dimension)")
        self.force_sd_axes = _axis_sds(force_sd, self.force.shape)
        self.acceleration_sd_axes = _axis_sds(acceleration_sd, self.force.shape)
        if not (np.all(self.force_sd_axes > 0) and np.all(self.acceleration_sd_axes > 0)):
            raise ValueError("supplied standard deviations must be positive")
        self.dimension = self.force.shape[2]
        if self.dimension not in DIMENSIONS:
            raise ValueError(f"dimension must be one of {DIMENSIONS}")
        # Per-series quantities an analyst stage computed and later stages reuse
        # (sliced along with the series by subset()).
        self.series_data = {}
        fa, aa = self.force_sd_axes, self.acceleration_sd_axes
        self.isotropic = bool(np.all(fa == fa[..., :1]) and np.all(aa == aa[..., :1]))
        if self.isotropic:
            self._force_sd = fa[..., 0].copy()
            self._acceleration_sd = aa[..., 0].copy()
            self._P, self._Q, self._D, self._s = standardize(self.force, self.acceleration,
                                                             self._force_sd, self._acceleration_sd)

    def _per_reading(self, name):
        if not self.isotropic:
            raise ValueError(PER_AXIS_ERROR)
        return getattr(self, name)

    force_sd = property(lambda self: self._per_reading("_force_sd"))
    acceleration_sd = property(lambda self: self._per_reading("_acceleration_sd"))
    P = property(lambda self: self._per_reading("_P"))
    Q = property(lambda self: self._per_reading("_Q"))
    D = property(lambda self: self._per_reading("_D"))
    s = property(lambda self: self._per_reading("_s"))

    def first_reading_log_ratios(self):
        """log(sigma_F,k / sigma_a,k) on each axis of the series' first reading, (B, d').

        After split_axes this is still the ORIGINAL first reading's d ratios (kept in
        series_data), so a prior centred on the first reading means the same reading
        whether or not the axes were split.
        """
        stored = self.series_data.get("first_reading_log_ratios")
        if stored is not None:
            return stored
        return np.log(self.force_sd_axes[:, 0, :] / self.acceleration_sd_axes[:, 0, :])

    def search_center(self):
        """Where the numerical search for the law starts, per series (B,): log s of the first
        reading; with per-axis noise, the mean of its axes' log ratios (of the original first
        reading after split_axes, so a vector series and its axes are searched alike). Numerics only."""
        ratios = self.first_reading_log_ratios()
        if self.isotropic and not np.any(np.ptp(ratios, axis=1) > 0):
            return np.log(self.s[:, 0])
        return np.mean(ratios, axis=1)

    def subset(self, index):
        """The series selected by `index` (a slice or index array)."""
        out = ReadingSet(self.force[index], self.acceleration[index],
                         self.force_sd_axes[index], self.acceleration_sd_axes[index])
        out.series_data = {k: v[index] for k, v in self.series_data.items()}
        return out

    def take_readings(self, index):
        """The same series, keeping only the readings selected by `index` (a slice)."""
        out = ReadingSet(self.force[:, index], self.acceleration[:, index],
                         self.force_sd_axes[:, index], self.acceleration_sd_axes[:, index])
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

        reference: 'flat' (21 Sept, df d(alpha) dOmega), or 'cartesian1' / 'cartesian2'
        (24 Sept, flat in the vector; two implementations, see the module docstring).
        Returns (log_like, log_ref, cond_alpha), each (B, N, G):
          log_like   log-likelihood of mass with this reading's nuisance integrated,
          log_ref    the mass factor the single-reading reference implies (density in u),
          cond_alpha E[alpha | m] for this reading.
        """
        if nuisance not in NUISANCES:
            raise ValueError(f"nuisance must be one of {NUISANCES}")
        if reference not in REFERENCES:
            raise ValueError(reference_error(reference))
        if not self.isotropic:
            return self._per_axis_terms(u, nuisance, reference)
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

    def _per_axis_terms(self, u, nuisance, reference):
        """ours26's law for per-axis noise: each reading's log-likelihood is the sum of its
        axes' one-dimensional cartesian log-likelihoods (axis splitting). Only the radius
        nuisance has this form, and the per-reading reference factor and E[|alpha_i| | m]
        are not built, so they come back as None."""
        if reference not in ("cartesian1", "cartesian2") or nuisance != "radius":
            raise ValueError(PER_AXIS_ERROR)
        b, n, d = self.force.shape
        like, _, _ = self.split_axes().reading_terms(u, "radius", reference=reference)
        return like.reshape(b, n, d, -1).sum(axis=2), None, None

    def stacked_cond_alpha(self, u, *, reference):
        """E[|a*| | m] for the whole series, with a* the N readings' true accelerations
        stacked into one vector of N d components, on log-mass grids u of shape (B, G).

        Under either cartesian reference the latent pushes are independent N(w_i, I_d)
        at fixed theta, so the stacked one is N(W, I_{Nd}) with |W|^2 = sum_i h_i^2:

            E[|a*| | m] = sigma_a cos(theta) E|N(W, I_{Nd})|.

        This is the acceleration weight of `ratio_of_means_rss` (ours25b), used when the
        series has one sigma_F and one sigma_a (one noise-unit frame for the whole series).
        Otherwise (ours26: SDs that differ between readings or between axes) each component
        c has its own weight sigma_a,c cos(theta_c) and the mean is the weighted one,
        `weighted_stacked_mean` (see `_weighted_cond_alpha`).
        """
        if reference not in ("cartesian1", "cartesian2"):
            raise ValueError("ratio_of_means_rss needs a cartesian reference (the stacked latent is "
                             "Gaussian only when the pushes are weighted flat as vectors)")
        u = np.asarray(u, dtype=float)
        if not (self.isotropic and np.all(self.force_sd == self.force_sd[:, :1])
                and np.all(self.acceleration_sd == self.acceleration_sd[:, :1])):
            return self._weighted_cond_alpha(u)
        z = u - np.log(self.s[:, :1])
        sn = np.exp(-0.5 * np.logaddexp(0.0, -2.0 * z))
        cs = np.exp(-0.5 * np.logaddexp(0.0, 2.0 * z))
        P, Q, D = (np.sum(v, axis=1)[:, None] for v in (self.P, self.Q, self.D))
        big_h = np.sqrt(np.maximum(0.0, P * sn * sn + Q * cs * cs + 2.0 * D * sn * cs))
        k = self.readings * self.dimension
        return self.acceleration_sd[:, :1] * cs * stacked_mean_radius(big_h, k)

    def _noise_groups(self):
        """The series' components grouped by their noise pair (sigma_F, sigma_a).

        Returns onehot (B, K, G) mapping the K = N d components to G groups (G the most
        groups any series has; unused groups have no components), dof (B, G), and the
        groups' sigma_F, sigma_a (B, G) (unused groups get 1)."""
        b, n, d = self.force.shape
        k = n * d
        pairs = np.stack([self.force_sd_axes.reshape(b, k), self.acceleration_sd_axes.reshape(b, k)], axis=-1)
        found = [np.unique(pairs[i], axis=0, return_inverse=True) for i in range(b)]
        g = max(len(values) for values, _ in found)
        onehot = np.zeros((b, k, g))
        sf, sa = np.ones((b, g)), np.ones((b, g))
        for i, (values, inverse) in enumerate(found):
            onehot[i, np.arange(k), inverse.ravel()] = 1.0
            sf[i, :len(values)], sa[i, :len(values)] = values[:, 0], values[:, 1]
        return onehot, onehot.sum(axis=1), sf, sa

    def _weighted_cond_alpha_at(self, u):
        """E|a*_stack| given m, evaluated directly at every point of u (B, S)."""
        b, n, d = self.force.shape
        onehot, dof, sf, sa = self._noise_groups()
        z = u[:, None, :] - np.log(sf / sa)[:, :, None]                          # (B, G, S)
        sn = np.exp(-0.5 * np.logaddexp(0.0, -2.0 * z))
        cs = np.exp(-0.5 * np.logaddexp(0.0, 2.0 * z))
        x = (self.force / self.force_sd_axes).reshape(b, n * d)
        y = (self.acceleration / self.acceleration_sd_axes).reshape(b, n * d)
        sn_c = np.einsum("bkg,bgs->bks", onehot, sn)
        cs_c = np.einsum("bkg,bgs->bks", onehot, cs)
        w2 = np.einsum("bkg,bks->bgs", onehot, (x[:, :, None] * sn_c + y[:, :, None] * cs_c) ** 2)
        c2 = (sa[:, :, None] * cs) ** 2
        mean = weighted_stacked_mean(np.moveaxis(c2, 1, -1), np.moveaxis(w2, 1, -1), dof[:, None, :])
        return mean

    def _weighted_cond_alpha(self, u):
        """ours26's stacked weight on grids u (B, G).

        On a uniform grid (as posterior_summaries uses) it is computed at points no more
        than WEIGHTED_SPACING apart in log m and carried to the grid by a cubic spline
        (the weight is smooth in log m; about 1e-9 relative); otherwise point by point."""
        b, g = u.shape
        span = u[:, -1] - u[:, 0]
        t = (u - u[:, :1]) / np.where(span > 0, span, 1.0)[:, None]
        common = np.linspace(0.0, 1.0, g)
        points = int(min(g, max(WEIGHTED_MIN_POINTS, math.ceil(float(np.max(span)) / WEIGHTED_SPACING) + 1)))
        if points >= g or not np.allclose(t, common[None, :], rtol=0.0, atol=1e-12):
            return self._weighted_cond_alpha_at(u)
        nodes = np.linspace(0.0, 1.0, points)
        values = self._weighted_cond_alpha_at(u[:, :1] + span[:, None] * nodes[None, :])
        return CubicSpline(nodes, values, axis=1)(common)

    def split_axes(self):
        """Each d-dimensional reading as d one-dimensional readings (axis splitting).

        The readings of a series become N d readings in 1D, in reading-major order;
        each keeps its axis's supplied SDs (its reading's, when they are per reading).
        The original first reading's axis ratios are kept in series_data, for a prior
        centred on the first reading.
        """
        b, n, d = self.force.shape
        out = ReadingSet(self.force.reshape(b, n * d, 1), self.acceleration.reshape(b, n * d, 1),
                         self.force_sd_axes.reshape(b, n * d), self.acceleration_sd_axes.reshape(b, n * d))
        out.series_data = dict(self.series_data)
        out.series_data.setdefault("first_reading_log_ratios", self.first_reading_log_ratios())
        return out

    def pooled(self):
        """Inverse-variance mean of each series' readings, as one reading per series.

        Exact sufficient reduction when every reading of a series measures the
        same latent vector pair with independent Gaussian errors (contribution 11,
        eq. 5). The channels are independent, so each is pooled separately; with
        per-axis noise, each axis of each channel.
        """
        if not self.isotropic:
            wf = 1.0 / self.force_sd_axes**2
            wa = 1.0 / self.acceleration_sd_axes**2
            force = np.sum(wf * self.force, axis=1) / np.sum(wf, axis=1)
            acceleration = np.sum(wa * self.acceleration, axis=1) / np.sum(wa, axis=1)
            return ReadingSet(force[:, None, :], acceleration[:, None, :],
                              (1.0 / np.sqrt(np.sum(wf, axis=1)))[:, None, :],
                              (1.0 / np.sqrt(np.sum(wa, axis=1)))[:, None, :])
        wf = 1.0 / self.force_sd**2
        wa = 1.0 / self.acceleration_sd**2
        force = np.sum(wf[:, :, None] * self.force, axis=1) / np.sum(wf, axis=1)[:, None]
        acceleration = np.sum(wa[:, :, None] * self.acceleration, axis=1) / np.sum(wa, axis=1)[:, None]
        return ReadingSet(force[:, None, :], acceleration[:, None, :],
                          (1.0 / np.sqrt(np.sum(wf, axis=1)))[:, None],
                          (1.0 / np.sqrt(np.sum(wa, axis=1)))[:, None])
