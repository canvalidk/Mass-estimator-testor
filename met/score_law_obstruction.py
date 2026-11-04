"""T56/T57 continuous scalar diagnostics, with a supplied half-Cauchy prior.

Full positive-mass quadrature uses m=c*tan(theta), theta in (0,pi/2).
Observed scatter eigenlines locate integration breakpoints, never select a
branch or remove a peak. Error estimates are numerical, not certificates.
Anchor schedules, their constants and any weight override are explicit
research choices, not production rules. The historical profile-oscillation
schedule remains the constructor default so recorded comparisons reproduce.
"""
import math

import numpy as np
from scipy import integrate, optimize


class ScalarScoreQuadrature:
    def __init__(self, scatter, covariance, *, prior_center, anchor_exponent,
                 epsabs, epsrel, limit, anchor_schedule="profile_oscillation",
                 anchor_weight_override=None):
        """Full-positive scalar law reference under a supplied half-Cauchy.

        The profile schedule uses the explicit anchor_exponent. The spectral
        comparison uses (1+largest whitened scatter eigenvalue)^(-1/2), with
        prefactor one as an unforced named research choice. A construction-time
        override is separately recorded and leaves cache keys unambiguous.
        No count, true signal, prior patch or branch selection enters.
        """
        scatter, covariance = np.asarray(scatter, dtype=float), np.asarray(covariance, dtype=float)
        if (scatter.shape != (2, 2) or covariance.shape != (2, 2)
                or np.any(~np.isfinite(scatter)) or np.any(~np.isfinite(covariance))
                or not np.allclose(scatter, scatter.T, rtol=1e-12, atol=0.)
                or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
            raise ValueError("finite symmetric two-channel scatter and covariance are required")
        self.center, self.alpha = float(prior_center), float(anchor_exponent)
        if not math.isfinite(self.center) or self.center <= 0 or not 0 < self.alpha < 1:
            raise ValueError("positive finite prior center and anchor exponent in (0,1) required")
        epsabs, epsrel = float(epsabs), float(epsrel)
        if (not math.isfinite(epsabs) or not math.isfinite(epsrel)
                or epsabs <= 0 or epsrel <= 0 or type(limit) is not int or limit < 1):
            raise ValueError("positive quadrature controls required")
        if anchor_schedule not in ("profile_oscillation", "spectral"):
            raise ValueError("anchor_schedule must be profile_oscillation or spectral")
        if anchor_weight_override is not None:
            anchor_weight_override = float(anchor_weight_override)
            if not math.isfinite(anchor_weight_override) or anchor_weight_override < 0:
                raise ValueError("an explicit anchor weight must be finite and nonnegative")
        chol = np.linalg.cholesky(covariance)
        self.whitening = np.linalg.solve(chol, np.eye(2))
        self.scatter = self.whitening@scatter@self.whitening.T
        self.scatter = .5*self.scatter+.5*self.scatter.T
        eigenvalues, eigenvectors = np.linalg.eigh(self.scatter)
        if eigenvalues[0] <= 0:
            raise ValueError("this continuous reference requires positive definite scatter")
        gap = eigenvalues[1]-eigenvalues[0]
        eigenvalues.setflags(write=False)
        eigenvectors.setflags(write=False)
        self.eigenvalues = self.scatter_eigenvalues = eigenvalues
        self._eigenvectors = eigenvectors
        self.smallest_eigenvalue, self.largest_eigenvalue = map(float, eigenvalues)
        self.gap = float(gap)
        candidates = [0., math.pi/2]
        low_squares = [self._coordinates(theta)[0]**2 for theta in candidates]
        roots = []
        for j in range(2):
            physical = chol@eigenvectors[:, j]
            if physical[0]*physical[1] > 0:
                theta = math.atan((physical[0]/physical[1])/self.center)
                candidates.append(theta)
                roots.append((theta, eigenvalues[j]))
                # At a recognized eigenline this coordinate is known exactly;
                # avoid a tiny rotation error in the profile baseline.
                low_squares.append(1. if j == 0 else 0.)
        self.minimum_low_squared = min(low_squares)
        self.maximum_low_squared = max(low_squares)
        self.maximum = self.largest_eigenvalue-gap*self.minimum_low_squared
        self.oscillation = max(0., gap*(self.maximum_low_squared-self.minimum_low_squared)/2)
        self.anchor_schedule = anchor_schedule
        self.profile_anchor_weight = (1+self.oscillation)**(-self.alpha)
        self.spectral_anchor_weight = math.exp(-.5*math.log1p(self.largest_eigenvalue))
        self.derived_anchor_weight = (self.profile_anchor_weight
            if anchor_schedule == "profile_oscillation" else self.spectral_anchor_weight)
        self.anchor_weight_override = anchor_weight_override
        self.anchor_weight = (self.derived_anchor_weight if anchor_weight_override is None
                              else anchor_weight_override)
        self.anchor_state = dict(schedule=self.anchor_schedule, profile_exponent=self.alpha,
            eigenvalues=[self.smallest_eigenvalue, self.largest_eigenvalue], gap=self.gap,
            profile_oscillation=self.oscillation, profile_anchor_weight=self.profile_anchor_weight,
            spectral_anchor_weight=self.spectral_anchor_weight, spectral_prefactor=1.,
            derived_anchor_weight=self.derived_anchor_weight,
            anchor_weight_override=self.anchor_weight_override, effective_anchor_weight=self.anchor_weight)
        points = []
        for theta, value in roots:
            vector = self.whitening@np.array([self.center*math.sin(theta), math.cos(theta)])
            speed = abs(np.linalg.det(self.whitening))*self.center/(vector@vector)
            width = min(math.sqrt(value)/gap, 1/math.sqrt(gap))/speed if gap else math.pi/2
            points.extend(theta+k*width for k in (-40, -12, -4, -1, 0, 1, 4, 12, 40))
        self.points = sorted(set(p for p in points if 0 < p < math.pi/2))
        self.controls = dict(epsabs=float(epsabs), epsrel=float(epsrel), limit=limit)
        self.roots = [self.center*math.tan(t) for t, _ in roots]
        self._totals = {}

    def _coordinates(self, theta):
        vector = self.whitening@np.array([self.center*math.sin(theta), math.cos(theta)])
        coefficients = self._eigenvectors.T@vector
        length = math.hypot(float(coefficients[0]), float(coefficients[1]))
        return float(coefficients[0]/length), float(coefficients[1]/length)

    def components(self, theta):
        low, high = self._coordinates(theta)
        projected = self.smallest_eigenvalue+self.gap*high*high
        standardized_cross = self.gap*low*high/math.sqrt(projected)
        return float(projected), float(standardized_cross*standardized_cross)

    def profile_log_factor(self, theta):
        """Maximum-shifted profile without subtracting O(N) projections."""
        low, _ = self._coordinates(theta)
        return -.5*self.gap*max(0., low*low-self.minimum_low_squared)

    def kernel(self, theta, method):
        projection, statistic = self.components(theta)
        profile = self.profile_log_factor(theta)
        if method == "raw":
            exponent = profile
        elif method == "score":
            exponent = -.5*statistic
        elif method == "anchored_score":
            exponent = -.5*statistic+self.anchor_weight*profile
        else:
            raise ValueError("method must be raw, score or anchored_score")
        return math.exp(min(0., exponent))

    def _integral(self, stop, method):
        if stop <= 0:
            return 0., 0.
        return integrate.quad(lambda t: self.kernel(t, method), 0., stop,
            points=[p for p in self.points if p < stop], **self.controls)

    def total(self, method):
        if method not in self._totals:
            self._totals[method] = self._integral(math.pi/2, method)
        return self._totals[method]

    def cdf(self, mass, method):
        if mass <= 0:
            return 0.
        total, _ = self.total(method)
        return self._integral(math.atan(mass/self.center), method)[0]/total

    def quantile(self, probability, method):
        if not 0 < probability < 1:
            raise ValueError("probability must lie in (0,1)")
        total, _ = self.total(method)
        root = optimize.brentq(lambda t: self._integral(t, method)[0]/total-probability,
                               0., math.pi/2, xtol=1e-11)
        return self.center*math.tan(root)


def witness_constants():
    """Proof quantities only; no truth is passed to the quadrature kernel."""
    true_mass, false_mass = 6/5, 4/15
    covariance = np.array([[1., .9], [.9, 1.]])
    mean = np.array([true_mass, 1.])
    population_scatter = covariance+np.outer(mean, mean)
    ratio = (1159/241)/math.sqrt(47/19)
    return dict(true_mass=true_mass, false_mass=false_mass, covariance=covariance,
                mean=mean, population_scatter=population_scatter,
                false_limit_weight=ratio/(1+ratio))
