"""T56/T57 continuous scalar diagnostics, with a supplied half-Cauchy prior.

Full positive-mass quadrature uses m=c*tan(theta), theta in (0,pi/2).
Observed scatter eigenlines locate integration breakpoints, never select a
branch or remove a peak. Error estimates are numerical, not certificates.
The anchor exponent is an explicit research choice, not a production rule.
"""
import math

import numpy as np
from scipy import integrate, optimize


class ScalarScoreQuadrature:
    def __init__(self, scatter, covariance, *, prior_center, anchor_exponent,
                 epsabs, epsrel, limit):
        scatter, covariance = np.asarray(scatter, dtype=float), np.asarray(covariance, dtype=float)
        if (scatter.shape != (2, 2) or covariance.shape != (2, 2)
                or np.any(~np.isfinite(scatter)) or np.any(~np.isfinite(covariance))
                or not np.allclose(scatter, scatter.T, rtol=1e-12, atol=0.)
                or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
            raise ValueError("finite symmetric two-channel scatter and covariance are required")
        self.center, self.alpha = float(prior_center), float(anchor_exponent)
        if not math.isfinite(self.center) or self.center <= 0 or not 0 < self.alpha < 1:
            raise ValueError("positive finite prior center and anchor exponent in (0,1) required")
        if epsabs <= 0 or epsrel <= 0 or type(limit) is not int or limit < 1:
            raise ValueError("positive quadrature controls required")
        chol = np.linalg.cholesky(covariance)
        self.whitening = np.linalg.solve(chol, np.eye(2))
        self.scatter = self.whitening@scatter@self.whitening.T
        self.scatter = .5*self.scatter+.5*self.scatter.T
        eigenvalues, eigenvectors = np.linalg.eigh(self.scatter)
        if eigenvalues[0] <= 0:
            raise ValueError("this continuous reference requires positive definite scatter")
        gap = eigenvalues[1]-eigenvalues[0]
        candidates = [0., math.pi/2]
        roots = []
        for j in range(2):
            physical = chol@eigenvectors[:, j]
            if physical[0]*physical[1] > 0:
                theta = math.atan((physical[0]/physical[1])/self.center)
                candidates.append(theta)
                roots.append((theta, eigenvalues[j]))
        extrema = [self.components(t)[0] for t in candidates]
        self.maximum = max(extrema)
        self.oscillation = max(0., (max(extrema)-min(extrema))/2)
        self.anchor_weight = (1+self.oscillation)**(-self.alpha)
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

    def components(self, theta):
        vector = self.whitening@np.array([self.center*math.sin(theta), math.cos(theta)])
        norm = vector@vector
        transverse = np.array([-vector[1], vector[0]])
        projected = (vector@self.scatter@vector)/norm
        cross = (transverse@self.scatter@vector)/norm
        return float(projected), float(cross*cross/projected)

    def kernel(self, theta, method):
        projection, statistic = self.components(theta)
        profile = .5*(projection-self.maximum)
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
