"""Continuous scalar T64 research comparison for heterogeneous scatters.

The supplied prior is a half-Cauchy conditioned on explicit finite bounds.
Those bounds are prior support, not a numerical truncation of a different
law. Scatters are sums of observed physical (force, acceleration) outer
products; no sample count, true signal or population is an evaluation input.
Quadrature errors and stationary-point searches are diagnostics, not proofs.
"""
import math

import numpy as np
from scipy import integrate, optimize


def prior_projector_variance(whitening, center, theta_bounds, nodes):
    """Scalar projection variance under the declared truncated half-Cauchy.

    tr(P)=1 implies E(P)-E(P)^2=v*I_2. Compute v through centered
    projectors, avoiding cancellation in a nearly deterministic projection.
    Legendre nodes approximate a uniform prior angle; nodes is explicit.
    """
    roots, weights = np.polynomial.legendre.leggauss(nodes)
    angles = theta_bounds[0]+(roots+1)*(theta_bounds[1]-theta_bounds[0])/2
    vectors = np.column_stack([center*np.sin(angles), np.cos(angles)])@whitening.T
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    projectors = np.einsum("ni,nj->nij", vectors, vectors)
    offsets = projectors-projectors[0]
    mean_offset = np.einsum("n,nij->ij", weights/2, offsets)
    centered = offsets-mean_offset
    variance = float(np.einsum("n,nij,nij->", weights/2, centered, centered)/2)
    return projectors[0]+mean_offset, variance


class ScalarMassContrastQuadrature:
    """Full supplied-prior support, retaining every detected score root.

    Each nonzero grouped scatter must be positive definite; exact zero
    groups are allowed. This bounded reference does not choose numerical
    rank gates for arbitrary singular scatter or certify all root searches.
    """
    def __init__(self, scatters, covariances, *, prior_center, prior_bounds,
                 anchor_prefactor, prior_nodes, epsabs, epsrel, limit):
        center, prefactor = float(prior_center), float(anchor_prefactor)
        bounds = np.asarray(prior_bounds, dtype=float)
        if (not math.isfinite(center) or center <= 0 or not math.isfinite(prefactor)
                or prefactor <= 0 or bounds.shape != (2,)
                or np.any(~np.isfinite(bounds)) or not 0 < bounds[0] < bounds[1]):
            raise ValueError("positive prior center/prefactor and finite positive prior bounds required")
        if type(prior_nodes) is not int or prior_nodes < 2 or type(limit) is not int or limit < 2:
            raise ValueError("explicit integer quadrature controls >=2 required")
        if not all(math.isfinite(v) and v > 0 for v in (epsabs, epsrel)):
            raise ValueError("positive finite integration controls required")
        scatters, covariances = list(scatters), list(covariances)
        if len(scatters) != len(covariances):
            raise ValueError("one declared covariance is required per scatter")
        self.center, self.bounds, self.prefactor = center, bounds, prefactor
        self.theta_bounds = tuple(np.arctan(bounds/center))
        self.controls = dict(epsabs=epsabs, epsrel=epsrel, limit=limit)
        self.groups, energies, contrast_energies = [], [], []
        reference = math.atan(math.sqrt(bounds[0]*bounds[1])/center)
        candidates = list(self.theta_bounds)
        for scatter, covariance in zip(scatters, covariances):
            scatter, covariance = np.asarray(scatter, float), np.asarray(covariance, float)
            if (scatter.shape != (2, 2) or covariance.shape != (2, 2)
                    or np.any(~np.isfinite(scatter)) or np.any(~np.isfinite(covariance))
                    or not np.allclose(scatter, scatter.T, rtol=1e-12, atol=0.)
                    or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
                raise ValueError("finite symmetric two-channel scatter and covariance required")
            chol = np.linalg.cholesky(covariance)
            whitening = np.linalg.solve(chol, np.eye(2))
            white = whitening@scatter@whitening.T
            white = white/2+white.T/2
            values, vectors = np.linalg.eigh(white)
            if values[0] <= 0 and np.any(scatter):
                raise ValueError("nonzero grouped scatter must be positive definite")
            mean_projector, variance = prior_projector_variance(
                whitening, center, self.theta_bounds, prior_nodes)
            energy = float(np.trace(white))
            contrast = variance*energy
            if not math.isfinite(energy) or not math.isfinite(contrast):
                raise ValueError("observed energy exceeds the floating-point range")
            energies.append(energy)
            contrast_energies.append(contrast)
            group = dict(whitening=whitening, determinant=float(np.linalg.det(whitening)),
                low=float(values[0]), gap=float(values[1]-values[0]), high_vector=vectors[:, 1],
                mean_projector=mean_projector, projector_variance=variance)
            group["reference_line"] = self._line(group, reference)[0]
            self.groups.append(group)
            for vector in vectors.T:
                physical = chol@vector
                if physical[0]*physical[1] > 0:
                    angle = math.atan((physical[0]/physical[1])/center)
                    if self.theta_bounds[0] < angle < self.theta_bounds[1]:
                        candidates.append(angle)
        self.energy = math.fsum(energies)
        self.contrast_energy = math.fsum(contrast_energies)
        self.weights = dict(trace=prefactor/math.sqrt(1+self.energy),
                            contrast=prefactor/math.sqrt(1+self.contrast_energy), score=0.)
        self.state = dict(observed_whitened_trace=self.energy,
            observed_mass_contrast_energy=self.contrast_energy, anchor_prefactor=prefactor,
            anchor_weights=dict(self.weights), prior_nodes=prior_nodes,
            prior_center=center, prior_bounds=bounds.tolist())
        mesh = sorted(set(np.linspace(*self.theta_bounds, 129).tolist()+candidates))
        scores = [self.fields(t)[0] for t in mesh]
        roots = []
        for left, right, sleft, sright in zip(mesh, mesh[1:], scores, scores[1:]):
            if sleft*sright < 0:
                roots.append(optimize.brentq(lambda t: self.fields(t)[0], left, right))
            elif sleft == 0:
                roots.append(left)
        self.roots = sorted(set(roots))
        points = candidates+self.roots
        for root in self.roots:
            distance = min(root-self.theta_bounds[0], self.theta_bounds[1]-root)
            step = min(1e-5, distance/4)
            if step > 0:
                slope = abs((self.fields(root+step)[0]-self.fields(root-step)[0])/(2*step))
                covariance = self.fields(root)[1]
                if slope > 0 and covariance > 0:
                    width = math.sqrt(covariance)/slope
                    points.extend(root+k*width for k in (-12, -4, -1, 1, 4, 12))
        self.points = sorted(set(t for t in points if self.theta_bounds[0] < t < self.theta_bounds[1]))
        self.maximum_profile = max(self.fields(t)[2] for t in candidates+self.roots)
        self._totals = {}

    def _line(self, group, theta):
        vector = group["whitening"]@np.array([self.center*math.sin(theta), math.cos(theta)])
        norm = math.hypot(*vector)
        unit = vector/norm
        speed = group["determinant"]*self.center/(norm*norm)
        return unit, speed

    def fields(self, theta):
        """Profile score/C in common prior-angle coordinates, shifted profile."""
        scores, covariances, profiles = [], [], []
        for group in self.groups:
            unit, speed = self._line(group, theta)
            tangent = np.array([unit[1], -unit[0]])
            high = group["high_vector"]@unit
            cross = group["high_vector"]@tangent
            projected = group["low"]+group["gap"]*high*high
            scores.append(speed*group["gap"]*high*cross)
            covariances.append(speed*speed*projected)
            difference = group["high_vector"]@(unit-group["reference_line"])
            sum_projection = group["high_vector"]@(unit+group["reference_line"])
            profiles.append(.5*group["gap"]*difference*sum_projection)
        return math.fsum(scores), math.fsum(covariances), math.fsum(profiles)

    def kernel(self, theta, method):
        if method not in self.weights:
            raise ValueError("method must be trace, contrast or score")
        score, covariance, profile = self.fields(theta)
        statistic = score*score/covariance if covariance > 0 else 0.
        exponent = -.5*statistic+self.weights[method]*(profile-self.maximum_profile)
        # A data-only shift need not equal the exact maximum. Do not clip
        # a positive exponent: that would change the normalized law.
        return math.exp(exponent)

    def _integral(self, theta, method):
        if theta <= self.theta_bounds[0]:
            return 0., 0.
        return integrate.quad(lambda t: self.kernel(t, method), self.theta_bounds[0], theta,
            points=[p for p in self.points if p < theta], **self.controls)

    def total(self, method):
        if method not in self._totals:
            self._totals[method] = self._integral(self.theta_bounds[1], method)
        return self._totals[method]

    def cdf(self, mass, method):
        if mass <= self.bounds[0]:
            return 0.
        if mass >= self.bounds[1]:
            return 1.
        return self._integral(math.atan(mass/self.center), method)[0]/self.total(method)[0]

    def quantile(self, probability, method):
        if not 0 < probability < 1:
            raise ValueError("probability must lie in (0,1)")
        total = self.total(method)[0]
        root = optimize.brentq(lambda t: self._integral(t, method)[0]/total-probability,
                               *self.theta_bounds, xtol=1e-11)
        return self.center*math.tan(root)
