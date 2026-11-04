"""T66/T67 research references: local projection energy and pairwise ordering.

The supplied prior and positive prefactor specify comparison constructions,
not production choices. Finite atoms check algebra, not continuous-law limits.
No candidate basin, true signal, or learned signal population is an input.
"""
import math

import numpy as np
from scipy import integrate
from scipy.special import logsumexp

from met.joint_mass import FiniteJointMassLaw, gaussian_joint_geometry
from met.joint_score_law import gaussian_joint_score_statistic
from met.score_law_obstruction import ScalarScoreQuadrature


def scalar_local_projection_energy(readings, mass, *, prior_density):
    """Sum ||dP/dq z||^2, q the supplied prior CDF, for scalar 2D events.

    prior_density is dPi/dm at this candidate. Under another parameter it
    must be transformed with its Jacobian; it is not a numerical density
    held fixed across coordinates. Atomic priors have no such derivative.
    """
    mass, density = float(mass), float(prior_density)
    if not math.isfinite(mass) or mass <= 0 or not math.isfinite(density) or density <= 0:
        raise ValueError("finite positive mass and prior density required")
    contributions = []
    for reading in readings:
        if reading.body_count != 1 or reading.dimension != 1 or reading.protocol != "common_acceleration":
            raise ValueError("full scalar two-channel events required")
        whitening = np.linalg.solve(reading._chol, np.eye(2))
        vector = whitening @ np.array([mass, 1.])
        speed = abs(float(np.linalg.det(whitening))) / float(vector @ vector)
        z = whitening @ reading.data
        contributions.append((speed/density)**2 * float(z @ z))
    value = math.fsum(contributions)
    if not math.isfinite(value):
        raise ValueError("local projection energy exceeds floating-point range")
    return value


def pairwise_profile_penalty(readings, prior):
    """Return A(m)=E_t [(L(t)-L(m))_+/sqrt(1+V(m,t))] on supplied atoms.

    V(m,t)=sum ||(P_i(t)-P_i(m))z_i||^2. Each row is a candidate m,
    each column its rival t. Profile contrasts use projection differences,
    so a common large residual constant never enters the subtraction.
    The same original joint prior weights every rival; no product coupling
    or empirical posterior is substituted. This bounded helper is O(J^2).
    """
    if not isinstance(prior, FiniteJointMassLaw):
        raise ValueError("prior must be a FiniteJointMassLaw")
    readings = tuple(readings)
    count = len(prior.atoms)
    contrast, variance = np.zeros((count, count)), np.zeros((count, count))
    for reading in readings:
        if (reading.protocol != "common_acceleration"
                or reading.body_indices != tuple(range(prior.mass_count))):
            raise ValueError("canonical full-body common-acceleration events required")
        projected = []
        for masses in prior.atoms:
            geometry = gaussian_joint_geometry(reading, masses)
            projected.append(geometry["whitened_design"] @ geometry["fitted_latent"])
        projected = np.asarray(projected)
        for i in range(count):
            difference = projected - projected[i]
            variance[i] += np.einsum("ij,ij->i", difference, difference)
            contrast[i] += .5*np.einsum("ij,ij->i", difference, projected+projected[i])
    if np.any(~np.isfinite(contrast)) or np.any(~np.isfinite(variance)):
        raise ValueError("pairwise projection contrasts exceed floating-point range")
    penalty = (np.maximum(contrast, 0.) / np.sqrt(1.+variance)) @ prior.probabilities
    return dict(penalty=penalty, profile_contrast=contrast, projection_difference=variance)


def joint_pairwise_anchor(readings, prior, *, anchor_prefactor, rank_rtol, rank_atol):
    """Proper finite joint score law with the explicit T67 ordering penalty.

    This operator uses the prior inside its kernel. It is not in T62's
    prior-independent sequential class, and no sequential closure is claimed.
    """
    prefactor = float(anchor_prefactor)
    if not math.isfinite(prefactor) or prefactor <= 0:
        raise ValueError("anchor_prefactor must be finite and positive")
    readings = tuple(readings)
    state = pairwise_profile_penalty(readings, prior)
    statistics = np.array([gaussian_joint_score_statistic(readings, masses,
        rank_rtol=rank_rtol, rank_atol=rank_atol)["statistic"] for masses in prior.atoms])
    # Remove the total contrast minimum before adding unequal prior weights.
    minimum_statistic = float(np.min(statistics))
    with np.errstate(over="ignore", invalid="ignore"):
        cost = (statistics-minimum_statistic)/2 + prefactor*state["penalty"]
    if np.any(~np.isfinite(cost)):
        raise ValueError("combined pairwise contrast exceeds floating-point range")
    minimum = float(np.min(cost))
    weights = np.log(prior.probabilities) - (cost-minimum)
    evidence = float(logsumexp(weights))
    log_evidence = evidence-minimum-minimum_statistic/2
    if not math.isfinite(log_evidence):
        raise ValueError("pairwise log evidence exceeds floating-point range")
    law = FiniteJointMassLaw(prior.atoms, np.exp(weights-evidence),
                             log_evidence=log_evidence)
    state.update(score_statistic=statistics, anchor_prefactor=prefactor,
                 scope="finite supplied joint prior; continuous limits require T67 conditions")
    return law, state


class ScalarPairwiseQuadrature(ScalarScoreQuadrature):
    """Continuous full-half-Cauchy reference, one common scalar covariance.

    All mass hypotheses remain in the outer and rival integrals. Integration
    points resolve both score roots; they do not select a root or prior basin.
    Parent routines supply only Gaussian score geometry and outer quadrature.
    """
    def __init__(self, scatter, covariance, *, prior_center, anchor_prefactor,
                 epsabs, epsrel, limit):
        self.prefactor = float(anchor_prefactor)
        if not math.isfinite(self.prefactor) or self.prefactor <= 0:
            raise ValueError("anchor_prefactor must be finite and positive")
        super().__init__(scatter, covariance, prior_center=prior_center,
                         anchor_exponent=.5, anchor_weight_override=0.,
                         epsabs=epsabs, epsrel=epsrel, limit=limit)

    def penalty(self, theta):
        low, high = self._coordinates(theta)
        def integrand(rival):
            other_low, other_high = self._coordinates(rival)
            gain = .5*self.gap*(other_high-high)*(other_high+high)
            sine = low*other_high-high*other_low
            variance = float(np.sum(self.eigenvalues))*sine*sine
            return max(gain, 0.)/math.sqrt(1.+variance)*2/math.pi
        # Equal-profile rivals are the candidate and its reflected eigenline.
        reflected = self._eigenvectors @ np.array([-low, high])
        physical = np.linalg.solve(self.whitening, reflected)
        points = list(self.points)+[theta]
        if physical[0]*physical[1] > 0:
            points.append(math.atan(physical[0]/physical[1]/self.center))
        points = sorted(set(p for p in points if 0 < p < math.pi/2))
        return integrate.quad(integrand, 0., math.pi/2, points=points, **self.controls)

    def kernel(self, theta, method):
        if method == "pairwise":
            _, statistic = self.components(theta)
            return math.exp(-statistic/2-self.prefactor*self.penalty(theta)[0])
        return super().kernel(theta, method)
