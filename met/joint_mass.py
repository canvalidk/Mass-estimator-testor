"""T54 research references for declared physical joint Gaussian experiments.

A shared acceleration is one observation with one full source covariance.
Separate accelerations are separately declared free latent vectors. Neither
protocol infers a signal population, independence from marginals, or a joint
prior. The finite-atom update below is a proper joint profile-weighted law;
it is not a certification of a continuous quadrature grid or a production
estimator. No multivariate calibration rule is selected.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy import special

from met.composition_calibration import FinitePositiveLaw


def _masses(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not values.size or np.any(~np.isfinite(values))
            or np.any(values <= 0)):
        raise ValueError("masses must be a finite nonempty positive vector")
    return values


def _readonly(values):
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _body_indices(indices, count):
    indices = tuple(range(count)) if indices is None else tuple(indices)
    if (len(indices) != count or len(set(indices)) != count
            or any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer))
                   or i < 0 for i in indices)):
        raise ValueError("body_indices must be distinct nonnegative integer labels")
    return tuple(int(i) for i in indices)


@dataclass(frozen=True)
class JointReading:
    """One immutable source event, constructed by an explicit protocol factory.

    Data order is all force blocks first, then the one common acceleration
    or the separately observed acceleration blocks in body_indices order.
    The covariance covers this unique observation vector in exactly that
    order. Readings are combined as independent events; within an event its
    full covariance is retained. Arrays are copied and made read-only.
    """
    data: np.ndarray
    covariance: np.ndarray
    dimension: int
    body_indices: tuple
    protocol: str
    _chol: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if (isinstance(self.dimension, (bool, np.bool_))
                or not isinstance(self.dimension, (int, np.integer)) or self.dimension < 1
                or self.protocol not in ("common_acceleration", "separate_accelerations")):
            raise ValueError("a positive dimension and a supported explicit protocol are required")
        indices = _body_indices(self.body_indices, len(self.body_indices))
        if not indices:
            raise ValueError("an event must declare at least one body")
        count, d = len(indices), int(self.dimension)
        size = (count+1)*d if self.protocol == "common_acceleration" else 2*count*d
        data, covariance = np.asarray(self.data, dtype=float), np.asarray(self.covariance, dtype=float)
        if (data.shape != (size,) or covariance.shape != (size, size)
                or np.any(~np.isfinite(data)) or np.any(~np.isfinite(covariance))
                or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.)):
            raise ValueError("finite data and matching symmetric full covariance are required")
        if not np.array_equal(covariance, covariance.T):
            covariance = covariance/2+covariance.T/2
        try:
            chol = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError("source covariance must be positive definite") from error
        object.__setattr__(self, "data", _readonly(data))
        object.__setattr__(self, "covariance", _readonly(covariance))
        object.__setattr__(self, "_chol", _readonly(chol))
        object.__setattr__(self, "dimension", d)
        object.__setattr__(self, "body_indices", indices)

    @property
    def body_count(self):
        return len(self.body_indices)

    @property
    def latent_dimension(self):
        return self.dimension if self.protocol == "common_acceleration" else self.body_count*self.dimension

    def _selected(self, masses):
        masses = _masses(masses)
        if max(self.body_indices) >= len(masses):
            raise ValueError("a declared body label lies outside the mass vector")
        return masses, masses[list(self.body_indices)]

    def design(self, masses):
        """Physical Gaussian mean design C(m), shape [observations, latents]."""
        _, selected = self._selected(masses)
        d = self.dimension
        if self.protocol == "common_acceleration":
            return np.vstack([mass*np.eye(d) for mass in selected]+[np.eye(d)])
        upper = np.diag(np.repeat(selected, d))
        return np.vstack([upper, np.eye(self.body_count*d)])

    def design_derivatives(self, masses):
        """dC/dm_a, shape [all masses, observations, latents].

        Unselected body labels have zero derivatives. Separate latent blocks
        follow body_indices, not the order of the full supplied mass vector.
        """
        masses, _ = self._selected(masses)
        d, q = self.dimension, self.latent_dimension
        derivatives = np.zeros((len(masses), len(self.data), q))
        for j, label in enumerate(self.body_indices):
            columns = slice(0, d) if self.protocol == "common_acceleration" else slice(j*d, (j+1)*d)
            derivatives[label, j*d:(j+1)*d, columns] = np.eye(d)
        return derivatives


def _forces(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or min(values.shape) < 1 or np.any(~np.isfinite(values)):
        raise ValueError("forces must be a finite nonempty body-by-dimension array")
    return values


def common_acceleration(forces, acceleration, covariance, *, body_indices=None):
    """Declare (F_1,...,F_K,a) with one shared free latent acceleration."""
    forces = _forces(forces)
    acceleration = np.asarray(acceleration, dtype=float)
    if acceleration.shape != (forces.shape[1],):
        raise ValueError("common acceleration must be one vector matching the force dimension")
    return JointReading(np.r_[forces.ravel(), acceleration], covariance, forces.shape[1],
                        _body_indices(body_indices, len(forces)), "common_acceleration")


def separate_accelerations(forces, accelerations, covariance, *, body_indices=None):
    """Declare (F_1,...,F_K,a_1,...,a_K) with separately free latent vectors.

    This declares the latent relation, not independent measurement errors:
    correlations within the supplied full covariance remain intact.
    """
    forces = _forces(forces)
    accelerations = np.asarray(accelerations, dtype=float)
    if accelerations.shape != forces.shape:
        raise ValueError("separate accelerations must match the full force array shape")
    return JointReading(np.r_[forces.ravel(), accelerations.ravel()], covariance, forces.shape[1],
                        _body_indices(body_indices, len(forces)), "separate_accelerations")


def gaussian_joint_geometry(reading, masses):
    """Reusable whitened QR geometry of a declared event at a fixed candidate.

    No true signals enter. The public physical design methods are also
    available for independent implementations and the T55 research reference.
    """
    if not isinstance(reading, JointReading):
        raise ValueError("reading must be a JointReading")
    design, derivatives = reading.design(masses), reading.design_derivatives(masses)
    basis = np.linalg.solve(reading._chol, design)
    whitened_derivatives = np.array([np.linalg.solve(reading._chol, item) for item in derivatives])
    whitened = np.linalg.solve(reading._chol, reading.data)
    rotation, triangular = np.linalg.qr(basis, mode="complete")
    q = reading.latent_dimension
    parallel, perpendicular = rotation[:, :q], rotation[:, q:]
    coefficients = parallel.T@whitened
    zero_projection = False
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            rhs = design.T@np.linalg.solve(reading.covariance, reading.data)
        if np.all(np.isfinite(rhs)) and np.all(rhs == 0.):
            coefficients = np.zeros_like(coefficients)
            zero_projection = True
    except np.linalg.LinAlgError:
        pass
    return dict(design=design, design_derivatives=derivatives,
                whitened_design=basis, whitened_derivatives=whitened_derivatives,
                whitened_data=whitened, parallel=parallel, perpendicular=perpendicular,
                triangular=triangular[:q], fitted_latent=np.linalg.solve(triangular[:q], coefficients),
                residual=perpendicular.T@whitened, zero_projection=zero_projection)


def gaussian_profile_cost(reading, masses):
    """Stable profiled squared residual, without a physical-residual determinant."""
    geometry = gaussian_joint_geometry(reading, masses)
    # If exact supplied floating normal equations recognize Pz=0, its cost
    # is ||z||². Re-rotating a huge orthogonal datum can otherwise introduce
    # false candidate-dependent residual-norm roundoff. This optional guard
    # is not a general exact-arithmetic cancellation detector.
    residual = geometry["whitened_data"] if geometry["zero_projection"] else geometry["residual"]
    cost = float(residual@residual)
    if not math.isfinite(cost):
        raise ValueError("profile cost exceeds the floating-point range")
    return cost


def efficient_mass_information(readings, masses, latent_signals):
    """Efficient Gaussian information at explicit theoretical latent signals.

    This returns sum D_i' R_i D_i after whitening, with the nuisance latents
    projected out. True vectors are required only by this reference; neither
    gaussian_profile_cost nor gaussian_joint_posterior receives them.
    It is not an observed covariance estimate or a joint calibration rule.
    """
    masses = _masses(masses)
    readings, latent_signals = tuple(readings), tuple(latent_signals)
    if len(readings) != len(latent_signals):
        raise ValueError("one explicit true latent vector is required per declared event")
    result = np.zeros((len(masses), len(masses)))
    for reading, signal in zip(readings, latent_signals):
        geometry = gaussian_joint_geometry(reading, masses)
        signal = np.asarray(signal, dtype=float)
        if signal.shape != (reading.latent_dimension,) or np.any(~np.isfinite(signal)):
            raise ValueError("true latent vectors must be finite and match the event latent dimension")
        tangent = np.einsum("anj,j->na", geometry["whitened_derivatives"], signal)
        transverse = geometry["perpendicular"].T@tangent
        result += transverse.T@transverse
    if np.any(~np.isfinite(result)):
        raise ValueError("efficient information exceeds the floating-point range")
    return result


@dataclass(frozen=True)
class FiniteJointMassLaw:
    """A supplied finite positive joint law, not a certified continuous grid.

    Nonnegative supplied weights are normalized; duplicate joint atoms are
    combined. Body labels are columns. A posterior's log_evidence records
    the bounded Gaussian profile update's normalizing constant.
    """
    atoms: np.ndarray
    probabilities: np.ndarray
    log_evidence: float = None

    def __post_init__(self):
        atoms, weights = np.asarray(self.atoms, dtype=float), np.asarray(self.probabilities, dtype=float)
        if (atoms.ndim != 2 or min(atoms.shape) < 1 or weights.shape != (len(atoms),)
                or np.any(~np.isfinite(atoms)) or np.any(atoms <= 0)
                or np.any(~np.isfinite(weights)) or np.any(weights < 0)
                or not np.any(weights > 0)):
            raise ValueError("positive finite joint atoms and matching nonnegative weights are required")
        if self.log_evidence is not None and not math.isfinite(float(self.log_evidence)):
            raise ValueError("log_evidence must be finite when supplied")
        if self.log_evidence is not None:
            object.__setattr__(self, "log_evidence", float(self.log_evidence))
        weights = weights/np.max(weights)
        total = math.fsum(float(x) for x in weights)
        unique, inverse = np.unique(atoms, axis=0, return_inverse=True)
        combined = np.array([math.fsum(float(x) for x in weights[inverse == j])/total
                             for j in range(len(unique))])
        positive = combined > 0
        object.__setattr__(self, "atoms", _readonly(unique[positive]))
        object.__setattr__(self, "probabilities", _readonly(combined[positive]))

    @property
    def mass_count(self):
        return self.atoms.shape[1]


def gaussian_joint_posterior(readings, prior):
    """Normalize a supplied joint prior times exp(-sum profile costs/2).

    The prior is supplied jointly, never inferred from its marginal laws.
    Independent declared events contribute additive costs. This is a
    profile-weighted raw research law, not integrated nuisance Bayes or a
    Gaussian calibration theorem for finite atoms.
    """
    if not isinstance(prior, FiniteJointMassLaw):
        raise ValueError("prior must be a FiniteJointMassLaw")
    readings = tuple(readings)
    try:
        costs = np.array([math.fsum(gaussian_profile_cost(reading, masses) for reading in readings)
                          for masses in prior.atoms])
    except OverflowError as error:
        raise ValueError("total profile cost exceeds the floating-point range") from error
    if np.any(~np.isfinite(costs)):
        raise ValueError("total profile cost exceeds the floating-point range")
    minimum = float(np.min(costs))
    # Remove the data-only common cost before adding log prior weights. At
    # huge equal costs, adding log(p) first would erase the supplied prior.
    log_weights = np.log(prior.probabilities)-(costs-minimum)/2
    log_shifted_evidence = float(special.logsumexp(log_weights))
    log_evidence = log_shifted_evidence-minimum/2
    weights = np.exp(log_weights-log_shifted_evidence)
    return FiniteJointMassLaw(prior.atoms, weights, log_evidence=log_evidence)


def group_joint_law(law, groups, *, operation):
    """Push the same joint law through explicitly disjoint body groups.

    Groups partition the listed bodies; unlisted bodies are marginalized
    out. If all bodies are grouped, this is a partition of the full label
    set. All mass coordinates must use the same physical units before
    composing. No product of marginals or new independence is inserted.
    """
    if not isinstance(law, FiniteJointMassLaw):
        raise ValueError("law must be a FiniteJointMassLaw")
    groups = tuple(tuple(group) for group in groups)
    flat = [index for group in groups for index in group]
    if (not groups or any(not group for group in groups)
            or len(set(flat)) != len(flat)
            or any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer))
                   or not 0 <= i < law.mass_count for i in flat)):
        raise ValueError("groups must be nonempty, disjoint and contain valid body labels")
    if operation not in ("addition", "parallel"):
        raise ValueError("operation must be explicitly 'addition' or 'parallel'")
    output = np.empty((len(law.atoms), len(groups)))
    for j, group in enumerate(groups):
        values = law.atoms[:, group]
        if operation == "addition":
            try:
                output[:, j] = [math.fsum(float(x) for x in row) for row in values]
            except OverflowError as error:
                raise ValueError("grouped mass exceeds the floating-point range") from error
        else:
            smallest = np.min(values, axis=1)
            output[:, j] = smallest/np.sum(smallest[:, None]/values, axis=1)
    return FiniteJointMassLaw(output, law.probabilities, log_evidence=law.log_evidence)


def scalar_composite_law(law, *, operation):
    """Sum or parallel-sum all bodies under the supplied joint coupling."""
    if not isinstance(law, FiniteJointMassLaw):
        raise ValueError("law must be a FiniteJointMassLaw")
    grouped = group_joint_law(law, [range(law.mass_count)], operation=operation)
    return FinitePositiveLaw(grouped.atoms[:, 0], grouped.probabilities)
