"""Finite Bayesian update for two masses under a known delivered total force.

For each independent push the observed channels are (F1, F2, A). At a mass
pair (x, y) their Gaussian mean is s * (x, y, 1) / (x + y), where s is the
calibrated signed load. See e12-mass-observables.md for the experiment.
"""
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp


@dataclass(frozen=True)
class MassPairPosterior:
    """Weights aligned with the input pairs, including zero-prior entries."""

    mass_pairs: np.ndarray
    probabilities: np.ndarray
    log_probabilities: np.ndarray

    def _centered_moments(self):
        # Center at the modal atom before summing; retain tiny contributions to
        # moments even when the probability itself underflows to zero.
        anchor = self.mass_pairs[np.argmax(self.log_probabilities)]
        offsets = self.mass_pairs - anchor
        root_weights = np.exp(0.5 * self.log_probabilities)
        mean_offset = root_weights @ (root_weights[:, None] * offsets)
        return anchor + mean_offset, root_weights[:, None] * (offsets - mean_offset)

    @property
    def mean(self):
        """Posterior arithmetic means of the two labelled masses."""
        return self._centered_moments()[0]

    @property
    def covariance(self):
        """Population covariance matrix of the two masses under the posterior."""
        _, centered = self._centered_moments()
        return centered.T @ centered


def controlled_force_posterior(mass_pairs, prior_weights, loads, readings, covariances, *,
                               log_prior=False):
    """Condition a finite positive mass-pair prior on independent scalar pushes.

    Supply one (m1, m2) per prior weight, one signed load per (F1, F2, A)
    reading, and either a common 3-by-3 covariance or one covariance per push.
    A single scalar load and reading triple are accepted. Empty data preserve
    the prior. Weights can be unnormalized. For sequential updates, pass the
    preceding log_probabilities with log_prior=True to retain underflowed atoms.

    Example: controlled_force_posterior([(1, 3), (2, 2)], [1, 1],
                                        4, [1, 3, 1], np.eye(3)).
    """
    pairs = np.array(mass_pairs, dtype=float, copy=True)
    if (pairs.ndim != 2 or pairs.shape[1] != 2 or not len(pairs)
            or not np.isfinite(pairs).all() or np.any(pairs <= 0)):
        raise ValueError("expected a nonempty array of finite positive mass pairs")
    weights = np.asarray(prior_weights, dtype=float)
    if weights.shape != (len(pairs),):
        raise ValueError("expected one prior weight per mass pair")
    if log_prior:
        if np.any(np.isnan(weights) | np.isposinf(weights)) or np.all(np.isneginf(weights)):
            raise ValueError("log weights must be finite or -inf with positive total weight")
        logs = weights.copy()
    else:
        if np.any(~np.isfinite(weights) | (weights < 0)) or not np.any(weights > 0):
            raise ValueError("expected finite nonnegative weights with positive total weight")
        with np.errstate(divide="ignore"):
            logs = np.log(weights)
    logs -= logs.max()
    active = np.isfinite(logs)

    loads = np.atleast_1d(np.asarray(loads, dtype=float))
    data = np.asarray(readings, dtype=float)
    if loads.size == 0 and data.shape == (0,):
        data = data.reshape(0, 3)
    if loads.size == 1 and data.shape == (3,):
        data = data[None, :]
    if (loads.ndim != 1 or data.shape != (len(loads), 3)
            or not np.isfinite(loads).all() or not np.isfinite(data).all()):
        raise ValueError("expected a finite scalar load and a finite reading triple per push")
    cov = np.asarray(covariances, dtype=float)
    if not len(loads) and cov.shape == (0,):
        cov = cov.reshape(0, 3, 3)
    if cov.shape == (3, 3):
        cov = cov[None, :, :]
    if (cov.shape not in ((1, 3, 3), (len(loads), 3, 3)) or not np.isfinite(cov).all()
            or not np.allclose(cov, cov.swapaxes(-1, -2), rtol=1e-12, atol=0)):
        raise ValueError("expected a symmetric finite 3-by-3 covariance common or per push")
    try:
        roots = np.broadcast_to(np.linalg.cholesky(cov), (len(loads), 3, 3))
    except np.linalg.LinAlgError as exc:
        raise ValueError("each covariance must be positive definite") from exc

    # Scale masses by a power of two so even an overflowing x+y can be used.
    masses = pairs[active]
    _, exponent = np.frexp(masses.max(axis=1))
    scaled = np.ldexp(masses, -exponent[:, None])
    totals = scaled.sum(axis=1)
    shares = scaled / totals[:, None]
    anchor = np.argmax(logs[active])
    likelihood = _relative_log_likelihood(anchor, shares, totals, exponent, loads, data, roots)
    # Recenter on a likely atom. Two near-tied candidates can each have log
    # odds of 1e20 against the first reference and differ by less than one.
    peak = np.argmax(likelihood)
    if peak != anchor:
        likelihood = _relative_log_likelihood(peak, shares, totals, exponent, loads, data, roots)
    scores = logs[active] + likelihood
    scores -= scores.max()
    logs[active] = scores - logsumexp(scores)
    probabilities = np.exp(logs)
    for array in (pairs, logs, probabilities):
        array.setflags(write=False)
    return MassPairPosterior(pairs, probabilities, logs)


def _relative_log_likelihood(anchor, shares, totals, exponent, loads, data, roots):
    # Either body specifies the change in allocation; use the smaller shares
    # to distinguish nearly coincident allocations near an endpoint.
    share_change = np.where(shares[:, 0] + shares[anchor, 0] < 1,
                            shares[:, 0] - shares[anchor, 0],
                            shares[anchor, 1] - shares[:, 1])
    directions = np.array([[1., 0.], [-1., 0.], [0., 1.]])
    scores = np.zeros(len(shares))
    for load, row, root in zip(loads, data, roots):
        if load == 0:
            continue
        mantissa, power = np.frexp(load)
        acceleration = np.ldexp(mantissa / totals, power - exponent)
        reference = np.r_[load * shares[anchor], acceleration[anchor]]
        delta = np.column_stack((load * share_change, acceleration - reference[2]))
        # Evaluate log odds to one supported atom. This cancels Gaussian
        # normalizers and candidate-independent energy before either is squared.
        basis = np.linalg.solve(root, directions)
        score = basis.T @ np.linalg.solve(root, row - reference)
        shift = delta @ basis.T
        scores += delta @ score - 0.5 * np.einsum("ij,ij->i", shift, shift)
    if not np.isfinite(scores).all():
        raise ValueError("Gaussian log odds exceed the floating-point range")
    return scores
