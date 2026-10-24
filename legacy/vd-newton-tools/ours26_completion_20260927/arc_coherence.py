"""Scratch checks of arc prior under Gaussian experiment equivalences."""
import json
import numpy as np


def terms(y, sigma, b, derivative):
    w = np.linalg.inv(sigma)
    k = np.linalg.inv(b.T @ w @ b)
    mean = k @ b.T @ w @ y
    r = w - w @ b @ k @ b.T @ w
    return float(y @ r @ y), float(np.trace(k @ derivative.T @ r @ derivative)), mean, k


def block_diag(blocks):
    n = sum(len(x) for x in blocks)
    result = np.zeros((n, n))
    i = 0
    for x in blocks:
        result[i:i+len(x), i:i+len(x)] = x
        i += len(x)
    return result


def main():
    rng = np.random.default_rng(2726)
    errors = {}
    covs = []
    for _ in range(4):
        c = rng.normal(size=(2, 2))
        covs.append(c @ c.T + np.eye(2))
    observations = rng.normal(size=(4, 2))
    precisions = [np.linalg.inv(c) for c in covs]
    effective_cov = np.linalg.inv(sum(precisions))
    effective_y = effective_cov @ sum(w @ y for w, y in zip(precisions, observations))
    constant_q = sum((y-effective_y) @ w @ (y-effective_y) for w, y in zip(precisions, observations))
    stable_errors = []
    transform_errors = []
    ancillary_errors = []
    for mass in np.geomspace(.03, 30, 21):
        b = np.array([[mass], [1.]])
        d = np.array([[1.], [0.]])
        original = terms(observations.ravel(), block_diag(covs), np.tile(b, (4,1)), np.tile(d, (4,1)))
        compressed = terms(effective_y, effective_cov, b, d)
        stable_errors.extend([original[0] - compressed[0] - constant_q, original[1] - compressed[1]])
        stable_errors.extend((original[2] - compressed[2]).ravel())
        stable_errors.extend((original[3] - compressed[3]).ravel())
        # Invertible processing of actual observations preserves Gaussian experiment.
        t = np.array([[2., .7], [-.3, 1.2]])
        transformed = terms(t @ effective_y, t @ effective_cov @ t.T, t @ b, t @ d)
        for a, c in zip(compressed, transformed):
            transform_errors.extend(np.atleast_1d(a-c).ravel())
        # Append an observed independent pure-noise channel, with zero mean design.
        augmented = terms(np.r_[effective_y, .8], block_diag([effective_cov, np.array([[2.]])]), np.vstack([b, [[0.]]]), np.vstack([d, [[0.]]]))
        ancillary_errors.extend([augmented[0] - compressed[0] - .8**2/2, augmented[1]-compressed[1]])
        ancillary_errors.extend((augmented[2]-compressed[2]).ravel())
        ancillary_errors.extend((augmented[3]-compressed[3]).ravel())
    errors['stable_repeat_compression_max_error'] = float(np.max(np.abs(stable_errors)))
    errors['invertible_experiment_transform_max_error'] = float(np.max(np.abs(transform_errors)))
    errors['ancillary_pure_noise_max_error'] = float(np.max(np.abs(ancillary_errors)))

    # Expanding reference update telescopes, including unequal scales and orders.
    masses = np.geomspace(1e-4, 1e4, 1001)
    speeds = [s/(s*s+masses*masses) for s in [0.3, 2., 9.]]
    likelihoods = [np.exp(-.5*(f-masses*a)**2/(s*s+masses*masses)) for s,f,a in [(.3,.8,.1),(2.,1.4,.7),(9.,-2.,.4)]]
    batch = np.sqrt(sum(g*g for g in speeds)) * np.prod(likelihoods, axis=0)
    sequential = []
    for order in [(0,1,2), (2,0,1), (1,2,0)]:
        first = order[0]
        g = speeds[first]**2
        density = speeds[first]*likelihoods[first]
        for index in order[1:]:
            new_g = g + speeds[index]**2
            density *= likelihoods[index] * np.sqrt(new_g/g)
            g = new_g
        sequential.append(np.max(np.abs(density-batch)))
    errors['expanding_reference_associativity_max_error'] = float(max(sequential))

    # Independent Fisher-information verification of the Gaussian regularization.
    matrix = rng.normal(size=(4, 4))
    sigma = matrix @ matrix.T + np.eye(4)
    e, u = np.linalg.eigh(sigma)
    whitening = (u * (1/np.sqrt(e))) @ u.T
    fisher_errors = []
    for mass in (.1, 1., 10.):
        b = np.vstack([mass*np.eye(2), np.eye(2)])
        d = np.vstack([np.eye(2), np.zeros((2,2))])
        g = terms(np.zeros(4), sigma, b, d)[1]
        bw, dw = whitening @ b, whitening @ d
        k = np.linalg.inv(bw.T @ bw)
        p = bw @ k @ bw.T
        term = (np.eye(4)-p) @ dw @ k @ bw.T
        p_derivative = term + term.T
        for tau in (.2, 1., 5.):
            cov = np.eye(4) + tau*tau*p
            cov_derivative = tau*tau*p_derivative
            fisher = .5*np.trace(np.linalg.solve(cov, cov_derivative) @ np.linalg.solve(cov, cov_derivative))
            fisher_errors.append(fisher - tau**4/(1+tau*tau)*g)
    errors['gaussian_regularization_fisher_max_error'] = float(np.max(np.abs(fisher_errors)))
    assert max(errors.values()) < 1e-10, errors
    print(json.dumps(errors, indent=2))


if __name__ == '__main__':
    main()
