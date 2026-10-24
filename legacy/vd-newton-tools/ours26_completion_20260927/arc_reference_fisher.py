"""Check arc prior as Jeffreys prior of an isotropic proper reference model.

This is a proposed mathematical completion, not a production change. Snapshot:
VD-docs@c877b1df131bc3a809ef25bc270160eb5b11b6d5.
"""
from pathlib import Path
import json
import numpy as np


def projector(m, whitening):
    d = whitening.shape[0] // 2
    b = whitening @ np.vstack((m*np.eye(d), np.eye(d)))
    e = whitening @ np.vstack((np.eye(d), np.zeros((d, d))))
    k = np.linalg.inv(b.T @ b)
    p = b @ k @ b.T
    t = (np.eye(2*d)-p) @ e @ k @ b.T
    dp = t + t.T
    g = float(.5*np.trace(dp @ dp))
    return p, dp, g


def main():
    rng = np.random.default_rng(94713)
    a = rng.normal(size=(6, 6))
    sigma = a @ a.T + np.eye(6)
    eigen, rotation = np.linalg.eigh(sigma)
    whitening = (rotation / np.sqrt(eigen)) @ rotation.T
    z = rng.normal(size=6)
    errors = {"projector_derivative": [], "fisher_relative": [], "likelihood_relative_log_constant": [], "conditional_mean": [], "conditional_covariance": []}
    for tau in (.3, 1., 5., 100.):
        beta = tau*tau/(1+tau*tau)
        constants = []
        for m in (.1, .5, 1., 3., 20.):
            p, dp, g = projector(m, whitening)
            h = 1e-5*max(1, m)
            fd = (projector(m+h, whitening)[0]-projector(m-h, whitening)[0])/(2*h)
            errors["projector_derivative"].append(float(np.max(np.abs(dp-fd))))
            c = np.eye(6) + tau*tau*p
            ci = np.linalg.inv(c)
            cp = tau*tau*dp
            fisher = .5*np.trace(ci @ cp @ ci @ cp)
            expected = tau**4/(1+tau*tau)*g
            errors["fisher_relative"].append(float(abs(fisher-expected)/expected))
            loglik = -.5*np.linalg.slogdet(c)[1] - .5*z @ ci @ z
            q = z @ (np.eye(6)-p) @ z
            constants.append(float(loglik + .5*beta*q))
            b = whitening @ np.vstack((m*np.eye(3), np.eye(3)))
            k = np.linalg.inv(b.T @ b)
            prior_cov = tau*tau*k
            posterior_cov = np.linalg.inv(np.linalg.inv(prior_cov)+b.T @ b)
            posterior_mean = posterior_cov @ b.T @ z
            errors["conditional_mean"].append(float(np.max(np.abs(posterior_mean-beta*k @ b.T @ z))))
            errors["conditional_covariance"].append(float(np.max(np.abs(posterior_cov-beta*k))))
        errors["likelihood_relative_log_constant"].append(float(np.ptp(constants)))
    maxima = {key: max(value) for key, value in errors.items()}
    assert maxima["projector_derivative"] < 2e-9
    assert maxima["fisher_relative"] < 2e-10
    assert maxima["likelihood_relative_log_constant"] < 2e-10
    assert maxima["conditional_mean"] < 2e-12
    assert maxima["conditional_covariance"] < 2e-12
    out = {"source_commit": "c877b1df131bc3a809ef25bc270160eb5b11b6d5", "dimension": 3, "noise_covariance": "random SPD including cross-channel and cross-component correlation", "checks_max_error": maxima}
    Path(__file__).with_name("arc_reference_fisher_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
