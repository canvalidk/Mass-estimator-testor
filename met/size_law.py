"""The shared size law's shrinkage, beside eta (O14, 4 Nov).

ours28b (docs/ours27/ours28_candidate_push_scale_as_nuisance_2026-09-28.md, section 4) declared every
push ~ N(0, omega^2 I) in noise units with one omega for the series. Marginalising the pushes gives a
genuine likelihood in beta = omega^2 / (1 + omega^2):

    L(theta, beta) = (1 - beta)^(M/2) exp(beta H^2(theta) / 2),     M = sum_i k_i,

with H^2(theta) the series' total explained energy on the line theta (equation.md, H_i summed). In
standard terms this is the structural errors-in-variables model (Lindley & El-Sayyad 1968) and beta is
its shrinkage factor. Profiling beta at the best line gives

    beta_hat = 1 - M / H^2(theta_hat) = Delta_hat / (Delta_hat + M),

the James-Stein / empirical-Bayes form: the share of the fitted explained energy that is signal
rather than the M units of noise the latent pushes absorb (one per latent component; M is the
Neyman-Scott count). ours28b integrated beta under the arcsine law instead; `beta_bar` is that
posterior mean. Both are what the declared size law supplies and L10 forbids declaring.

eta (T31, `met.full.bartlett_eta`) estimates the same quantity from the identity's failure instead
of from a declaration: on long series eta -> Delta / (Delta + M) (T18's inverse; O14, 24 Oct). The
tool `tools/eta_vs_size_law.py` compares the two on seeded worlds; `tests/test_eta_target.py` pins
the identities. Nothing here is the equation or a candidate for it: beta_hat fails Z (an exact-zero
reading adds k to M and nothing to H^2), which is the shared-omega row of the 28 Sept table.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from met.full import _h_terms, bartlett_eta


def explained_energy(readings, c, phi):
    """H^2(phi): the series' total explained energy on the line phi, in noise units."""
    return float(sum(_h_terms(r.pre(), c, np.array([phi]))[0] for r in readings))


def latent_count(readings):
    """M = sum_i k_i: the number of latent push components (incidental parameters)."""
    return sum(r.k for r in readings)


def flat_eta(readings, c):
    """eta of the flat law (g = 1, every reading): (eta, curvature, J, phi_hat)."""
    return bartlett_eta(readings, c, [1.0] * len(readings))


def beta_hat(readings, c, phi_hat=None):
    """The profile shrinkage of the shared size law at the series' best line, clipped at 0."""
    if phi_hat is None:
        phi_hat = flat_eta(readings, c)[3]
    h2 = explained_energy(readings, c, phi_hat)
    if h2 <= 0:
        return 0.0
    return max(0.0, 1.0 - latent_count(readings) / h2)


def beta_profile(h2, m):
    """argmax over beta of (m/2) log(1 - beta) + beta h2 / 2, by search (checks the closed form)."""
    if h2 <= m:
        return 0.0
    res = minimize_scalar(lambda b: -(0.5 * m * math.log1p(-b) + 0.5 * b * h2),
                          bounds=(0.0, 1.0 - 1e-12), method="bounded", options={"xatol": 1e-13})
    return float(res.x)


def beta_bar(h2, m):
    """The posterior mean of beta under ours28b's arcsine law, Beta(1/2, 1/2), at explained energy h2."""
    def weight(b):
        return b ** -0.5 * (1.0 - b) ** (-0.5 + 0.5 * m) * math.exp(0.5 * b * h2 - 0.5 * h2)
    num = quad(lambda b: b * weight(b), 0.0, 1.0, limit=200)[0]
    den = quad(weight, 0.0, 1.0, limit=200)[0]
    return num / den


def oracle(delta, m):
    """Delta / (Delta + M) with the true whitened signal energy Delta."""
    return delta / (delta + m)
