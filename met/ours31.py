"""ours31: a probe built on T37 (the sigmas believed, each reading's noise size its shared state).

NOT the definition. ours30 is the working equation (R-2026-10-29-1); this is the fourteenth probe,
built to the premise list as it stands on 30 Oct with the candidates L24 and the G18 sort, and
built to break (GOAL.md). Its holes are stated in problems/O18.md before any result.

Per reading i (full type: any Sigma_i, declared subspace N_i, k_i = dim N_i, nu_i = k_i - 1):

  H_i(m)   explained energy, as in met.full (the projection onto the mass subspace, T36)
  Q+_i     the misfit at the best positive mass (R-2026-09-07-1, T32), E_i - max_{m>0} H_i
  A_i(m)   = max_{m>0} H_i - H_i(m) >= 0 on m > 0, the along-line misfit at m (T37 (4))

  the noise size v_i ~ InvGamma(nu_eff_i/2, nu_eff_i/2): the stated sigma estimated from nu_eff_i
  Gaussian degrees of freedom, precision right on average (T37, G20; GUM Annex G). Under L21's cut
  (the size is read from the misfit at the best mass) the reading's factor is the exact average of
  its Gaussian factor over the size's posterior:

  f_i(m)   = (1 + A_i(m) / (nu_eff_i + Q+_i)) ^ (-(nu_eff_i + nu_i) / 2)

  series:  no eta. Each reading's size is its own; nothing is shared unless declared (L24, B12).
  law:     p(m) ~ pi(m) prod_i f_i(m) on m > 0, pi the half-Cauchy with centre c (the input).
  point:   the law's median (as ours30).

nu_eff is an input with no default: the GUM's fallback for an unstated nu is infinity, which is the
fixed-size Gaussian (f = exp(-A/2), ours27's raw law), and that fails T30 (O18, hole 1).

Known holes, stated before results (problems/O18.md): the size posterior uses the strong-signal
count Q | v ~ v chi2_nu, also for Q+ (the sign condition is not counted); weak signal is not modelled
in the size posterior; the power form of L21 is replaced by the exact cut factor (G19, a ruling).
"""
import math

import numpy as np

from met.full import _h_terms, _log_prior, _quantiles, reading_state

COARSE_HALF = 25.0
COARSE_NODES = 4001
FINE_NODES = 2001
SUPPORT_NATS = 40.0


def _per_reading(readings, c, nu_eff):
    states = []
    for r, ne in zip(readings, nu_eff):
        st = reading_state(r, c)
        h_pos_max = r.pre()["energy"] - st["Q_plus"]
        states.append((r, ne, st, h_pos_max))
    return states


def _log_factors(states, c, lm):
    phi = np.arctan(np.exp(lm - math.log(c)))
    total = np.zeros_like(lm)
    for r, ne, st, h_pos_max in states:
        along = np.maximum(h_pos_max - _h_terms(r.pre(), c, phi), 0.0)
        nu = r.k - 1
        if math.isinf(ne):
            total -= along / 2
        else:
            total -= 0.5 * (ne + nu) * np.log1p(along / (ne + st["Q_plus"]))
    return total


def ours31(readings, c, nu_eff, *, grid=None):
    """The ours31 probe on one series. readings: a list of met.full.Reading. c: pi's centre (input).
    nu_eff: the effective degrees of freedom of each reading's stated uncertainty, a positive number
    or a list (one per reading); math.inf gives the fixed-size Gaussian (no T30).

    Returns a dict: point (the law's median), median, interval (central 95%), law, lm, Q_plus, nu_eff."""
    if not c > 0:
        raise ValueError("the starting law's centre c must be positive")
    readings = list(readings)
    if np.ndim(nu_eff) == 0:
        nu_eff = [float(nu_eff)] * len(readings)
    nu_eff = [float(x) for x in nu_eff]
    if len(nu_eff) != len(readings) or not all(x > 0 for x in nu_eff):
        raise ValueError("nu_eff: one positive value per reading (no default: O18 hole 1)")
    states = _per_reading(readings, c, nu_eff)

    def log_law(lm):
        return _log_prior(lm, c) + (_log_factors(states, c, lm) if states else 0.0)

    if grid is not None:
        half, nodes = grid
        lm = np.linspace(math.log(c) - half, math.log(c) + half, int(nodes))
    else:
        coarse = np.linspace(math.log(c) - COARSE_HALF, math.log(c) + COARSE_HALF, COARSE_NODES)
        lp = log_law(coarse)
        keep = np.nonzero(lp > lp.max() - SUPPORT_NATS)[0]
        step = coarse[1] - coarse[0]
        lm = np.linspace(coarse[keep[0]] - 2 * step, coarse[keep[-1]] + 2 * step, FINE_NODES)
    lp = log_law(lm)
    (lo, med, hi), p = _quantiles(lm, lp, [0.025, 0.5, 0.975])
    return dict(point=float(med), median=float(med), interval=(float(lo), float(hi)), law=p, lm=lm,
                Q_plus=[s[2]["Q_plus"] for s in states], nu_eff=nu_eff)
