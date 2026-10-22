"""Worlds the tester cannot generate, run on the ours29 reference (law only; no readout).

(a) rescue: 20 readings of signal 3 plus one reading misaligned far beyond its noise (2D);
(b) one strong reading (signal 6) among 100 pure-noise readings (3D);
(c) every reading its own instrument (sigma_F, sigma_a each 10^U(-2, 2)) and its own signal
    size (0, 0, 0.5, 1, 3 or 10 noise units), 3D, 40 readings;
(d) pure noise with the planted mass drawn from the starting law itself, 2D and 3D, 100 readings.
Coverage of the 95% interval of ours27, ours29's law (misfit weight only) and ours29 (quoted).
Run:  python ours29_reference_worlds.py   (seed 30; a few minutes)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ours29_reference import _span_quantities, reading_state, misfit_weight

rng = np.random.default_rng(30)
LM = np.linspace(-12, 12, 2401); M = np.exp(LM)
PHI = np.linspace(-np.pi/2, np.pi/2, 4001, endpoint=False)


def variants(readings, c):
    """Coverage-relevant CDFs for ours27, ours29_law, ours29 from one pass over the readings."""
    phi_law = np.arctan(M/c)
    Hl = np.array([_span_quantities(r[0], r[1], None, c, phi_law)[0] for r in readings])      # (N, G)
    Hf = np.array([_span_quantities(r[0], r[1], None, c, PHI)[0] for r in readings])          # (N, K)
    st = [reading_state(r[0], r[1], None, c) for r in readings]
    g = np.array([misfit_weight(s["Q"], s["nu"]) for s in st])
    logpi = np.log(c/(c*c + M*M)) + LM
    tot = (0.5*g[:, None]*Hf).sum(0); j = int(np.argmax(tot)); h = PHI[1] - PHI[0]
    ell = 0.5*g[:, None]*Hf[:, [j-1, j, (j+1) % len(PHI)]]
    psi = (ell[:, 2] - ell[:, 0])/(2*h); curv = -np.sum(ell[:, 2] - 2*ell[:, 1] + ell[:, 0])/h**2
    J = np.sum(psi**2); eta = 1.0 if J <= 1e-300 or curv <= 0 else min(1.0, curv/J)
    out = {}
    for name, gg, ee in (("ours27", np.ones_like(g), 1.0), ("ours29_law", g, 1.0), ("ours29", g, eta)):
        lp = logpi + ee*(0.5*gg[:, None]*Hl).sum(0); p = np.exp(lp - lp.max()); p /= p.sum()
        out[name] = np.cumsum(p)
    return out, eta


def covered(cdf, m0):
    k = np.searchsorted(M, m0); return 0.025 < cdf[min(k, len(M)-1)] < 0.975


def reading(m0, snr, d, sf=1.0, sa=1.0):
    u = rng.normal(size=d); u /= np.linalg.norm(u)
    a_star = snr*sa*u
    return (np.r_[m0*a_star + sf*rng.normal(size=d), a_star + sa*rng.normal(size=d)], np.diag([sf**2]*d + [sa**2]*d))


def run(name, make, series):
    hits = {k: [] for k in ("ours27", "ours29_law", "ours29")}; etas = []
    for _ in range(series):
        R, m0, c = make()
        cdfs, eta = variants(R, c); etas.append(eta)
        for k in hits: hits[k].append(covered(cdfs[k], m0))
    row = "  ".join(f"{k} {100*np.mean(v):5.1f}" for k, v in hits.items())
    print(f"{name:<52} {row}   median eta {np.median(etas):.2f}   ({series} series)")


def rescue():
    m0 = np.exp(rng.uniform(-1, 1))
    R = [reading(m0, 3.0, 2) for _ in range(20)]
    R.append((np.r_[[2.0, 0.0], [0.0, 1.0]], 0.01*np.eye(4)))               # sigma = 0.1: Q = 100
    return R, m0, 1.0


def strong_among_nulls():
    m0 = np.exp(rng.uniform(-1, 1))
    return [reading(m0, 6.0, 3)] + [reading(m0, 0.0, 3) for _ in range(100)], m0, 1.0


def own_setups():
    m0 = np.exp(rng.uniform(-1, 1))
    R = [reading(m0, rng.choice([0, 0, 0.5, 1, 3, 10]), 3, 10**rng.uniform(-2, 2), 10**rng.uniform(-2, 2)) for _ in range(40)]
    return R, m0, 1.0


def null_from_prior(d):
    def make():
        m0 = np.tan(rng.uniform(0, np.pi/2))                                   # half-Cauchy, c = 1
        return [reading(m0, 0.0, d) for _ in range(100)], m0, 1.0
    return make


print("Coverage of the 95% interval (percent)")
run("(a) rescue: 20 x signal 3 + one impossible reading, 2D", rescue, 150)
run("(b) one signal-6 reading among 100 nulls, 3D", strong_among_nulls, 100)
run("(c) every reading its own instrument and signal, 3D", own_setups, 100)
run("(d) pure noise, mass drawn from pi, 2D, N = 100", null_from_prior(2), 100)
run("(d) pure noise, mass drawn from pi, 3D, N = 100", null_from_prior(3), 100)
