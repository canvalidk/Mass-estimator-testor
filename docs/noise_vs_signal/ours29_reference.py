"""ours29: reference implementation in the full-covariance form, with the readout.

ours29 = the flat per-reading factors (ours27), each tempered by its own off-line misfit (L21, T30),
the product raised to the Bartlett power eta (30 Sept), the prior counted once.

A reading i declares: y_i = (F_i; a_i) in R^{2d}; an error covariance Sigma_i (2d x 2d, positive
definite, any correlation including between the channels); and a candidate subspace N_i (d x k_i,
orthonormal columns; default the identity, k_i = d). Starting law pi: half-Cauchy c/(c^2+m^2).

Per reading, on the whole projective line of candidate slopes (m = c tan(phi), phi in (-pi/2, pi/2],
negative and infinite m included, so every reading has a best line):
   C(phi) = (c sin(phi) N_i; cos(phi) N_i),  G = C^T Sigma^-1 C,  b = C^T Sigma^-1 y,
   H_i(phi) = b^T G^-1 b                                  explained energy (depends only on the span)
   E_i = energy of the whitened data inside the span of all candidates, (N_i 0; 0 N_i)
   Q_i = E_i - max_phi H_i                                 off-line misfit at the best line (inside N_i)
   g_i = min(1, nu_i / Q_i),  nu_i = k_i - 1               "evidence, just less" (kink form; O12 open)
   l_i(phi) = g_i H_i(phi) / 2                             the reading's tempered log factor
Series:
   phi_hat = argmax_phi sum_i l_i                           on the whole line (never at a boundary)
   curv = -sum_i l_i''(phi_hat),  J = sum_i l_i'(phi_hat)^2  (raw: Z forbids counting readings)
   eta = min(1, curv / J)   (eta = 1 when J = 0: one reading, or nothing that disagrees)
Law:        p(m) ~ pi(m) exp(eta sum_i l_i(m))              on m > 0
Candidates: alpha_i | m ~ N(G^-1 b, G^-1), a*_i = N_i alpha_i, F*_i = m a*_i: each reading's own flat
            cloud. eta and g_i change how much a reading counts, never the size of its candidate:
            shrinking candidates by eta would infer signal size from the series (L19).
Readout:    m_hat = int m A(m) p(m) dm / int A(m) p(m) dm,  A(m) = E[ |stacked a*| | m ]   (T9, T10)

calibrate=False gives the raw law (eta = 1); misfit=False gives g_i = 1; both off is ours27.
"""
import numpy as np


def _span_quantities(y, S, Nsub, c, phi):
    """H(phi) on a grid of phi, and the in-span energy E. Also MU, LAM, ROT (for the readout) when asked."""
    y = np.asarray(y, float); S = np.asarray(S, float); d = len(y)//2
    Nsub = np.eye(d) if Nsub is None else np.asarray(Nsub, float)
    Si = np.linalg.inv(S); k = Nsub.shape[1]
    sp, cp = np.sin(phi), np.cos(phi)
    C = np.concatenate([c*sp[:, None, None]*Nsub[None], cp[:, None, None]*Nsub[None]], 1)      # (n, 2d, k)
    G = np.einsum('nai,ab,nbj->nij', C, Si, C); b = np.einsum('nai,ab,b->ni', C, Si, y)
    H = np.einsum('ni,ni->n', b, np.linalg.solve(G, b[..., None])[..., 0])
    P = np.block([[Nsub, np.zeros((d, k))], [np.zeros((d, k)), Nsub]])                        # all candidates' span
    Gp = P.T@Si@P; bp = P.T@Si@y; E = float(bp@np.linalg.solve(Gp, bp))
    return H, E, k


def reading_state(y, S, Nsub, c):
    """The reading's best line on the whole projective line, its misfit Q and nu."""
    phi = np.linspace(-np.pi/2, np.pi/2, 2001, endpoint=False)
    H, E, k = _span_quantities(y, S, Nsub, c, phi)
    j = int(np.argmax(H)); lo, hi = phi[j]-np.pi/2000, phi[j]+np.pi/2000
    for _ in range(60):                                                  # golden-section refine
        a1, a2 = lo+0.382*(hi-lo), lo+0.618*(hi-lo)
        h1, h2 = _span_quantities(y, S, Nsub, c, np.array([a1, a2]))[0]
        lo, hi = (a1, hi) if h2 > h1 else (lo, a2)
    Hmax = float(_span_quantities(y, S, Nsub, c, np.array([(lo+hi)/2]))[0][0])
    return dict(Q=max(E - Hmax, 0.0), E=E, k=k, nu=k-1)


def misfit_weight(Q, nu, form="kink"):
    if nu <= 0 or Q <= 0: return 1.0
    if form == "kink": return min(1.0, nu/Q)
    if form == "smooth": return 1.0/(1.0 + Q/nu)
    raise ValueError(form)


def _series_log(readings, c, g, phi):
    tot = np.zeros_like(phi); per = []
    for r, gi in zip(readings, g):
        y, S = r[0], r[1]; Nsub = r[2] if len(r) > 2 else None
        l = 0.5*gi*_span_quantities(y, S, Nsub, c, phi)[0]; per.append(l); tot += l
    return tot, per


def bartlett_eta(readings, c, g, h=1e-4):
    """eta = min(1, curv/J) at the series' best line on the whole projective line. Returns (eta, curv, J, phi_hat)."""
    if len(readings) == 0: return 1.0, 0.0, 0.0, 0.0
    phi = np.linspace(-np.pi/2, np.pi/2, 4001, endpoint=False)
    tot, _ = _series_log(readings, c, g, phi)
    j = int(np.argmax(tot)); lo, hi = phi[j]-np.pi/4000, phi[j]+np.pi/4000
    for _ in range(60):
        a1, a2 = lo+0.382*(hi-lo), lo+0.618*(hi-lo)
        t = _series_log(readings, c, g, np.array([a1, a2]))[0]
        lo, hi = (a1, hi) if t[1] > t[0] else (lo, a2)
    p0 = (lo+hi)/2
    _, per = _series_log(readings, c, g, np.array([p0-h, p0, p0+h]))
    per = np.array(per)                                                   # (N, 3)
    psi = (per[:, 2]-per[:, 0])/(2*h); curv_i = -(per[:, 2]-2*per[:, 1]+per[:, 0])/h**2
    curv, J = float(curv_i.sum()), float((psi**2).sum())
    if J <= 1e-300: return 1.0, curv, J, p0
    return float(min(1.0, max(curv, 0.0)/J)) if curv > 0 else 1.0, curv, J, p0


def lse(a, axis):
    mx = np.max(a, axis=axis, keepdims=True); return np.squeeze(mx + np.log(np.sum(np.exp(a-mx), axis=axis, keepdims=True)), axis)


def _readout_terms(y, S, Nsub, m):
    """MU = G^-1 b and the eigen-decomposition of G^-1 at each m (the reading's flat candidate cloud)."""
    y = np.asarray(y, float); S = np.asarray(S, float); d = len(y)//2
    Nsub = np.eye(d) if Nsub is None else np.asarray(Nsub, float)
    Si = np.linalg.inv(S)
    C = np.concatenate([m[:, None, None]*Nsub[None], np.broadcast_to(Nsub, (len(m),) + Nsub.shape)], 1)
    G = np.einsum('mai,ab,mbj->mij', C, Si, C); b = np.einsum('mai,ab,b->mi', C, Si, y)
    MU = np.linalg.solve(G, b[..., None])[..., 0]
    LAM, ROT = np.linalg.eigh(np.linalg.inv(G))
    return MU, LAM, ROT


def ours29(readings, c=1.0, misfit=True, calibrate=True, g_form="kink", lm_half=9.0, nodes=721,
           readout=True, tnodes=241):
    """readings: list of (y, Sigma) or (y, Sigma, N). Returns readout, median, 95% interval, law, grid, eta, g."""
    states = [reading_state(r[0], r[1], r[2] if len(r) > 2 else None, c) for r in readings]
    g = [misfit_weight(s["Q"], s["nu"], g_form) if misfit else 1.0 for s in states]
    eta, curv, J, phi_hat = bartlett_eta(readings, c, g) if calibrate else (1.0, np.nan, np.nan, np.nan)
    lm = np.linspace(np.log(c)-lm_half, np.log(c)+lm_half, nodes); m = np.exp(lm)
    logpi = np.log(c/(c*c+m*m)) + lm
    phi_m = np.arctan(m/c)
    tot, _ = _series_log(readings, c, g, phi_m) if readings else (np.zeros_like(m), [])
    lp = logpi + eta*tot
    p = np.exp(lp - lp.max()); p /= p.sum()
    cdf = np.cumsum(p) - 0.5*p; q = np.exp(np.interp([0.025, 0.5, 0.975], cdf, lm))
    out = dict(median=q[1], lo=q[0], hi=q[2], law=p, lm=lm, eta=eta, g=g, Q=[s["Q"] for s in states],
               curv=curv, J=J)
    if not readout or not readings:
        return out
    per = [_readout_terms(r[0], r[1], r[2] if len(r) > 2 else None, m) for r in readings]
    A = np.zeros_like(m); lt = np.linspace(-20, 20, tnodes); keep = p > 1e-14*p.max()
    for j in np.nonzero(keep)[0]:
        EX = sum(((ROT[j].T@MU[j])**2).sum() + LAM[j].sum() for MU, LAM, ROT in per)
        tt = np.exp(lt)/max(EX, 1e-300); logL = np.zeros_like(tt)
        for MU, LAM, ROT in per:
            mt2 = (ROT[j].T@MU[j])**2
            logL += (-0.5*np.log1p(2*tt[:, None]*LAM[j][None, :])
                     - tt[:, None]*mt2[None, :]/(1+2*tt[:, None]*LAM[j][None, :])).sum(1)
        A[j] = np.trapezoid(-np.expm1(logL)*tt**(-0.5), lt)/(2*np.sqrt(np.pi))
    out["readout"] = (m*A*p).sum()/(A*p).sum(); out["A"] = A
    return out
