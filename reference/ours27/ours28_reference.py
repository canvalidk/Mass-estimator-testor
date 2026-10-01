"""ours28: reference implementation in the full-covariance form, with the readout.

A reading i declares: y_i = (F_i; a_i) in R^{2d}; an error covariance Sigma_i (2d x 2d, positive definite, any
correlation including between the channels); and a candidate subspace N_i (d x k, orthonormal columns;
default the identity, k = d: nothing known about the push's direction).
Starting law pi: half-Cauchy c/(c^2+m^2). Per-reading scale law G on beta = w^2/(1+w^2): arcsine by default.

Per reading, at trial mass m, with C = B_m N_i, B_m = (m I; I):
   Gm = C^T Sigma^-1 C (k x k),  b = C^T Sigma^-1 y,  H_i(m) = b^T Gm^-1 b   (explained energy)
Law:        p(m) ~ pi(m) prod_i int G(dbeta) (1-beta)^{k_i/2} exp(beta H_i(m)/2)     [arcsine: 1F1(1/2; 1+k/2; H/2)]
Candidates: beta_i | m ~ G(dbeta)(1-beta)^{k/2} e^{beta H/2} (normalised);  alpha_i | m, beta ~ N(beta Gm^-1 b, beta Gm^-1);
            a*_i = N_i alpha_i,  F*_i = m a*_i.   Nothing is shared between readings.
Readout:    m_hat = int m A(m) p(m) dm / int A(m) p(m) dm,   A(m) = E[ |stacked a*| | m ]   (T9, T10 unchanged)
flat=True puts G at beta = 1: ours27."""
import numpy as np

LB = np.linspace(-30, 30, 481); BETA = 1/(1+np.exp(-LB)); L1B = -np.log1p(np.exp(LB))
def lse(a, axis):
    mx = np.max(a, axis=axis, keepdims=True); return np.squeeze(mx + np.log(np.sum(np.exp(a-mx), axis=axis, keepdims=True)), axis)
def log_arcsine():            # arcsine density times d beta / d logit
    return 0.5*np.log(BETA) + 0.5*L1B - np.log(np.pi)

def reading_quantities(y, S, Nsub, m):
    """For one reading on a grid of m: H (nm,), mu (nm,k) = Gm^-1 b, eigen-decomposition of Gm^-1."""
    y = np.asarray(y, float); S = np.asarray(S, float); d = len(y)//2
    Nsub = np.eye(d) if Nsub is None else np.asarray(Nsub, float)
    Si = np.linalg.inv(S); k = Nsub.shape[1]
    C = np.concatenate([m[:, None, None]*Nsub[None], np.broadcast_to(Nsub, (len(m),) + Nsub.shape)], 1)   # (nm, 2d, k)
    Gm = np.einsum('mai,ab,mbj->mij', C, Si, C); b = np.einsum('mai,ab,b->mi', C, Si, y)
    MU = np.linalg.solve(Gm, b[..., None])[..., 0]; H = np.einsum('mi,mi->m', b, MU)
    LAM, ROT = np.linalg.eigh(np.linalg.inv(Gm))
    return H, MU, LAM, ROT, k

def ours28(readings, c=1.0, flat=False, log_scale_law=None, lm_half=9.0, nodes=721, readout=True, bins=8, tnodes=241):
    """readings: list of (y, Sigma) or (y, Sigma, N). Returns readout, median, 95% interval, law, grid."""
    lm = np.linspace(np.log(c)-lm_half, np.log(c)+lm_half, nodes); m = np.exp(lm)
    logpi = np.log(c/(c*c+m*m)) + lm
    lg = log_arcsine() if log_scale_law is None else log_scale_law(BETA, LB)
    lp = logpi.copy(); per = []
    for r in readings:
        y, S = r[0], r[1]; Nsub = r[2] if len(r) > 2 else None
        H, MU, LAM, ROT, k = reading_quantities(y, S, Nsub, m)
        if flat:
            lp += 0.5*H; w = None
        else:
            lj = lg[None, :] + 0.5*k*L1B[None, :] + 0.5*BETA[None, :]*H[:, None]
            lev = lse(lj, 1); lp += lev; w = np.exp(lj - lev[:, None])        # beta | m weights, per m
        per.append((H, MU, LAM, ROT, k, w))
    p = np.exp(lp - lp.max()); p /= p.sum()
    cdf = np.cumsum(p) - 0.5*p; q = np.exp(np.interp([0.025, 0.5, 0.975], cdf, lm))
    out = dict(median=q[1], lo=q[0], hi=q[2], law=p, lm=lm)
    if not readout:
        return out
    # A(m) = E|a*| over the stacked candidates, by the Laplace transform of X = |a*|^2 = sum_i |alpha_i|^2
    A = np.zeros_like(m); lt = np.linspace(-20, 20, tnodes); keep = p > 1e-14*p.max()
    nb = len(BETA)//bins; BB = BETA[:nb*bins].reshape(nb, bins)
    for j in np.nonzero(keep)[0]:
        EX = 0.0
        for H, MU, LAM, ROT, k, w in per:
            mt2 = (ROT[j].T@MU[j])**2; EX += (mt2.sum() + LAM[j].sum())      # E X at beta = 1 (upper scale)
        tt = np.exp(lt)/max(EX, 1e-300); logL = np.zeros_like(tt)
        for H, MU, LAM, ROT, k, w in per:
            mt2 = (ROT[j].T@MU[j])**2                                          # squared means in the eigenbasis of Gm^-1
            if w is None:
                logL += (-0.5*np.log1p(2*tt[:, None]*LAM[j][None, :]) - tt[:, None]*mt2[None, :]/(1+2*tt[:, None]*LAM[j][None, :])).sum(1)
            else:
                wb = w[j][:nb*bins].reshape(nb, bins); ws = wb.sum(1)             # coarse-grain the beta weights
                sel = ws > 1e-14; bb = (wb[sel]*BB[sel]).sum(1)/ws[sel]; lw = np.log(ws[sel])
                den = 1 + 2*tt[:, None, None]*bb[None, :, None]*LAM[j][None, None, :]
                terms = (-0.5*np.log(den) - tt[:, None, None]*(bb[None, :, None]**2)*mt2[None, None, :]/den).sum(2)
                logL += lse(terms + lw[None, :], 1)
        A[j] = np.trapezoid(-np.expm1(logL)*tt**(-0.5), lt)/(2*np.sqrt(np.pi))
    out["readout"] = (m*A*p).sum()/(A*p).sum(); out["A"] = A
    return out
