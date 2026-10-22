"""ours28b: reference implementation (candidate, 2026-09-28).

Inputs per reading i: measured force F_i, acceleration a_i in R^d, and a declared error covariance
Sigma_i = blockdiag(S_FF, S_aa) whose two blocks commute (per-axis instruments seen in any frame; one isotropic
pair is the special case). Starting law pi (half-Cauchy at c). Scale law G on beta (arcsine = Beta(1/2,1/2)).

Construction:
 1. Rotate each reading into the common eigenbasis of its S_FF and S_aa (covariant; own frames).
    Each component k then has a declared noise pair (sF, sa).
 2. Classes g = components with the same declared pair, pooled over the whole series.
 3. Noise units x = F/sF, y = a/sa, s_g = sF/sa, theta_g(m) = arctan(m/s_g), A_c = x_c sin(theta) + y_c cos(theta),
    H_g(m) = sum_{c in g} A_c^2, M_g = |g|.
 4. Law: p(m) ~ pi(m) prod_g int G(dbeta) (1-beta)^{M_g/2} exp(beta H_g/2).
 5. Candidates given (m, beta_g): v_c ~ N(beta A_c, beta), a*_c = sa_c cos(theta_g) v_c, F*_c = m a*_c.
 6. Readout: m_hat = int m A(m) p(m) / int A(m) p(m), A(m) = E[ |stacked a*| | m ].
Flat limit (G = point mass at beta = 1): ours27 on this input type."""
import numpy as np

def logsumexp(a, axis):
    mx = np.max(a, axis=axis, keepdims=True)
    return np.squeeze(mx + np.log(np.sum(np.exp(a - mx), axis=axis, keepdims=True)), axis)

def to_components(readings, tol=1e-9):
    comps = []   # (x-value F_k, a_k, sF, sa)
    for F, a, S in readings:
        F = np.atleast_1d(np.asarray(F, float)); a = np.atleast_1d(np.asarray(a, float)); d = len(F)
        S = np.asarray(S, float); SFF, SFa, Saa = S[:d, :d], S[:d, d:], S[d:, d:]
        assert np.abs(SFa).max() < tol, "cross-channel correlation: not in ours28b's input type (open)"
        assert np.abs(SFF@Saa - Saa@SFF).max() < tol*max(1, np.abs(SFF).max()*np.abs(Saa).max()), \
            "S_FF and S_aa do not commute: not in ours28b's input type (open)"
        _, U = np.linalg.eigh(SFF + np.pi*Saa)          # generic combination separates the common eigenbasis
        for k in range(d):
            u = U[:, k]
            comps.append((u@F, u@a, np.sqrt(u@SFF@u), np.sqrt(u@Saa@u)))
    return comps

def classes(comps, digits=10):
    out = {}
    for F, a, sF, sa in comps:
        key = (round(sF, digits), round(sa, digits))
        g = out.setdefault(key, {"x": [], "y": [], "sF": sF, "sa": sa})
        g["x"].append(F/sF); g["y"].append(a/sa)
    for g in out.values():
        x, y = np.array(g["x"]), np.array(g["y"])
        g.update(P=(x*x).sum(), Q=(y*y).sum(), D=(x*y).sum(), M=len(x), s=g["sF"]/g["sa"])
    return list(out.values())        # carried state per class: (P, Q, D, M) and the declared pair

# scale-law grid in logit(beta): d beta = beta(1-beta) dl
L = np.linspace(-40, 40, 1601); BETA = 1/(1+np.exp(-L)); LOG1MB = -np.log1p(np.exp(L))
def log_prior_arcsine():                      # arcsine density beta^-1/2 (1-beta)^-1/2 / pi, times beta(1-beta)
    return 0.5*np.log(BETA) + 0.5*LOG1MB - np.log(np.pi)

def ours28b(readings, c=1.0, scale_law="arcsine", lm_half=9.0, nodes=801, flat=False):
    cls = classes(to_components(readings))
    lm = np.linspace(np.log(c)-lm_half, np.log(c)+lm_half, nodes); m = np.exp(lm)
    logpi = np.log(c/(c*c+m*m)) + lm                          # half-Cauchy, density in log m
    lp = logpi.copy()
    lg = log_prior_arcsine() if scale_law == "arcsine" else scale_law(BETA, L)
    post_w = []                                               # beta posterior per class, per m node
    for g in cls:
        th = np.arctan(m/g["s"]); S_, C_ = np.sin(th), np.cos(th)
        H = g["P"]*S_**2 + g["Q"]*C_**2 + 2*g["D"]*S_*C_
        g["th"], g["H"] = th, H
        if flat:
            lp += 0.5*H; post_w.append(None); continue
        lj = lg[None, :] + 0.5*g["M"]*LOG1MB[None, :] + 0.5*BETA[None, :]*H[:, None]
        lev = logsumexp(lj, 1); lp += lev
        post_w.append(np.exp(lj - lev[:, None]))              # normalised weights over the beta grid
    p = np.exp(lp - lp.max()); p /= p.sum()
    # A(m) = E|a*| via the Laplace transform of X = |a*|^2 (independent classes given m)
    A = np.zeros_like(m)
    X0 = sum(g["sa"]**2*np.cos(g["th"])**2*(g["H"] + g["M"]) for g in cls) + 1e-300
    lt = np.linspace(-18, 18, 721)
    for j in range(len(m)):
        tt = np.exp(lt)/X0[j]; logL = np.zeros_like(tt)
        for g, w in zip(cls, post_w):
            k = g["sa"]**2*np.cos(g["th"][j])**2                 # scale of a*^2 per unit v^2
            if w is None:   # flat: v ~ N(A, 1)
                logL += -0.5*g["M"]*np.log1p(2*tt*k) - tt*k*g["H"][j]/(1+2*tt*k)
            else:           # mixture over beta: v ~ N(beta A, beta)
                kb = k*BETA[None, :]; tb = tt[:, None]
                terms = -0.5*g["M"]*np.log1p(2*tb*kb) - tb*k*BETA[None, :]**2*g["H"][j]/(1+2*tb*kb)
                wj = w[j][None, :]; keep = wj[0] > 1e-300
                logL += logsumexp(terms[:, keep] + np.log(wj[:, keep]), 1)
        integrand = -np.expm1(logL)*tt**(-0.5)                   # (1 - L(t)) t^-1/2 per d log t
        A[j] = np.trapezoid(integrand, lt)/(2*np.sqrt(np.pi))
    mhat = (m*A*p).sum()/(A*p).sum()
    cdf = np.cumsum(p) - 0.5*p; q = np.exp(np.interp([0.025, 0.5, 0.975], cdf, lm))
    return dict(readout=mhat, median=q[1], lo=q[0], hi=q[2], law=p, lm=lm, A=A, classes=cls)
