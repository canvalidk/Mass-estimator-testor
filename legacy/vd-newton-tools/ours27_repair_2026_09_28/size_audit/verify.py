"""Independent analytic / numerical size and width audit; no production writes."""
import json
import numpy as np

rng = np.random.default_rng(271028)
d, m = 3, 1.7
C = rng.normal(size=(2*d, 2*d))
Sigma = C @ C.T / 3 + np.eye(2*d) * 0.3
a = np.array([1.0, -0.4, 0.8])
S = np.linalg.inv(Sigma)
ew, ev = np.linalg.eigh(S)
W = (ev * np.sqrt(ew)) @ ev.T

def projector(x):
    B = np.vstack([x*np.eye(d), np.eye(d)])
    H = W @ B
    K = np.linalg.inv(H.T @ H)
    return H @ K @ H.T, K, H

P, K, H = projector(m)
D = W @ np.vstack([np.eye(d), np.zeros((d,d))])
Pp = (np.eye(2*d)-P) @ D @ K @ H.T
Pp = Pp + Pp.T
V = Sigma[:d,:d] - m*(Sigma[:d,d:]+Sigma[d:,:d]) + m*m*Sigma[d:,d:]
J = float(a @ np.linalg.solve(V,a))
G = float(np.trace(np.linalg.solve(V,K)))
mu = H @ a
step = 1e-4
Ppp = (projector(m+step)[0] - 2*P + projector(m-step)[0]) / step**2
z = mu + rng.normal(size=(500000,2*d))
score = .5*np.einsum('ni,ij,nj->n',z,Pp,z)
negative_curvature = -.5*np.einsum('ni,ij,nj->n',z,Ppp,z)
result = {'full_covariance': {
    'J': J, 'G': G,
    'projector_J': float(mu @ Pp @ Pp @ mu),
    'projector_G': float(.5*np.trace(Pp @ Pp)),
    'score_mean_simulated': float(score.mean()),
    'score_variance_expected': J+G,
    'score_variance_simulated': float(score.var()),
    'negative_curvature_expected': J,
    'negative_curvature_simulated': float(negative_curvature.mean()),
    'flat_width_inflation': 1+G/J,
}}

# A common learned scale hears unweighted size; mass curvature hears geometric weights.
m=1.0
sf=np.array([1.,1.]); sa=np.array([.2,1.25])
V=sf**2+m*m*sa**2
K=sf**2*sa**2/V
g=K/V
result['anisotropic_scale_examples']=[]
for v in [np.array([3.,.3]), np.array([.3,3.])]:
    J=float(g@(v*v)); G=float(g.sum())
    beta_iso=float(v@v/(v@v+len(v)))
    result['anisotropic_scale_examples'].append({
        'v':v.tolist(), 'geometric_weights':g.tolist(),
        'correct_scalar_beta':J/(J+G), 'isotropic_learned_beta':beta_iso,
        'sampling_variance_over_reported_variance':beta_iso*(1+G/J),
    })

# Ordered reference scale: q(beta) proportional 1/(1-beta).
# Relative mass evidence odds for h^2=4 versus h^2=0, then append one null.
def evidence(M,t):
    x,w=np.polynomial.legendre.leggauss(160)
    u=(1+x)/2
    return float(w@(u**(M-1)*np.exp((1-u*u)*t/2)))
result['reference_scale_null_counterexample']={
    'before_odds':evidence(1,4)/evidence(1,0),
    'after_one_null_odds':evidence(2,4)/evidence(2,0),
    'after_nine_nulls_odds':evidence(10,4)/evidence(10,0),
}
print(json.dumps(result,indent=2))
