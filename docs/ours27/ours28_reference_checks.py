"""Checks on the ours28 reference implementation (full-covariance form, declared subspace, readout)."""
import numpy as np
from scipy.special import hyp1f1
from ours28_reference import ours28, reading_quantities, BETA, L1B, log_arcsine, lse
rng = np.random.default_rng(7)
def show(name, r): print(f"  {name:<52} readout {r['readout']:.6f}  median {r['median']:.6f}  95% [{r['lo']:.4f}, {r['hi']:.4f}]")
def rot3(): return np.linalg.qr(rng.normal(size=(3, 3)))[0]
def swap(rd):
    y, S = rd[0], rd[1]; d = len(y)//2; P = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])
    return (P@y, P@S@P.T) + tuple(rd[2:])
def rotate(rd, R):
    y, S = rd[0], rd[1]; B = np.kron(np.eye(2), R); out = (B@y, B@S@B.T)
    return out + ((R@rd[2],) if len(rd) > 2 else ())
def units(rd, lam):
    y, S = rd[0], rd[1]; d = len(y)//2; D = np.diag(np.r_[lam*np.ones(d), np.ones(d)]); return (D@y, D@S@D) + tuple(rd[2:])

S1 = np.diag([1.0]*3 + [0.25]*3)
R1 = [(np.r_[rng.normal(size=3) + [1.0, 0.4, -0.3], rng.normal(size=3)*0.5 + [0.5, 0.2, -0.15]], S1) for _ in range(6)]
L = rng.normal(size=(6, 6)); S2 = L@L.T/6 + 0.3*np.eye(6)          # full covariance, cross-channel correlation included
R2 = [(np.r_[[2.0, 1.0, 0.5], [0.9, 0.6, 0.1]], S2), (np.r_[[-1.0, 2.0, 1.2], [-0.3, 1.1, 0.7]], S2),
      (np.r_[[0.4, -0.8, 1.9], [0.3, -0.2, 1.1]], S2, np.array([[0.2], [-0.35], [0.915]])/np.linalg.norm([0.2, -0.35, 0.915]))]

print("1. EVIDENCE against the closed form 1F1(1/2; 1+k/2; H/2) (arcsine)")
for H, k in [(0.5, 3), (7.0, 3), (40.0, 1), (200.0, 2)]:
    lev = lse(log_arcsine() + 0.5*k*L1B + 0.5*BETA*H, 0) - lse(log_arcsine() + 0.5*k*L1B, 0)
    print(f"  H={H:<6} k={k}: grid {lev:.8f}  closed form {np.log(hyp1f1(0.5, 1+k/2, H/2)):.8f}")

print("\n2. ONE PAIR, six weak 3D readings, c = 2")
b1 = ours28(R1, c=2.0); show("ours28", b1); show("ours27 (flat)", ours28(R1, c=2.0, flat=True))
s1 = ours28([swap(r) for r in R1], c=0.5); print(f"  swap: readout product {b1['readout']*s1['readout']:.6f}, median product {b1['median']*s1['median']:.6f}")
z1 = ours28(R1 + [(np.zeros(6), S1)], c=2.0)
print(f"  + one exact 3D zero reading: law max change {np.abs(z1['law']-b1['law']).max()/b1['law'].max():.1e}; readout {b1['readout']:.6f} -> {z1['readout']:.6f}")
z0 = ours28([(np.zeros(6), S1)]*3, c=2.0); print(f"  zero experiment alone: median {z0['median']:.6f} (c = 2), readout {z0['readout']:.6f} (s = 2)")
split = [(np.r_[r[0][k], r[0][3+k]], np.diag([S1[k, k], S1[3+k, 3+k]])) for r in R1 for k in range(3)]
sp = ours28(split, c=2.0); show("axis split into 18 separate 1D readings (scale NOT carried)", sp)

print("\n3. FULL COVARIANCE with cross-channel correlation; third reading declares a 1D subspace")
b2 = ours28(R2, c=2.0); show("ours28", b2); show("ours27 (flat)", ours28(R2, c=2.0, flat=True))
rr = ours28([rotate(r, rot3()) for r in R2], c=2.0); print(f"  each reading rotated with its covariance and subspace: readout diff {rr['readout']-b2['readout']:+.1e}, median diff {rr['median']-b2['median']:+.1e}")
s2 = ours28([swap(r) for r in R2], c=0.5); print(f"  swap: readout product {b2['readout']*s2['readout']:.6f}, median product {b2['median']*s2['median']:.6f}")
u2 = ours28([units(r, 3.0) for r in R2], c=6.0); print(f"  units (force x3): readout ratio {u2['readout']/b2['readout']:.6f}, median ratio {u2['median']/b2['median']:.6f}")
z2 = ours28(R2 + [(np.zeros(6), S2)], c=2.0); print(f"  + exact zero reading: law max change {np.abs(z2['law']-b2['law']).max()/b2['law'].max():.1e}")

print("\n4. DECLARED SUBSPACE: one 3D reading, push along x (F = 4, a = 2, unit noise), transverse axes read 0")
Su = np.eye(6); y = np.r_[[4.0, 0, 0], [2.0, 0, 0]]
for label, N in [("k = 3 (direction undeclared)", None), ("k = 1 along x (declared)", np.array([[1.0], [0], [0]]))]:
    r = ours28([(y, Su) if N is None else (y, Su, N)], c=1.0); show(label, r)
r1d = ours28([(np.array([4.0, 2.0]), np.eye(2))], c=1.0); show("a 1D reading with the same numbers", r1d)

print("\n5. TOTALITY: awkward single readings (unit noise, c = 1)")
for label, F, a in [("perpendicular", [2.0, 0, 0], [0, 1.0, 0]), ("anti-parallel", [2.0, 0, 0], [-1.0, 0, 0]),
                    ("force only", [3.0, 0, 0], [0, 0, 0]), ("acceleration only", [0, 0, 0], [3.0, 0, 0])]:
    r = ours28([(np.r_[F, a], np.eye(6))], c=1.0); print(f"  {label:<20} readout {r['readout']:.4f}  median {r['median']:.4f}  finite: {np.isfinite(r['readout'])}")
