import numpy as np
from ours28b import ours28b
rng = np.random.default_rng(7)
def diag(sF, sa): return np.diag(np.concatenate([np.atleast_1d(sF), np.atleast_1d(sa)]).astype(float)**2)
def rot3(): return np.linalg.qr(rng.normal(size=(3, 3)))[0]
def rotated(F, a, S, R):
    B = np.kron(np.eye(2), R); return R@np.asarray(F), R@np.asarray(a), B@S@B.T
def swapped(R_):
    out = []
    for F, a, S in R_:
        d = len(F); P = np.block([[np.zeros((d, d)), np.eye(d)], [np.eye(d), np.zeros((d, d))]])
        out.append((a, F, P@S@P.T))
    return out
def show(name, r): print(f"  {name:<44} readout {r['readout']:.6f}  median {r['median']:.6f}  95% [{r['lo']:.4f}, {r['hi']:.4f}]")

# weak-ish one-pair data, 3D, several readings
S1 = diag([1, 1, 1], [0.5, 0.5, 0.5])
R1 = [(rng.normal(size=3)*1.0 + np.array([1.0, 0.4, -0.3]), rng.normal(size=3)*0.5 + np.array([0.5, 0.2, -0.15]), S1) for _ in range(6)]
print("ONE PAIR, six 3D readings, c = 2")
base = ours28b(R1, c=2.0); show("ours28b", base); show("ours27 (flat limit)", ours28b(R1, c=2.0, flat=True))
sw = ours28b(swapped(R1), c=0.5); print(f"  swap: readout product {base['readout']*sw['readout']:.6f}, median product {base['median']*sw['median']:.6f}")
un = ours28b([(3*np.asarray(F), a, np.diag(np.r_[9*np.ones(3), np.ones(3)])@S) for F, a, S in R1], c=6.0)
print(f"  units (force x3): readout ratio {un['readout']/base['readout']:.6f}")
rt = ours28b([rotated(F, a, S, rot3()) for F, a, S in R1], c=2.0)
print(f"  each reading rotated in its own frame: readout {rt['readout']:.6f} (diff {rt['readout']-base['readout']:+.1e})")
sp = ours28b([([F[k]], [a[k]], np.diag([S[k, k], S[3+k, 3+k]])) for F, a, S in R1 for k in range(3)], c=2.0)
print(f"  axis splitting (6 x 3D as 18 x 1D): readout {sp['readout']:.6f} (diff {sp['readout']-base['readout']:+.1e})")
z = ours28b(R1 + [(np.zeros(3), np.zeros(3), S1)], c=2.0)
print(f"  strong zero-reading property: + one exact 3D zero reading -> readout {z['readout']:.6f}, median {z['median']:.6f}, 95% [{z['lo']:.4f}, {z['hi']:.4f}]")
z0 = ours28b([(np.zeros(3), np.zeros(3), S1)]*3, c=2.0)
print(f"  weak zero-reading property: zero experiment alone -> median {z0['median']:.6f} (starting law median 2), readout {z0['readout']:.6f} (s = 2)")

print("\nPER-AXIS, three 3D readings, sigma_F = 1, sigma_a = (0.2, 0.5, 1.25), c = 2")
S2 = diag([1, 1, 1], [0.2, 0.5, 1.25])
R2 = [(np.array([2.0, 1.0, 0.5]), np.array([0.9, 0.6, 0.1]), S2), (np.array([-1.0, 2.0, 1.2]), np.array([-0.3, 1.1, 0.7]), S2),
      (np.array([0.4, -0.8, 1.9]), np.array([0.3, -0.2, 1.1]), S2)]
b2 = ours28b(R2, c=2.0); show("ours28b", b2); show("ours27 (flat limit)", ours28b(R2, c=2.0, flat=True))
rt2 = ours28b([rotated(F, a, S, rot3()) for F, a, S in R2], c=2.0)
print(f"  each reading rotated with its covariance: readout {rt2['readout']:.6f} (diff {rt2['readout']-b2['readout']:+.1e})")
sw2 = ours28b(swapped(R2), c=0.5); print(f"  swap: readout product {b2['readout']*sw2['readout']:.6f}")
print(f"  classes: {[(round(g['sF'],3), round(g['sa'],3), g['M']) for g in b2['classes']]}")

print("\nPRECISE NULLS -> KNOWN EMPTY: active 1D reading F=3, a=1 (unit noise), plus two 0/0 axes with noise eps (both channels)")
one = ours28b([([3.0], [1.0], np.diag([1.0, 1.0]))], c=1.0); show("1D declared", one)
for eps in [1.0, 0.1, 0.01, 0.001]:
    S3 = np.diag([1.0, eps**2, eps**2, 1.0, eps**2, eps**2])
    r3 = ours28b([([3.0, 0.0, 0.0], [1.0, 0.0, 0.0], S3)], c=1.0); show(f"3D, empty axes noise {eps}", r3)
