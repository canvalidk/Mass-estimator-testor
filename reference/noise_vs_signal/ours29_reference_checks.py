"""Checks on the ours29 reference implementation (full covariance, declared subspace, misfit weight, Bartlett power).

Run from this folder:  python ours29_reference_checks.py
Needs ../ours27/ours28_reference.py for check 1 (ours29 with both tempering steps off is ours27).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ours27"))
from ours29_reference import ours29, bartlett_eta
rng = np.random.default_rng(7)


def show(name, r):
    print(f"  {name:<50} readout {r.get('readout', float('nan')):.6f}  median {r['median']:.6f}  "
          f"95% [{r['lo']:.4f}, {r['hi']:.4f}]  eta {r['eta']:.4f}")


def rot3():
    return np.linalg.qr(rng.normal(size=(3, 3)))[0]


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
# a weak series where both tempering steps act: 3D pure noise under the full covariance S2
Lc = np.linalg.cholesky(S2)
R3 = [(Lc@rng.normal(size=6), S2) for _ in range(40)]

print("1. BOTH STEPS OFF IS ours27 (against the ours28 reference's flat branch)")
from ours28_reference import ours28
for name, R in (("R1", R1), ("R2", R2)):
    a = ours28(R, c=2.0, flat=True); b = ours29(R, c=2.0, misfit=False, calibrate=False)
    print(f"  {name}: ours28_reference(flat) readout {a['readout']:.10f} median {a['median']:.10f} | ours29(off) {b['readout']:.10f} {b['median']:.10f}")

print("\n2. 1D: eta against the carried-state closed form (raw J, cap at 1)")
F, a = rng.normal(size=40) + .2, rng.normal(size=40) + .1
Ci = (a**2 - F**2) + 2j*F*a; C = Ci.sum(); th = np.angle(C)/2
Jc = (np.sum(abs(Ci)**2) - np.real(np.exp(4j*th)*np.conj(np.sum(Ci**2))))/8
e = bartlett_eta([(np.array([F[i], a[i]]), np.eye(2)) for i in range(40)], 1.0, [1.0]*40)[0]
print(f"  numeric {e:.10f}  closed form {min(1, abs(C)/Jc):.10f}")

print("\n3. THE WEAK 3D SERIES (40 readings of pure noise, full covariance with cross-channel correlation)")
b3 = ours29(R3, c=2.0); show("ours29", b3)
show("ours29, eta = 1 (the raw law)", ours29(R3, c=2.0, calibrate=False))
show("ours27", ours29(R3, c=2.0, misfit=False, calibrate=False))
print(f"  misfit weights: min {min(b3['g']):.3f}, {sum(g < 1 for g in b3['g'])} of 40 below 1")
rr = ours29([rotate(r, rot3()) for r in R3], c=2.0)
print(f"  each reading rotated with its covariance: eta diff {rr['eta']-b3['eta']:+.1e}, readout diff {rr['readout']-b3['readout']:+.1e}, median diff {rr['median']-b3['median']:+.1e}")
s3 = ours29([swap(r) for r in R3], c=0.5)
print(f"  swap: eta diff {s3['eta']-b3['eta']:+.1e}, readout product {b3['readout']*s3['readout']:.6f}, median product {b3['median']*s3['median']:.6f}")
u3 = ours29([units(r, 3.0) for r in R3], c=6.0)
print(f"  units (force x3): readout ratio {u3['readout']/b3['readout']:.6f}, median ratio {u3['median']/b3['median']:.6f}")
z3 = ours29(R3 + [(np.zeros(6), S2)]*10, c=2.0)
print(f"  + ten exact zero readings: eta diff {z3['eta']-b3['eta']:+.1e}, law max change {np.abs(z3['law']-b3['law']).max()/b3['law'].max():.1e}, "
      f"median {b3['median']:.6f} -> {z3['median']:.6f}, readout {b3['readout']:.6f} -> {z3['readout']:.6f} (T29: the cloud readout counts them)")

print("\n4. DECLARED SUBSPACE: one 3D reading along x (F = 4, a = 2, unit noise), transverse axes read 0")
Su = np.eye(6); y = np.r_[[4.0, 0, 0], [2.0, 0, 0]]
r_k3 = ours29([(y, Su)], c=1.0); show("k = 3", r_k3)
r_k1 = ours29([(y, Su, np.array([[1.0], [0], [0]]))], c=1.0); show("k = 1 along x", r_k1)
show("a 1D reading with the same numbers", ours29([(np.array([4.0, 2.0]), np.eye(2))], c=1.0))
yq = np.r_[[4.0, 1.5, 0], [2.0, -1.0, 0]]
print(f"  a misaligned 3D reading: Q = {ours29([(yq, Su)], c=1.0, readout=False)['Q'][0]:.4f}, "
      f"g = {ours29([(yq, Su)], c=1.0, readout=False)['g'][0]:.4f}; declared k = 1 along x: "
      f"Q = {ours29([(yq, Su, np.array([[1.0], [0], [0]]))], c=1.0, readout=False)['Q'][0]:.2e}, g = 1")

print("\n5. ONE READING: eta = 1 whatever it is")
for yy in (np.r_[[1.0, 0, 0], [0, 1.0, 0]], np.r_[[0.3, -0.2, 0.1], [0.1, 0.2, -0.4]], np.zeros(6)):
    print(f"  y = {np.round(yy, 2)}: eta {ours29([(yy, S2)], c=2.0, readout=False)['eta']}")

print("\n6. THE IMPOSSIBLE PAIR (2D, F = (2, 0), a = (0, 1), sigma -> 0): tilt of the log law, kink form; limit nu a/(1-a) = 1.5")
for sg in (1.0, 0.5, 0.1, 0.005):
    r = ours29([(np.r_[[2.0, 0], [0, 1.0]], sg**2*np.eye(4))], c=2.0, readout=False, lm_half=14, nodes=4001)
    lpi = np.log(2.0/(4.0+np.exp(2*r['lm']))) + r['lm']
    lr = np.log(r['law']) - lpi
    print(f"  sigma {sg:<6} Q {r['Q'][0]:>10.1f}  g {r['g'][0]:.5f}  tilt {np.ptp(lr):.3f}")

print("\n7. TOTALITY: awkward single readings (unit noise, c = 1), and the all-zero experiment")
for label, F_, a_ in [("perpendicular", [2.0, 0, 0], [0, 1.0, 0]), ("anti-parallel", [2.0, 0, 0], [-1.0, 0, 0]),
                      ("force only", [3.0, 0, 0], [0, 0, 0]), ("acceleration only", [0, 0, 0], [3.0, 0, 0])]:
    r = ours29([(np.r_[F_, a_], np.eye(6))], c=1.0); print(f"  {label:<20} readout {r['readout']:.4f}  median {r['median']:.4f}  g {r['g'][0]:.3f}  finite: {np.isfinite(r['readout'])}")
r0 = ours29([(np.zeros(6), np.eye(6))]*3, c=2.0); print(f"  all-zero experiment (3 readings): median {r0['median']:.6f} (c = 2), readout {r0['readout']:.6f}, eta {r0['eta']}")

print("\n8. AXIS SPLITTING: a 3D series and its axes read as separate 1D readings (intended to differ)")
Rw = [(np.r_[rng.normal(size=3)*1.0 + [0.3, 0, 0], rng.normal(size=3)*1.0 + [0.2, 0, 0]], np.eye(6)) for _ in range(30)]
sp = [(np.r_[r[0][k], r[0][3+k]], np.eye(2)) for r in Rw for k in range(3)]
a3, a1 = ours29(Rw, c=1.0, readout=False), ours29(sp, c=1.0, readout=False)
print(f"  30 x 3D: eta {a3['eta']:.4f}, median {a3['median']:.4f}, 95% [{a3['lo']:.3f}, {a3['hi']:.3f}]")
print(f"  90 x 1D: eta {a1['eta']:.4f}, median {a1['median']:.4f}, 95% [{a1['lo']:.3f}, {a1['hi']:.3f}]")
