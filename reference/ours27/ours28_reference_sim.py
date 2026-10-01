import numpy as np, sys, time
from ours28_reference import ours28
def world(N, omega, rng, m0=2.0, sF=1.0, sa=0.5):
    s = sF/sa; th = np.arctan(m0/s); S = np.diag([sF**2]*3 + [sa**2]*3); out = []
    for _ in range(N):
        v = rng.normal(size=3)*omega                           # push in noise units
        x = v*np.sin(th) + rng.normal(size=3); y = v*np.cos(th) + rng.normal(size=3)
        out.append((np.r_[x*sF, y*sa], S))
    return out
if __name__ == "__main__":
    N, omega, reps, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(seed); m0 = 2.0; res = []; t0 = time.time()
    for _ in range(reps):
        r = ours28(world(N, omega, rng), c=2.0, lm_half=6.0, nodes=481)
        res.append([r['readout'], r['median'], r['lo'], r['hi']])
    res = np.array(res); lr = np.log(res[:, 0]/m0); lmed = np.log(res[:, 1]/m0)
    cover = np.mean((res[:, 2] <= m0) & (m0 <= res[:, 3])); law_sd = np.median(np.log(res[:, 3]/res[:, 2]))/3.92
    print(f"N={N:<4} omega={omega}: law coverage {cover:.3f} | readout/m0 median {np.exp(np.median(lr)):.3f}, "
          f"median/m0 {np.exp(np.median(lmed)):.3f} | SD(log readout) {lr.std():.3f} vs law SD {law_sd:.3f} | {reps} series, {time.time()-t0:.0f}s", flush=True)
