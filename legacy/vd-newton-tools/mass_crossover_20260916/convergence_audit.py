"""Fixed physical pair / vanishing noise check for the 3D mass readout.

Data remain F*=8e1, a*=3e1 (mass 8/3). Both standard uncertainties are 1/k.
The exact angular moment integrals are split at the observed likelihood peak,
not at a true parameter. Common draws provide paired error comparisons.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "03_trace_and_evaluator/mass_estimator/investigation/crossover_vector_checks.py"
spec = importlib.util.spec_from_file_location("crossover", SOURCE)
cross = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cross)


def point(F, A, order=96):
    p = np.sum(F*F, axis=1)
    q = np.sum(A*A, axis=1)
    d = np.sum(F*A, axis=1)
    theta = np.where(d > 0, .5*np.arctan2(2*d, q-p),
                     np.where(p >= q, math.pi/2, 0.))
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x+1)/2, w/2
    t = np.concatenate((theta[:, None]*x,
                        theta[:, None]+(math.pi/2-theta[:, None])*x), axis=1)
    weights = np.concatenate((theta[:, None]*w,
                              (math.pi/2-theta[:, None])*w), axis=1)
    sn, cs = np.sin(t), np.cos(t)
    h2 = p[:, None]*sn*sn + q[:, None]*cs*cs + 2*d[:, None]*sn*cs
    peak = np.where(d >= 0, (p+q+np.hypot(p-q,2*d))/2, np.maximum(p,q))
    density = np.exp((h2-peak[:, None])/2)*weights
    return np.sum(sn*density,axis=1)/np.sum(cs*density,axis=1)


def summarize(err):
    return {"log_mse": float(np.mean(err*err)),
            "mean_log_error": float(np.mean(err)),
            "log_variance": float(np.var(err))}


N = 65536
rng = np.random.default_rng(2026091601)
errors = rng.standard_normal((N,6))
truth = 8/3
results = []
for k in [1,2,4,8,16]:
    F = errors[:,:3] + [8*k,0,0]
    A = errors[:,3:] + [3*k,0,0]
    ours = np.concatenate([point(F[i:i+512],A[i:i+512]) for i in range(0,N,512)])
    norm = np.linalg.norm(F,axis=1)/np.linalg.norm(A,axis=1)
    ref = point(F[:1024], A[:1024], order=192)
    refinement = float(np.max(np.abs(ref/ours[:1024]-1)))
    original, _ = cross.infer(F[:256], A[:256], order=192)
    independent = float(np.max(np.abs(original[:,0]/ours[:256]-1)))
    assert refinement < 1e-7, (k,refinement)
    assert independent < 1e-7, (k,independent)
    eo = np.log(ours/truth)
    en = np.log(norm/truth)
    difference = eo*eo-en*en
    first_order = (errors[:,0]/8-errors[:,3]/3)/k
    delta = np.log(ours/norm)
    direction_cosine = np.sum(F*A,axis=1)/(np.linalg.norm(F,axis=1)*np.linalg.norm(A,axis=1))
    p, q = np.sum(F*F,axis=1), np.sum(A*A,axis=1)
    angular_expansion = (p-q)/(p+q)*(1-direction_cosine)
    results.append({"precision_multiplier":k,"physical_sigmas":1/k,
                    "force_snr":8*k,"acceleration_snr":3*k,
                    "ours":summarize(eo), "norm":summarize(en),
                    "risk_ratio":float(np.mean(eo*eo)/np.mean(en*en)),
                    "paired_risk_gap":float(np.mean(difference)),
                    "paired_gap_mcse":float(np.std(difference,ddof=1)/math.sqrt(N)),
                    "rms_log_estimator_difference":float(np.sqrt(np.mean(delta*delta))),
                    "rms_ours_minus_first_order":float(np.sqrt(np.mean((eo-first_order)**2))),
                    "rms_norm_minus_first_order":float(np.sqrt(np.mean((en-first_order)**2))),
                    "rms_angular_expansion_residual":float(np.sqrt(np.mean((delta-angular_expansion)**2))),
                    "quadrature_refinement":refinement,"original_integral_check":independent})
    print(json.dumps(results[-1]),flush=True)

report={"fixed_force":[8,0,0],"fixed_acceleration":[3,0,0],"fixed_mass":truth,
        "sample_count":N,"seed":2026091601,
        "assumptions":"Unknown common direction; independent isotropic Gaussian errors; both sigmas=1/k. No pooling across trials.",
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "results":results}
Path(__file__).with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
