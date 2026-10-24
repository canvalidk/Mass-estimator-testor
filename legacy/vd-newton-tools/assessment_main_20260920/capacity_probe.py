"""Bounded runtime/memory probe, not an estimator-performance study."""
import argparse
import ctypes
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "03_trace_and_evaluator/mass_estimator/investigation"))
from compare_estimators import measurement_noise, concatenate_inference, refinement_check
from comparison_diagonal_posterior import infer_batch
from comparison_posterior import infer_batch as isotropic_infer


class MemoryCounters(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + [
        (name, ctypes.c_size_t) for name in ("PeakWorkingSetSize", "WorkingSetSize",
        "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
        "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]


def peak_memory():
    counters = MemoryCounters()
    counters.cb = ctypes.sizeof(counters)
    current_process = ctypes.windll.kernel32.GetCurrentProcess
    current_process.restype = ctypes.c_void_p
    memory_info = ctypes.windll.psapi.GetProcessMemoryInfo
    memory_info.argtypes = [ctypes.c_void_p, ctypes.POINTER(MemoryCounters), ctypes.c_ulong]
    if not memory_info(current_process(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError()
    return counters.PeakWorkingSetSize


parser = argparse.ArgumentParser()
parser.add_argument("--samples", type=int, default=10000)
parser.add_argument("--batch-size", type=int, default=256)
parser.add_argument("--order", type=int, default=96)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
if args.output.exists():
    raise FileExistsError(args.output)

config = {"measurement": {"repeats": 4, "calibration": {"mode": "estimated", "samples": 30, "seed": 2026092011}}}
force_sd, acceleration_sd = np.array([.5, 1., 2.]), np.array([2., 1., .5])
sf, sa, metadata = measurement_noise(config, args.samples, force_sd, acceleration_sd)
noise = np.random.default_rng(2026092012).standard_normal((args.samples, 6))
force = noise[:, :3] * force_sd / 2 + [8., 0., 0.]
acceleration = noise[:, 3:] * acceleration_sd / 2 + [2., 0., 0.]
began = time.perf_counter()
chunks = []
for start in range(0, args.samples, args.batch_size):
    end = start + args.batch_size
    chunks.append(infer_batch(force[start:end], acceleration[start:end], 4.,
                             force_sd=sf[start:end], acceleration_sd=sa[start:end], order=args.order))
    if len(chunks) % 10 == 0:
        print(json.dumps({"completed": min(end, args.samples), "seconds": time.perf_counter()-began}), flush=True)
inferred = concatenate_inference(chunks)
primary_seconds = time.perf_counter()-began
refined = infer_batch(force[:32], acceleration[:32], 4., force_sd=sf[:32], acceleration_sd=sa[:32], order=2*args.order)
refinement = refinement_check(inferred, refined)
primary_peak = peak_memory()
began_rms = time.perf_counter()
rms_f = np.hypot.reduce(sf, axis=1) / np.sqrt(3)
rms_a = np.hypot.reduce(sa, axis=1) / np.sqrt(3)
rms_points = []
for start in range(0, args.samples, args.batch_size):
    end = start + args.batch_size
    result = isotropic_infer(force[start:end]/rms_f[start:end, None], acceleration[start:end]/rms_a[start:end, None],
                            4./(rms_f[start:end]/rms_a[start:end]), order=args.order)
    rms_points.append(result["points"]["flat_joint"]*(rms_f[start:end]/rms_a[start:end]))
report = {"purpose": "Resource capacity probe only; all inference settings unchanged", "samples": args.samples,
          "batch_size": args.batch_size, "order": args.order, "numpy": np.__version__,
          "measurement": config["measurement"], "primary_seconds": primary_seconds,
          "batched_rms_seconds": time.perf_counter()-began_rms,
          "primary_peak_working_set_bytes": primary_peak, "overall_peak_working_set_bytes": peak_memory(),
          "all_primary_points_finite_positive": all(np.all(np.isfinite(v)&(v>0)) for v in inferred["points"].values()),
          "all_rms_points_finite_positive": bool(np.all(np.isfinite(np.concatenate(rms_points)) & (np.concatenate(rms_points)>0))),
          "refinement": refinement}
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("x", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, allow_nan=False)
print(json.dumps(report), flush=True)
