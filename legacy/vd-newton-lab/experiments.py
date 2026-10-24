"""Reproduce the case table, sweeps and seeded sampling checks. No network."""
from pathlib import Path
import csv
import json
import math
import time
import numpy as np
from estimator import estimate, isotropic_covariance, isotropic_qmin, profile, compatibility_intervals, gate

ROOT = Path(__file__).resolve().parent


def run_case(name, force, acceleration, sf=1., sa=1.):
    cov = isotropic_covariance(2, sf, sa)
    e = estimate(force, acceleration, cov, direction_order=256, ratio_order=2048)
    coarse = estimate(force, acceleration, cov, direction_order=128, ratio_order=1024)
    qmin = isotropic_qmin(force, acceleration, sf, sa)
    intervals = compatibility_intervals(force, acceleration, sf, sa, 9)
    x = np.linspace(-8, 8, 321)
    theta = np.arctan(np.exp(x))
    density = np.interp(theta, e.ratio_angle, e.angle_density) * np.sin(theta) * np.cos(theta)
    profile_values = [profile(force, acceleration, cov, e.mass_scale * math.exp(v)) for v in x]
    return dict(name=name, force=force, acceleration=acceleration, sigma_force=sf, sigma_acceleration=sa,
                mass=e.mass, mean_force=e.mean_force, mean_acceleration=e.mean_acceleration,
                lower=e.quantile(.025), median=e.quantile(.5), upper=e.quantile(.975),
                mass_scale=e.mass_scale, qmin=qmin, q_at_point=profile(force, acceleration, cov, e.mass),
                intervals=[[lo, hi if math.isfinite(hi) else None] for lo, hi in intervals],
                point_in_compatibility_set=any(lo <= e.mass <= hi for lo, hi in intervals),
                gate=gate(e.mass, qmin, 9)["status"],
                refinement_mass_relative=abs(coarse.mass / e.mass - 1),
                refinement_upper_relative=abs(coarse.quantile(.975) / e.quantile(.975) - 1),
                log_ratio=x.tolist(), log_density=density.tolist(), profile=profile_values,
                displayed_probability=float(np.trapezoid(density, x)))


def batch_qmin(f, a):
    p, q, dot = np.sum(f * f, axis=1), np.sum(a * a, axis=1), np.sum(f * a, axis=1)
    root = np.hypot(p - q, 2 * dot)
    answer = np.minimum(p, q)
    positive = dot > 0
    answer[positive] = np.maximum(0, 2 * (p[positive] * q[positive] - dot[positive] ** 2) / (p[positive] + q[positive] + root[positive]))
    return answer


def main():
    started = time.monotonic()
    p = math.sqrt(6)
    specs = [("Measured zero-zero", [0, 0], [0, 0], 1, 1),
             ("Zero-zero, force uncertainty doubled", [0, 0], [0, 0], 2, 1),
             ("Force zero, acceleration resolved", [0, 0], [10, 0], 1, 1),
             ("Force resolved, acceleration zero", [10, 0], [0, 0], 1, 1),
             ("One-zero but unresolved", [0, 0], [2, 0], 1, 1),
             ("Strong aligned", [20, 0], [10, 0], 1, 1),
             ("Weak aligned", [1, 0], [1, 0], 1, 1),
             ("Balanced perpendicular", [0, p], [p, 0], 1, 1),
             ("Anti-aligned, rejected", [-10, 0], [10, 0], 2, 2),
             ("Anti-aligned, two compatible ranges", [-p, 0], [p, 0], 1, 1),
             ("Anti-aligned, force less resolved", [-2, 0], [5, 0], 1, 1),
             ("Moderate angular mismatch", [4, 2], [2, 0], 1, 1)]
    cases = [run_case(*spec) for spec in specs]
    angles = []
    for degrees in range(0, 181, 5):
        t = math.radians(degrees)
        angles.append(run_case(f"Angle {degrees} degrees", [5 * math.cos(t), 5 * math.sin(t)], [5, 0]))
        angles[-1]["angle"] = degrees
    excitation = []
    for snr in (0, .1, .25, .5, 1, 2, 3, 5, 10, 20):
        excitation.append(run_case(f"Equal-channel SNR {snr:g}", [snr, 0], [snr, 0]))
        excitation[-1]["snr"] = snr
    rng = np.random.default_rng(20260906)
    calibration = []
    n = 100000
    for dimension in (1, 2, 3):
        for signal in (0., 1., 10.):
            f, a = rng.normal(size=(2, n, dimension))
            f[:, 0] += signal
            a[:, 0] += signal
            scores = batch_qmin(f, a)
            rejected = int(np.count_nonzero(scores > 9))
            calibration.append(dict(dimension=dimension, true_signal_per_channel=signal, n=n,
                                    q95=float(np.quantile(scores, .95)), q99=float(np.quantile(scores, .99)),
                                    rejection_fraction_at_9=rejected / n,
                                    rejection_standard_error=math.sqrt((rejected / n) * (1 - rejected / n) / n)))
    # Repeated-sampling distribution: rerun the entire two-dimensional estimator
    # on independent trials at a stated true latent pair. No bootstrap ambiguity.
    sampling = []
    for signal in (0., 1., 5.):
        n = 400
        readings = rng.normal(size=(n, 4)) + [2 * signal, 0, signal, 0]
        masses, lower, upper, rejected = [], [], [], 0
        for f1, f2, a1, a2 in readings:
            e = estimate([f1, f2], [a1, a2], np.eye(4), direction_order=64, ratio_order=512)
            masses.append(e.mass)
            lower.append(e.quantile(.025))
            upper.append(e.quantile(.975))
            rejected += isotropic_qmin([f1, f2], [a1, a2], 1, 1) > 9
        masses = np.asarray(masses)
        covered = np.count_nonzero((np.asarray(lower) <= 2) & (np.asarray(upper) >= 2))
        sampling.append(dict(true_mass=2, true_acceleration=signal, n=n,
                             mean_candidate=float(masses.mean()), median_candidate=float(np.median(masses)),
                             candidate_sd=float(masses.std(ddof=1)), mean_candidate_mc_se=float(masses.std(ddof=1) / math.sqrt(n)),
                             induced_interval_coverage=covered / n,
                             coverage_mc_se=math.sqrt((covered / n) * (1 - covered / n) / n),
                             rejection_fraction=rejected / n,
                             note="Candidates include rejected trials; induced intervals are not asserted to have 95% sampling coverage."))
    result = dict(date="2026-09-06", seed=20260906,
                  managed_commit="743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6",
                  model="Gaussian; independent isotropic channels; 2D unknown common direction; flat positive magnitudes",
                  threshold=9, threshold_status="Illustrative absolute discrepancy threshold, not a universal 3-sigma rule",
                  cases=cases, angles=angles, excitation=excitation, calibration=calibration, sampling=sampling,
                  runtime_seconds=time.monotonic() - started)
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    (out / "experiments.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    columns = ["name", "mass", "lower", "median", "upper", "qmin", "q_at_point", "gate", "point_in_compatibility_set", "refinement_mass_relative", "refinement_upper_relative"]
    with (out / "cases.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(cases)
    compact = [{k: case[k] for k in columns} for case in cases]
    print(json.dumps(dict(cases=compact, calibration=calibration, sampling=sampling,
                          max_mass_refinement=max(c["refinement_mass_relative"] for c in cases + angles + excitation),
                          max_upper_refinement=max(c["refinement_upper_relative"] for c in cases + angles + excitation),
                          runtime_seconds=result["runtime_seconds"]), indent=2))


if __name__ == "__main__":
    main()
