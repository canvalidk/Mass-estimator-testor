"""Reanalyse published air-track timings and compare identical-input readouts.

The scalar system coordinate is string travel, not the Cartesian COM of the
glider plus hanging load. This tests the readout after a mechanical reduction.
Uncertainty scenarios are explicitly supplied by this analysis, not the paper.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from estimator import estimate, isotropic_covariance

HERE = Path(__file__).resolve().parent
SOURCE = json.loads((HERE / "air_track_source.json").read_text())


def acceleration(intervals, spacing):
    """OLS slope of average speed against interval midpoint, batched on last axis."""
    time = np.cumsum(intervals, axis=-1) - intervals / 2
    speed = spacing / intervals
    tc = time - time.mean(axis=-1, keepdims=True)
    vc = speed - speed.mean(axis=-1, keepdims=True)
    return (tc * vc).sum(axis=-1) / (tc * tc).sum(axis=-1)


def exact_1d_point(force, accel, sf, sa):
    """Analytic independent-Gaussian S^0 version of the full latent integral.

    Direction +/- has equal measure. D=Phi(F/sf)+Phi(a/sa)-1;
    normalized moment denominators cancel, leaving the expression below.
    Used here only for positive readings; no general negative-tail claims.
    """
    zf, za = force / sf, accel / sa
    cdf = lambda z: .5 * math.erfc(-z / math.sqrt(2))
    phi = lambda z: math.exp(-z*z/2) / math.sqrt(2*math.pi)
    d = cdf(zf) + cdf(za) - 1
    return (force*d + sf*phi(zf)) / (accel*d + sa*phi(za))


def run():
    intervals = np.array(SOURCE["intervals_ms"], float) / 1000
    spacing = SOURCE["post_spacing_m"]
    gravity = SOURCE["g_m_s2"]
    loads = np.array(SOURCE["hanging_load_g"]) / 1000
    forces = loads * gravity
    accels = acceleration(intervals, spacing)
    # Reproduce Table 2 as a transcription check. Do NOT round speeds before fit.
    assert np.array_equal(np.round(spacing / intervals[0] * 100, 1),
                          SOURCE["first_run_table_2_speed_cm_s"])
    # Weights implement the slope per 50g step and remove constant resistance.
    # sum(w)=0, w @ [0,1,2,3,4]=1. Same contrast applied to F and a.
    w = np.array([-2., -1., 0., 1., 2.]) / 10
    f_effective, a_effective = float(w @ forces), float(w @ accels)
    baseline_mass = f_effective / a_effective
    intercept = float(accels.mean() - a_effective * 3)
    resistance = -baseline_mass * intercept
    residuals = accels - (intercept + a_effective * np.arange(1, 6))

    # Analyst timing scenario: six independent photogate crossing-time errors,
    # uniform +/- half a clock period. Differencing gives correlated intervals.
    # Shared 0.1mm spacing quantization is also assumed, not a calibration budget.
    rng = np.random.default_rng(20260906)
    count = 100000
    crossing_error = rng.uniform(-.0005, .0005, (count, 5, 6))
    sampled_intervals = intervals + np.diff(crossing_error, axis=-1)
    sampled_spacing = spacing + rng.uniform(-.00005, .00005, (count, 1, 1))
    sampled_accels = acceleration(sampled_intervals, sampled_spacing)
    sampled_contrast = sampled_accels @ w
    sigma_accels = sampled_accels.std(axis=0, ddof=1)
    sigma_a = float(sampled_contrast.std(ddof=1))
    # Separate sensitivity scenario: estimate slope error from between-run scatter.
    sigma_a_scatter = float(np.sqrt((residuals @ residuals) / 3 / 10))

    checks, scenarios = [], []
    for uncertainty_source, sa in [("timing_and_spacing_scenario", sigma_a),
                                   ("force_line_residual_scenario", sigma_a_scatter)]:
        for fraction in [.001, .01, .05]:
            sf = fraction * f_effective
            analytic = exact_1d_point(f_effective, a_effective, sf, sa)
            coarse = estimate([f_effective], [a_effective],
                              isotropic_covariance(1, sf, sa), ratio_order=8192)
            ratio_order = 16384
            fine = estimate([f_effective], [a_effective],
                            isotropic_covariance(1, sf, sa), ratio_order=ratio_order)
            quantile_change = max(abs(coarse.quantile(p)-fine.quantile(p))
                                  for p in [.025, .975])
            while quantile_change > .000005 and ratio_order < 262144:
                coarse = fine
                ratio_order *= 2
                fine = estimate([f_effective], [a_effective],
                                isotropic_covariance(1, sf, sa), ratio_order=ratio_order)
                quantile_change = max(abs(coarse.quantile(p)-fine.quantile(p))
                                      for p in [.025, .975])
            rel_error = abs(fine.mass / analytic - 1)
            assert rel_error < 1e-8, (uncertainty_source, fraction, rel_error)
            assert abs(coarse.mass / fine.mass - 1) < 1e-8
            assert quantile_change < .000005
            scenarios.append({"acceleration_uncertainty_model": uncertainty_source,
                "assumed_relative_force_sigma": fraction,
                "sigma_force_N": sf, "sigma_acceleration_m_s2": sa,
                "ordinary_ratio_kg": baseline_mass, "integral_mass_kg": fine.mass,
                "analytic_mass_kg": analytic,
                "integral_minus_ratio_g": 1000 * (fine.mass - baseline_mass),
                "ratio_distribution_equal_tail_95_percent_kg":
                    [fine.quantile(.025), fine.quantile(.975)]})
            checks.append({"model": uncertainty_source, "force_fraction": fraction,
                           "relative_analytic_error": rel_error,
                           "relative_refinement_change": abs(coarse.mass/fine.mass-1),
                           "final_ratio_order": ratio_order,
                           "max_quantile_refinement_change_kg": quantile_change})

    # Per-run diagnostic: mg is driving force, NOT corrected generalized net force.
    # These apparent masses isolate preprocessing/drag from estimator effects.
    rows = []
    for i in range(5):
        sf, sa = .01 * forces[i], float(sigma_accels[i])
        result = estimate([forces[i]], [accels[i]],
                          isotropic_covariance(1, sf, sa), ratio_order=16384)
        baseline = forces[i] / accels[i]
        assert abs(result.mass / baseline - 1) < 1e-8
        rows.append({"hanging_load_g": float(loads[i]*1000),
                     "driving_force_N": float(forces[i]),
                     "acceleration_m_s2": float(accels[i]),
                     "assumed_sigma_acceleration_m_s2": sa,
                     "uncorrected_ratio_mass_g": float(baseline*1000),
                     "uncorrected_integral_mass_g": result.mass*1000})

    # Deliberately enlarged error bars on the same central contrast, NOT real data.
    sensitivity = []
    for af in [.05, .1, .2, .5, 1.]:
        mass = exact_1d_point(f_effective, a_effective, .01*f_effective, af*a_effective)
        sensitivity.append({"hypothetical_relative_acceleration_sigma": af,
                            "mass_kg": mass,
                            "change_percent": (mass/baseline_mass-1)*100})

    out = {"source": SOURCE, "analysis_date": "2026-09-06",
           "model": "1D string-travel generalized force/acceleration; flat positive magnitudes and uniform S^0 direction; independent Gaussian contrast inputs",
           "assumptions": ["Massless inextensible string, negligible pulley inertia, common string-travel acceleration; constant resistance across the five runs",
                           "Force-error scenarios are a common multiplicative scale error independent of timing, not independent errors on all hanging loads",
                           "Independent uniform crossing-time errors and a shared uniform spacing error are analyst assumptions; Gaussian summary used only at the estimator input",
                           "No complete systematic uncertainty budget, no calibrated coverage claim, and no 3D vector admission claim"],
           "contrast_weights": w.tolist(), "contrast_force_N": f_effective,
           "contrast_acceleration_m_s2": a_effective,
           "reconstructed_slope_cm_s2_per_g": a_effective*2,
           "reconstructed_intercept_cm_s2": intercept*100,
           "reconstructed_mass_kg": baseline_mass,
           "mass_from_published_rounded_slope_kg": gravity/(SOURCE["reported_acceleration_fit_cm_s2"]["slope_per_g_hanging_load"]*10),
           "fitted_resistance_N": resistance,
           "residual_accelerations_m_s2": residuals.tolist(),
           "timing_mc_draws": count, "timing_mc_seed": 20260906,
           "timing_sigma_contrast_m_s2": sigma_a,
           "scatter_sigma_contrast_m_s2": sigma_a_scatter,
           "scenarios": scenarios, "per_run_diagnostics": rows,
           "hypothetical_error_sensitivity_not_measured": sensitivity,
           "numerical_checks": checks,
           "sha256": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ["air_track_source.json", "air_track_experiment.py", "estimator.py"]}}
    results = HERE / "results"
    results.mkdir(exist_ok=True)
    (results / "air_track.json").write_text(json.dumps(out, indent=2)+"\n")
    with (results / "air_track_per_run.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({key: out[key] for key in ["reconstructed_mass_kg", "reconstructed_slope_cm_s2_per_g", "reconstructed_intercept_cm_s2", "fitted_resistance_N", "timing_sigma_contrast_m_s2", "scatter_sigma_contrast_m_s2", "scenarios", "hypothetical_error_sensitivity_not_measured"]}, indent=2))


if __name__ == "__main__":
    run()
