"""Compare mass points and uncertainty laws against independent simulation truth.

This is the independent, isotropic, unknown-direction 3D experiment. It does
not change the full-covariance production estimator or discard misaligned
measurements. Requires Python and NumPy. See comparison_guide.md.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np

from comparison_estimators import METHOD_METADATA, baseline_points
from comparison_posterior import infer_batch
from comparison_scores import interval_score, paired_summary, point_losses, score_point

HERE = Path(__file__).resolve().parent
SOURCE_COMMIT = "8782a293297330f7a354d87ebe62792018fa8866"
TRAIN_SEED = 2026091601
VALIDATION_SEED = 2026091702
FACTORS = np.geomspace(1., 20., 401)
LEVELS = (50, 80, 95)


def scenarios(selection="all"):
    """Old changing-mass path, and fixed-mass paths with changing excitation."""
    rows = []
    if selection in ("all", "original"):
        pairs = [(8., a) for a in (.25, .5, 1., 2., 3., 5., 8.)]
        pairs += [(4., a) for a in (.5, 1., 2., 3.)]
        pairs += [(a, a) for a in (.25, .5, 1., 2., 3.)]
        for f, a in pairs:
            rows.append(dict(family="original", true_force_snr=f,
                             true_acceleration_snr=a, true_mass=f/a,
                             total_signal_snr=math.hypot(f, a)))
    if selection in ("all", "fixed_mass"):
        for mass in (.25, 1., 4., 16., 32.):
            for strength in (.5, 1., 2., 4., 8., 16.):
                a = strength/math.hypot(mass, 1.)
                rows.append(dict(family="fixed_mass", true_force_snr=mass*a,
                                 true_acceleration_snr=a, true_mass=mass,
                                 total_signal_snr=strength))
    if not rows:
        raise ValueError("Unknown scenario selection")
    for i, row in enumerate(rows):
        row["id"] = f"{i:02d}_{row['family']}"
    return rows


def sample_summary(values):
    """All-trial finite-sample summary; never condition silently on validity."""
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        return {"mean": None, "empirical_mcse": None,
                "status": "nonfinite_loss_present", "nonfinite_count": int(np.sum(~np.isfinite(values)))}
    return {"mean": float(np.mean(values)),
            "empirical_mcse": float(np.std(values, ddof=1)/np.sqrt(values.size)),
            "status": "finite_sample"}


def predictive_losses(point, truth, acceleration, force_noise, acceleration_noise):
    """Independent held-out reading, at twice the training's true excitation.

    Input acceleration is either supplied exactly or measured with unit
    Gaussian noise per component. Force is newly measured with unit noise.
    These are specified prediction tasks, not additional mass accuracy scores.
    """
    valid = (point > 0) & np.isfinite(point)
    a_true = np.array([2*acceleration, 0., 0.])
    f_true = truth*a_true
    measured_force = f_true + force_noise
    losses = {}
    for label, a_new in (("known_acceleration", np.broadcast_to(a_true, measured_force.shape)),
                         ("noisy_acceleration", a_true+acceleration_noise)):
        residual = np.zeros_like(measured_force)
        residual[valid] = point[valid, None]*a_new[valid]-measured_force[valid]
        loss = np.full(point.shape, np.inf)
        loss[valid] = np.sum(residual[valid]**2, axis=1)/np.dot(f_true, f_true)
        losses[f"heldout_{label}_relative_squared_force_error"] = loss
    return losses


def concatenate_inference(chunks):
    first = chunks[0]
    return {"points": {name: np.concatenate([c["points"][name] for c in chunks])
                       for name in first["points"]},
            "distributions": {
                name: {key: np.concatenate([c["distributions"][name][key] for c in chunks])
                       for key in distribution}
                for name, distribution in first["distributions"].items()}}


def distribution_summary(distribution, truth):
    y = math.log(truth)
    scores = {key: distribution[key] for key in ("log_crps", "log_density_score")}
    intervals = {}
    for level in LEVELS:
        lower, upper = distribution[f"log_lower_{level}"], distribution[f"log_upper_{level}"]
        alpha = 1-level/100
        scores[f"interval_score_{level}"] = interval_score(lower, upper, y, alpha)
        covered = ((y >= lower) & (y <= upper)).astype(float)
        intervals[str(level)] = {
            "coverage": sample_summary(covered),
            "mean_log_width": sample_summary(upper-lower),
            "interval_score": sample_summary(scores[f"interval_score_{level}"]),
            "nominal_coverage": level/100}
    return {"scores": {key: sample_summary(v) for key, v in scores.items()},
            "intervals": intervals}, scores


def refinement_check(coarse, fine):
    result = {"maximum_relative_point_difference": 0.,
              "maximum_absolute_distribution_difference": {}}
    for name, estimate in fine["points"].items():
        old = coarse["points"][name][:len(estimate)]
        if not np.all(np.isfinite(estimate)) or not np.all(estimate > 0) or not np.all(np.isfinite(old)):
            raise FloatingPointError("Nonfinite or nonpositive posterior point in refinement")
        result["maximum_relative_point_difference"] = max(
            result["maximum_relative_point_difference"], float(np.max(np.abs(old/estimate-1))))
    for name, distribution in fine["distributions"].items():
        for key, values in distribution.items():
            if not np.all(np.isfinite(values)):
                raise FloatingPointError(f"Nonfinite posterior quantity: {name}/{key}")
            error = float(np.max(np.abs(values-coarse["distributions"][name][key][:len(values)])))
            if not math.isfinite(error):
                raise FloatingPointError(f"Nonfinite refinement difference: {name}/{key}")
            result["maximum_absolute_distribution_difference"][key] = max(
                error, result["maximum_absolute_distribution_difference"].get(key, 0.))
    limits = {key: (5e-5 if key.startswith(("log_lower_", "log_upper_")) or key == "log_median" else 2e-6)
              for key in result["maximum_absolute_distribution_difference"]}
    result["relative_point_limit"] = 2e-6
    result["absolute_distribution_limits"] = limits
    result["passed"] = (result["maximum_relative_point_difference"] <= result["relative_point_limit"]
                        and all(value <= limits[key] for key, value in result["maximum_absolute_distribution_difference"].items()))
    return result


def run(samples, selection, order, batch_size, output, save_trials=False, limit=None):
    if samples < 2 or batch_size < 1:
        raise ValueError("At least two samples and a positive batch size are required")
    errors = np.random.default_rng(TRAIN_SEED).standard_normal((samples, 6))
    heldout = np.random.default_rng(VALIDATION_SEED).standard_normal((samples, 6))
    output.mkdir(parents=True, exist_ok=True)
    if save_trials:
        (output/"trials").mkdir(exist_ok=True)
    rows = []
    start_time = time.perf_counter()
    for row in scenarios(selection)[:limit]:
        f, a, truth = row["true_force_snr"], row["true_acceleration_snr"], row["true_mass"]
        row["heldout_reference_values"] = {
            "known_acceleration_expected_loss_at_true_mass": 3/(4*f*f),
            "noisy_acceleration_expected_loss_at_true_mass": 3*(truth*truth+1)/(4*f*f),
            "noisy_acceleration_optimal_fixed_coefficient": truth*(4*a*a)/(4*a*a+3)}
        force = errors[:, :3]+[f, 0., 0.]
        acceleration = errors[:, 3:]+[a, 0., 0.]
        chunks = [infer_batch(force[i:i+batch_size], acceleration[i:i+batch_size], truth, order=order)
                  for i in range(0, samples, batch_size)]
        inferred = concatenate_inference(chunks)
        check_n = min(32, samples)
        refined = infer_batch(force[:check_n], acceleration[:check_n], truth, order=2*order)
        row["refinement"] = refinement_check(inferred, refined)
        if not row["refinement"]["passed"]:
            raise FloatingPointError(f"Numerical refinement failed in {row['id']}: {row['refinement']}")
        methods = dict(inferred["points"], **baseline_points(force, acceleration))
        row["methods"], losses = {}, {}
        arrays = {}
        for name, values in methods.items():
            summary = score_point(values, truth)
            method_losses = point_losses(values, truth)
            extra = predictive_losses(values, truth, a, heldout[:, :3], heldout[:, 3:])
            # The original joint loss is sqrt(f*a) times normalized root loss.
            extra["joint_root_loss"] = math.sqrt(f*a)*method_losses["reciprocal_root"]
            summary["additional_tasks"] = {key: sample_summary(v) for key, v in extra.items()}
            valid = (values > 0) & np.isfinite(values)
            abslog = np.full(samples, np.inf)
            abslog[valid] = np.abs(np.log(values[valid])-math.log(truth))
            summary["tolerance_curve"] = (np.searchsorted(np.sort(abslog), np.log(FACTORS), side="right")/samples).tolist()
            if name in ("positive_profile", "forward_ols", "reverse_ols", "noise_corrected_ols"):
                summary["population_unbounded_risk"] = "infinite_under_declared_failure_policy"
                summary["population_risk_note"] = (
                    "Gaussian noise gives a positive probability of invalid mass output; "
                    "under the declared infinite failure penalty all unbounded losses have infinite risk, "
                    "even when this finite sample happens to contain no failures.")
            if name == "norm_ratio":
                summary["population_risk_note"] = (
                    "In 3D the squared mass loss has finite mean but infinite variance "
                    "(inverse fourth moment of a Gaussian vector norm diverges). "
                    "Empirical MCSEs for squared mass/relative and force prediction losses "
                    "are descriptive and cannot support ordinary normal confidence intervals.")
            losses[name] = dict(method_losses, **extra)
            row["methods"][name] = summary
            arrays[f"point__{name}"] = values
        row["paired_flat_minus_other"] = {
            name: {key: paired_summary(losses["flat_joint"][key], values)
                   for key, values in method_losses.items()}
            for name, method_losses in losses.items() if name != "flat_joint"}
        for key in ("squared_mass", "squared_relative",
                    "heldout_known_acceleration_relative_squared_force_error",
                    "heldout_noisy_acceleration_relative_squared_force_error"):
            row["paired_flat_minus_other"]["norm_ratio"][key]["mcse_interpretation"] = (
                "Descriptive only: this paired loss has infinite population variance.")
        row["uncertainty"], distribution_losses = {}, {}
        for name, distribution in inferred["distributions"].items():
            row["uncertainty"][name], distribution_losses[name] = distribution_summary(distribution, truth)
            arrays.update({f"posterior__{name}__{key}": values for key, values in distribution.items()})
        row["paired_uncertainty_flat_minus_tube"] = {
            key: paired_summary(values, distribution_losses["tube_joint"][key])
            for key, values in distribution_losses["flat_joint"].items()}
        if save_trials:
            np.savez_compressed(output/"trials"/(row["id"]+".npz"), **arrays)
        rows.append(row)
        print(json.dumps({"completed": row["id"], "true_mass": truth,
                          "seconds": round(time.perf_counter()-start_time, 1),
                          "maximum_point_refinement": row["refinement"]["maximum_relative_point_difference"],
                          "maximum_interval_refinement": max(v for k, v in row["refinement"]["maximum_absolute_distribution_difference"].items()
                                                               if k.startswith(("log_lower_", "log_upper_")))}), flush=True)
    report = {
        "managed_commit": SOURCE_COMMIT,
        "python": platform.python_version(), "numpy": np.__version__,
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [Path(__file__), HERE/"comparison_scores.py",
                                    HERE/"comparison_posterior.py", HERE/"comparison_estimators.py",
                                    HERE/"crossover_vector_checks.py"]},
        "baseline_metadata": METHOD_METADATA,
        "samples_per_scenario": samples, "scenario_selection": selection,
        "training_seed": TRAIN_SEED, "heldout_seed": VALIDATION_SEED,
        "quadrature_order": order, "batch_size": batch_size,
        "tolerance_factors": FACTORS.tolist(), "elapsed_seconds": time.perf_counter()-start_time,
        "assumptions": (
            "One noisy codirectional latent 3D pair per trial; true vectors (f,0,0),(a,0,0). "
            "Independent isotropic unit component Gaussian SD in each channel; direction not supplied "
            "to any fitting method. All methods use identical training observations. "
            "True mass and true magnitudes enter only simulation and scoring. No rejection gate. "
            "The same noise draws are reused across scenarios: do not pool them as independent trials. "
            "Masses are in sigma_force/sigma_acceleration units. No global mixture-weighted winner."),
        "failure_policy": (
            "Nonpositive or nonfinite masses fail every tolerance, receive the maximum capped loss, "
            "and receive infinite loss for unbounded criteria. Null summaries never mean zero loss. "
            "Finite empirical means do not establish finite population risk."),
        "paired_difference_convention": "flat_joint minus comparator; negative means lower loss for flat_joint",
        "uncertainty_scope": (
            "Scores apply only to the flat and tube probability laws. Flat geometric/median/root "
            "points share the flat law; deterministic baselines are not assigned invented intervals. "
            "CRPS, density scores and central interval scores are evaluated on log mass."),
        "rows": rows}
    (output/"results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=16384)
    parser.add_argument("--scenario-set", choices=("all", "original", "fixed_mass"), default="all")
    parser.add_argument("--order", type=int, default=96)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--save-trials", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.samples, args.scenario_set, args.order, args.batch_size,
                 args.output, args.save_trials, args.limit)
    print(json.dumps({"output": str(args.output.resolve()/"results.json"),
                      "scenarios": len(report["rows"]),
                      "samples_per_scenario": args.samples,
                      "seconds": round(report["elapsed_seconds"], 1)}))


if __name__ == "__main__":
    main()
