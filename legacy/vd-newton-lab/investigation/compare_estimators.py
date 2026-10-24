"""Compare mass points and uncertainty laws against independent simulation truth.

Supports isotropic and independent diagonal-noise unknown-direction 3D experiments. It does
not change the full-covariance production estimator or discard misaligned
measurements. Requires Python and NumPy. See comparison_guide.md.
"""
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import time

import numpy as np

from comparison_estimators import METHOD_METADATA, baseline_points
from comparison_posterior import infer_batch
from comparison_scores import interval_score, paired_summary, point_losses, score_point
from trial_archive import ARCHIVE_SCHEMA_VERSION, write_trial_archive

HERE = Path(__file__).resolve().parent
SOURCE_COMMIT = "8782a293297330f7a354d87ebe62792018fa8866"
TRAIN_SEED = 2026091601
VALIDATION_SEED = 2026091702
FACTORS = np.geomspace(1., 20., 401)
LEVELS = (50, 80, 95)
ENGINE_SOURCES = ("compare_estimators.py", "comparison_scores.py", "comparison_posterior.py",
                  "comparison_estimators.py", "crossover_vector_checks.py", "experiment_config.py",
                  "comparison_diagonal_posterior.py", "comparison_diagonal_baselines.py", "calibration_posterior.py", "trial_archive.py")
PREDICTION_METRICS = tuple(f"heldout_{kind}_relative_squared_force_error"
                           for kind in ("known_acceleration", "noisy_acceleration"))


def inapplicable_summary(reason, paired=False):
    """No score is defined; this is neither zero loss nor estimator failure."""
    if paired:
        return {"mean_difference": None, "mcse": None, "status": "not_applicable", "reason": reason}
    return {"mean": None, "empirical_mcse": None, "status": "not_applicable", "reason": reason}


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


def predictive_losses(point, truth, acceleration, force_noise, acceleration_noise, direction=None):
    """Independent held-out reading, at twice the training's true excitation.

    Input acceleration is either supplied exactly or measured with unit
    Gaussian noise per component. Force is newly measured with unit noise.
    These are specified prediction tasks, not additional mass accuracy scores.
    """
    valid = (point > 0) & np.isfinite(point)
    direction = np.array([1., 0., 0.]) if direction is None else np.asarray(direction, float)
    a_true = 2*acceleration*direction
    f_true = truth*a_true
    measured_force = f_true + force_noise
    if acceleration == 0:
        return {key: None for key in PREDICTION_METRICS}
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
    paired_metrics = {}
    for level in LEVELS:
        lower, upper = distribution[f"log_lower_{level}"], distribution[f"log_upper_{level}"]
        alpha = 1-level/100
        scores[f"interval_score_{level}"] = interval_score(lower, upper, y, alpha)
        covered = ((y >= lower) & (y <= upper)).astype(float)
        paired_metrics.update({f"coverage_{level}": covered,
                               f"below_interval_{level}": (y < lower).astype(float),
                               f"above_interval_{level}": (y > upper).astype(float),
                               f"log_width_{level}": upper-lower})
        intervals[str(level)] = {
            "coverage": sample_summary(covered),
            "below_interval": sample_summary((y < lower).astype(float)),
            "above_interval": sample_summary((y > upper).astype(float)),
            "mean_log_width": sample_summary(upper-lower),
            "interval_score": sample_summary(scores[f"interval_score_{level}"]),
            "nominal_coverage": level/100}
    return {"scores": {key: sample_summary(v) for key, v in scores.items()},
            "intervals": intervals}, dict(scores, **paired_metrics)


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


def measurement_noise(config, samples, force_sd, acceleration_sd):
    """Generate a separate calibration experiment, never using training data.

    Separate coordinate variances use k-1 degrees of freedom. Pooled isotropic
    channel variances use 3(k-1), after fitting the three coordinate means.
    These are exact sufficient-statistic simulations of independent external
    Gaussian calibration experiments; repeated-reading scatter is not used.
    """
    from experiment_config import MAX_SD_RATIO, MAX_CALIBRATED_SD_RATIO
    specification = (config or {}).get("measurement", {"repeats": 1, "calibration": {"mode": "known"}})
    repeats = specification["repeats"]
    calibration = specification["calibration"]
    true_sd = np.concatenate((force_sd, acceleration_sd))
    mode = calibration["mode"]
    ratio_limit = MAX_CALIBRATED_SD_RATIO if mode == "estimated" else MAX_SD_RATIO
    if mode == "estimated":
        degrees = calibration["samples"] - 1
        variance_factors = np.random.default_rng(calibration["seed"]).chisquare(degrees, size=(samples, 6)) / degrees
        supplied = true_sd * np.sqrt(variance_factors)
    elif mode == "pooled_isotropic":
        if any(np.any(values != values[0]) for values in (force_sd, acceleration_sd)):
            raise ValueError("pooled_isotropic calibration requires isotropic true noise in each channel")
        degrees = 3 * (calibration["samples"] - 1)
        variance_factors = np.random.default_rng(calibration["seed"]).chisquare(degrees, size=(samples, 2)) / degrees
        supplied = true_sd * np.repeat(np.sqrt(variance_factors), 3, axis=1)
    elif mode == "scaled":
        supplied = true_sd * np.array(calibration["force_scale"] + calibration["acceleration_scale"])
    else:
        supplied = true_sd.copy()
    root_n = math.sqrt(repeats)
    effective = supplied / root_n
    if not np.all(np.isfinite(effective)) or not np.all(effective > 0):
        raise FloatingPointError("Calibration produced a nonfinite or nonpositive supplied mean SD; the entire run is stopped without clipping or dropping trials")
    for start, channel in ((0, "force"), (3, "acceleration")):
        values = effective[..., start:start+3]
        ratios = np.max(values, axis=-1) / np.min(values, axis=-1)
        if np.any(ratios > ratio_limit):
            maximum = float(np.max(ratios))
            raise FloatingPointError(f"Supplied {channel} calibration exceeds the validated coordinate SD ratio {ratio_limit} (realized maximum {maximum:.6g}); the entire run is stopped without clipping or dropping trials")
    rms = [np.hypot.reduce(effective[..., start:start+3], axis=-1) / math.sqrt(3) for start in (0, 3)]
    ratio = rms[0] / rms[1]
    if not np.all(np.isfinite(ratio)) or not np.all(ratio > 0):
        raise FloatingPointError("Calibration produced an unsupported force / acceleration RMS noise ratio")

    def summary(values):
        values = np.broadcast_to(values, (samples, 3))
        return {"mean": np.mean(values, axis=0).tolist(), "min": np.min(values, axis=0).tolist(),
                "max": np.max(values, axis=0).tolist()}

    metadata = {**specification,
                "protocol": "Independent repeated readings of the same latent force and acceleration vectors; exact Gaussian sufficient mean with SD divided by sqrt(repeats). No varying excitation within a trial.",
                "calibration_scope": ("Independent per-trial calibration, shared by all noise-aware estimators. Calibration draws are reused across scenario cells. Plug-in SDs treat calibration estimates as fixed; calibration uncertainty is not integrated."
                                      if mode == "estimated" else "The same supplied calibration is shared by all noise-aware estimators; it need not equal the data-generating noise."),
                "true_force_sd": force_sd.tolist(), "true_acceleration_sd": acceleration_sd.tolist(),
                "true_mean_force_sd": (force_sd/root_n).tolist(), "true_mean_acceleration_sd": (acceleration_sd/root_n).tolist(),
                "supplied_force_sd_summary": summary(supplied[..., :3]),
                "supplied_acceleration_sd_summary": summary(supplied[..., 3:]),
                "supplied_mean_force_sd_summary": summary(effective[..., :3]),
                "supplied_mean_acceleration_sd_summary": summary(effective[..., 3:])}
    metadata["maximum_supported_supplied_coordinate_sd_ratio"] = ratio_limit
    if mode == "pooled_isotropic":
        metadata.update(
            calibration_degrees=degrees,
            calibration_statistic="Independent pooled channel variances: sum of squared residuals over k external Gaussian 3-vectors about their fitted vector mean, divided by 3(k-1).",
            inference_evidence="Mean of n stable observations and external calibration sufficient statistics only; no repeated-reading residual scatter is supplied to inference.",
            calibration_scope="Fresh independent calibration per replicate and channel, shared by plug-in and integrated estimators and reused across scenario cells. Plug-in methods hold the pooled SD fixed; calibrated methods integrate calibration uncertainty. Named oracle methods alone receive the true channel SDs.")
    return effective[..., :3], effective[..., 3:], metadata


def calibration_refinement_indices(force, acceleration, force_sd, acceleration_sd, *, include_scale_extremes=False):
    """Include calibration extremes alongside the ordinary first 32 trials."""
    samples = len(force)
    indices = set(range(min(32, samples)))
    rms = []
    for observation, supplied in ((force, force_sd), (acceleration, acceleration_sd)):
        supplied = np.broadcast_to(supplied, observation.shape)
        contrast = np.max(supplied, axis=1) / np.min(supplied, axis=1)
        standardized_strength = np.linalg.norm(observation / supplied, axis=1)
        indices.add(int(np.argmax(contrast)))
        indices.add(int(np.argmax(standardized_strength)))
        rms.append(np.hypot.reduce(supplied, axis=1) / math.sqrt(3))
        if include_scale_extremes:
            indices.add(int(np.argmin(rms[-1])))
            indices.add(int(np.argmax(rms[-1])))
    mass_scale = rms[0]/rms[1]
    indices.add(int(np.argmin(mass_scale)))
    indices.add(int(np.argmax(mass_scale)))
    return np.array(sorted(indices), dtype=int)


def inference_subset(inferred, indices):
    return {"points": {key: values[indices] for key, values in inferred["points"].items()},
            "distributions": {name: {key: values[indices] for key, values in distribution.items()}
                              for name, distribution in inferred["distributions"].items()}}


def pooled_calibration_inference(force, acceleration, truth, supplied_sf, supplied_sa,
                                 true_sf, true_sa, degrees, order, calibration_order, batch_size, indices):
    """Fit shared calibration evidence and a separately named true-noise oracle."""
    from calibration_posterior import infer_batch as calibrated_infer
    from comparison_diagonal_posterior import infer_batch as diagonal_infer

    def candidate(f, a, sf, sa, angular_order, scale_order):
        return calibrated_infer(f, a, truth, force_sd=sf, acceleration_sd=sa, degrees=degrees,
                                order=angular_order, calibration_order=scale_order)

    def oracle(f, a, angular_order):
        fitted = diagonal_infer(f, a, truth, force_sd=true_sf, acceleration_sd=true_sa, order=angular_order)
        return {"points": {"oracle_" + name.removeprefix("flat_"): fitted["points"][name]
                           for name in ("flat_joint", "flat_geometric", "flat_median")},
                "distributions": {"oracle_joint": fitted["distributions"]["flat_joint"]}}

    calibrated = concatenate_inference([
        candidate(force[i:i+batch_size], acceleration[i:i+batch_size], supplied_sf[i:i+batch_size, 0],
                  supplied_sa[i:i+batch_size, 0], order, calibration_order)
        for i in range(0, len(force), batch_size)])
    oracle_fit = concatenate_inference([oracle(force[i:i+batch_size], acceleration[i:i+batch_size], order)
                                        for i in range(0, len(force), batch_size)])
    coarse = inference_subset(calibrated, indices)
    refinement = {}
    for name, angular_order, scale_order in (("integrated_angular", 2*order, calibration_order),
                                            ("integrated_calibration", order, 2*calibration_order)):
        fine = candidate(force[indices], acceleration[indices], supplied_sf[indices, 0],
                         supplied_sa[indices, 0], angular_order, scale_order)
        refinement[name] = refinement_check(coarse, fine)
    refinement["oracle_angular"] = refinement_check(inference_subset(oracle_fit, indices),
                                                    oracle(force[indices], acceleration[indices], 2*order))
    calibrated["points"].update(oracle_fit["points"])
    calibrated["distributions"].update(oracle_fit["distributions"])
    return calibrated, refinement


def combine_refinement(components, indices):
    """Retain both integration checks while exposing existing progress fields."""
    first = next(iter(components.values()))
    result = {"components": components,
              "maximum_relative_point_difference": max(value["maximum_relative_point_difference"] for value in components.values()),
              "maximum_absolute_distribution_difference": {},
              "relative_point_limit": first["relative_point_limit"],
              "absolute_distribution_limits": first["absolute_distribution_limits"],
              "passed": all(value["passed"] for value in components.values()),
              "checked_trial_indices": indices.tolist(),
              "selection": "First 32 trials plus each channel's supplied SD extrema, coordinate contrast, standardized strength and RMS mass-scale extremes"}
    for value in components.values():
        for key, error in value["maximum_absolute_distribution_difference"].items():
            result["maximum_absolute_distribution_difference"][key] = max(error, result["maximum_absolute_distribution_difference"].get(key, 0.))
    return result


def run(samples, selection, order, batch_size, output, save_trials=False, limit=None,
        *, scenario_rows=None, training_seed=TRAIN_SEED, heldout_seed=VALIDATION_SEED,
        progress=None, config=None):
    """Run a comparison, retaining the original positional command interface.

    ``experiment_config.run_config`` is the validated entry point for saved
    configurations and the local interface. A progress callback receives one
    completed-scenario event at a time; stdout retains the same JSON stream.
    """
    if samples < 2 or batch_size < 1:
        raise ValueError("At least two samples and a positive batch size are required")
    source_bytes = {name: (HERE/name).read_bytes() for name in ENGINE_SOURCES}
    source_hashes = {name: hashlib.sha256(content).hexdigest() for name, content in source_bytes.items()}
    from experiment_config import uses_measurement_engine
    true_diagonal = config is not None and config.get("noise_model") == "diagonal_gaussian_3d"
    diagonal = true_diagonal or uses_measurement_engine(config or {})
    sf = np.array(config["noise"]["force_sd"] if true_diagonal else [1., 1., 1.])
    sa = np.array(config["noise"]["acceleration_sd"] if true_diagonal else [1., 1., 1.])
    direction = np.array(config["direction"] if true_diagonal else [1., 0., 0.])
    rms_f, rms_a = math.hypot(*sf)/math.sqrt(3), math.hypot(*sa)/math.sqrt(3)
    supplied_sf, supplied_sa, measurement = measurement_noise(config, samples, sf, sa)
    repeats = measurement["repeats"]
    root_n = math.sqrt(repeats)
    pooled_calibration = measurement["calibration"]["mode"] == "pooled_isotropic"
    method_metadata = METHOD_METADATA
    infer = infer_batch
    baseline = baseline_points
    if diagonal:
        from comparison_diagonal_posterior import infer_batch as diagonal_infer
        from comparison_diagonal_baselines import baseline_points as diagonal_baseline, METHOD_METADATA as DIAGONAL_METADATA
        def infer(force, acceleration, truth, order, start=0):
            sf_batch = supplied_sf[start:start+len(force)] if supplied_sf.ndim == 2 else supplied_sf
            sa_batch = supplied_sa[start:start+len(force)] if supplied_sa.ndim == 2 else supplied_sa
            return diagonal_infer(force, acceleration, truth, force_sd=sf_batch, acceleration_sd=sa_batch, order=order)
        def baseline(force, acceleration):
            return diagonal_baseline(force, acceleration, supplied_sf, supplied_sa)
        method_metadata = deepcopy(DIAGONAL_METADATA)
        flat_assumption = ("Flat physical df d(alpha) dOmega measure; common unknown direction is integrated out. "
                           "Uses the supplied independent coordinate noise SDs in both channels.")
        method_metadata.update({
            "flat_joint": {"label": "Our ratio of means · coordinate noise", "assumptions": flat_assumption},
            "flat_geometric": {"label": "Flat-law geometric", "assumptions": flat_assumption},
            "flat_median": {"label": "Flat-law median", "assumptions": flat_assumption},
            "flat_reciprocal_root": {"label": "Flat-law root-loss point", "assumptions": flat_assumption},
            "tube_joint": {"label": "Tube ratio of means · coordinate noise",
                           "assumptions": "Supplied diagonal noise; reference measure weighted by (f/RMS_force_SD)^2 + (alpha/RMS_acceleration_SD)^2. Direction remains unknown."},
            "flat_joint_rms": {"label": "Our equation · RMS noise approximation",
                               "formula": "E_RMS[f] / E_RMS[alpha]",
                               "assumptions": "Control: replaces the supplied coordinate variances with each channel's mean variance before applying the old isotropic equation. This deliberately discards anisotropy."}})
        if measurement["calibration"]["mode"] != "known":
            for metadata in method_metadata.values():
                if "assumptions" in metadata:
                    metadata["assumptions"] = re.sub(r"\bknown\b", "supplied", metadata["assumptions"], flags=re.IGNORECASE)
                metadata["noise_information"] = "Noise-aware formulas receive the supplied calibration SDs divided by sqrt(repeats), never the true data-generating SDs. Estimated calibration is used as a plug-in input."
    if pooled_calibration:
        for prefix, label, assumptions in (
                ("calibrated", "Calibration integrated", "Flat physical df d(alpha) dOmega measure; independent channel variance priors d(sigma^2)/sigma^2. Integrates the likelihood of the supplied stable-reading mean and independent pooled calibration with 3(k-1) degrees of freedom. No repeated-reading scatter is used."),
                ("oracle", "Known-noise oracle", "Flat physical df d(alpha) dOmega measure with the true generating channel SDs divided by sqrt(repeats); an explicitly privileged oracle, using the same mean observations.")):
            for suffix, functional in (("joint", "ratio of means"), ("geometric", "geometric point"), ("median", "median")):
                method_metadata[prefix + "_" + suffix] = {"label": label + " · " + functional,
                                                          "assumptions": assumptions,
                                                          "noise_information": "True channel SDs" if prefix == "oracle" else "External pooled calibration sufficient statistics only"}
    selected_rows = [dict(row) for row in (scenarios(selection) if scenario_rows is None else scenario_rows)[:limit]]
    training_rng, heldout_rng = np.random.default_rng(training_seed), np.random.default_rng(heldout_seed)
    errors = training_rng.standard_normal((samples, 6))
    heldout = heldout_rng.standard_normal((samples, 6))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if save_trials:
        (output/"trials").mkdir(exist_ok=True)
        (output/"source").mkdir(exist_ok=True)
        for name, content in source_bytes.items():
            (output/"source"/name).write_bytes(content)
    archive_entries = []
    rows = []
    start_time = time.perf_counter()
    for row in selected_rows:
        f, a, truth = row["true_force_snr"], row["true_acceleration_snr"], row["true_mass"]
        physical_f, physical_a = row.get("true_force_magnitude", f), row.get("true_acceleration_magnitude", a)
        if config is not None and "measurement" in config:
            row.update(repeats=repeats, effective_force_snr=f*root_n,
                       effective_acceleration_snr=a*root_n,
                       effective_total_signal_snr=row["total_signal_snr"]*root_n)
        if pooled_calibration:
            row["calibration_evidence"] = {
                "mode": "pooled_isotropic", "samples": measurement["calibration"]["samples"],
                "degrees": measurement["calibration_degrees"], "seed": measurement["calibration"]["seed"],
                "input": measurement["inference_evidence"], "variance_prior": "d(sigma^2)/sigma^2 per channel"}
        row["metric_applicability"] = {}
        if physical_a == 0:
            row["metric_applicability"] = {
                key: {"status": "not_applicable", "reason": "Relative force error is undefined at zero true force."}
                for key in PREDICTION_METRICS}
            row["metric_applicability"]["joint_root_loss"] = {
                "status": "not_applicable", "reason": "The joint magnitude loss is identically zero at zero excitation and cannot assess mass accuracy."}
        row["heldout_reference_values"] = {
            "known_acceleration_expected_loss_at_true_mass": None if physical_a == 0 else float(np.sum(sf*sf))/(4*physical_f*physical_f),
            "noisy_acceleration_expected_loss_at_true_mass": None if physical_a == 0 else float(np.sum(sf*sf)+truth*truth*np.sum(sa*sa))/(4*physical_f*physical_f),
            "noisy_acceleration_optimal_fixed_coefficient": None if physical_a == 0 else truth*(4*physical_a*physical_a)/(4*physical_a*physical_a+float(np.sum(sa*sa)))}
        force = errors[:, :3]*(sf/root_n)+physical_f*direction
        acceleration = errors[:, 3:]*(sa/root_n)+physical_a*direction
        chunks = [infer(force[i:i+batch_size], acceleration[i:i+batch_size], truth, order=order,
                        **({"start": i} if diagonal else {}))
                  for i in range(0, samples, batch_size)]
        inferred = concatenate_inference(chunks)
        if measurement["calibration"]["mode"] in ("estimated", "pooled_isotropic"):
            indices = calibration_refinement_indices(force, acceleration, supplied_sf, supplied_sa,
                                                       include_scale_extremes=pooled_calibration)
            refined = diagonal_infer(force[indices], acceleration[indices], truth,
                                     force_sd=supplied_sf[indices], acceleration_sd=supplied_sa[indices], order=2*order)
            row["refinement"] = refinement_check(inference_subset(inferred, indices), refined)
            row["refinement"]["checked_trial_indices"] = indices.tolist()
            row["refinement"]["selection"] = "First 32 trials plus calibration contrast, standardized strength and RMS mass-scale extremes"
        else:
            check_n = min(32, samples)
            refined = infer(force[:check_n], acceleration[:check_n], truth, order=2*order)
            row["refinement"] = refinement_check(inferred, refined)
        if pooled_calibration:
            extra_inference, checks = pooled_calibration_inference(
                force, acceleration, truth, supplied_sf, supplied_sa, sf/root_n, sa/root_n,
                measurement["calibration_degrees"], order, measurement["calibration"]["calibration_order"], batch_size, indices)
            checks = {"plugin_angular": row["refinement"], **checks}
            row["refinement"] = combine_refinement(checks, indices)
            for key in ("points", "distributions"):
                inferred[key].update(extra_inference[key])
        if not row["refinement"]["passed"]:
            raise FloatingPointError(f"Numerical refinement failed in {row['id']}: {row['refinement']}")
        methods = dict(inferred["points"], **baseline(force, acceleration))
        if diagonal:
            supplied_rms_f = (np.array([math.hypot(*value)/math.sqrt(3) for value in supplied_sf])
                              if supplied_sf.ndim == 2 else math.hypot(*supplied_sf)/math.sqrt(3))
            supplied_rms_a = (np.array([math.hypot(*value)/math.sqrt(3) for value in supplied_sa])
                              if supplied_sa.ndim == 2 else math.hypot(*supplied_sa)/math.sqrt(3))
            mass_scale = supplied_rms_f/supplied_rms_a
            f_divisor = supplied_rms_f[:, None] if supplied_sf.ndim == 2 else supplied_rms_f
            a_divisor = supplied_rms_a[:, None] if supplied_sa.ndim == 2 else supplied_rms_a
            # Keep this control within the same memory bound as the primary
            # inference. Estimated calibration gives a distinct mass scale for
            # every replicate, so slice that scale with the observations.
            standardized_f, standardized_a = force/f_divisor, acceleration/a_divisor
            approximation_points = []
            for i in range(0, samples, batch_size):
                scale = mass_scale[i:i+batch_size] if np.ndim(mass_scale) else mass_scale
                approximation = infer_batch(standardized_f[i:i+batch_size],
                                            standardized_a[i:i+batch_size],
                                            truth/scale, order=order)
                approximation_points.append(approximation["points"]["flat_joint"]*scale)
            methods["flat_joint_rms"] = np.concatenate(approximation_points)
        row["methods"], losses = {}, {}
        arrays = {"observation__force": force, "observation__acceleration": acceleration,
                  "heldout__true_force": np.broadcast_to(2*physical_f*direction, (samples, 3)),
                  "heldout__true_acceleration": np.broadcast_to(2*physical_a*direction, (samples, 3)),
                  "heldout__force": 2*physical_f*direction+heldout[:, :3]*sf,
                  "heldout__acceleration": 2*physical_a*direction+heldout[:, 3:]*sa}
        arrays.update(noise__true_force_sd=np.broadcast_to(sf, (samples, 3)),
                          noise__true_acceleration_sd=np.broadcast_to(sa, (samples, 3)),
                          noise__supplied_force_sd=np.broadcast_to(supplied_sf*root_n, (samples, 3)),
                          noise__supplied_acceleration_sd=np.broadcast_to(supplied_sa*root_n, (samples, 3)),
                          noise__supplied_mean_force_sd=np.broadcast_to(supplied_sf, (samples, 3)),
                          noise__supplied_mean_acceleration_sd=np.broadcast_to(supplied_sa, (samples, 3)))
        for name, values in methods.items():
            summary = score_point(values, truth)
            method_losses = point_losses(values, truth)
            extra = predictive_losses(values, truth, physical_a, heldout[:, :3]*sf, heldout[:, 3:]*sa, direction)
            # The original joint loss is sqrt(f*a) times normalized root loss.
            extra["joint_root_loss"] = (None if physical_a == 0 else
                                        math.sqrt(physical_f*physical_a)*method_losses["reciprocal_root"])
            summary["additional_tasks"] = {
                key: inapplicable_summary(row["metric_applicability"][key]["reason"]) if v is None else sample_summary(v)
                for key, v in extra.items()}
            valid = (values > 0) & np.isfinite(values)
            abslog = np.full(samples, np.inf)
            abslog[valid] = np.abs(np.log(values[valid])-math.log(truth))
            summary["tolerance_curve"] = (np.searchsorted(np.sort(abslog), np.log(FACTORS), side="right")/samples).tolist()
            if name in ("positive_profile", "forward_ols", "reverse_ols", "noise_corrected_ols",
                        "covariance_profile", "weighted_forward_ols", "weighted_reverse_ols", "weighted_corrected_ols"):
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
            name: {key: (inapplicable_summary(row["metric_applicability"][key]["reason"], paired=True)
                         if values is None else paired_summary(losses["flat_joint"][key], values))
                   for key, values in method_losses.items()}
            for name, method_losses in losses.items() if name != "flat_joint"}
        for key in ("squared_mass", "squared_relative",
                    "heldout_known_acceleration_relative_squared_force_error",
                    "heldout_noisy_acceleration_relative_squared_force_error"):
            if key not in row["metric_applicability"]:
                row["paired_flat_minus_other"]["norm_ratio"][key]["mcse_interpretation"] = (
                    "Descriptive only: this paired loss has infinite population variance.")
        row["uncertainty"], distribution_losses = {}, {}
        for name, distribution in inferred["distributions"].items():
            row["uncertainty"][name], distribution_losses[name] = distribution_summary(distribution, truth)
            arrays.update({f"posterior__{name}__{key}": values for key, values in distribution.items()})
        row["paired_uncertainty_flat_minus_tube"] = {
            key: paired_summary(values, distribution_losses["tube_joint"][key])
            for key, values in distribution_losses["flat_joint"].items()
            if key in row["uncertainty"]["flat_joint"]["scores"]}
        row["paired_uncertainty_flat_minus_other"] = {
            name: {key: paired_summary(values, comparator[key])
                   for key, values in distribution_losses["flat_joint"].items()}
            for name, comparator in distribution_losses.items() if name != "flat_joint"}
        if save_trials:
            archive_entries.append(write_trial_archive(output, row, arrays))
        rows.append(row)
        event = {"completed": row["id"], "completed_count": len(rows), "total": len(selected_rows),
                 "true_mass": truth, "seconds": round(time.perf_counter()-start_time, 1),
                 "maximum_point_refinement": row["refinement"]["maximum_relative_point_difference"],
                 "maximum_interval_refinement": max(v for k, v in row["refinement"]["maximum_absolute_distribution_difference"].items()
                                                    if k.startswith(("log_lower_", "log_upper_")))}
        print(json.dumps(event), flush=True)
        if progress is not None:
            progress(dict(event))
    report = {
        "managed_commit": SOURCE_COMMIT,
        "python": platform.python_version(), "numpy": np.__version__,
        "source_sha256": source_hashes,
        "random_generators": {"training": type(training_rng.bit_generator).__name__,
                              "heldout": type(heldout_rng.bit_generator).__name__,
                              "replicate_identity": "Zero-based row index; the same index shares noise across scenarios."},
        "baseline_metadata": method_metadata,
        "samples_per_scenario": samples, "scenario_selection": selection,
        "training_seed": training_seed, "heldout_seed": heldout_seed,
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
    if true_diagonal:
        report.update(noise_model="diagonal_gaussian_3d", noise=config["noise"], direction=config["direction"],
                      rms_force_sd=rms_f, rms_acceleration_sd=rms_a,
                      mass_units="Physical force units divided by physical acceleration units",
                      snr_convention="Latent channel magnitude divided by that channel's RMS component SD; directional information is not fixed by this number alone.",
                      assumptions=("One noisy codirectional latent 3D pair per trial. True vectors have the configured common direction; "
                                   "coordinate errors are independent Gaussian with the supplied known SDs. Direction is supplied only "
                                   "to data generation, never to an estimator. All methods use the same observations; covariance-aware "
                                   "methods use the full diagonal noise information, while named RMS approximations deliberately discard it. "
                                   "Truth enters only simulation and scoring. No rejection gate. Noise draws are shared across scenarios; "
                                   "cells must not be pooled as independent trials. Mass estimates are in physical F/a units. "
                                   "No global mixture-weighted winner."))
    if config is not None and "measurement" in config:
        report["measurement"] = measurement
        report["noise_model"] = config["noise_model"]
        report["inference_noise_model"] = "diagonal_gaussian_3d" if diagonal else "isotropic_unit_3d"
        report["snr_convention"] = "Configured SNR is per reading, using true channel RMS component SD. Effective mean SNR is sqrt(repeats) times per-reading SNR; supplied calibration does not redefine truth."
        report["assumptions"] = (
            measurement["protocol"] + " " + measurement["calibration_scope"]
            + " True component errors are independent Gaussian. Latent direction enters only data generation, never fitting. "
            "All methods use the same mean observations. All noise-aware methods receive the same supplied mean SDs. "
            "Fitting receives only the declared supplied scales; true noise SDs are kept separately for data generation and scoring. "
            "The held-out task remains one independent reading at twice the latent training excitation, with original true per-reading noise. "
            "Training and calibration draws are reused across scenario cells: do not pool cells as independent trials. "
            "No rejection gate or global mixture-weighted winner.")
    if pooled_calibration:
        report["calibration_quadrature_order"] = measurement["calibration"]["calibration_order"]
        report["uncertainty_scope"] = (
            "Separate log-mass probability laws: plug-in flat and tube, calibration-integrated flat, and known-SD oracle flat. "
            "Each law's joint/geometric/median points share its intervals; deterministic baselines have no invented intervals. "
            "Scores and coverage remain fixed-truth repeated-experiment assessments, not prior-predictive calibration.")
        report["assumptions"] = (
            measurement["protocol"] + " " + measurement["inference_evidence"] + " " + measurement["calibration_scope"]
            + " True component errors are independent Gaussian and isotropic within each channel. "
            "All methods use the same mean observations; plug-in and calibrated methods share external pooled calibration. "
            "Only named oracle methods receive true noise SDs; true mass and latent direction never enter fitting. "
            "The held-out task remains one independent reading at twice the latent training excitation with original per-reading noise. "
            "Training and calibration draws are reused across scenario cells; do not pool them as independent trials. "
            "No rejection gate or global mixture-weighted winner.")
    if config is not None:
        report["experiment_config"] = config
    if any(hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest for name, digest in source_hashes.items()):
        raise RuntimeError("Engine source changed during calculation; results were not published")
    if save_trials:
        manifest = {"schema_version": ARCHIVE_SCHEMA_VERSION, "entries": archive_entries,
                    "provenance": {key: report[key] for key in ("source_sha256", "python", "numpy", "managed_commit", "random_generators", "training_seed", "heldout_seed")}}
        (output/"trial_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        report["trial_archive"] = {"schema_version": ARCHIVE_SCHEMA_VERSION, "manifest": "trial_manifest.json",
                                   "scope": "Observations, held-out readings, supplied noise, points and selected posterior summaries; full posterior curves are not saved."}
    (output/"results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--scenario-set", choices=("all", "original", "fixed_mass"))
    parser.add_argument("--order", type=int)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--save-trials", action="store_true", default=None)
    parser.add_argument("--config", type=Path, help="Saved JSON experiment configuration; replaces other experiment flags")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.config is not None:
        if any(value is not None for value in (args.samples, args.scenario_set, args.order,
                                               args.batch_size, args.limit, args.save_trials)):
            parser.error("--config cannot be combined with other experiment settings; edit the saved configuration instead")
        from experiment_config import run_config
        try:
            config = json.loads(args.config.read_text(encoding="utf-8-sig"))
            report = run_config(config, args.output)
        except (OSError, ValueError) as error:
            parser.error(str(error))
    else:
        report = run(16384 if args.samples is None else args.samples,
                     args.scenario_set or "all", 96 if args.order is None else args.order,
                     256 if args.batch_size is None else args.batch_size,
                     args.output, bool(args.save_trials), args.limit)
    print(json.dumps({"output": str(args.output.resolve()/"results.json"),
                      "scenarios": len(report["rows"]),
                      "samples_per_scenario": report["samples_per_scenario"],
                      "seconds": round(report["elapsed_seconds"], 1)}))


if __name__ == "__main__":
    main()
