"""Validated, reproducible configurations for the shared comparison engine.

The engine supports independent isotropic or diagonal Gaussian noise in 3D.
Signal-to-noise values use each channel's RMS component standard deviation;
configured masses are physical masses, including for fixed-mass designs.
Configurations contain settings, never executable code or output paths.
"""

import json
import math
from pathlib import Path


SCHEMA_VERSION = 1
MAX_SCENARIOS = 64
MAX_TOTAL_TRIALS = 50_000_000
# Per-cell and whole-run work guards are separate. Quadrature is batched;
# 10,000 supports the fixed-replication assessment without raising total work.
MAX_DIAGONAL_SAMPLES = 10_000
MAX_DIAGONAL_TOTAL_TRIALS = 262_144
MAX_SD_RATIO = 16
MAX_CALIBRATED_SD_RATIO = 32
DEFAULTS = {
    "schema_version": SCHEMA_VERSION,
    "name": "Untitled experiment",
    "samples": 1024,
    "order": 96,
    "batch_size": 256,
    "training_seed": 2026091601,
    "heldout_seed": 2026091702,
    "save_trials": False,
    "noise_model": "isotropic_unit_3d",
}


def _fields(value, allowed, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    unknown = set(value) - set(allowed)
    if unknown:
        raise ValueError(f"Unknown {label} field(s): {', '.join(sorted(map(str, unknown)))}")


def _integer(value, low, high, label):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{label} must be an integer from {low} to {high}")
    return value


def _positive(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite positive number")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{label} must be a finite positive number") from error
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{label} must be a finite positive number")
    return result


def _positive_list(value, label):
    if not isinstance(value, list) or not 1 <= len(value) <= MAX_SCENARIOS:
        raise ValueError(f"{label} must contain 1 to {MAX_SCENARIOS} numbers")
    return [_positive(item, f"{label}[{index}]") for index, item in enumerate(value)]


def _nonnegative_list(value, label):
    if not isinstance(value, list) or not 1 <= len(value) <= MAX_SCENARIOS:
        raise ValueError(f"{label} must contain 1 to {MAX_SCENARIOS} numbers")
    result = []
    for index, item in enumerate(value):
        # Only an explicit numeric zero extends the existing positive domain.
        # In particular False, strings, and positive values that underflow in
        # later physical conversions must never create a null experiment.
        if not isinstance(item, bool) and isinstance(item, (int, float)) and item == 0:
            result.append(0.)
        else:
            result.append(_positive(item, f"{label}[{index}] (zero is also allowed)"))
    return result


def _noise_vector(value, label):
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{label} must contain exactly three positive standard deviations")
    result = [_positive(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if max(result) / min(result) > MAX_SD_RATIO:
        raise ValueError(f"{label}: largest / smallest standard deviation must be at most {MAX_SD_RATIO} for the validated diagonal engine")
    return result


def _direction(value):
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError("direction must contain exactly three finite numbers")
    result = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError("direction must contain exactly three finite numbers")
        try:
            number = float(item)
        except OverflowError as error:
            raise ValueError("direction must contain exactly three finite numbers") from error
        if not math.isfinite(number):
            raise ValueError("direction must contain exactly three finite numbers")
        result.append(number)
    scale = max(map(abs, result))
    if scale == 0:
        raise ValueError("direction must be nonzero")
    # Retain an already normalized vector so saved configurations are exactly
    # idempotent when validated again by the worker and configuration runner.
    if math.isclose(math.hypot(*result), 1., rel_tol=5e-15, abs_tol=0.):
        return result
    scaled = [number / scale for number in result]
    norm = math.hypot(*scaled)
    return [number / norm for number in scaled]


def rms_sd(values):
    """RMS component standard deviation, without squaring large inputs."""
    scale = max(values)
    return scale * (math.hypot(*(number / scale for number in values)) / math.sqrt(3.))


def uses_measurement_engine(config):
    """Whether repeated readings or uncertain calibration need diagonal fitting."""
    measurement = config.get("measurement", {})
    return (measurement.get("repeats", 1) != 1
            or measurement.get("calibration", {}).get("mode", "known") != "known")


def _normalize_measurement(value, config):
    _fields(value, ("repeats", "calibration"), "measurement")
    repeats = _integer(value.get("repeats", 1), 1, 1000, "measurement.repeats")
    calibration = value.get("calibration", {"mode": "known"})
    if not isinstance(calibration, dict):
        raise ValueError("measurement.calibration must be an object")
    mode = calibration.get("mode", "known")
    if mode == "known":
        _fields(calibration, ("mode",), "measurement.calibration")
        normalized = {"mode": mode}
    elif mode == "scaled":
        _fields(calibration, ("mode", "force_scale", "acceleration_scale"), "measurement.calibration")
        normalized = {"mode": mode}
        for channel in ("force", "acceleration"):
            key = channel + "_scale"
            factors = calibration.get(key)
            if not isinstance(factors, list) or len(factors) != 3:
                raise ValueError(f"measurement.calibration.{key} must contain exactly three positive factors")
            normalized[key] = [_positive(item, f"measurement.calibration.{key}[{index}]")
                               for index, item in enumerate(factors)]
    elif mode in ("estimated", "pooled_isotropic"):
        allowed = ("mode", "samples", "seed", "calibration_order") if mode == "pooled_isotropic" else ("mode", "samples", "seed")
        _fields(calibration, allowed, "measurement.calibration")
        seed = _integer(calibration.get("seed"), 0, 2**32 - 1, "measurement.calibration.seed")
        if seed in (config["training_seed"], config["heldout_seed"]):
            raise ValueError("The calibration seed must differ from training_seed and heldout_seed")
        normalized = {"mode": mode, "samples": _integer(calibration.get("samples"), 4, 100000,
                                                         "measurement.calibration.samples"), "seed": seed}
        if mode == "pooled_isotropic":
            normalized["calibration_order"] = _integer(calibration.get("calibration_order", 24), 8, 96,
                                                        "measurement.calibration.calibration_order")
    else:
        raise ValueError("measurement.calibration.mode must be known, scaled, estimated, or pooled_isotropic")
    # Deterministic supplied scales can be checked before creating output. For
    # estimated calibration each realized draw receives the same check at run
    # time; an unsupported draw fails the entire run, never a selected trial.
    noise = config.get("noise", {"force_sd": [1., 1., 1.], "acceleration_sd": [1., 1., 1.]})
    if mode == "pooled_isotropic":
        for channel in ("force", "acceleration"):
            values = noise[channel + "_sd"]
            if any(value != values[0] for value in values[1:]):
                raise ValueError(f"pooled_isotropic calibration requires isotropic true {channel} noise")
    supplied_rms = []
    for channel in ("force", "acceleration"):
        scales = normalized.get(channel + "_scale", [1., 1., 1.])
        supplied = [sd * factor / math.sqrt(repeats) for sd, factor in zip(noise[channel + "_sd"], scales)]
        supplied = _noise_vector(supplied, f"Supplied mean {channel} noise")
        supplied_rms.append(rms_sd(supplied))
    ratio = supplied_rms[0] / supplied_rms[1]
    if not math.isfinite(ratio) or ratio <= 0:
        raise ValueError("The supplied force / acceleration RMS noise ratio must remain finite and positive")
    return {"repeats": repeats, "calibration": normalized}


def _rows_for_config(config):
    if config["scenario"]["type"] == "mass_excitation":
        return _mass_excitation_rows(config)
    if config["noise_model"] == "isotropic_unit_3d":
        return _rows_for_scenario(config["scenario"])
    force_sd = rms_sd(config["noise"]["force_sd"])
    acceleration_sd = rms_sd(config["noise"]["acceleration_sd"])
    noise_scale = force_sd / acceleration_sd
    if not math.isfinite(noise_scale) or noise_scale <= 0:
        raise ValueError("The force / acceleration RMS noise ratio must remain finite and positive")
    scenario = config["scenario"]
    if scenario["type"] == "fixed_mass":
        scenario = {**scenario, "masses": [mass / noise_scale for mass in scenario["masses"]]}
        if not all(math.isfinite(mass) and mass > 0 for mass in scenario["masses"]):
            raise ValueError("Fixed masses divided by the RMS noise ratio must remain finite and positive")
    rows = _rows_for_scenario(scenario)
    for row in rows:
        # Preset fixed-mass cells also specify physical masses. Original preset
        # cells and explicit pairs instead specify the two RMS-based SNRs.
        if scenario["type"] == "preset" and row["family"] == "fixed_mass":
            normalized_mass = row["true_mass"] / noise_scale
            acceleration = row["total_signal_snr"] / math.hypot(normalized_mass, 1.)
            row["true_acceleration_snr"] = acceleration
            row["true_force_snr"] = normalized_mass * acceleration
        force = row["true_force_snr"] * force_sd
        acceleration = row["true_acceleration_snr"] * acceleration_sd
        if not all(math.isfinite(number) and number > 0 for number in (force, acceleration)):
            raise ValueError("Physical force and acceleration magnitudes must remain finite and positive")
        mass = force / acceleration
        if not math.isfinite(mass) or mass <= 0:
            raise ValueError("Physical mass must remain finite and positive")
        row.update(true_mass=mass, true_force_magnitude=force,
                   true_acceleration_magnitude=acceleration, direction=list(config["direction"]))
    return rows


def _mass_excitation_rows(config):
    """Keep mass independent of excitation, including the true-zero boundary."""
    scenario = config["scenario"]
    if len(scenario["masses"]) * len(scenario["acceleration_snrs"]) > MAX_SCENARIOS:
        raise ValueError(f"An experiment may contain at most {MAX_SCENARIOS} scenarios")
    diagonal = config["noise_model"] == "diagonal_gaussian_3d"
    force_sd = rms_sd(config["noise"]["force_sd"]) if diagonal else 1.
    acceleration_sd = rms_sd(config["noise"]["acceleration_sd"]) if diagonal else 1.
    noise_scale = force_sd / acceleration_sd
    if not math.isfinite(noise_scale) or noise_scale <= 0:
        raise ValueError("The force / acceleration RMS noise ratio must remain finite and positive")
    rows = []
    for mass in scenario["masses"]:
        for acceleration_snr in scenario["acceleration_snrs"]:
            acceleration = acceleration_snr * acceleration_sd
            force = mass * acceleration
            force_snr = force / force_sd
            strength = math.hypot(force_snr, acceleration_snr)
            derived = (acceleration, force, force_snr, strength)
            if (not all(math.isfinite(number) for number in derived)
                    or (acceleration_snr > 0 and not all(number > 0 for number in derived))):
                raise ValueError("Positive excitation must give finite positive physical magnitudes and signal-to-noise values")
            row = {"id": f"{len(rows):02d}_mass_excitation", "family": "mass_excitation",
                   "true_force_snr": force_snr, "true_acceleration_snr": acceleration_snr,
                   "true_mass": mass, "total_signal_snr": strength,
                   "true_force_magnitude": force, "true_acceleration_magnitude": acceleration}
            if diagonal:
                row["direction"] = list(config["direction"])
            rows.append(row)
    return rows


def _normalize_scenario(value):
    if not isinstance(value, dict):
        raise ValueError("scenario must be an object")
    kind = value.get("type", "preset")
    if kind == "preset":
        _fields(value, ("type", "selection"), "scenario")
        selection = value.get("selection", "original")
        if selection not in ("all", "original", "fixed_mass"):
            raise ValueError("scenario.selection must be all, original, or fixed_mass")
        return {"type": kind, "selection": selection}
    if kind == "force_path":
        _fields(value, ("type", "force_snr", "acceleration_snrs"), "scenario")
        return {"type": kind, "force_snr": _positive(value.get("force_snr"), "scenario.force_snr"),
                "acceleration_snrs": _positive_list(value.get("acceleration_snrs"), "scenario.acceleration_snrs")}
    if kind == "fixed_mass":
        _fields(value, ("type", "masses", "signal_snrs"), "scenario")
        return {"type": kind, "masses": _positive_list(value.get("masses"), "scenario.masses"),
                "signal_snrs": _positive_list(value.get("signal_snrs"), "scenario.signal_snrs")}
    if kind == "mass_excitation":
        _fields(value, ("type", "masses", "acceleration_snrs"), "scenario")
        return {"type": kind, "masses": _positive_list(value.get("masses"), "scenario.masses"),
                "acceleration_snrs": _nonnegative_list(value.get("acceleration_snrs"), "scenario.acceleration_snrs")}
    if kind == "pairs":
        _fields(value, ("type", "pairs"), "scenario")
        pairs = value.get("pairs")
        if not isinstance(pairs, list) or not 1 <= len(pairs) <= MAX_SCENARIOS:
            raise ValueError(f"scenario.pairs must contain 1 to {MAX_SCENARIOS} pairs")
        normalized = []
        for index, pair in enumerate(pairs):
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError(f"scenario.pairs[{index}] must be [force_snr, acceleration_snr]")
            normalized.append([_positive(item, f"scenario.pairs[{index}][{column}]")
                               for column, item in enumerate(pair)])
        return {"type": kind, "pairs": normalized}
    raise ValueError("scenario.type must be preset, force_path, fixed_mass, mass_excitation, or pairs")


def _rows_for_scenario(scenario):
    kind = scenario["type"]
    if kind == "preset":
        from compare_estimators import scenarios
        return scenarios(scenario["selection"])
    if kind == "force_path":
        pairs = [(scenario["force_snr"], acceleration) for acceleration in scenario["acceleration_snrs"]]
    elif kind == "fixed_mass":
        if len(scenario["masses"]) * len(scenario["signal_snrs"]) > MAX_SCENARIOS:
            raise ValueError(f"An experiment may contain at most {MAX_SCENARIOS} scenarios")
        pairs = []
        for mass in scenario["masses"]:
            for strength in scenario["signal_snrs"]:
                acceleration = strength / math.hypot(mass, 1.)
                pairs.append((mass * acceleration, acceleration))
    else:
        pairs = scenario["pairs"]
    rows = []
    for index, (force, acceleration) in enumerate(pairs):
        # The inputs can be finite while their derived mass overflows, or a
        # fixed-mass acceleration underflows. Reject these before any files
        # are written or random draws are generated.
        if force <= 0 or acceleration <= 0:
            raise ValueError("Scenario magnitudes must remain positive at machine precision")
        mass, strength = force / acceleration, math.hypot(force, acceleration)
        if not all(math.isfinite(number) and number > 0 for number in (force, acceleration, mass, strength)):
            raise ValueError("Scenario magnitudes, mass, and total signal must remain finite and positive")
        rows.append({"id": f"{index:02d}_{kind}", "family": kind,
                     "true_force_snr": force, "true_acceleration_snr": acceleration,
                     "true_mass": mass, "total_signal_snr": strength})
    return rows


def normalize_config(data):
    """Return a canonical JSON-safe configuration, or raise ValueError.

    Validation is performed in full before the runner creates its output
    directory. Defaults are copied; the caller's object is never modified.
    """
    _fields(data, (*DEFAULTS, "scenario", "noise", "direction", "measurement"), "configuration")
    config = {key: data.get(key, default) for key, default in DEFAULTS.items()}
    _integer(config["schema_version"], SCHEMA_VERSION, SCHEMA_VERSION, "schema_version")
    name = config["name"]
    if not isinstance(name, str) or not name.strip() or len(name) > 100:
        raise ValueError("name must be a nonempty string of at most 100 characters")
    config["name"] = name.strip()
    _integer(config["samples"], 2, 1_048_576, "samples")
    _integer(config["order"], 48, 384, "order")
    _integer(config["batch_size"], 1, 4096, "batch_size")
    for key in ("training_seed", "heldout_seed"):
        _integer(config[key], 0, 2**32 - 1, key)
    if config["training_seed"] == config["heldout_seed"]:
        raise ValueError("training_seed and heldout_seed must be different")
    if type(config["save_trials"]) is not bool:
        raise ValueError("save_trials must be true or false")
    if config["noise_model"] not in ("isotropic_unit_3d", "diagonal_gaussian_3d"):
        raise ValueError("noise_model must be isotropic_unit_3d or diagonal_gaussian_3d")
    if config["noise_model"] == "diagonal_gaussian_3d":
        noise = data.get("noise", {})
        _fields(noise, ("force_sd", "acceleration_sd"), "noise")
        config["noise"] = {key: _noise_vector(noise.get(key, [1., 1., 1.]), f"noise.{key}")
                           for key in ("force_sd", "acceleration_sd")}
        config["direction"] = _direction(data.get("direction", [1., 0., 0.]))
        _integer(config["samples"], 2, MAX_DIAGONAL_SAMPLES, "samples for diagonal noise")
    elif "noise" in data or "direction" in data:
        raise ValueError("noise and direction settings require noise_model diagonal_gaussian_3d")
    if "measurement" in data:
        config["measurement"] = _normalize_measurement(data["measurement"], config)
    diagonal_engine = config["noise_model"] == "diagonal_gaussian_3d" or uses_measurement_engine(config)
    if diagonal_engine:
        _integer(config["samples"], 2, MAX_DIAGONAL_SAMPLES, "samples for diagonal measurement engine")
    config["scenario"] = _normalize_scenario(data.get("scenario", {"type": "preset", "selection": "original"}))
    rows = _rows_for_config(config)
    if len(rows) > MAX_SCENARIOS:
        raise ValueError(f"An experiment may contain at most {MAX_SCENARIOS} scenarios")
    total_limit = MAX_DIAGONAL_TOTAL_TRIALS if diagonal_engine else MAX_TOTAL_TRIALS
    if config["samples"] * len(rows) > total_limit:
        raise ValueError(f"An experiment may contain at most {total_limit:,} total trials for this noise model")
    return config


def scenario_rows(config):
    """Expand a configuration into the comparison engine's scenario rows."""
    return _rows_for_config(normalize_config(config))


def run_config(config, output: Path, progress=None):
    """Save validated settings and execute them through the existing engine."""
    from compare_estimators import run
    normalized = normalize_config(config)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "config.json").write_text(json.dumps(normalized, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    scenario = normalized["scenario"]
    selection = scenario["selection"] if scenario["type"] == "preset" else scenario["type"]
    return run(normalized["samples"], selection, normalized["order"], normalized["batch_size"],
               output, normalized["save_trials"], scenario_rows=_rows_for_config(normalized),
               training_seed=normalized["training_seed"], heldout_seed=normalized["heldout_seed"],
               progress=progress, config=normalized)
