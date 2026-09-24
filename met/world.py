"""The world: stage 1 of the pipeline. Generates series of readings.

The world knows the true mass, how each reading was excited, the true noise,
and which readings share a latent pair. The analyst sees none of this; it
receives only the readings and the supplied standard deviations.

A cell of a study is one fully specified world. Its settings:

  dimension         1, 2 or 3
  design            same_pair       one latent vector pair per series, measured N times
                    new_excitation  each reading draws its own latent pair (same mass)
  mass              true mass, physical units (force unit / acceleration unit)
  acceleration_snr  true |a*| / sigma_a per latent pair: a number, or
                    {uniform: [lo, hi]} drawn per latent pair. 0 means no excitation.
  direction         fixed (the first axis) or random (uniform on the sphere),
                    drawn per latent pair
  readings          N readings per series
  noise             {force_sd, acceleration_sd}: true per-coordinate SDs
  supplied_noise    exact (the analyst is told the true SDs) or
                    {force_scale, acceleration_scale}: supplied = scale x true
"""

import hashlib
import json

import numpy as np

from .law import DIMENSIONS, ReadingSet

DESIGNS = ("same_pair", "new_excitation")
DIRECTIONS = ("fixed", "random")
WORLD_KEYS = ("dimension", "design", "mass", "acceleration_snr", "direction", "readings",
              "noise", "supplied_noise")


def _positive(value, label, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a number")
    if not (value > 0 or (allow_zero and value == 0)) or not np.isfinite(value):
        raise ValueError(f"{label} must be {'nonnegative' if allow_zero else 'positive'} and finite")
    return float(value)


def validate_cell(cell):
    """Check one fully expanded world (no grid lists left)."""
    missing = set(WORLD_KEYS) - set(cell)
    if missing:
        raise ValueError(f"world is missing required setting(s): {sorted(missing)}")
    extra = set(cell) - set(WORLD_KEYS)
    if extra:
        raise ValueError(f"unknown world setting(s): {sorted(extra)}")
    if type(cell["dimension"]) is not int or cell["dimension"] not in DIMENSIONS:
        raise ValueError(f"world.dimension must be one of {DIMENSIONS}")
    if cell["design"] not in DESIGNS:
        raise ValueError(f"world.design must be one of {DESIGNS}")
    if cell["direction"] not in DIRECTIONS:
        raise ValueError(f"world.direction must be one of {DIRECTIONS}")
    _positive(cell["mass"], "world.mass")
    snr = cell["acceleration_snr"]
    if isinstance(snr, dict):
        if set(snr) != {"uniform"} or not isinstance(snr["uniform"], list) or len(snr["uniform"]) != 2:
            raise ValueError("world.acceleration_snr as a distribution must be {uniform: [lo, hi]}")
        lo, hi = (_positive(v, "world.acceleration_snr.uniform", allow_zero=True) for v in snr["uniform"])
        if hi < lo:
            raise ValueError("world.acceleration_snr.uniform needs lo <= hi")
    else:
        _positive(snr, "world.acceleration_snr", allow_zero=True)
    if type(cell["readings"]) is not int or cell["readings"] < 1:
        raise ValueError("world.readings must be a positive integer")
    if (cell["design"] == "new_excitation" and cell["direction"] == "fixed"
            and not isinstance(snr, dict) and cell["readings"] > 1):
        raise ValueError("world.design new_excitation with a fixed direction and a single acceleration_snr "
                         "gives every reading the same latent pair; say design: same_pair, or vary the "
                         "direction or the excitation")
    noise = cell["noise"]
    if not isinstance(noise, dict) or set(noise) != {"force_sd", "acceleration_sd"}:
        raise ValueError("world.noise needs exactly force_sd and acceleration_sd")
    for key in noise:
        _positive(noise[key], f"world.noise.{key}")
    supplied = cell["supplied_noise"]
    if supplied != "exact":
        if not isinstance(supplied, dict) or set(supplied) != {"force_scale", "acceleration_scale"}:
            raise ValueError("world.supplied_noise must be 'exact' or {force_scale, acceleration_scale}")
        for key in supplied:
            _positive(supplied[key], f"world.supplied_noise.{key}")


def _directions(rng, count, d, rule):
    if rule == "fixed":
        out = np.zeros((count, d))
        out[:, 0] = 1.0
        return out
    if d == 1:
        return rng.choice([-1.0, 1.0], size=(count, 1))
    v = rng.standard_normal((count, d))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _snrs(rng, count, spec):
    if isinstance(spec, dict):
        lo, hi = spec["uniform"]
        return rng.uniform(lo, hi, size=count)
    return np.full(count, float(spec))


def cell_key(cell):
    """A stable integer from the cell's settings, so its data do not depend on
    where the cell sits in the grid."""
    canonical = json.dumps(cell, sort_keys=True, default=float)
    return int(hashlib.sha256(canonical.encode()).hexdigest()[:15], 16)


def series_rngs(seed, cell, count):
    """One independent random stream per series: series i of a cell is the same
    whether the run asks for 50 series or 5000."""
    key = cell_key(cell)
    return [np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(key, i))) for i in range(count)]


def generate(cell, rngs):
    """Draw one series of readings per random stream in `rngs`.

    Returns (ReadingSet as the analyst sees it, truth dict).
    """
    d, n = cell["dimension"], cell["readings"]
    sf, sa = cell["noise"]["force_sd"], cell["noise"]["acceleration_sd"]
    mass = float(cell["mass"])
    same = cell["design"] == "same_pair"
    series = len(rngs)
    true_a = np.empty((series, n, d))
    force = np.empty((series, n, d))
    acceleration = np.empty((series, n, d))
    for i, rng in enumerate(rngs):
        pairs = 1 if same else n
        direction = _directions(rng, pairs, d, cell["direction"])
        a = (_snrs(rng, pairs, cell["acceleration_snr"]) * sa)[:, None] * direction
        true_a[i] = np.repeat(a, n, axis=0) if same else a
        force[i] = mass * true_a[i] + sf * rng.standard_normal((n, d))
        acceleration[i] = true_a[i] + sa * rng.standard_normal((n, d))
    supplied = cell["supplied_noise"]
    ksf, ksa = (1.0, 1.0) if supplied == "exact" else (supplied["force_scale"], supplied["acceleration_scale"])
    readings = ReadingSet(force, acceleration, np.full((series, n), sf * ksf), np.full((series, n), sa * ksa))
    return readings, {"mass": mass, "true_force": mass * true_a, "true_acceleration": true_a}
