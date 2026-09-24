"""Presets: a recorded study is one YAML file.

Required top-level settings (no defaults):

  name        short identifier
  question    the one question this study answers
  world       see world.py; any setting given as a list becomes a grid axis
              (noise.force_sd, noise.acceleration_sd and supplied_noise scales too)
  estimators  list of estimators (see analyst.py)
  scores      {within_factor: [K, ...]}
  replicates  series per cell
  seed        integer

Optional: `notes` (free text), `numerics` ({batch_elements}).
"""

import copy
import itertools
from pathlib import Path

import yaml

from .analyst import Estimator
from .world import validate_cell

TOP_REQUIRED = ("name", "question", "world", "estimators", "scores", "replicates", "seed")
TOP_OPTIONAL = ("notes", "numerics")
GRIDDABLE = ("dimension", "design", "mass", "acceleration_snr", "direction", "readings")
DEFAULT_RUN_NUMERICS = {"batch_elements": 3_000_000}


def load(path):
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return validate(data)


def _axes(world):
    """Grid axes: (dotted key, list of values) for every list-valued setting."""
    axes = []
    for key in GRIDDABLE:
        if isinstance(world.get(key), list):
            axes.append((key, world[key]))
    for group in ("noise", "supplied_noise"):
        if isinstance(world.get(group), dict):
            for key, value in world[group].items():
                if isinstance(value, list):
                    axes.append((f"{group}.{key}", value))
    return axes


def _set(world, dotted, value):
    if "." in dotted:
        group, key = dotted.split(".", 1)
        world[group][key] = value
    else:
        world[dotted] = value


def expand(world):
    """All cells of the grid, each with an id and its axis values."""
    axes = _axes(world)
    for key, values in axes:
        if not values:
            raise ValueError(f"grid axis world.{key} is empty")
    cells = []
    for index, combo in enumerate(itertools.product(*[values for _, values in axes])):
        cell = copy.deepcopy(world)
        label = {}
        for (key, _), value in zip(axes, combo):
            _set(cell, key, value)
            label[key] = value
        validate_cell(cell)
        cells.append({"id": f"c{index:03d}", "axes": label, "world": cell})
    return cells


def validate(data):
    if not isinstance(data, dict):
        raise ValueError("a preset must be a mapping")
    missing = [k for k in TOP_REQUIRED if k not in data]
    if missing:
        raise ValueError(f"preset is missing required setting(s): {missing}")
    extra = set(data) - set(TOP_REQUIRED) - set(TOP_OPTIONAL)
    if extra:
        raise ValueError(f"unknown preset setting(s): {sorted(extra)}")
    for key in ("name", "question"):
        if not isinstance(data[key], str) or not data[key].strip():
            raise ValueError(f"preset.{key} must be a nonempty string")
    if type(data["replicates"]) is not int or data["replicates"] < 2:
        raise ValueError("preset.replicates must be an integer >= 2")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("preset.seed must be a nonnegative integer")
    scores = data["scores"]
    if not isinstance(scores, dict) or set(scores) != {"within_factor"}:
        raise ValueError("preset.scores needs exactly within_factor: [K, ...]")
    factors = scores["within_factor"]
    if not isinstance(factors, list) or not factors or not all(
            isinstance(k, (int, float)) and k > 1 for k in factors):
        raise ValueError("scores.within_factor must be a nonempty list of numbers > 1")
    if not isinstance(data["estimators"], list) or not data["estimators"]:
        raise ValueError("preset.estimators must be a nonempty list")
    estimators = [Estimator(spec) for spec in data["estimators"]]
    names = [e.name for e in estimators]
    if len(set(names)) != len(names):
        raise ValueError("estimator names must be unique")
    cells = expand(data["world"])
    numerics = dict(DEFAULT_RUN_NUMERICS, **(data.get("numerics") or {}))
    if set(numerics) != set(DEFAULT_RUN_NUMERICS):
        raise ValueError(f"preset.numerics accepts only {sorted(DEFAULT_RUN_NUMERICS)}")
    return {"raw": data, "cells": cells, "estimators": estimators, "numerics": numerics}


def list_presets(folder):
    out = []
    for path in sorted(Path(folder).glob("*.yaml")):
        try:
            with open(path, encoding="utf-8") as handle:
                data = yaml.safe_load(handle)
            out.append((path.name, data.get("question", "")))
        except (OSError, yaml.YAMLError) as error:
            out.append((path.name, f"<unreadable: {error}>"))
    return out
