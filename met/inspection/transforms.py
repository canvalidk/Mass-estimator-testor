"""Small transformations reused by checks and future searches; declarations travel too."""
from copy import deepcopy
import math

import numpy as np

from .case import Inputs, Observation


def orthogonal_frames(inputs):
    readings = []
    for index, r in enumerate(inputs.readings):
        d = r.force.size
        rotation = np.eye(d)
        if d == 1:
            rotation[0, 0] = -1
        else:
            angle = 0.37 * (index + 1)
            rotation[:2, :2] = [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
        transform = np.kron(np.eye(2), rotation)
        readings.append(Observation(r.id, rotation @ r.force, rotation @ r.acceleration,
            transform @ r.covariance @ transform.T, rotation @ r.subspace, r.window, r.metadata, r.extras))
    return Inputs(tuple(readings), deepcopy(inputs.prior))


def channel_units(inputs, force_scale=3.0, acceleration_scale=0.7):
    if not all(math.isfinite(x) and x > 0 for x in (force_scale, acceleration_scale)):
        raise ValueError("unit scales must be finite and positive")
    if inputs.prior.get("family") != "half_cauchy" or set(inputs.prior) != {"family", "center"}:
        raise ValueError("the built-in unit transform supports the half_cauchy center; other prior parameters need their own transform")
    readings = []
    for r in inputs.readings:
        transform = np.diag([force_scale] * r.force.size + [acceleration_scale] * r.force.size)
        readings.append(Observation(r.id, force_scale * r.force, acceleration_scale * r.acceleration,
            transform @ r.covariance @ transform.T, r.subspace, r.window, r.metadata, r.extras))
    prior = deepcopy(inputs.prior)
    prior["center"] *= force_scale / acceleration_scale
    return Inputs(tuple(readings), prior), force_scale / acceleration_scale


def reversed_readings(inputs):
    return Inputs(tuple(reversed(inputs.readings)), deepcopy(inputs.prior))


def append_null(inputs, count=5):
    if type(count) is not int or count < 1:
        raise ValueError("null padding count must be a positive integer")
    reference = inputs.readings[-1]
    readings = list(inputs.readings)
    identifiers = {r.id for r in readings}
    for index in range(count):
        identifier = f"null_{index}"
        while identifier in identifiers:
            identifier += "_"
        identifiers.add(identifier)
        readings.append(Observation(identifier, np.zeros_like(reference.force), np.zeros_like(reference.acceleration),
            reference.covariance, reference.subspace, reference.window, {"purpose": "L20 exact-null check"}))
    return Inputs(tuple(readings), deepcopy(inputs.prior))
