"""Small, replaceable case transformations for traversing declared problems.

A path is an ordinary iterable of cases. Neither the inspector nor its reports
need a new schema for each parameter someone wants to explore.
"""
from copy import deepcopy
import math

from .case import Case


def scale_noise(case, scale):
    """Change declared SDs by a common factor, holding measured data and prior fixed.

    Correlations and subspaces are preserved. Optional investigator truth is copied,
    not regenerated: this changes declarations, not a synthetic world's actual errors.
    """
    if isinstance(scale, bool) or not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 0:
        raise ValueError("noise scale must be finite and positive")
    data = case.to_dict(include_truth=True)
    for original, reading in zip(case.inputs.readings, data["readings"]):
        reading["covariance"] = (original.covariance * scale ** 2).tolist()
    return Case.from_dict(data)


def expand_paths(cases, values, *, transform=scale_noise, parameter="noise_scale"):
    """Apply a supplied function independently to each source and coordinate.

    The built-in noise coordinate is positive. Callers exploring other parameter
    domains can construct Case iterables directly, or supply another transform here.
    """
    values, cases = list(values), list(cases)
    if not values or not cases:
        raise ValueError("a path needs cases and coordinate values")
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in values):
        raise ValueError("path coordinates must be finite numbers")
    if len(set(values)) != len(values):
        raise ValueError("path coordinates must be distinct")
    out = []
    for case in cases:
        for index, value in enumerate(values):
            changed = transform(Case.from_dict(case.to_dict(include_truth=True)), value)
            changed.name = f"{case.name}--{parameter}-{index:03d}-{value:g}"
            changed.metadata = deepcopy(changed.metadata)
            changed.metadata["path"] = dict(source=case.name, parameter=parameter, value=value, index=index)
            out.append(changed)
    if len({c.name for c in out}) != len(out):
        raise ValueError("path case names collide; use distinguishable coordinates or source names")
    return out
