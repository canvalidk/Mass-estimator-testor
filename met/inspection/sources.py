"""Optional source adapters: reuse existing generators without coupling case storage to them."""
from copy import deepcopy

import numpy as np

from .case import Case, Inputs, Observation


def from_reading_set(readings, *, prior, name, index=0, question="", truth=None, metadata=None):
    """Select one generated series; use supplied SDs and keep actual errors separate.

    Pass the ReadingSet and truth returned by met.world.generate. A prior is always
    explicit; it is never derived from truth['mass']. World/seed annotations can be
    carried in metadata. This adapter does not generate data or calibrate a check.
    """
    if type(index) is not int or not 0 <= index < len(readings.force):
        raise ValueError("index must select an existing generated series")
    f, a = readings.force[index], readings.acceleration[index]
    observations = tuple(Observation(f"r{i:03d}", fi, ai,
        np.diag(np.r_[sf ** 2, sa ** 2])) for i, (fi, ai, sf, sa) in enumerate(zip(
            f, a, readings.force_sd_axes[index], readings.acceleration_sd_axes[index])))
    investigator = None
    if truth is not None:
        investigator = dict(readings=[dict(id=r.id) for r in observations])
        if "mass" in truth:
            investigator["mass"] = float(truth["mass"])
        for channel, key, observed in (("force", "true_force", f), ("acceleration", "true_acceleration", a)):
            if key not in truth:
                continue
            latent = np.asarray(truth[key])[index]
            if latent.shape != observed.shape or not np.all(np.isfinite(latent)):
                raise ValueError(f"{key} must match the selected observed series")
            for row, true_value, error in zip(investigator["readings"], latent, observed - latent):
                row[channel] = true_value.tolist()
                row[channel + "_error"] = error.tolist()
    return Case(name, question, Inputs(observations, deepcopy(prior)), investigator, deepcopy(metadata or {}))
