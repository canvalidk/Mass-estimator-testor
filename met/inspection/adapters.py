"""Translate an estimator into inspectable curves; this is not an independent oracle.

The inspector accepts another callable with the same small signature. Prior
families and estimator-specific assumptions are handled here, not by case storage.
"""
import math

import numpy as np

from ..full import Reading, _h_terms, _stacked_mean_length, ours29, reading_state


def ours29_adapter(inputs, masses, *, variants=("unweighted", "untempered", "quoted"), include_size=True):
    """Only Inputs crosses this boundary; it contains no investigator truth.

    Native grids supply the reported law/quantiles. log_kernel and explained-energy
    curves use the requested diagnostic grid; they do not pretend that window
    contains the full probability law. size_weight is the global stacked cloud
    mean, never the sum of individual reading lengths.
    """
    prior = inputs.prior
    if prior.get("family") != "half_cauchy" or set(prior) - {"family", "center"}:
        raise ValueError("ours29 adapter currently supports only an explicit half_cauchy prior with center")
    c = prior.get("center")
    if isinstance(c, bool) or not isinstance(c, (int, float)) or not math.isfinite(c) or c <= 0:
        raise ValueError("ours29 prior center must be finite and positive")
    readings = []
    for r in inputs.readings:
        if r.force.size not in (1, 2, 3):
            raise ValueError("ours29 adapter supports dimensions 1 to 3; independent geometry supports higher dimensions")
        readings.append(Reading(r.force.copy(), r.acceleration.copy(), r.covariance.copy(),
                                r.subspace.copy(), r.window))
    masses = np.asarray(masses, dtype=float)
    lm = np.log(masses)
    phi = np.arctan(masses / c)
    prior_log = math.log(2 / math.pi) + math.log(c) + lm - np.logaddexp(2 * math.log(c), 2 * lm)
    energy = np.array([_h_terms(r.pre(), c, phi) for r in readings])
    quoted = ours29(readings, c)
    states = [reading_state(r, c) for r in readings]
    out = dict(adapter="ours29", mass=masses.copy(), prior_log_density=prior_log, laws={}, per_reading=[])
    for index, (r, state, g) in enumerate(zip(inputs.readings, states, quoted["g"])):
        out["per_reading"].append(dict(id=r.id, Q=state["Q"], g=g, nu=state["nu"],
            projective_phi=state["phi"], explained_energy=energy[index],
            weighted_log_evidence=0.5 * g * energy[index]))
    settings = {"unweighted": (False, False), "untempered": (True, False), "quoted": (True, True)}
    for name in variants:
        if name not in settings:
            raise ValueError(f"unknown ours29 diagnostic variant {name}")
        result = quoted if name == "quoted" else ours29(readings, c, misfit=settings[name][0], calibrate=settings[name][1])
        density = result["law"] / np.trapezoid(result["law"], result["lm"])
        edge_fraction = float(max(density[0], density[-1]) / density.max())
        log_kernel = prior_log + result["eta"] * np.sum(0.5 * np.asarray(result["g"])[:, None] * energy, axis=0)
        out["laws"][name] = dict(log_mass=result["lm"], density=density, point=result["readout"],
            median=result["median"], interval=list(result["interval"]), log_kernel=log_kernel,
            interval_probability=0.95, interval_type="equal-tailed",
            support_edge_density_fraction=edge_fraction, support_edge_warning=edge_fraction > 1e-7,
            log_grid_bounds=[float(result["lm"][0]), float(result["lm"][-1])],
            eta=result["eta"], g=result["g"], Q=result["Q"], curvature=result["curv"],
            J=result["J"], projective_phi=result["phi_hat"])
    if include_size:
        out["size_weight"] = _stacked_mean_length(masses, readings)
    return out
