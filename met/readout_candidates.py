"""O2's explicit size candidates on one unchanged ours29 probability law.

These are pressure diagnostics, not selectable production estimators. Conditional
RMS size and the global RMS quotient are deliberately distinct constructions.
The median is a law-only control, not a choice of series size. No candidate uses
g or eta to shrink its Gaussian conditional cloud.
"""
import math

import numpy as np

from .full import _stacked_mean_length, ours29


CANDIDATE_IDS = ("cloud_mean", "mean_vector", "conditional_rms", "plain_sum_cloud",
                 "global_rms", "law_median")


def candidate_weights(readings, masses, cloud=None):
    """Gaussian conditional series sizes at supplied positive candidate masses.

    The stacked length uses each reading's orthonormal subspace coordinates, so
    squared lengths add even when the readings have their own physical frames.
    ``cloud`` may supply production's native-grid A, which contains support-cutoff
    zeros. Omitting it evaluates the expected cloud norm across the whole grid.
    The per-reading sum is an explicitly labelled alternative to T9's join.
    """
    readings = list(readings)
    if not readings:
        raise ValueError("candidate comparison requires at least one measured reading")
    masses = np.asarray(masses, dtype=float)
    if (masses.ndim != 1 or not masses.size
            or not np.all(np.isfinite(masses) & (masses > 0))):
        raise ValueError("masses must be a nonempty vector of finite positive values")
    mean_squared = np.zeros_like(masses)
    second_moment = np.zeros_like(masses)
    plain_sum = np.zeros_like(masses)
    mean_norm_cache = {}
    for reading in readings:
        pre = reading.pre()
        g = (masses[:, None, None] ** 2 * pre["k_ff"][None]
             + masses[:, None, None] * pre["k_s"][None] + pre["k_aa"][None])
        b = masses[:, None] * pre["u"][None] + pre["v"][None]
        mu = np.linalg.solve(g, b[..., None])[..., 0]
        mean2 = np.einsum("mi,mi->m", mu, mu)
        variance_trace = np.trace(np.linalg.inv(g), axis1=1, axis2=2)
        mean_squared += mean2
        second_moment += mean2 + variance_trace
        # Repeated exact-null declarations often have identical clouds. Cache
        # only this local quadrature, not a result across cases or mass grids.
        key = (reading.force.tobytes(), reading.acceleration.tobytes(),
               reading.cov.tobytes(), reading.subspace.tobytes())
        if key not in mean_norm_cache:
            mean_norm_cache[key] = _stacked_mean_length(masses, [reading])
        plain_sum += mean_norm_cache[key]
    if cloud is None:
        cloud = (_stacked_mean_length(masses, readings) if len(readings) > 1
                 else next(iter(mean_norm_cache.values())).copy())
    cloud = np.asarray(cloud, dtype=float)
    if cloud.shape != masses.shape or not np.all(np.isfinite(cloud) & (cloud >= 0)):
        raise ValueError("supplied cloud weight must match the masses and be finite and nonnegative")
    return dict(weights={"cloud_mean": cloud.copy(), "mean_vector": np.sqrt(mean_squared),
                         "conditional_rms": np.sqrt(second_moment), "plain_sum_cloud": plain_sum},
                second_moment=second_moment)


def _ratio_point(masses, probability, weight):
    denominator = float(np.sum(probability * weight))
    numerator = float(np.sum(probability * masses * weight))
    if denominator <= 0 or not (math.isfinite(denominator) and math.isfinite(numerator)):
        return dict(status="undefined", point=None)
    point = numerator / denominator
    if not math.isfinite(point) or point <= 0:
        return dict(status="undefined", point=None)
    return dict(status="defined", point=float(point))


def evaluate(readings, c):
    """Compare six explicitly defined readouts on the current law and native grid.

    The production equation is called once. Its probability, g, eta and intervals
    are returned without change. The current cloud readout preserves production's
    exact support mask and value; the other candidates do not alter that law.
    All-zero mean-vector size is undefined rather than silently patched to c.
    """
    readings = list(readings)
    if not readings:
        raise ValueError("candidate comparison requires at least one measured reading")
    if isinstance(c, (bool, np.bool_)) or not np.isscalar(c) or not math.isfinite(c) or c <= 0:
        raise ValueError("prior centre must be finite and positive")
    native = ours29(readings, float(c), readout=True)
    masses, probability = np.exp(native["lm"]), native["law"]
    sizes = candidate_weights(readings, masses, cloud=native["A"])
    candidates = {name: _ratio_point(masses, probability, weight)
                  for name, weight in sizes["weights"].items()}
    # The cloud row is the unchanged production output, including its cutoff.
    candidates["cloud_mean"] = dict(status="defined", point=float(native["readout"]))
    denominator = float(np.sum(probability * sizes["second_moment"]))
    numerator = float(np.sum(probability * masses ** 2 * sizes["second_moment"]))
    candidates["global_rms"] = (dict(status="defined", point=math.sqrt(numerator / denominator))
                                if denominator > 0 and math.isfinite(numerator) and math.isfinite(denominator)
                                else dict(status="undefined", point=None))
    candidates["law_median"] = dict(status="defined", point=float(native["median"]))
    return dict(
        law=dict(log_mass=native["lm"], probability=probability, median=native["median"],
                 interval=list(native["interval"]), raw_interval=list(native["raw_interval"]),
                 eta=native["eta"], g=native["g"], Q=native["Q"]),
        candidates=candidates, weights=sizes["weights"], second_moment=sizes["second_moment"],
        assumptions=[
            "Every row uses the same current ours29 law and Gaussian conditional clouds.",
            "Conditional RMS integrates sqrt(E[length squared|m]); global RMS takes the square root after integrating both squared channels.",
            "The cloud weight retains production support-cutoff zeros; other sampled sizes are not a numerical error bound.",
            "Undefined mean-vector weights receive no fallback value.",
            "The perpendicular prior-return case is a diagnostic, not a declaration that it is null under Z.",
        ])
