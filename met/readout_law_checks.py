"""O2 checks on one supplied finite candidate law, without re-fitting evidence.

The composition theorem's domain includes arbitrary positive finite-mean laws and
common-channel couplings. These atomic examples are in that domain; Gaussian
measurement-adapter reachability is not asserted. Undefined expressions return
None, rather than introducing a fallback point or a median tie convention.
"""
import math

import numpy as np


MAGNITUDE_CANDIDATES = ("cloud_mean", "conditional_rms", "global_rms", "law_median")


def _probabilities(probabilities, count):
    if probabilities is None:
        return np.full(count, 1 / count)
    p = np.asarray(probabilities, dtype=float)
    if (p.shape != (count,) or not np.all(np.isfinite(p))
            or np.any(p < 0) or not np.any(p > 0)):
        raise ValueError("probabilities must be finite nonnegative atom weights with positive total")
    p = p / p.max()
    return p / p.sum()


def _magnitude_law(x, y, probabilities):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if (x.ndim != 1 or not x.size or y.shape != x.shape
            or not np.all(np.isfinite(x) & (x > 0))
            or not np.all(np.isfinite(y) & (y > 0))):
        raise ValueError("magnitudes must be nonempty equal-length finite positive arrays")
    p = _probabilities(probabilities, len(x))
    active = p > 0
    return x[active], y[active], p[active]


def _ratio(numerator, denominator):
    if denominator == 0:
        return None
    value = float(numerator / denominator)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("readout exceeds the finite positive numerical range")
    return value


def magnitude_readout(candidate, x, y, probabilities=None):
    """Read positive magnitude atoms, preserving their original coupling.

    Conditional candidates group exactly equal floating-point M=X/Y values;
    no near-mass bins or data-dependent tolerance are introduced. Median examples
    use unique medians; this function returns None on a nonunique median interval.
    """
    x, y, p = _magnitude_law(x, y, probabilities)
    if candidate == "cloud_mean":
        return _ratio(p @ x, p @ y)
    if candidate == "global_rms":
        ratio = _ratio(p @ (x * x), p @ (y * y))
        return math.sqrt(ratio) if ratio is not None else None
    mass = x / y
    if not np.all(np.isfinite(mass) & (mass > 0)):
        raise ValueError("mass ratios exceed the finite positive numerical range")
    values, group = np.unique(mass, return_inverse=True)
    weights = np.bincount(group, weights=p)
    if candidate == "conditional_rms":
        size = np.sqrt(np.bincount(group, weights=p * y * y) / weights)
        return _ratio(np.sum(weights * values * size), weights @ size)
    if candidate == "law_median":
        cumulative = np.cumsum(weights)
        index = int(np.searchsorted(cumulative, .5, side="left"))
        if cumulative[index] == .5 and index + 1 < len(values):
            return None
        return float(values[index])
    raise ValueError(f"unknown magnitude candidate {candidate!r}")


def vector_readout(candidate, force, acceleration, probabilities=None):
    """Read compatible vector atoms; the mean-vector variant conditions on mass.

    Every atom has a positive codirectional pair, F=M a. Direction cancellation
    may nevertheless make its conditional-mean size zero. No latent direction is
    added to a magnitude-only candidate.
    """
    f, a = np.asarray(force, dtype=float), np.asarray(acceleration, dtype=float)
    if (f.ndim != 2 or not len(f) or not f.shape[1] or a.shape != f.shape
            or not np.all(np.isfinite(f)) or not np.all(np.isfinite(a))):
        raise ValueError("force and acceleration must be equal finite atom-by-dimension arrays")
    x, y = np.linalg.norm(f, axis=1), np.linalg.norm(a, axis=1)
    x, y, p = _magnitude_law(x, y, probabilities)
    original_p = _probabilities(probabilities, len(f))
    f, a = f[original_p > 0], a[original_p > 0]
    mass = x / y
    if np.any(np.linalg.norm(f - mass[:, None] * a, axis=1) > 1e-12 * x):
        raise ValueError("vector atoms must be positive codirectional pairs")
    if candidate != "mean_vector":
        return magnitude_readout(candidate, x, y, p)
    values, group = np.unique(mass, return_inverse=True)
    weights = np.bincount(group, weights=p)
    size = np.array([np.linalg.norm(np.sum(p[group == j, None] * a[group == j], axis=0) / weights[j])
                     for j in range(len(values))])
    return _ratio(np.sum(weights * values * size), weights @ size)


def _composition_record(candidate, x1, x2, y, *, vector=False, exact_left, exact_right):
    x1, x2, y = np.asarray(x1, float), np.asarray(x2, float), np.asarray(y, float)
    p = np.full(len(y), 1 / len(y))
    read = vector_readout if vector else magnitude_readout
    left = read(candidate, x1 + x2, y, p)
    right = read(candidate, x1, y, p) + read(candidate, x2, y, p)
    inverse_left = 1 / read(candidate, y, x1 + x2, p)
    inverse_right = 1 / read(candidate, y, x1, p) + 1 / read(candidate, y, x2, p)
    return dict(candidate=candidate, status="defined", probabilities=p.tolist(),
                first_force=x1.tolist(), second_force=x2.tolist(), common_acceleration=y.tolist(),
                vector_atoms=vector, L5_left=left, L5_right=right,
                L14_inverse_left=inverse_left, L14_inverse_right=inverse_right,
                L5_absolute_error=abs(left - right), L14_absolute_error=abs(inverse_left - inverse_right),
                exact_left=exact_left, exact_right=exact_right,
                domain="One propagated finite joint candidate law with its common channel; no fresh fit or reweighting.",
                scope="Full composition-law domain; reachability through the Gaussian measurement adapter is not asserted.")


def counterexamples():
    """Return JSON-ready composition identities or counterexamples for the table."""
    common = dict(x1=[1, 1, 4], x2=[1, 4, 1], y=[1, 1, 1])
    return [
        _composition_record("cloud_mean", **common, exact_left="4", exact_right="4"),
        _composition_record("law_median", **common, exact_left="5", exact_right="2"),
        _composition_record("conditional_rms", [1, 3, 2, 6], [1, 6, 1, 6], [1, 3, 1, 3],
                            exact_left="1+sqrt(5)", exact_right="13/4"),
        _composition_record("global_rms", [1, 3], [3, 1], [1, 1],
                            exact_left="4", exact_right="2*sqrt(5)"),
        _composition_record("mean_vector", [[1], [-1], [2]], [[1], [-2], [1]], [[1], [-1], [1]],
                            vector=True, exact_left="2", exact_right="10/3"),
    ]


def direction_counterexample():
    """Same original radial law, different conditional direction cancellation.

    This diagnoses dependence beyond that radial law. A proposal redefining the
    joined-size variables needs a separate argument; direct composition examples
    above do not rely on resolving that terminology.
    """
    force_a, acceleration_a, p_a = [[1], [2]], [[1], [1]], [.5, .5]
    force_b, acceleration_b, p_b = [[1], [-1], [2]], [[1], [-1], [1]], [.25, .25, .5]
    return dict(candidate="mean_vector", requirement="original magnitude-law dependence",
                first=dict(force=force_a, acceleration=acceleration_a, probabilities=p_a),
                second=dict(force=force_b, acceleration=acceleration_b, probabilities=p_b),
                first_readout=vector_readout("mean_vector", force_a, acceleration_a, p_a),
                second_readout=vector_readout("mean_vector", force_b, acceleration_b, p_b),
                common_magnitude_law=[dict(force=1, acceleration=1, probability=.5),
                                      dict(force=2, acceleration=1, probability=.5)],
                scope="Original candidate magnitudes on the full finite-law domain; not a Gaussian reachability claim.")


def cancellation_counterexample():
    force, acceleration, probabilities = [[2], [-2]], [[1], [-1]], [.5, .5]
    result = vector_readout("mean_vector", force, acceleration, probabilities)
    return dict(candidate="mean_vector", status="undefined" if result is None else "defined", value=result,
                reason="Conditional mean directions cancel even though the candidate magnitudes are surely (2,1).",
                force=force, acceleration=acceleration, probabilities=probabilities,
                cloud_mean=vector_readout("cloud_mean", force, acceleration, probabilities),
                scope="Full finite candidate-law domain, not a claim about Gaussian measurement-adapter reachability.")
