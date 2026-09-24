"""Stage 6: scores. How each estimator's readouts did against the true mass.

A readout that is NaN (improper law, off-grid, or an invalid direct rule)
counts as a failure for every score and is reported in `invalid`.
"""

import math

import numpy as np


def _fraction(hits):
    hits = np.asarray(hits, dtype=float)
    p = float(np.mean(hits))
    return {"value": p, "mcse": math.sqrt(max(p * (1 - p), 0.0) / len(hits))}


def score_point(values, truth, factors):
    values = np.asarray(values, dtype=float)
    ratio = values / truth
    finite = np.isfinite(ratio) & (ratio > 0)
    log_err = np.where(finite, np.log(np.where(finite, ratio, 1.0)), np.inf)
    cap = math.log(2.0) ** 2
    result = {"invalid": int(np.sum(~finite)),
              "within_factor": {str(k): _fraction(np.abs(log_err) < math.log(k)) for k in factors},
              "capped_squared_log_error": float(np.mean(np.minimum(log_err**2, cap)))}
    if np.any(finite):
        q = np.quantile(ratio[finite], [0.25, 0.5, 0.75])
        result["ratio_quartiles"] = [float(v) for v in q]
    else:
        result["ratio_quartiles"] = [math.nan] * 3
    return result


def score_interval(bounds, truth):
    bounds = np.asarray(bounds, dtype=float)
    lo, hi = bounds[:, 0], bounds[:, 1]
    finite = np.isfinite(lo) & np.isfinite(hi)
    with np.errstate(invalid="ignore", divide="ignore"):
        width = np.log(hi[finite] / lo[finite])
    return {"invalid": int(np.sum(~finite)),
            "coverage": _fraction(finite & (lo <= truth) & (truth <= hi)),
            "truth_below": _fraction(finite & (truth < lo)),
            "truth_above": _fraction(finite & (truth > hi)),
            "median_log_width": float(np.median(width)) if width.size else math.nan}


def add_regret(cell_scores, factors):
    """Per cell and factor K: best success among all point readouts, and each one's regret."""
    for k in map(str, factors):
        entries = [(name, s["within_factor"][k]["value"]) for name, s in cell_scores["points"].items()]
        best = max(v for _, v in entries)
        cell_scores.setdefault("best", {})[k] = {"value": best, "by": [n for n, v in entries if v == best]}
        for name, value in entries:
            cell_scores["points"][name].setdefault("regret", {})[k] = best - value


def worst_regret(cells, factors):
    """Across cells: each point readout's largest regret, per K, and where it happened."""
    names = list(cells[0]["scores"]["points"])
    out = {}
    for name in names:
        out[name] = {}
        for k in map(str, factors):
            values = [c["scores"]["points"][name]["regret"][k] for c in cells]
            worst = max(values)
            out[name][k] = {"value": worst, "cell": cells[values.index(worst)]["id"],
                            "ties": sum(v == worst for v in values)}
    return out


def agreement(values, reference):
    """How far one estimator's readout sits from a reference estimator's, series by series.

    For points: |log(value / reference)|. For intervals: the larger of the two
    endpoints' |log ratio|. Reports quantiles, and the fraction within 1% and 5%.
    """
    values, reference = np.asarray(values, dtype=float), np.asarray(reference, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        dist = np.abs(np.log(values / reference))
    if dist.ndim == 2:
        dist = np.max(dist, axis=1)
    finite = np.isfinite(dist)
    d = dist[finite]
    if d.size == 0:
        return {"invalid": int(np.sum(~finite))}
    q = np.quantile(d, [0.5, 0.9, 0.99])
    return {"invalid": int(np.sum(~finite)), "median": float(q[0]), "p90": float(q[1]), "p99": float(q[2]),
            "max": float(d.max()), "within_1pct": float(np.mean(d < math.log(1.01))),
            "within_5pct": float(np.mean(d < math.log(1.05)))}
