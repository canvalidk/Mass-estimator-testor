"""Replay B18's weak direct measurements with the propagated half-Cauchy prior.

Uses the three original presets and their full seed streams. Only the direct
laws are fitted, in one process. The parts' priors are propagated independently
as in B18. For ours33, replace the raw law's prior before its own CDF warp.
"""

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from met.composition_coverage import _series, jobs
from met.composition_laws import LogMassLaw, analyse
from met.composition_probe import validate
from met.full import _log_prior


PROBABILITIES = (.025, .5, .975)
OPERATIONS = ("L5", "L14")


def lumped_log_density(lm, c, operation):
    """Normalized log density with respect to d(log mass), for equal centres."""
    y = np.asarray(lm, dtype=float) - math.log(c)
    if operation == "L14":
        y = -y
    elif operation != "L5":
        raise ValueError("operation must be L5 or L14")
    # At |y| > 350 the leading tail expressions are accurate to float precision.
    central = np.clip(y, -350., 350.)
    z = np.exp(central)
    numerator = np.logaddexp(np.log(np.log1p(z*z)), central + np.log(np.arctan(z)))
    value = math.log(8 / math.pi**2) + numerator - np.logaddexp(2*central, math.log(4))
    return np.where(y < -350., math.log(4 / math.pi**2) + 2*y,
                    np.where(y > 350., math.log(4 / math.pi) - y, value))


def _sum_cdf(z):
    # In the angle coordinate the first half-Cauchy is uniform. Conditioning
    # on it leaves the second half-Cauchy's CDF at the available mass.
    return (2 / math.pi)**2 * quad(lambda angle: math.atan(z - math.tan(angle)),
                                  0., math.atan(z), epsabs=1e-12, epsrel=1e-12)[0]


def lumped_quantiles(c, operation, probabilities=PROBABILITIES):
    """Quantiles by one-dimensional quadrature and reciprocal transformation."""
    if operation not in OPERATIONS:
        raise ValueError("operation must be L5 or L14")
    quantiles = []
    for probability in probabilities:
        if not 0 < probability < 1:
            raise ValueError("probabilities must be strictly between zero and one")
        rank = probability if operation == "L5" else 1 - probability
        upper = 1.
        while _sum_cdf(upper) < rank:
            upper *= 2
        quantile = brentq(lambda x: _sum_cdf(x) - rank, 0., upper, xtol=1e-13)
        quantiles.append(c * quantile if operation == "L5" else c / quantile)
    return np.asarray(quantiles)


def reweight(own, c, operation, *, refine=False):
    """Replace the prior factor of the direct law and retain its fitted eta."""
    nodes = own.lm
    prior = _log_prior(nodes, c)
    if refine:
        # Interpolate the log data factor, then evaluate the new prior exactly.
        midpoints = (nodes[:-1] + nodes[1:]) / 2
        nodes = np.sort(np.r_[nodes, midpoints])
        density = np.interp(nodes, own.lm, own.log_density - prior)
    else:
        density = own.log_density - prior
    replacement = lumped_log_density(nodes, c, operation)
    change = replacement[[0, -1]] - prior[[0, -1]]
    # Match integrated endpoint tail masses. The sum has left log-mass
    # rate 2 and right rate 1, reversed for the parallel sum. LogMassLaw
    # stores unit-rate tails, so divide the rate-2 tail mass in half.
    return LogMassLaw(
        nodes, density + replacement, eta=own.eta, warp=own.warp,
        left_log_constant=own.left_log_constant + change[0] - (operation == "L5")*math.log(2),
        right_log_constant=own.right_log_constant + change[1] - (operation == "L14")*math.log(2),
    )


def _direct_readings(seed, cell):
    rng = np.random.default_rng(seed)
    size = cell["size"]
    _series(rng, cell["m1"], size)
    _series(rng, cell["m2"], size)
    if cell["truth_case"] == "historical" and cell["composition"] == "L14":
        _series(rng, 3., size)
    return _series(rng, cell["truth"], size)


def replay(spec):
    raw = spec["raw"]
    cells = [cell for cell in spec["cells"] if cell["size"] == .1]
    seed_stream = jobs(raw["seed"], raw["replicates"], raw["sizes"])
    c, equation = raw["c"], raw["equation"]
    direct_interval = c * np.tan(np.pi/2 * np.array([.025, .975]))
    intervals = {op: lumped_quantiles(c, op)[[0, 2]] for op in OPERATIONS}
    rows = []
    for cell in cells:
        operation, truth = cell["composition"], cell["truth"]
        rank = _sum_cdf(truth/c) if operation == "L5" else 1 - _sum_cdf(c/truth)
        rows.append(dict(
            cell=cell, equation=equation,
            zero_parts=[bool(direct_interval[0] <= cell[key] <= direct_interval[1])
                        for key in ("m1", "m2")],
            zero_direct=bool(direct_interval[0] <= truth <= direct_interval[1]),
            zero_computed=bool(intervals[operation][0] <= truth <= intervals[operation][1]),
            zero_computed_cdf=rank, seed=[], direct=[], replacement=[], grid_rechecks=0,
        ))

    for replicate in range(raw["replicates"]):
        fitted = {}
        for row in rows:
            cell = row["cell"]
            seed = seed_stream[cell["size_index"]*raw["replicates"] + replicate][0]
            operation, truth = cell["composition"], cell["truth"]
            # Slots one and two consume the same RNG draws regardless of mass.
            # Thus the c and 4c direct data can be shared across both operations.
            slot = 4 if cell["truth_case"] == "historical" and operation == "L14" else 3
            key = (seed, slot, truth)
            if key not in fitted:
                readings = _direct_readings(seed, cell)
                fitted[key] = analyse(readings, c, equation, raw.get("nu_eff"))
            own = fitted[key]
            log_truth = math.log(truth)
            lo, hi = own.log_interval()
            changed = reweight(own, c, operation)
            smallest_tail = math.exp(min(changed.logcdf(log_truth)))
            covered = smallest_tail >= .025
            if abs(smallest_tail - .025) < .0001:
                refined = reweight(own, c, operation, refine=True)
                if covered != (min(refined.logcdf(log_truth)) >= math.log(.025)):
                    raise RuntimeError(f"grid refinement changed coverage in {equation} {cell['id']} seed {seed}")
                row["grid_rechecks"] += 1
            row["seed"].append(seed)
            row["direct"].append(lo <= log_truth <= hi)
            row["replacement"].append(covered)

    for row in rows:
        delta = np.subtract(row["replacement"], row["direct"], dtype=float)
        row.update(
            direct_coverage=float(np.mean(row["direct"])),
            replacement_coverage=float(np.mean(row["replacement"])),
            gain=float(delta.mean()),
            paired_mcse=float(delta.std(ddof=1) / math.sqrt(len(delta))),
            gained=int(np.count_nonzero(delta > 0)), lost=int(np.count_nonzero(delta < 0)),
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/b26-coverage.json")
    args = parser.parse_args()
    results = []
    for equation in ("ours30", "ours32", "ours33"):
        preset = ROOT / "presets" / f"composition_coverage_{equation}.yaml"
        spec = validate(yaml.safe_load(preset.read_bytes()))
        results.extend(replay(spec))
        print(f"{equation}: six weak cells complete", flush=True)
    quantiles = {op: lumped_quantiles(1., op).tolist() for op in OPERATIONS}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(unit_quantiles=quantiles, rows=results), indent=2) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
