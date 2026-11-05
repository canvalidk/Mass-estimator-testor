"""B18's explicitly selected extension of the locked B03 study.

The independent-measurement acquisition and product propagation are the
historical benchmark, not a declaration of physical shared-force or
shared-acceleration noise. Warping order and post-warp eta choices are
reported sensitivity procedures; none is adopted as an equation.
"""
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
from pathlib import Path

import numpy as np
from scipy.special import log_ndtr, ndtri
import yaml

from .composition_coverage import FIXED, S, _series, jobs


PROCEDURES = ("warp_parts", "warp_after_min", "warp_after_geometric")
LOG_ALPHA = math.log(.025)
Z = float(ndtri(.975))


def validate(data):
    required = {"kind", "protocol", "name", "question", "replicates", "seed",
                "sizes", "truth_cases", *FIXED}
    if required - data.keys():
        raise ValueError(f"B18 composition_coverage needs {sorted(required)}")
    if set(data) - required - {"notes", "nu_eff"}:
        raise ValueError("unknown B18 composition_coverage setting")
    if data["kind"] != "composition_coverage" or data["protocol"] != "b18":
        raise ValueError("expected composition_coverage protocol b18")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..") or
            any(ch not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for ch in name)):
        raise ValueError("composition_coverage name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("composition_coverage needs a nonempty question")
    for key, minimum in (("replicates", 2), ("seed", 0)):
        if type(data[key]) is not int or data[key] < minimum:
            raise ValueError(f"composition_coverage {key} must be an integer >= {minimum}")
    sizes = data["sizes"]
    if (not isinstance(sizes, list) or not sizes or any(
            isinstance(v, bool) or not isinstance(v, (int, float)) or
            not math.isfinite(v) or v <= 0 for v in sizes) or len(set(sizes)) != len(sizes)):
        raise ValueError("sizes must be unique finite positive numbers")
    cases = data["truth_cases"]
    if (not isinstance(cases, list) or not cases or any(
            value not in ("historical", "c", "4c") for value in cases) or len(set(cases)) != len(cases)):
        raise ValueError("truth_cases must select unique historical, c or 4c cases")
    equation = data["equation"]
    if equation not in ("ours30", "ours32", "ours33"):
        raise ValueError("B18 equation must be ours30, ours32 or ours33")
    for key, value in FIXED.items():
        if key == "equation":
            continue
        if isinstance(data[key], bool) or data[key] != value:
            raise ValueError(f"B18 preserves acquisition {key}={value!r}")
        if key in ("dimension", "n_readings", "draws") and type(data[key]) is not int:
            raise ValueError(f"composition_coverage {key} must be an integer")
    if equation in ("ours32", "ours33") and "nu_eff" not in data:
        raise ValueError("ours32/ours33 need explicit nu_eff")
    if "nu_eff" in data and (isinstance(data["nu_eff"], bool) or
            not isinstance(data["nu_eff"], (int, float)) or
            not math.isfinite(data["nu_eff"]) or data["nu_eff"] <= 0):
        raise ValueError("nu_eff must be finite and positive")
    raw = deepcopy(data)
    raw["sizes"] = [float(value) for value in sizes]
    cells = []
    for si, size in enumerate(raw["sizes"]):
        for case in cases:
            for composition in ("L5", "L14"):
                m1, m2, truth = masses(case, composition)
                cells.append(dict(id=f"c{len(cells):03d}", size=size, size_index=si,
                                  truth_case=case, composition=composition,
                                  m1=m1, m2=m2, truth=truth, target_multiple=truth/S))
    return dict(kind="composition_coverage", raw=raw,
                settings={**deepcopy(FIXED), "equation": equation, "protocol": "b18"}, cells=cells)


def masses(case, composition):
    if composition not in ("L5", "L14"):
        raise ValueError("composition must be L5 or L14")
    if case == "historical":
        return (1., 2., 3. if composition == "L5" else 2./3.)
    if case not in ("c", "4c"):
        raise ValueError("unknown truth case")
    truth = S * (1 if case == "c" else 4)
    return ((truth/3, 2*truth/3, truth) if composition == "L5" else
            (1.5*truth, 3*truth, truth))


def post_eta(eta1, eta2, procedure):
    if procedure == "warp_after_min":
        return min(eta1, eta2)
    if procedure == "warp_after_geometric":
        return math.sqrt(eta1)*math.sqrt(eta2)
    raise ValueError("unknown post-warp procedure")


def _decision(log_probability, log_error, log_threshold):
    """Nominal, lower and upper threshold decisions under an absolute error."""
    lower = (-math.inf if log_error >= log_probability else
             log_probability + math.log(-math.expm1(log_error-log_probability)))
    upper = float(np.logaddexp(log_probability, log_error))
    return (log_probability >= log_threshold, lower >= log_threshold, upper >= log_threshold)


def _one(arguments):
    from .composition_laws import analyse, composed_logcdf
    seed, cell, raw = arguments
    rng = np.random.default_rng(seed)
    equation, nu = raw["equation"], raw.get("nu_eff")
    m1, m2, truth = cell["m1"], cell["m2"], cell["truth"]
    rd1, rd2 = _series(rng, m1, cell["size"]), _series(rng, m2, cell["size"])
    if cell["truth_case"] == "historical":
        # Preserve the original acquisition RNG positions: its third series
        # is the lump and its fourth the parallel sum, before any propagation.
        lump, parallel = _series(rng, 3., cell["size"]), _series(rng, 2./3., cell["size"])
        own_readings = lump if cell["composition"] == "L5" else parallel
    else:
        own_readings = _series(rng, truth, cell["size"])
    p1, p2, own = [analyse(rd, S, equation, nu) for rd in (rd1, rd2, own_readings)]
    lo, hi = own.log_interval()
    log_truth = math.log(truth)
    own_covers = lo <= log_truth <= hi
    parts_cover = [a <= math.log(m) <= b for (a, b), m in
                   zip((p1.log_interval(), p2.log_interval()), (m1, m2))]
    procedures = PROCEDURES if equation == "ours33" else ("quoted_parts",)
    # Refine only when the monotone integration bounds can change a verdict.
    # Starting every unambiguous decision with 20,000 nodes wastes work;
    # extreme post-warp ranks still use exactly the same full-tail CDF.
    cache, results, nodes = {}, {}, min(256, raw["draws"])
    maximum_nodes = 16*raw["draws"]

    def cdf(value, quoted):
        key = value, quoted, nodes
        if key not in cache:
            cache[key] = composed_logcdf(p1, p2, value, cell["composition"], quoted=quoted, nodes=nodes)
        return cache[key]

    while True:
        results = {}
        for procedure in procedures:
            quoted = procedure in ("warp_parts", "quoted_parts")
            threshold = (LOG_ALPHA if quoted else
                         float(log_ndtr(-Z/math.sqrt(post_eta(p1.eta, p2.eta, procedure)))))
            tc, ts, te = cdf(log_truth, quoted)
            hc, _, he = cdf(hi, quoted)
            _, ls, le = cdf(lo, quoted)
            coverage_checks = [_decision(tc, te, threshold), _decision(ts, te, threshold)]
            overlap_checks = [_decision(hc, he, threshold), _decision(ls, le, threshold)]
            cv = tuple(all(check[i] for check in coverage_checks) for i in range(3))
            ov = tuple(all(check[i] for check in overlap_checks) for i in range(3))
            results[procedure] = dict(computed_covers=cv[0], computed_lower=cv[1], computed_upper=cv[2],
                                      own_covers=own_covers, overlap=ov[0], overlap_lower=ov[1],
                                      overlap_upper=ov[2], unresolved=cv[1] != cv[2] or ov[1] != ov[2])
        if not any(value["unresolved"] for value in results.values()) or nodes >= maximum_nodes:
            break
        nodes = min(nodes*4, maximum_nodes)
    return dict(seed=seed, cell_id=cell["id"], parts_cover=parts_cover,
                eta_parts=[p1.eta, p2.eta], eta_own=own.eta, nodes=nodes, procedures=results)


def _mean(values):
    values = np.asarray(values, dtype=float)
    return float(values.mean()), float(values.std(ddof=1)/math.sqrt(len(values)))


def summarize(results, cells):
    rows = []
    for cell in cells:
        block = [value for value in results if value["cell_id"] == cell["id"]]
        for procedure in block[0]["procedures"]:
            row = {key: cell[key] for key in ("size", "truth_case", "composition", "m1", "m2", "truth", "target_multiple")}
            row.update(cell=cell["id"], procedure=procedure, replicates=len(block))
            for key in block[0]["procedures"][procedure]:
                row[key], row[key+"_mcse"] = _mean([v["procedures"][procedure][key] for v in block])
            row["parts_coverage"], row["parts_coverage_mcse"] = _mean([np.mean(v["parts_cover"]) for v in block])
            row["parts_cap"], row["parts_cap_mcse"] = _mean([np.mean(np.array(v["eta_parts"]) == 1.) for v in block])
            row["own_cap"], row["own_cap_mcse"] = _mean([v["eta_own"] == 1. for v in block])
            row["max_nodes"] = max(v["nodes"] for v in block)
            rows.append(row)
    return rows


def run(path, out_root="results", replicates=None, progress=print, workers=None):
    from .runner import code_version
    path = Path(path)
    source = path.read_bytes()
    raw = validate(yaml.safe_load(source))["raw"]
    if replicates is not None:
        if type(replicates) is not int or replicates < 2:
            raise ValueError("replicates must be an integer >= 2")
        raw["replicates"] = replicates
    workers = 1 if workers is None else workers
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    spec = validate(raw)
    work_seeds = jobs(raw["seed"], raw["replicates"], raw["sizes"])
    work = [(work_seeds[cell["size_index"]*raw["replicates"]+r][0], cell, raw)
            for cell in spec["cells"] for r in range(raw["replicates"])]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    base = Path(out_root)/raw["name"]
    out, suffix = base/stamp, 1
    while out.exists():
        suffix += 1
        out = base/f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    (out/"source_preset.yaml").write_bytes(source)
    (out/"preset.yaml").write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    code = code_version()
    code["source_preset_sha256"] = hashlib.sha256(source).hexdigest()
    code["effective_preset_sha256"] = hashlib.sha256((out/"preset.yaml").read_bytes()).hexdigest()
    for name in ("composition_probe.py", "composition_laws.py", "composition_coverage.py"):
        code[name+"_sha256"] = hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
    results = []
    executor = None if workers == 1 else ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"))
    iterator = map(_one, work) if executor is None else executor.map(_one, work, chunksize=4)
    try:
        for i, result in enumerate(iterator):
            results.append(result)
            if (i+1) % raw["replicates"] == 0:
                (out/"replicates.json").write_text(json.dumps(results, allow_nan=False)+"\n", encoding="utf-8")
                if progress:
                    progress(f"{raw['equation']} {result['cell_id']}: {raw['replicates']} replicates complete")
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    summary = summarize(results, spec["cells"])
    record = dict(kind="composition_coverage", protocol="b18", equation=raw["equation"],
                  nu_eff=raw.get("nu_eff"), seed=raw["seed"], replicates=raw["replicates"],
                  sizes=raw["sizes"], truth_cases=raw["truth_cases"], cells=spec["cells"],
                  draws=raw["draws"], c=S, workers=workers, code=code, summary=summary,
                  sampling="independent parts and independent direct composite; historical acquisition",
                  propagation="independent product benchmark; midpoint probability quadrature with log tails",
                  integration=dict(initial_nodes=min(256, raw["draws"]), maximum_nodes=16*raw["draws"],
                                   refinement="multiply by four while midpoint bounds straddle a decision"),
                  eta_after="min(eta1,eta2) and sqrt(eta1*eta2), unforced sensitivity procedures",
                  uncertainty="midpoint bounds conditional on marginal CDF grids; no certified total numerical bound",
                  seed_layout="historical size-major uint32 seeds; same seed reused across equations and truth cases")
    (out/"run.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    with (out/"summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    arrays = dict(job_seeds=np.array([v["seed"] for v in results], dtype=np.uint32),
                  seed=np.array([v["seed"] for v in results], dtype=np.uint32),
                  size=np.array([next(c["size"] for c in spec["cells"] if c["id"] == v["cell_id"]) for v in results]),
                  cell_id=np.array([v["cell_id"] for v in results]),
                  parts_cover=np.array([v["parts_cover"] for v in results], dtype=bool),
                  eta_parts=np.array([v["eta_parts"] for v in results]),
                  eta_own=np.array([v["eta_own"] for v in results]),
                  nodes=np.array([v["nodes"] for v in results]))
    for procedure in results[0]["procedures"]:
        for key in results[0]["procedures"][procedure]:
            arrays[procedure+"_"+key] = np.array([v["procedures"][procedure][key] for v in results], dtype=bool)
    np.savez_compressed(out/"values.npz", **arrays)
    if progress:
        progress(f"saved {out}")
    return out, record
