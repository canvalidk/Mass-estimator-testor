"""Reproducible full-positive-law T56/T57 Gaussian witness study.

Population-scatter calculations are labeled separately from simulated data.
The preset's anchor exponent is a research choice, not a selected rule.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import stats
import yaml

from met.score_law_obstruction import ScalarScoreQuadrature


def write_summary(records, path):
    lines = ["# Full positive scalar score-law sanity check", "",
        "The true mass and covariance are recorded in preset.yaml; run.json records the false root.",
        "Coverage is descriptive at the recorded small world count; see run.json for scope.", "",
        "| Readings | Method | Worlds | Median of medians | Medians on false side | Coverage | 95% Monte Carlo interval |",
        "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in records:
        lines.append(f"| {row['count']} | {row['method']} | {row['replicates']} | "
            f"{float(row['median_of_medians']):.6f} | {row['medians_on_false_side']} | "
            f"{float(row['coverage']):.2%} | {float(row['mc_low']):.2%}–{float(row['mc_high']):.2%} |")
    lines += ["", "Population-scatter calculations are in population_scatter.csv; they are not simulated coverage.",
        "The anchored comparison uses the preset's unforced research exponent. No law is selected by this table."]
    Path(path).write_text("\n".join(lines)+"\n", encoding="utf-8")


def validate(config):
    required = {"kind", "name", "question", "seed", "replicates", "mass", "acceleration",
        "correlation", "prior_center", "anchor_exponent", "counts", "population_counts",
        "probability", "epsabs", "epsrel", "limit"}
    if not isinstance(config, dict) or set(config) != required or config["kind"] != "score_law_sanity":
        raise ValueError("score_law_sanity requires its explicit settings")
    if config["name"] != "score-law-sanity":
        raise ValueError("this bounded study uses the name score-law-sanity")
    for key in ("seed", "replicates", "limit"):
        if type(config[key]) is not int or config[key] < (0 if key == "seed" else 2):
            raise ValueError("invalid integer control")
    for key in ("mass", "acceleration", "prior_center", "epsabs", "epsrel"):
        if not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError("positive finite model and integration controls required")
    if not -1 < config["correlation"] < 1 or not 0 < config["anchor_exponent"] < 1 or not 0 < config["probability"] < 1:
        raise ValueError("correlation/probability/anchor exponent out of range")
    for key in ("counts", "population_counts"):
        if not config[key] or any(type(n) is not int or n < 2 for n in config[key]):
            raise ValueError("nonempty integer count lists required")
        if len(set(config[key])) != len(config[key]):
            raise ValueError("count lists must contain distinct values for replay")
    mass, rho = config["mass"], config["correlation"]
    if mass == rho or not 0 < (1-rho*mass)/(rho-mass) < mass:
        raise ValueError("this bounded study requires a lower positive orthogonal branch")
    return dict(kind="score_law_sanity", raw=config, cells=config["counts"])


def run(path, out_root="results", *, replicates=None, progress=print):
    config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if replicates is not None:
        config["replicates"] = replicates
    config = validate(config)["raw"]
    output = Path(out_root)/config["name"]/datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")
    output.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for name in ("score_law_study.py", "score_law_obstruction.py"):
        source = Path(__file__).with_name(name).read_bytes()
        (output/name).write_bytes(source)
        hashes[name] = hashlib.sha256(source).hexdigest()
    rng = np.random.default_rng(config["seed"])
    mass, rho, center = config["mass"], config["correlation"], config["prior_center"]
    covariance = np.array([[1., rho], [rho, 1.]])
    mean = config["acceleration"]*np.array([mass, 1.])
    false = (1-rho*mass)/(rho-mass)
    if false <= 0 or false >= mass:
        raise ValueError("this study requires its lower positive orthogonal branch")
    split = np.sqrt(mass*false)  # Investigator scoring cut; never a kernel input.
    ranks = [(1-config["probability"])/2, .5, (1+config["probability"])/2]
    rows, population, arrays = [], [], {}
    methods = ("raw", "score", "anchored_score")

    def evaluate(scatter):
        q = ScalarScoreQuadrature(scatter, covariance, prior_center=center,
            anchor_exponent=config["anchor_exponent"], epsabs=config["epsabs"],
            epsrel=config["epsrel"], limit=config["limit"])
        return [dict(method=method, false_band_probability=q.cdf(split, method),
            lower=q.quantile(ranks[0], method), median=q.quantile(.5, method),
            upper=q.quantile(ranks[2], method), integral=q.total(method)[0],
            reported_integral_error=q.total(method)[1], anchor_weight=q.anchor_weight)
            for method in methods]

    for count in config["population_counts"]:
        for result in evaluate(count*(covariance+np.outer(mean, mean))):
            population.append(dict(count=count, **result))
    for count in config["counts"]:
        observations = mean+rng.multivariate_normal(np.zeros(2), covariance,
                                                   size=(config["replicates"], count))
        arrays[f"observations_{count}"] = observations
        for rep, data in enumerate(observations):
            for result in evaluate(data.T@data):
                rows.append(dict(count=count, replicate=rep, **result))
        if progress:
            progress(f"score-law count {count}: {config['replicates']} full-law worlds complete")
    summary = []
    for count in config["counts"]:
        for method in methods:
            selected = [r for r in rows if r["count"] == count and r["method"] == method]
            covered = sum(r["lower"] <= mass <= r["upper"] for r in selected)
            n = len(selected)
            summary.append(dict(count=count, method=method, replicates=n, coverage=covered/n,
                mc_low=float(stats.beta.ppf(.025, covered, n-covered+1)) if covered else 0.,
                mc_high=float(stats.beta.ppf(.975, covered+1, n-covered)) if covered < n else 1.,
                median_of_medians=float(np.median([r["median"] for r in selected])),
                medians_on_false_side=sum(r["median"] < split for r in selected),
                mean_false_band_probability=float(np.mean([r["false_band_probability"] for r in selected]))))
    for name, records in (("summary", summary), ("population_scatter", population), ("worlds", rows)):
        with (output/(name+".csv")).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    np.savez_compressed(output/"values.npz", **arrays)
    write_summary(summary, output/"summary.md")
    (output/"preset.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    metadata = dict(preset=config, source_sha256=hashes, false_root=float(false),
        scoring_split=float(split), limitations=[
        "Population-scatter rows are deterministic, not sample coverage.",
        "Adaptive full-domain integrals have ordinary numerical error estimates, not certified tails.",
        "Anchor exponent and rule are unforced research choices; no adoption by score.",
        "32-world default coverage is descriptive; the analytic counterexample is decisive."])
    (output/"run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if progress:
        progress(str(output))
    return output, dict(summary=summary, **metadata)
