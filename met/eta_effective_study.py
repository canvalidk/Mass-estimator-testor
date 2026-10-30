"""B04 repeated-world check of fitted-score bias; eta-prime is study-only.

The regression identity has fixed nonnegative curvatures. Here the fitted
curvatures are random and can be signed, so the report separates means of
ratios, ratios of means and their paired discrepancies. Neither coverage
column selects a replacement for the production Bartlett power.
"""
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy
import yaml

from .eta_effective import intervals, statistics
from .world import generate, series_rngs


def validate(data):
    from .preset import expand

    required = {"kind", "name", "question", "world", "replicates", "seed"}
    if not isinstance(data, dict) or required - data.keys():
        raise ValueError(f"eta_effective needs {sorted(required)}")
    if set(data) - required - {"notes", "numerics"}:
        raise ValueError("unknown eta_effective setting")
    if data["kind"] != "eta_effective":
        raise ValueError("expected kind: eta_effective")
    name = data["name"]
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_" for c in name)):
        raise ValueError("eta_effective name must be a simple identifier")
    if not isinstance(data["question"], str) or not data["question"].strip():
        raise ValueError("eta_effective needs a nonempty question")
    if type(data["seed"]) is not int or data["seed"] < 0:
        raise ValueError("eta_effective seed must be a nonnegative integer")
    if type(data["replicates"]) is not int or data["replicates"] < 2:
        raise ValueError("eta_effective replicates must be an integer >= 2")
    settings = {"grid_points": 2001, "batch_size": 128}
    supplied = data.get("numerics", {})
    if not isinstance(supplied, dict) or set(supplied) - settings.keys():
        raise ValueError("eta_effective numerics accepts grid_points and batch_size only")
    settings.update(supplied)
    if type(settings["grid_points"]) is not int or settings["grid_points"] < 201:
        raise ValueError("eta_effective grid_points must be an integer >= 201")
    if type(settings["batch_size"]) is not int or settings["batch_size"] < 1:
        raise ValueError("eta_effective batch_size must be a positive integer")
    if not isinstance(data["world"], dict):
        raise ValueError("eta_effective world must be a mapping")
    cells = expand(data["world"])
    for cell in cells:
        world = cell["world"]
        if world["dimension"] not in (1, 2, 3):
            raise ValueError("eta_effective requires dimension 1, 2 or 3")
        if isinstance(world["acceleration_snr"], (bool, dict)):
            raise ValueError("eta_effective needs a fixed scalar acceleration_snr")
        if any(isinstance(world["noise"][key], dict) for key in ("force_sd", "acceleration_sd")):
            raise ValueError("eta_effective requires isotropic channels with one SD each")
        if world["supplied_noise"] != "exact":
            raise ValueError("eta_effective compares against exact supplied covariance")
    return {"kind": "eta_effective", "raw": data, "cells": cells, "numerics": settings}


def _summary(values):
    """Unconditional summaries are undefined if any replicate is undefined."""
    values = np.asarray(values, dtype=float)
    finite_count = int(np.isfinite(values).sum())
    valid = finite_count == len(values)
    sd = float(np.std(values, ddof=1)) if valid else None
    quantiles = np.quantile(values, [0.05, 0.5, 0.95]) if valid else [None] * 3
    return dict(mean=float(values.mean()) if valid else None,
                sd=sd, mc_se=sd / math.sqrt(len(values)) if valid else None,
                q05=float(quantiles[0]) if valid else None,
                q50=float(quantiles[1]) if valid else None,
                q95=float(quantiles[2]) if valid else None,
                finite_count=finite_count, undefined_count=len(values) - finite_count)


def _coverage(covered):
    values = np.asarray(covered, dtype=float)
    valid = bool(np.isfinite(values).all())
    p = float(values.mean()) if valid else None
    return dict(coverage=p, mc_se=math.sqrt(p * (1 - p) / len(values)) if valid else None,
                count=len(values), undefined_count=int((~np.isfinite(values)).sum()))


def summarize(values):
    """Finite-sample summaries; some empirical SEs have no finite-variance CLT.

    In particular, summary ratio/linear_target/paired_gap mc_se values and
    expectation_ratio.paired_gap_mc_se are empirical diagnostics only. See
    run.json's moment_status and O14's B04 cancellation argument.
    """
    scalar = {name: _summary(value) for name, value in values.items()
              if name not in ("covered", "prime_covered")}
    c, j, target = values["C"], values["J"], values["linear_target"]
    ratio = weighted_target = gap = se = gap_se = None
    finite = np.isfinite(c) & np.isfinite(j) & np.isfinite(target)
    denominator_positive = bool(np.isfinite(c).all() and np.mean(c) > 0)
    valid = bool(finite.all() and denominator_positive)
    if valid:
        mean_c = float(c.mean())
        ratio = float(j.mean() / mean_c)
        weighted_target = float(np.mean(c * target) / mean_c)
        gap = ratio - weighted_target
        se = float(np.std((j - ratio * c) / mean_c, ddof=1) / math.sqrt(len(c)))
        influence = (j - c * target - gap * c) / mean_c
        gap_se = float(np.std(influence, ddof=1) / math.sqrt(len(c)))
    return dict(summary=scalar,
                expectation_ratio=dict(value=ratio, mc_se=se,
                                       curvature_weighted_target=weighted_target,
                                       paired_gap=gap, paired_gap_mc_se=gap_se,
                                       undefined_count=int((~finite).sum()),
                                       denominator_positive=denominator_positive),
                coverage=dict(ours30=_coverage(values["covered"]),
                              eta_prime=_coverage(values["prime_covered"]),
                              paired_change=_summary(values["prime_covered"] - values["covered"])))


def _cell_values(readings, truth, numerics):
    parts = []
    for start in range(0, readings.force.shape[0], numerics["batch_size"]):
        part = readings.subset(slice(start, start + numerics["batch_size"]))
        stats = statistics(part)
        bounds = intervals(part, stats, grid_points=numerics["grid_points"])
        # Retain every scalar statistic and interval, not large per-reading arrays.
        values = {name: np.asarray(value).copy() for name, value in stats.items()
                  if np.asarray(value).ndim == 1}
        curvatures = np.asarray(stats["curvatures"])
        values["negative_curvature_fraction"] = np.mean(curvatures < 0, axis=1)
        values["any_negative_curvature"] = np.any(curvatures < 0, axis=1).astype(float)
        values["neff_le_one"] = (np.asarray(stats["n_eff"]) <= 1).astype(float)
        values["degenerate_curvature"] = (np.asarray(stats["C"]) <= 0).astype(float)
        values["paired_gap"] = values["ratio"] - values["linear_target"]
        values.update({name: np.asarray(value).copy() for name, value in bounds.items()})
        for prefix, out_name in (("interval", "covered"), ("prime", "prime_covered")):
            lo, hi = values[f"{prefix}_lo"], values[f"{prefix}_hi"]
            valid = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi >= lo)
            values[out_name] = np.where(valid, ((lo <= truth) & (truth <= hi)).astype(float), np.nan)
            values[f"{prefix}_log_width"] = np.where(valid, np.log(hi) - np.log(lo), np.nan)
        parts.append(values)
    return {name: np.concatenate([part[name] for part in parts]) for name in parts[0]}


def _fmt(value):
    return "undefined" if value is None else f"{value:.6g}"


def _markdown(record):
    lines = ["# Effective readings and fitted-score bias", "", record["question"], "",
             f"Seed {record['seed']}; {record['replicates']} paired series per cell. All runs are serial.",
             "C is random here. The ratio-of-means estimate E[J]/E[C] and the sample mean of J/C are reported separately.",
             "Targets 1-1/n_eff use signed fitted curvatures. Negative targets and n_eff<=1 are retained.",
             "For N>=2 at every finite signal strength, cancellation near C=0 gives target expectation -infinity, gap expectation +infinity and infinite second moment of J/C (O14, B04).",
             "SEs for sample J/C, target and gap, including the C-weighted gap, are empirical diagnostics with no finite-variance CLT justification. The ratio-of-means E[J]/E[C] delta-method SE and binomial coverage SEs remain ordinary Monte Carlo SEs.",
             "Any undefined replicate makes its unconditional summary undefined; counts are reported in run.json.",
             "Eta-prime is a study candidate, not adopted. Coverage is descriptive, not a selection rule.", "",
             "| Cell | k | SNR | N | E[J]/E[C] estimate (MC SE) | sample C-weighted target | J/C sample mean (empirical SE) | target sample mean | gap sample mean (empirical SE) | J/C sample 5%,50%,95% | n_eff<=1 | negative c_i fraction |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in record["cells"]:
        w, s, e = row["world"], row["summary"], row["expectation_ratio"]
        lines.append(f"| {row['id']} | {w['dimension']} | {w['acceleration_snr']:g} | {w['readings']} | "
                     f"{_fmt(e['value'])} ({_fmt(e['mc_se'])}) | {_fmt(e['curvature_weighted_target'])} | "
                     f"{_fmt(s['ratio']['mean'])} ({_fmt(s['ratio']['mc_se'])}) | {_fmt(s['linear_target']['mean'])} | "
                     f"{_fmt(s['paired_gap']['mean'])} ({_fmt(s['paired_gap']['mc_se'])}) | "
                     f"{_fmt(s['ratio']['q05'])}, {_fmt(s['ratio']['q50'])}, {_fmt(s['ratio']['q95'])} | "
                     f"{_fmt(s['neff_le_one']['mean'])} | {_fmt(s['negative_curvature_fraction']['mean'])} |")
    lines += ["", "| Cell | k | SNR | N | eta mean / SD | eta-prime mean / SD | ours30 coverage (MC SE) | eta-prime coverage (MC SE) | paired change (MC SE) | log-width / prime |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in record["cells"]:
        w, s, c = row["world"], row["summary"], row["coverage"]
        lines.append(f"| {row['id']} | {w['dimension']} | {w['acceleration_snr']:g} | {w['readings']} | "
                     f"{_fmt(s['eta']['mean'])} / {_fmt(s['eta']['sd'])} | "
                     f"{_fmt(s['eta_prime']['mean'])} / {_fmt(s['eta_prime']['sd'])} | "
                     f"{_fmt(c['ours30']['coverage'])} ({_fmt(c['ours30']['mc_se'])}) | "
                     f"{_fmt(c['eta_prime']['coverage'])} ({_fmt(c['eta_prime']['mc_se'])}) | "
                     f"{_fmt(c['paired_change']['mean'])} ({_fmt(c['paired_change']['mc_se'])}) | "
                     f"{_fmt(s['interval_log_width']['mean'])} / {_fmt(s['prime_log_width']['mean'])} |")
    return "\n".join(lines) + "\n"


def run(path, out_root="results", replicates=None, progress=print):
    """Run the saved paired experiment and retain every replicate for replay."""
    from .runner import code_version

    source = Path(path).read_bytes()
    spec = validate(yaml.safe_load(source))
    reps = spec["raw"]["replicates"] if replicates is None else replicates
    if type(reps) is not int or reps < 2:
        raise ValueError("eta_effective replicates must be an integer >= 2")
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base = Path(out_root) / spec["raw"]["name"]
    out, suffix = base / stamp, 1
    while out.exists():
        suffix += 1
        out = base / f"{stamp}-{suffix}"
    out.mkdir(parents=True)
    if replicates is None:
        (out / "preset.yaml").write_bytes(source)
    else:
        effective = dict(spec["raw"], replicates=reps)
        (out / "source_preset.yaml").write_bytes(source)
        (out / "preset.yaml").write_text(yaml.safe_dump(effective, sort_keys=False), encoding="utf-8")
    start, arrays, results = time.perf_counter(), {}, []
    for cell in spec["cells"]:
        if progress:
            progress(f"cell {cell['id']} {cell['axes']}")
        world = cell["world"]
        readings, truth = generate(world, series_rngs(spec["raw"]["seed"], world, reps))
        values = _cell_values(readings, truth["mass"], spec["numerics"])
        arrays.update({f"{cell['id']}/{name}": value for name, value in values.items()})
        results.append({**cell, **summarize(values)})
    code = code_version()
    code.update(scipy=scipy.__version__, preset_sha256=hashlib.sha256(source).hexdigest(),
                effective_preset_sha256=hashlib.sha256((out / "preset.yaml").read_bytes()).hexdigest())
    record = dict(kind="eta_effective", preset=spec["raw"]["name"], question=spec["raw"]["question"],
                  seed=spec["raw"]["seed"], replicates=reps, replicates_override=replicates is not None,
                  workers=1, numerics=spec["numerics"], code=code, started=stamp,
                  seconds=time.perf_counter() - start, cells=results,
                  moment_status={
                      "scope": "N>=2, k=1,2,3, every finite SNR under the full-support Gaussian world; see O14's B04 cancellation proof.",
                      "linear_target_expectation": "-infinity",
                      "paired_gap_expectation": "+infinity",
                      "ratio_second_moment": "infinite",
                      "curvature_weighted_target_numerator_second_moment": "infinite",
                      "curvature_weighted_target_expectation": "not established by the local cancellation argument",
                      "empirical_se_only": ["summary.ratio.mc_se", "summary.linear_target.mc_se",
                                            "summary.paired_gap.mc_se", "expectation_ratio.paired_gap_mc_se"],
                      "ordinary_mc_se": ["expectation_ratio.mc_se", "coverage.ours30.mc_se",
                                         "coverage.eta_prime.mc_se", "coverage.paired_change.mc_se"]},
                  scope="Isotropic independent Gaussian channels; ours30 g; whole-line fitted eta; half-Cauchy c=s. Eta-prime is study-only.")
    np.savez_compressed(out / "values.npz", **arrays)
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False), encoding="utf-8")
    (out / "summary.md").write_text(_markdown(record), encoding="utf-8")
    fields = ["cell", "dimension", "mass", "snr", "readings", "J_over_C_sample_mean", "J_over_C_empirical_se", "J_over_C_sample_sd",
              "target_sample_mean", "gap_sample_mean", "gap_empirical_se", "ratio_of_means", "ratio_of_means_mc_se", "sample_weighted_target",
              "ratio_gap_sample", "ratio_gap_empirical_se", "neff_le_one_rate", "negative_curvature_fraction",
              "eta_mean", "eta_prime_mean", "coverage", "coverage_mc_se", "prime_coverage", "prime_coverage_mc_se",
              "coverage_change", "coverage_change_mc_se", "ratio_undefined_count"]
    with (out / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in results:
            w, s, e, c = row["world"], row["summary"], row["expectation_ratio"], row["coverage"]
            writer.writerow(dict(cell=row["id"], dimension=w["dimension"], mass=w["mass"], snr=w["acceleration_snr"],
                readings=w["readings"], J_over_C_sample_mean=s["ratio"]["mean"], J_over_C_empirical_se=s["ratio"]["mc_se"],
                J_over_C_sample_sd=s["ratio"]["sd"], target_sample_mean=s["linear_target"]["mean"], gap_sample_mean=s["paired_gap"]["mean"],
                gap_empirical_se=s["paired_gap"]["mc_se"], ratio_of_means=e["value"], ratio_of_means_mc_se=e["mc_se"],
                sample_weighted_target=e["curvature_weighted_target"], ratio_gap_sample=e["paired_gap"], ratio_gap_empirical_se=e["paired_gap_mc_se"],
                neff_le_one_rate=s["neff_le_one"]["mean"], negative_curvature_fraction=s["negative_curvature_fraction"]["mean"],
                eta_mean=s["eta"]["mean"], eta_prime_mean=s["eta_prime"]["mean"], coverage=c["ours30"]["coverage"],
                coverage_mc_se=c["ours30"]["mc_se"], prime_coverage=c["eta_prime"]["coverage"],
                prime_coverage_mc_se=c["eta_prime"]["mc_se"], coverage_change=c["paired_change"]["mean"],
                coverage_change_mc_se=c["paired_change"]["mc_se"], ratio_undefined_count=s["ratio"]["undefined_count"]))
    if progress:
        progress(f"wrote {out}")
    return out, record
