"""T64 scalar Gaussian sanity study with exact heterogeneous scatter draws.

This is a triangular-array counterpart of T63, not its individual fixed
prefix history: M informative events and M^4 central-noise events with
covariance diag(M^8,1). Both groups are sampled exactly, never omitted.
Counts/truth enter generation and scoring only, not the law constructor.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import stats
import yaml

from met.mass_contrast_scalar import ScalarMassContrastQuadrature
from met.weak_anchor_study import sample_scatter


def validate(config):
    required = {"kind", "name", "question", "seed", "replicates", "informative_counts",
        "mass", "acceleration", "correlation", "prior_center", "prior_bounds",
        "anchor_prefactor", "prior_nodes", "epsabs", "epsrel", "limit", "probability"}
    if not isinstance(config, dict) or set(config) != required or config["kind"] != "mass_contrast_sanity":
        raise ValueError("mass_contrast_sanity requires all explicit study declarations")
    if config["name"] != "mass-contrast-sanity":
        raise ValueError("the bounded study uses the name mass-contrast-sanity")
    for key, minimum in (("seed", 0), ("replicates", 2), ("prior_nodes", 2), ("limit", 2)):
        if type(config[key]) is not int or config[key] < minimum:
            raise ValueError("invalid integer control")
    for key in ("mass", "acceleration", "prior_center", "anchor_prefactor", "epsabs", "epsrel"):
        if not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError("positive finite model and numerical controls required")
    if not -1 < config["correlation"] < 1 or not 0 < config["probability"] < 1:
        raise ValueError("invalid correlation or probability")
    counts = config["informative_counts"]
    if (not counts or any(type(n) is not int or not 3 <= n <= 100000 for n in counts)
            or len(set(counts)) != len(counts)):
        raise ValueError("distinct informative counts in [3,100000] required")
    bounds = np.asarray(config["prior_bounds"], dtype=float)
    if bounds.shape != (2,) or np.any(~np.isfinite(bounds)) or not 0 < bounds[0] < bounds[1]:
        raise ValueError("the supplied prior needs finite positive bounds")
    mass, rho = config["mass"], config["correlation"]
    false = (rho*mass-1)/(mass-rho) if mass != rho else np.nan
    if not bounds[0] < false < mass < bounds[1]:
        raise ValueError("this witness requires both positive roots inside the supplied prior")
    return dict(kind="mass_contrast_sanity", raw=config, cells=counts)


def central_physical_scatter(rng, count, covariance):
    """Exact two-channel central Gaussian sufficient statistic."""
    bartlett = np.array([[np.sqrt(rng.chisquare(count)), 0.],
                         [rng.normal(), np.sqrt(rng.chisquare(count-1))]])
    factor = np.linalg.cholesky(covariance)@bartlett
    return factor@factor.T


def run(path, out_root="results", *, replicates=None, progress=print):
    config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if replicates is not None:
        config["replicates"] = replicates
    config = validate(config)["raw"]
    output = Path(out_root)/config["name"]/datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")
    output.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for name in ("mass_contrast_study.py", "mass_contrast_scalar.py", "weak_anchor_study.py"):
        source = Path(__file__).with_name(name).read_bytes()
        (output/name).write_bytes(source)
        hashes[name] = hashlib.sha256(source).hexdigest()
    mass, rho = config["mass"], config["correlation"]
    covariance = np.array([[1., rho], [rho, 1.]])
    mean = config["acceleration"]*np.array([mass, 1.])
    signal_energy = mean@np.linalg.solve(covariance, mean)
    false = (rho*mass-1)/(mass-rho)
    split = np.sqrt(mass*false)  # Investigator scoring cut, never passed to the law.
    ranks = [(1-config["probability"])/2, .5, (1+config["probability"])/2]
    controls = {key: config[key] for key in ("prior_center", "prior_bounds", "anchor_prefactor",
                                             "prior_nodes", "epsabs", "epsrel", "limit")}
    rng, rows, population, arrays = np.random.default_rng(config["seed"]), [], [], {}

    def evaluate(scatters, covariances):
        reference = ScalarMassContrastQuadrature(scatters, covariances, **controls)
        records = []
        for method in ("score", "trace", "contrast"):
            lower, median, upper = [reference.quantile(p, method) for p in ranks]
            records.append(dict(method=method, lower=lower, median=median, upper=upper,
                false_band_probability=reference.cdf(split, method),
                trace=reference.energy, contrast_energy=reference.contrast_energy,
                anchor_weight=reference.weights[method],
                shifted_normalizer=reference.total(method)[0],
                reported_integral_error=reference.total(method)[1]))
        return records

    for count in config["informative_counts"]:
        null_count = count**4
        null_covariance = np.diag([float(null_count)**2, 1.])
        covariances = [covariance, null_covariance]
        expected = [count*(covariance+np.outer(mean, mean)), null_count*null_covariance]
        for result in evaluate(expected, covariances):
            population.append(dict(informative_count=count, null_count=null_count, **result))
        scatters = []
        for rep in range(config["replicates"]):
            observed = [sample_scatter(rng, count, count*signal_energy, mass, covariance),
                        central_physical_scatter(rng, null_count, null_covariance)]
            scatters.append(observed)
            for result in evaluate(observed, covariances):
                rows.append(dict(informative_count=count, null_count=null_count, replicate=rep, **result))
        arrays[f"scatters_{count}"] = np.array(scatters)
        if progress:
            progress(f"mass contrast M={count}, null events={null_count}: {config['replicates']} exact-scatter worlds")
    summary = []
    for count in config["informative_counts"]:
        for method in ("score", "trace", "contrast"):
            selected = [r for r in rows if r["informative_count"] == count and r["method"] == method]
            n = len(selected)
            covered = sum(r["lower"] <= mass <= r["upper"] for r in selected)
            summary.append(dict(informative_count=count, null_count=count**4, method=method,
                replicates=n, median_of_medians=float(np.median([r["median"] for r in selected])),
                median_false_band_probability=float(np.median([r["false_band_probability"] for r in selected])),
                medians_on_false_side=sum(r["median"] < split for r in selected), coverage=covered/n,
                mc_low=float(stats.beta.ppf(.025, covered, n-covered+1)) if covered else 0.,
                mc_high=float(stats.beta.ppf(.975, covered+1, n-covered)) if covered < n else 1.))
    for name, records in (("summary", summary), ("population_scatter", population), ("worlds", rows)):
        with (output/(name+".csv")).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    np.savez_compressed(output/"values.npz", **arrays)
    (output/"preset.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    metadata = dict(preset=config, source_sha256=hashes,
        generator="Exact independent noncentral and central Gaussian scatter draws; M informative events, M^4 null events with force SD=M^4.",
        limitations=["Triangular arrays, not the individual prefixes of T63's fixed history.",
            "Both covariance groups are actually sampled; no limiting field is substituted.",
            "The prior is explicitly conditioned on its declared bounds, retaining both roots.",
            "Counts, true signals and the investigator cut enter generation/scoring only.",
            "Small-world coverage is descriptive; quadrature errors are numerical diagnostics.",
            "Neither this rule nor its finite prefactor is selected by study scores."])
    (output/"run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    lines = ["# Mass-contrast anchor sanity", "",
        "Both Gaussian covariance groups are sampled exactly. Prior support is explicitly declared.",
        "Small-world coverage is descriptive; see run.json for the triangular-array scope.", "",
        "| Informative events | Null events | Method | Median of medians | Median false-band probability | False-side medians | Coverage |",
        "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in summary:
        lines.append(f"| {r['informative_count']} | {r['null_count']} | {r['method']} | {r['median_of_medians']:.6g} | "
            f"{r['median_false_band_probability']:.6g} | {r['medians_on_false_side']}/{r['replicates']} | {r['coverage']:.2%} |")
    (output/"summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    if progress:
        progress(str(output))
    return output, dict(summary=summary, **metadata)
