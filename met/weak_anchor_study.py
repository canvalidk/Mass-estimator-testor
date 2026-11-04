"""T58/T59 fixed weak-history study, sampling exact sufficient scatters.

The deterministic standardized history has r_i^2=i^beta-(i-1)^beta,
so its prefix energy is N^beta. Orthogonally rotating this fixed vector
shows scatter = outer(sqrt(energy)*u+z0)+Wishart(N-1,I). Bartlett sampling
avoids constructing huge observation arrays; it is not a signal population
or a nuisance prior. Only scatter/covariance/prior enter law evaluation.
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


def scatter_from_components(energy, direction, gaussian_pair, chi_squares, off_diagonal):
    direction = np.asarray(direction, dtype=float)
    gaussian_pair = np.asarray(gaussian_pair, dtype=float)
    chi_squares = np.asarray(chi_squares, dtype=float)
    if (direction.shape != (2,) or gaussian_pair.shape != (2,) or chi_squares.shape != (2,)
            or np.any(~np.isfinite(direction)) or np.any(~np.isfinite(gaussian_pair))
            or np.any(~np.isfinite(chi_squares)) or np.any(chi_squares < 0)
            or not np.isfinite(energy) or energy < 0 or not np.isfinite(off_diagonal)
            or not np.isclose(direction@direction, 1., rtol=1e-12, atol=0.)):
        raise ValueError("finite Bartlett components, nonnegative energy and unit direction required")
    shifted = np.sqrt(energy)*direction+gaussian_pair
    bartlett = np.array([[np.sqrt(chi_squares[0]), 0.],
                         [off_diagonal, np.sqrt(chi_squares[1])]])
    return np.outer(shifted, shifted)+bartlett@bartlett.T


def sample_scatter(rng, count, energy, mass, covariance):
    chol = np.linalg.cholesky(covariance)
    direction = np.linalg.solve(chol, np.array([mass, 1.]))
    direction /= np.linalg.norm(direction)
    white = scatter_from_components(energy, direction, rng.normal(size=2),
        [rng.chisquare(count-1), rng.chisquare(count-2)], rng.normal())
    return chol@white@chol.T


def validate(config):
    required = {"kind", "name", "question", "seed", "replicates", "beta", "counts",
        "probability", "prior_center", "profile_anchor_exponent", "cases", "generator",
        "epsabs", "epsrel", "limit"}
    if not isinstance(config, dict) or set(config) != required or config["kind"] != "weak_anchor_sanity":
        raise ValueError("weak_anchor_sanity requires all explicit settings")
    if config["name"] != "weak-anchor-sanity" or config["generator"] != "exact_noncentral_scatter":
        raise ValueError("this study uses the declared name and exact scatter generator")
    for key, minimum in (("seed", 0), ("replicates", 2), ("limit", 2)):
        if type(config[key]) is not int or config[key] < minimum:
            raise ValueError("invalid integer control")
    if not .5 < config["beta"] < 1 or not 0 < config["profile_anchor_exponent"] < 1 or not 0 < config["probability"] < 1:
        raise ValueError("this benchmark requires an identifying weak history and valid exponents")
    for key in ("prior_center", "epsabs", "epsrel"):
        if not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError("positive finite prior/integration controls required")
    if (not config["counts"] or any(type(n) is not int or n < 3 for n in config["counts"])
            or len(set(config["counts"])) != len(config["counts"])):
        raise ValueError("distinct integer counts >=3 required")
    if not config["cases"]:
        raise ValueError("at least one declared case required")
    names = []
    for case in config["cases"]:
        if (set(case) != {"name", "mass", "correlation"} or not isinstance(case["name"], str)
                or not case["name"] or not np.isfinite(case["mass"]) or case["mass"] <= 0
                or not -1 < case["correlation"] < 1):
            raise ValueError("cases require names, positive masses and valid correlations")
        names.append(case["name"])
    if len(set(names)) != len(names):
        raise ValueError("case names must be unique")
    return dict(kind="weak_anchor_sanity", raw=config,
                cells=[(case["name"], n) for case in config["cases"] for n in config["counts"]])


def run(path, out_root="results", *, replicates=None, progress=print):
    config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if replicates is not None:
        config["replicates"] = replicates
    config = validate(config)["raw"]
    output = Path(out_root)/config["name"]/datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")
    output.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for name in ("weak_anchor_study.py", "score_law_obstruction.py"):
        source = Path(__file__).with_name(name).read_bytes()
        (output/name).write_bytes(source)
        hashes[name] = hashlib.sha256(source).hexdigest()
    rng = np.random.default_rng(config["seed"])
    ranks = [(1-config["probability"])/2, .5, (1+config["probability"])/2]
    rows, population, arrays = [], [], {}

    def evaluate(scatter, covariance):
        results = []
        for schedule in ("profile_oscillation", "spectral"):
            q = ScalarScoreQuadrature(scatter, covariance, prior_center=config["prior_center"],
                anchor_exponent=config["profile_anchor_exponent"], anchor_schedule=schedule,
                epsabs=config["epsabs"], epsrel=config["epsrel"], limit=config["limit"])
            lower, median, upper = [q.quantile(rank, "anchored_score") for rank in ranks]
            results.append(dict(schedule=schedule, lower=lower, median=median, upper=upper,
                anchor_weight=q.anchor_weight, gap=q.gap, largest_eigenvalue=q.largest_eigenvalue,
                local_anchor_ratio=q.anchor_weight*q.largest_eigenvalue/q.gap,
                separation_strength=q.anchor_weight*q.gap,
                normalizer=q.total("anchored_score")[0],
                reported_integral_error=q.total("anchored_score")[1]))
        return results

    for case in config["cases"]:
        mass, rho = case["mass"], case["correlation"]
        covariance = np.array([[1., rho], [rho, 1.]])
        c = np.array([mass, 1.])
        normalized = c/np.sqrt(c@np.linalg.solve(covariance, c))
        for count in config["counts"]:
            energy = count**config["beta"]
            expected = count*covariance+energy*np.outer(normalized, normalized)
            for result in evaluate(expected, covariance):
                population.append(dict(case=case["name"], count=count, energy=energy, **result))
            samples = []
            for rep in range(config["replicates"]):
                scatter = sample_scatter(rng, count, energy, mass, covariance)
                samples.append(scatter)
                for result in evaluate(scatter, covariance):
                    rows.append(dict(case=case["name"], count=count, replicate=rep, **result))
            arrays[f"scatter_{case['name']}_{count}"] = np.array(samples)
            if progress:
                progress(f"weak anchor {case['name']} N={count}: {config['replicates']} scatter worlds complete")
    summary = []
    for case in config["cases"]:
        for count in config["counts"]:
            for schedule in ("profile_oscillation", "spectral"):
                selected = [r for r in rows if r["case"] == case["name"] and r["count"] == count and r["schedule"] == schedule]
                n = len(selected)
                covered = sum(r["lower"] <= case["mass"] <= r["upper"] for r in selected)
                summary.append(dict(case=case["name"], count=count, schedule=schedule,
                    replicates=n, coverage=covered/n,
                    mc_low=float(stats.beta.ppf(.025, covered, n-covered+1)) if covered else 0.,
                    mc_high=float(stats.beta.ppf(.975, covered+1, n-covered)) if covered < n else 1.,
                    median_of_medians=float(np.median([r["median"] for r in selected])),
                    median_interval_width=float(np.median([r["upper"]-r["lower"] for r in selected])),
                    median_anchor_ratio=float(np.median([r["local_anchor_ratio"] for r in selected]))))
    for name, records in (("summary", summary), ("population_scatter", population), ("worlds", rows)):
        with (output/(name+".csv")).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    np.savez_compressed(output/"values.npz", **arrays)
    (output/"preset.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    metadata = dict(preset=config, source_sha256=hashes,
        fixed_history="r_i^2=i^beta-(i-1)^beta; normalized physical acceleration r_i/sqrt(C'cov^-1 C)",
        limitations=["Scatters are exact Gaussian sufficient-statistic draws, not individually stored readings.",
            "No truth/energy/count is supplied to the law; those enter generation and scoring only.",
            "Different counts use independent prefix experiments, not a shared nested noise realization.",
            "Adaptive integration errors are ordinary numerical diagnostics, not certificates.",
            "32-world default coverage is descriptive; theorem rate results are separate.",
            "Neither anchor family nor its finite constants are adopted by score."])
    (output/"run.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    lines = ["# Weak fixed-history anchor sanity", "", "Nominal probability is declared in preset.yaml. See run.json for generator and scope.", "",
        "| Case | N | Schedule | Coverage | 95% Monte Carlo interval | Median interval width | Median anchor/score precision ratio |",
        "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in summary:
        lines.append(f"| {r['case']} | {r['count']} | {r['schedule']} | {r['coverage']:.2%} | "
            f"{r['mc_low']:.2%}–{r['mc_high']:.2%} | {r['median_interval_width']:.6g} | {r['median_anchor_ratio']:.6g} |")
    (output/"summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    if progress:
        progress(str(output))
    return output, dict(summary=summary, **metadata)
