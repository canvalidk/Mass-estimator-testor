"""B15: fitted-score eta versus the T45/T48 conditional-variance candidate.

    python tools/eta_conditional_study.py --workers 2

Configuration: tools/eta_conditional_study.yaml. Outputs retain every replicate,
the effective configuration and summaries. The grid reproduces O14's Gaussian
random-amplitude world and stream convention; the investigator uses its true
energy only for the comparison target. Both probes declare nu_eff=10, while
the flat reference is T45's unit-Gaussian model. The ours33 C is exact for its
specified integrated-nuisance residual model, not for arbitrary fixed-signal
sampling conditionals under inverse-gamma noise.

All examples have c=s=true mass=1. Isotropic sufficient statistics implement
the same continuous factors and derivatives as the core functions. Whole-line
fits refine every scan maximum, adding own-line seeds for sub-grid peaks.
Coverage is equivalent central-interval membership computed from continuous
angular CDF ranks; no estimate of the median or local-normal approximation is
used. Its quadrature errors are numerical estimates, not certified bounds.
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import hashlib
import math
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
from scipy import integrate, optimize, special, stats
import scipy
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from met.full import isotropic  # noqa: E402
from met.gaussian_information import gaussian_angle_cdf  # noqa: E402

HALF_PI = math.pi / 2


def _eta(curvature, denominator):
    if denominator <= 1e-300 or curvature <= 0:
        return 1.0
    return float(min(1.0, curvature / denominator))


def gram_terms(x, y, nu, *, factors="ours33"):
    """Independent isotropic factor implementation used by this study only."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 2:
        raise ValueError("x and y must have equal (readings, axes) shapes")
    aa, bb, ab = np.sum(x*x, axis=1), np.sum(y*y, axis=1), np.sum(x*y, axis=1)
    exponent = (nu + x.shape[1]) / 2

    def terms(theta, order=0):
        theta = np.atleast_1d(theta)
        sine, cosine = np.sin(2*theta)[None, :], np.cos(2*theta)[None, :]
        residual = np.maximum((aa+bb)[:, None]/2 + (aa-bb)[:, None]*cosine/2-ab[:, None]*sine, 0.)
        if factors == "flat":
            value = -residual/2
        else:
            value = -exponent*np.log1p(residual/nu)
        if order == 0:
            return value
        h1 = (aa-bb)[:, None]*sine+2*ab[:, None]*cosine
        h2 = 2*(aa-bb)[:, None]*cosine-4*ab[:, None]*sine
        if factors == "flat":
            return value, h1/2, h2/2
        denominator = nu+residual
        return value, exponent*h1/denominator, exponent*(h2/denominator+(h1/denominator)**2)

    return terms, (aa, bb, ab)


def fit_isotropic(x, y, nu=10., *, factors="ours33", scan_nodes=512):
    """Return the two powers and their common fit, with scan extrema for CDFs."""
    if factors not in ("flat", "ours33"):
        raise ValueError("factors must be flat or ours33")
    terms, (aa, bb, ab) = gram_terms(x, y, nu, factors=factors)
    k = x.shape[1]
    if factors == "flat":
        theta = 0.5*math.atan2(2*float(ab.sum()), float(bb.sum()-aa.sum()))
        candidates = [theta, ((theta+HALF_PI+HALF_PI) % math.pi)-HALF_PI]
    else:
        step = math.pi/scan_nodes
        grid = -HALF_PI+step*(np.arange(scan_nodes)+0.5)
        narrow = (aa+bb)*step**2 > nu/4
        own = 0.5*np.arctan2(2*ab[narrow], bb[narrow]-aa[narrow])
        grid = np.unique(np.r_[grid, own])
        values = terms(grid).sum(axis=0)
        maxima = np.flatnonzero((values >= np.roll(values, 1)) & (values >= np.roll(values, -1)))
        candidates = []
        for index in maxima:
            low = grid[index-1] if index else grid[-1]-math.pi
            high = grid[index+1] if index+1 < len(grid) else grid[0]+math.pi
            seed = float(grid[index])
            optimum = optimize.minimize_scalar(
                lambda angle: -float(terms(angle).sum()), bounds=(low, high), method="bounded",
                options={"xatol": 1e-14})
            candidates.extend((seed, float(optimum.x)))
            score = lambda angle: float(terms(angle, 2)[1].sum())
            if score(low) > 0 and score(high) < 0:
                candidates.append(float(optimize.brentq(score, low, high, xtol=1e-14, rtol=1e-14)))
        if not candidates:
            candidates = [float(grid[np.argmax(values)])]
        theta = max(candidates, key=lambda angle: float(terms(angle).sum()))
        for _ in range(8):
            _, first, second = terms(theta, 2)
            first, second = float(first.sum()), float(second.sum())
            if second >= 0:
                break
            move = -first/second
            if abs(move) > step:
                break
            theta += move
            if abs(move) < 1e-15:
                break
    theta = (theta+HALF_PI) % math.pi-HALF_PI
    # Direct vector residuals avoid cancellation in the final derivatives.
    off = x*math.cos(theta)-y*math.sin(theta)
    along = x*math.sin(theta)+y*math.cos(theta)
    residual = np.sum(off*off, axis=1)
    h1, h2 = 2*np.sum(off*along, axis=1), 2*(residual-np.sum(along*along, axis=1))
    if factors == "flat":
        first, second = h1/2, h2/2
        variance = float(np.sum(along*along))
    else:
        exponent = (nu+k)/2
        first = exponent*h1/(nu+residual)
        second = exponent*(h2/(nu+residual)+(h1/(nu+residual))**2)
        variance = float((nu+k)/(nu+k+2)*np.sum(along*along))
    curvature, j = float(-second.sum()), float(np.sum(first*first))
    if abs(float(first.sum())) > 2e-7*(1+math.sqrt(j)+math.sqrt(abs(curvature))):
        raise ArithmeticError("whole-line fit failed its stationarity check")
    return dict(eta=_eta(curvature, j), eta_c=_eta(curvature, variance),
                H=curvature, J=j, C=variance, phi=theta,
                candidates=[float((a+HALF_PI) % math.pi-HALF_PI) for a in candidates])


def student_log_ranks(x, y, fit, *, nu, truth_angle=math.pi/4):
    """Log lower/upper raw-law CDF ranks, with separate tail scaling."""
    terms, _ = gram_terms(x, y, nu)
    log_law = lambda theta: float(terms(theta).sum())
    angles = [0., HALF_PI, truth_angle] + [a for a in fit["candidates"] if 0 < a < HALF_PI]

    def piece(low, high):
        choices = [a for a in angles if low <= a <= high]
        best = max(choices, key=log_law)
        baseline = log_law(best)
        curvature = max(1., abs(float(terms(best, 2)[2].sum())))
        width = 1/math.sqrt(curvature)
        points = sorted(set(a for a in choices + [best+j*width for j in (-10, -1, 1, 10)]
                            if low < a < high))
        value, error = integrate.quad(lambda theta: math.exp(log_law(theta)-baseline), low, high,
                                      points=points, epsabs=2e-12, epsrel=2e-10, limit=300)
        if value <= 0 or not math.isfinite(value) or error/value > 1e-6:
            raise ArithmeticError("raw-law tail quadrature did not resolve its relative mass")
        return baseline+math.log(value), error/value

    left, le = piece(0., truth_angle)
    right, re = piece(truth_angle, HALF_PI)
    total = float(np.logaddexp(left, right))
    return left-total, right-total, le+re


def _mean_summary(values):
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.all(np.isfinite(values)):
        return None, None
    return float(values.mean()), float(values.std(ddof=1)/math.sqrt(len(values))) if len(values)>1 else None


def _summarize(rows, *, coverage=False):
    answer = dict(replicates=len(rows), failures=sum(bool(r.get("error")) for r in rows))
    for key in ("eta", "eta_c", "target", "eta_c_minus_eta", "eta_cap", "eta_c_cap"):
        answer[key+"_mean"], answer[key+"_mcse"] = _mean_summary([r.get(key, np.nan) for r in rows])
    if coverage:
        for key in ("covered", "covered_c", "covered_c_minus_covered"):
            answer[key+"_mean"], answer[key+"_mcse"] = _mean_summary([r.get(key, np.nan) for r in rows])
        answer["max_cdf_error"] = max((r.get("cdf_error", 0.) for r in rows), default=0.)
    return answer


def _observed(rng, n, k, snr):
    acceleration = rng.normal(scale=snr, size=(n, k))/math.sqrt(2.)
    noise = rng.normal(size=(n, 2, k))
    return acceleration+noise[:, 0], acceleration+noise[:, 1], 2*float(np.sum(acceleration**2))


def run_grid_cell(job):
    config, index, k, n, snr = job
    rng = np.random.default_rng([config["seed"], index])
    rows = {factor: [] for factor in ("flat", "ours33")}
    for replicate in range(config["series"]):
        x, y, delta = _observed(rng, n, k, snr)
        for factor in rows:
            try:
                fit = fit_isotropic(x, y, config["nu_eff"], factors=factor,
                                    scan_nodes=config["scan_nodes"])
                rows[factor].append(dict(replicate=replicate, eta=fit["eta"], eta_c=fit["eta_c"],
                                         target=delta/(delta+n*k), eta_c_minus_eta=fit["eta_c"]-fit["eta"],
                                         eta_cap=float(fit["eta"] == 1.), eta_c_cap=float(fit["eta_c"] == 1.)))
            except (ValueError, ArithmeticError, FloatingPointError) as error:
                rows[factor].append(dict(replicate=replicate, error=str(error)))
    summary = [dict(kind="grid", factors=factor, k=k, n=n, snr=snr, **_summarize(values))
               for factor, values in rows.items()]
    return summary, {f"grid_{index}_{factor}": values for factor, values in rows.items()}


def run_leverage_cell(job):
    config, index, n = job
    rng = np.random.default_rng([config["seed"], 100+index])
    rows = {factor: [] for factor in ("flat", "ours33")}
    probability = config["probability"]
    z = float(stats.norm.ppf((1+probability)/2))
    for replicate in range(config["leverage_series"]):
        noise = rng.normal(size=(n, 2, 1))
        x, y = noise[:, 0].copy(), noise[:, 1].copy()
        signal = math.sqrt(config["leverage"]["amplitude_squared_per_n"]*n/2)
        x[0, 0] += signal
        y[0, 0] += signal
        for factor in rows:
            try:
                fit = fit_isotropic(x, y, config["nu_eff"], factors=factor,
                                    scan_nodes=config["scan_nodes"])
                if factor == "ours33":
                    lower, upper, error = student_log_ranks(x, y, fit, nu=config["nu_eff"])
                    log_rank = min(lower, upper)
                    threshold = float(stats.norm.logcdf(-z/math.sqrt(fit["eta"])))
                    threshold_c = float(stats.norm.logcdf(-z/math.sqrt(fit["eta_c"])))
                    covered, covered_c = log_rank >= threshold, log_rank >= threshold_c
                    if min(abs(log_rank-threshold), abs(log_rank-threshold_c)) < 5*error:
                        raise ArithmeticError("quadrature uncertainty straddles coverage boundary")
                else:
                    ranks = [gaussian_angle_cdf(math.pi/4, fitted_angle=fit["phi"],
                                                precision=eta*fit["H"])
                             for eta in (fit["eta"], fit["eta_c"])]
                    error = max(e for _, e in ranks)
                    low = (1-probability)/2
                    covered, covered_c = [low <= rank <= 1-low for rank, _ in ranks]
                    if any(min(abs(rank-low), abs(rank-(1-low))) < 5*err for rank, err in ranks):
                        raise ArithmeticError("quadrature uncertainty straddles coverage boundary")
                rows[factor].append(dict(replicate=replicate, eta=fit["eta"], eta_c=fit["eta_c"],
                                         target=config["leverage"]["amplitude_squared_per_n"]
                                         /(config["leverage"]["amplitude_squared_per_n"]+1),
                                         eta_c_minus_eta=fit["eta_c"]-fit["eta"],
                                         eta_cap=float(fit["eta"] == 1.), eta_c_cap=float(fit["eta_c"] == 1.),
                                         covered=float(covered), covered_c=float(covered_c),
                                         covered_c_minus_covered=float(covered_c)-float(covered),
                                         cdf_error=error, phi=fit["phi"]))
            except (ValueError, ArithmeticError, FloatingPointError) as error:
                rows[factor].append(dict(replicate=replicate, error=str(error)))
    summary = [dict(kind="leverage", factors=factor, k=1, n=n, snr=None,
                    **_summarize(values, coverage=True)) for factor, values in rows.items()]
    return summary, {f"leverage_{index}_{factor}": values for factor, values in rows.items()}


def core_checks(config):
    """Check sufficient-statistic fits against the general core implementation."""
    from met.eta_conditional import eta_c
    from met.full import bartlett_eta as gaussian_bartlett
    from met.ours31 import _per_reading
    from met.ours32 import bartlett_eta
    from met.ours33 import _l_terms, ours33
    from met.analyst import _cdf

    checks, rank_checks = [], []
    rng = np.random.default_rng([config["seed"], 999])
    examples = [(k, *_observed(rng, n, k, snr)[:2])
                for k, n, snr in ((1, 5, .5), (3, 20, 2.), (1, 100, 5.))]
    for n in (50, 1000):
        x, y = rng.normal(size=(2, n, 1))
        x[0, 0] += math.sqrt(2*n)
        y[0, 0] += math.sqrt(2*n)
        examples.append((1, x, y))
    for k, x, y in examples:
        n = len(x)
        readings = [isotropic(f, a, 1., 1.) for f, a in zip(x, y)]
        for factor in ("flat", "ours33"):
            fit = fit_isotropic(x, y, config["nu_eff"], factors=factor, scan_nodes=config["scan_nodes"])
            exact = eta_c(readings, 1., config["nu_eff"], factors=factor)
            values = np.array([fit["eta_c"], fit["H"], fit["C"]])
            expected = np.array(exact[:3])
            if not np.allclose(values, expected, rtol=3e-8, atol=2e-8):
                raise ArithmeticError(f"study/core fit mismatch: {factor}, k={k}, n={n}")
            fitted = (gaussian_bartlett(readings, 1., [1.]*n) if factor == "flat" else
                      bartlett_eta(_per_reading(readings, 1., [config["nu_eff"]]*n), 1., terms=_l_terms))
            if not math.isclose(fit["eta"], fitted[0], rel_tol=3e-8, abs_tol=2e-8):
                raise ArithmeticError(f"study/core fitted eta mismatch: {factor}, k={k}, n={n}")
            checks.append(dict(factors=factor, k=k, n=n, max_abs_error=float(np.max(abs(values-expected))),
                               fitted_eta_abs_error=abs(fit["eta"]-fitted[0])))
            if factor == "ours33" and n <= 20:
                lower, upper, error = student_log_ranks(x, y, fit, nu=config["nu_eff"])
                native = ours33(readings, 1., config["nu_eff"], calibrate=False)
                h = native["lm"][1]-native["lm"][0]
                cdf = _cdf(native["law"][None, :], np.array([h]))[0]
                rank = float(np.interp(0., native["lm"], cdf))
                difference = abs(math.exp(lower)-rank)
                if difference > 3e-4:
                    raise ArithmeticError("continuous rank disagrees with the native mass-law grid")
                rank_checks.append(dict(k=k, n=n, native_grid_abs_error=difference,
                                        quadrature_relative_error=error,
                                        log_lower_rank=lower, log_upper_rank=upper))
    angle, radius, d = math.atan(2.), 1e3, 3.
    u, v = np.array([math.sin(angle), math.cos(angle)]), np.array([math.cos(angle), -math.sin(angle)])
    points = [radius*u+d*v, radius*u-d*v]
    readings = [isotropic([p[0]], [p[1]], 1., 1.) for p in points]
    conditional = eta_c(readings, 1., config["nu_eff"], factors="ours33")
    fitted = bartlett_eta(_per_reading(readings, 1., [config["nu_eff"]]*2), 1., terms=_l_terms)
    padded = eta_c(readings+[isotropic([0.], [0.], 1., 1.)]*4, 1., config["nu_eff"], factors="ours33")
    if padded != conditional:
        raise ArithmeticError("exact-zero padding changed eta_C or its fit")
    return dict(core_agreement=checks, rank_agreement=rank_checks, opposed_pair=dict(radius=radius, error_size=d,
                eta=fitted[0], eta_c=conditional[0], H=conditional[1], C=conditional[2],
                phi=conditional[3], zero_padding_bitwise_equal=True))


def _write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _provenance(preset):
    """Capture source bytes at invocation, including uncommitted study code."""
    root = Path(__file__).resolve().parent.parent
    paths = [root/"tools/eta_conditional_study.py", root/"met/eta_conditional.py",
             root/"met/ours33.py", preset.resolve()]
    hashes = {str(path.relative_to(root)) if path.is_relative_to(root) else str(path):
              hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, cwd=root).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], text=True, cwd=root))
    except (OSError, subprocess.SubprocessError):
        commit, dirty = None, None
    return dict(git_head=commit, git_dirty=dirty, run_with_uncommitted_changes=dirty,
                source_provenance=dict(captured_utc=datetime.now(timezone.utc).isoformat(),
                    timing="Captured at invocation before core checks and simulations.", source_sha256=hashes))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", type=Path, default=Path("tools/eta_conditional_study.yaml"))
    parser.add_argument("--series", type=int)
    parser.add_argument("--leverage-series", type=int)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--part", choices=("all", "grid", "leverage", "checks"), default="all")
    parser.add_argument("--output", type=Path, default=Path("results/b15_eta_conditional"))
    args = parser.parse_args(argv)
    config = yaml.safe_load(args.preset.read_text(encoding="utf-8"))
    for key, value in (("series", args.series), ("leverage_series", args.leverage_series)):
        if value is not None:
            config[key] = value
    if args.workers < 1 or any(config[key] < 1 for key in ("series", "leverage_series")):
        parser.error("workers and replicate counts must be positive")
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    provenance = _provenance(args.preset)
    checks = core_checks(config)
    jobs = []
    if args.part in ("all", "grid"):
        jobs += [("grid", (config, i, k, n, snr)) for i, (k, n, snr) in enumerate(
            (k, n, snr) for k in config["grid"]["k"] for n in config["grid"]["n"] for snr in config["grid"]["snr"])]
    if args.part in ("all", "leverage"):
        jobs += [("leverage", (config, i, n)) for i, n in enumerate(config["leverage"]["n"])]
    summaries, details = [], {}
    metadata = dict(config=config, part=args.part, workers=args.workers, **provenance,
                    started_utc=started, python=platform.python_version(),
                    numpy=np.__version__, scipy=scipy.__version__, checks=checks)
    (output/"run.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    if args.workers == 1:
        iterator = map(_dispatch, jobs)
        executor = None
    else:
        executor = ProcessPoolExecutor(args.workers)
        iterator = executor.map(_dispatch, jobs)
    try:
        for summary, values in iterator:
            summaries.extend(summary)
            details.update(values)
            _write_csv(output/"summary.csv", summaries)
            (output/"replicates.json").write_text(json.dumps(details, allow_nan=False)+"\n", encoding="utf-8")
            for row in summary:
                print(json.dumps(row, allow_nan=False), flush=True)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["failures"] = sum(row["failures"] for row in summaries)
    (output/"run.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(checks, allow_nan=False), flush=True)
    return 1 if metadata["failures"] else 0


def _dispatch(job):
    kind, arguments = job
    return run_grid_cell(arguments) if kind == "grid" else run_leverage_cell(arguments)


if __name__ == "__main__":
    raise SystemExit(main())
