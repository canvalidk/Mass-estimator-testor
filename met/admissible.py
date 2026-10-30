"""B12 study: strong-signal, whole-line cut-precision families.

This module is not imported by the production estimator. Complete monotonicity,
existence of a proper noise-size prior, and its two normalizations are separate
checks. At nu=0, Q is physically zero almost surely; statements about the whole
Q curve concern its exponential-tilt continuation off that support.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc, gammaln


FORMS = ("gaussian", "ours29", "ours30", "smooth", "floor_mle",
         "floor_map", "floor_exact", "invgamma")


def _finite(name, value):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite number") from None
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _parameters(form, nu, alpha, exponent, normalization):
    if form not in FORMS:
        raise ValueError(f"unknown form {form!r}")
    nu = _finite("nu", nu)
    alpha = _finite("alpha", alpha)
    exponent = _finite("exponent", exponent)
    if nu < 0 or alpha <= 0:
        raise ValueError("nu must be nonnegative and alpha positive")
    if normalization not in ("covariance", "precision"):
        raise ValueError("normalization must be covariance or precision")
    if form == "invgamma" and normalization == "covariance" and alpha <= 1:
        raise ValueError("covariance-normalized inverse gamma requires alpha > 1")
    return nu, alpha, exponent


def _nonnegative_q(q):
    q = _finite("q", q)
    if q < 0:
        raise ValueError("q must be nonnegative")
    return q


def _log1p_ratio(q, denominator):
    if q <= denominator:
        return math.log1p(q / denominator)
    return math.log(q) - math.log(denominator) + math.log1p(denominator / q)


def _kink_numerator(form, nu, exponent):
    if form in ("ours29", "floor_mle"):
        return nu
    if form == "ours30":
        return nu + 0.5
    return nu + 2 * exponent  # Posterior-density-kernel mode in v, not in 1/v.


def _floor_pair(q, r):
    """Return normalized L and g for mu(du) proportional to u^(r-1), 0<u<=1."""
    if r <= 0:
        raise ValueError("floor cut integral is infinite when exponent - 1 + nu/2 <= 0")
    if q == 0:
        return 1.0, r / (r + 1)
    t = q / 2
    if t <= 1:
        # Algebraic-weight quadrature handles the integrable singularity r<1.
        integral = lambda s: quad(lambda u: math.exp(-t * u), 0, 1,
                                 weight="alg", wvar=(s - 1, 0),
                                 epsabs=1e-12, epsrel=1e-12)[0]
        denominator, numerator = integral(r), integral(r + 1)
        return r * denominator, numerator / denominator

    def log_integral(s):
        probability = float(gammainc(s, t))
        if probability > 0:
            return float(gammaln(s)) + math.log(probability) - s * math.log(t)
        # When s is much larger than t, incomplete gamma can underflow.
        # Subtract the maximum log integrand before integrating instead.
        mode = min(1.0, (s - 1) / t)
        peak = (s - 1) * math.log(mode) - t * mode
        def scaled(u):
            return 0.0 if u == 0 else math.exp((s - 1) * math.log(u) - t * u - peak)
        value = quad(scaled, 0, 1, epsabs=1e-12, epsrel=1e-12, limit=250)[0]
        if value <= 0:
            raise ValueError("floor integral could not be resolved numerically")
        return peak + math.log(value)

    log_denominator = log_integral(r)
    log_numerator = log_integral(r + 1)
    return math.exp(math.log(r) + log_denominator), math.exp(log_numerator - log_denominator)


def weight(form, q, nu, *, alpha=3.0, exponent=3.0, normalization="covariance"):
    """Return a named g(Q); floor_map is reciprocal posterior-kernel mode in v.

    The nu=0 floor MLE has no unique value even at Q=0, and raises. The
    ours29/smooth nu=0 cases reproduce the estimator's explicit g=1 guard.
    """
    nu, alpha, exponent = _parameters(form, nu, alpha, exponent, normalization)
    q = _nonnegative_q(q)
    if form == "floor_mle" and nu == 0:
        raise ValueError("nu=0 floor MLE is not uniquely defined")
    if form == "gaussian" or (form in ("ours29", "smooth") and nu == 0):
        return 1.0
    if form == "smooth":
        return nu / (nu + q)
    if form == "invgamma":
        rate = alpha - 1 if normalization == "covariance" else alpha
        return (alpha + nu / 2) / (rate + q / 2)
    if form == "floor_exact":
        return _floor_pair(q, exponent - 1 + nu / 2)[1]
    numerator = _kink_numerator(form, nu, exponent)
    if numerator <= 0:
        raise ValueError("this plug-in has no finite positive variance mode")
    return 1.0 if q <= numerator else numerator / q


def laplace(form, q, nu, *, alpha=3.0, exponent=3.0, normalization="covariance"):
    """Return the L(0)=1 solution of -2 L'/L=g, where that curve is defined."""
    nu, alpha, exponent = _parameters(form, nu, alpha, exponent, normalization)
    q = _nonnegative_q(q)
    if form == "floor_mle" and nu == 0:
        raise ValueError("nu=0 floor MLE has no unique full curve")
    if form == "gaussian" or (form in ("ours29", "smooth") and nu == 0):
        return math.exp(-q / 2)
    if form == "smooth":
        return math.exp(-nu / 2 * _log1p_ratio(q, nu))
    if form == "invgamma":
        rate = alpha - 1 if normalization == "covariance" else alpha
        return math.exp(-(alpha + nu / 2) * _log1p_ratio(q / 2, rate))
    if form == "floor_exact":
        return _floor_pair(q, exponent - 1 + nu / 2)[0]
    numerator = _kink_numerator(form, nu, exponent)
    if numerator <= 0:
        raise ValueError("this plug-in has no finite positive variance mode")
    if q <= numerator:
        return math.exp(-q / 2)
    return math.exp(-numerator / 2 * (1 + math.log(q) - math.log(numerator)))


def check_form(form, nu, *, alpha=3.0, exponent=3.0, normalization="covariance"):
    """Analytic verdicts; nonexistent P has moments None, rather than infinity.

    covariance and precision test E[v]=1 and E[1/v]=1 respectively. The
    normalization keyword selects the inverse-gamma parameterization only.
    """
    nu, alpha, exponent = _parameters(form, nu, alpha, exponent, normalization)
    out = dict(form=form, nu=nu, cm="fail", proper="fail", covariance="fail",
               precision="fail", mean_v=None, mean_precision=None, reason="", analytic=True)
    out["scope"] = ("nu=0: only Q=0 is observed; positive-Q results are a formal exponential-tilt continuation, not an identifiable physical cut."
                    if nu == 0 else "Strong-signal, whole-line chi-square cut.")
    if form == "floor_mle" and nu == 0:
        for key in ("cm", "proper", "covariance", "precision"):
            out[key] = "unknown"
        out["reason"] = "Q=0 almost surely; every v>=1 is an MLE there, so no unique full curve exists."
        return out
    if form == "gaussian" or (form in ("ours29", "smooth") and nu == 0):
        out.update(cm="pass", proper="pass", covariance="pass", precision="pass",
                   mean_v=1.0, mean_precision=1.0,
                   reason="L=exp(-Q/2) is the Laplace transform of delta_u=1; P=delta_v=1. This fixed Gaussian fails T30.")
        if nu == 0:
            out["reason"] += " The nu=0 guard uses g=1; only Q=0 is physically observed."
        return out
    if form == "smooth":
        out.update(cm="pass", reason="L=(1+Q/nu)^(-nu/2) is a gamma Laplace transform, but inverse tilting gives P_u proportional to u^-1 exp(-nu*u/2), which is not normalizable at zero. No proper P or its moments exist.")
        return out
    if form == "invgamma":
        rate = alpha - 1 if normalization == "covariance" else alpha
        out.update(cm="pass", proper="pass",
                   covariance="pass" if normalization == "covariance" else "fail",
                   precision="pass" if normalization == "precision" else "fail",
                   mean_v=rate / (alpha - 1) if alpha > 1 else math.inf,
                   mean_precision=alpha / rate,
                   reason="Gamma precision prior: L=(1+Q/(2*rate))^(-(alpha+nu/2)); inverse tilting recovers the proper inverse-gamma P. Both normalization verdicts use its prior moments.")
    elif form == "floor_exact":
        r = exponent - 1 + nu / 2
        if r <= 0:
            out["reason"] = "The floor cut integral diverges at u=0: no finite L exists (exponent-1+nu/2<=0)."
        elif exponent <= 1:
            out.update(cm="pass", reason="L is a positive Laplace integral on 0<u<=1, but inverse tilting gives improper P_v proportional to v^-exponent on v>=1 (exponent<=1).")
        else:
            out.update(cm="pass", proper="pass",
                       mean_v=(exponent - 1) / (exponent - 2) if exponent > 2 else math.inf,
                       mean_precision=(exponent - 1) / exponent,
                       reason="Exact floor cut: L proportional to integral_0^1 u^(exponent-2+nu/2) exp(-Q*u/2) du. Proper P requires exponent>1; its E[v]>1 (possibly infinite) and E[1/v]<1 fail both mean-one normalizations.")
    else:
        numerator = _kink_numerator(form, nu, exponent)
        if numerator <= 0:
            out["reason"] = "The plug-in has no finite positive variance mode, so no full g or finite L is defined."
        else:
            out["reason"] = "The capped reciprocal has L=C1 but not C2 at Q=" + str(numerator) + "; every finite Laplace transform is C-infinity for Q>0. Hence Bernstein fails, and no proper noise-size P exists."
            if form == "floor_map":
                out["reason"] += " This plug-in is reciprocal posterior-density-kernel mode in v (MAP when proper), with numerator nu+2*exponent."
    if nu == 0:
        out["reason"] += " This full-curve verdict is an off-support continuation: nu=0 observes only Q=0."
    return out


def numerical_check(g, nu, *, q_min=0.05, q_max=20.0, nodes=129,
                    max_order=4, tolerance=1e-8):
    """Finite-grid diagnostic for an arbitrary callable, never an analytic verdict.

    Reconstruct L anchored at Q=1. On the uniform grid, CM requires every
    alternating forward difference to be nonnegative. Finite orders/regions,
    quadrature error and rounding prevent either sampled result being a proof.
    This does not attempt Laplace inversion or infer a mixing probability.
    """
    nu = _finite("nu", nu)
    q_min, q_max = _finite("q_min", q_min), _finite("q_max", q_max)
    tolerance = _finite("tolerance", tolerance)
    if nu < 0 or not 0 < q_min < q_max or tolerance <= 0:
        raise ValueError("require nu>=0, 0<q_min<q_max and tolerance>0")
    if (isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer))
            or isinstance(max_order, bool) or not isinstance(max_order, (int, np.integer))
            or not 1 <= max_order < nodes):
        raise ValueError("require integer nodes>max_order>=1")
    grid = np.linspace(q_min, q_max, nodes)
    integrals = [quad(lambda q: _finite("g(Q)", g(q)), 1.0, float(q),
                      epsabs=1e-10, epsrel=1e-10, limit=200) for q in grid]
    log_l = -0.5 * np.array([value for value, error in integrals])
    if not np.isfinite(log_l).all():
        raise ValueError("the anchored integral of g must be finite on the grid")
    # A common positive factor does not affect any sign; this prevents overflow.
    differences = np.exp(log_l - log_l.max())
    minima = [float(differences.min())]
    for order in range(1, max_order + 1):
        differences = np.diff(differences)
        minima.append(float(((-1) ** order * differences).min()))
    sampled = "pass" if min(minima) >= -tolerance else "fail"
    return dict(form="callable", nu=nu, cm="unknown", proper="unknown",
                covariance="unknown", precision="unknown", mean_v=None,
                mean_precision=None, analytic=False, sampled_cm=sampled,
                reason="A finite alternating-difference diagnostic cannot certify global complete monotonicity, recovery of P, or its moments.",
                diagnostics=dict(q_min=q_min, q_max=q_max, nodes=int(nodes),
                                 max_order=int(max_order), tolerance=tolerance,
                                 signed_difference_minima=minima,
                                 quadrature_error_max=max(error for value, error in integrals),
                                 anchor=1.0))


def linear_gaussian(x, c):
    """B04's exact linear limit: residual statistic S and Bartlett H/J."""
    x, c = np.asarray(x, dtype=float), np.asarray(c, dtype=float)
    if (x.ndim != 1 or c.shape != x.shape or x.size == 0
            or not np.isfinite(x).all() or not np.isfinite(c).all() or np.any(c <= 0)):
        raise ValueError("x and c must be equal nonempty finite vectors with c>0")
    h = float(c.sum())
    if not math.isfinite(h):
        raise ValueError("the total precision must be finite")
    theta = float(np.dot(c / h, x))
    residual = x - theta
    s = float(np.dot(c, residual ** 2))
    j = float(np.dot(c * residual, c * residual))
    if not all(math.isfinite(value) for value in (theta, s, j)):
        raise ValueError("the residual statistics must be finite")
    return dict(theta=theta, H=h, J=j, S=s, nu=int(x.size - 1),
                eta=1.0 if j == 0 else min(1.0, h / j))


def check_eta(n):
    """Analytic obstruction for N>=2, including the degenerate N=1 control."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")
    if n == 1:
        out = check_form("gaussian", 0)
        out.update(form="eta", reason="N=1 has S=J=0 and eta=1, compatible with delta_v=1; the sole observed cut value does not identify a full mixing law.")
        return out
    return dict(form="eta", nu=float(n - 1), cm="fail", proper="fail",
                covariance="fail", precision="fail", mean_v=None,
                mean_precision=None, analytic=True,
                reason="Equal precisions give eta(S)=min(1,N/S), whose positive plateau and nonconstant tail cannot be a posterior precision. Heterogeneous precisions additionally make J depend on residual direction beyond sufficient S.",
                scope="Linear-Gaussian shared-noise-size cut, not a statement about weak signal.")


def run_preset(path):
    """Evaluate the isolated deterministic B12 grid; no production runner imports it."""
    from pathlib import Path
    import yaml

    with Path(path).open(encoding="utf-8") as source:
        preset = yaml.safe_load(source)
    if not isinstance(preset, dict) or preset.get("kind") != "admissible_tempering":
        raise ValueError("expected an admissible_tempering preset")
    seed = preset.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("the deterministic preset must record an integer seed")
    dimensions = preset["dimensions"]
    if (not isinstance(dimensions, list) or not dimensions
            or any(isinstance(k, bool) or not isinstance(k, int) or k not in (1, 2, 3)
                   for k in dimensions)):
        raise ValueError("dimensions must list k in 1,2,3")
    rows = []
    for k in dimensions:
        nu = k - 1
        cases = [(form, {}) for form in ("gaussian", "ours29", "ours30", "smooth", "floor_mle")]
        for exponent in preset["floor_exponents"]:
            cases.extend((form, dict(exponent=exponent)) for form in ("floor_map", "floor_exact"))
        for alpha in preset["inverse_gamma_shapes"]:
            cases.extend(("invgamma", dict(alpha=alpha, normalization=normalization))
                         for normalization in ("covariance", "precision"))
        for form, kwargs in cases:
            result = check_form(form, nu, **kwargs)
            result.update(k=k, parameters=kwargs)
            rows.append(result)
    return dict(name=preset.get("name", "admissible-tempering"), seed=seed,
                sampling=False, rows=rows, eta=[dict(n=n, **check_eta(n)) for n in (1, 2, 5)])


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("preset", help="path to the isolated deterministic study preset")
    args = parser.parse_args(argv)
    result = run_preset(args.preset)
    print(f"{result['name']}; seed={result['seed']}; deterministic (no sampling)")
    print("k=1: Q=0 almost surely; positive-Q classifications are formal continuations.")
    print("| k | form and parameters | CM | proper P | E[v] | E[1/v] | covariance | precision |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    def moment(value):
        return "undefined" if value is None else ("infinity" if math.isinf(value) else f"{value:.8g}")
    for row in result["rows"]:
        label = row["form"] + " " + ", ".join(f"{key}={value}" for key, value in row["parameters"].items())
        print(f"| {row['k']} | {label.strip()} | {row['cm']} | {row['proper']} | "
              f"{moment(row['mean_v'])} | {moment(row['mean_precision'])} | "
              f"{row['covariance']} | {row['precision']} |")
    for row in result["eta"]:
        print(f"eta, N={row['n']}: {row['proper']}; {row['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
