"""Dependence, prediction and information paths for the mass estimator.

The null sample is an exact sampler of the declared flat-magnitude Gaussian
law, not a bootstrap of measured central zeros. Non-null paths use estimator.py.
"""
from pathlib import Path
import json
import math
import numpy as np
from estimator import estimate, isotropic_covariance, isotropic_qmin, profile, compatibility_intervals, gauss

ROOT = Path(__file__).resolve().parent


def null_draws(n, seed=20260907, sigma_force=1., sigma_acceleration=1.):
    if n < 1 or sigma_force <= 0 or sigma_acceleration <= 0:
        raise ValueError("Positive sample count and uncertainty scales required")
    rng = np.random.default_rng(seed)
    f = sigma_force * np.abs(rng.standard_normal(n))
    a = sigma_acceleration * np.abs(rng.standard_normal(n))
    # A zero is a finite-generator event, not an atom of the continuous law.
    # Resample it explicitly instead of introducing an epsilon mass policy.
    for values, scale in ((f, sigma_force), (a, sigma_acceleration)):
        while np.any(values == 0):
            mask = values == 0
            values[mask] = scale * np.abs(rng.standard_normal(mask.sum()))
    return f, a, f / a


def half_cauchy_quantile(p, scale=1.):
    if not 0 < p < 1 or not math.isfinite(scale) or scale <= 0:
        raise ValueError("Interior probability and positive finite scale required")
    return scale * math.tan(math.pi * p / 2)


def positive_interval_product(left, right):
    """Range product, not a 95% joint probability statement."""
    return [left[0] * right[0], left[1] * right[1]]


def detached_force_quantile(p):
    """Independent oracle for M*A', M half-Cauchy(1), A' half-normal(1).

    Survival at z is E[(2/pi)*atan(A'/z)]. The Gaussian integral is on
    [0,12]; its omitted probability is below 4e-33. This is the deliberately
    detached prediction, not the original joint posterior reconstruction.
    """
    if not 0 < p < 1:
        raise ValueError("Interior probability required")
    x, w = gauss(256)
    v, weights = 6 * (x + 1), 6 * w
    weights = weights * math.sqrt(2 / math.pi) * np.exp(-v * v / 2)
    def cdf(value):
        return 1 - float(weights @ ((2 / math.pi) * np.arctan(v / value)))
    lo, hi = 0., 1.
    while cdf(hi) < p:
        hi *= 2
    for _ in range(80):
        mid = (lo + hi) / 2
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def summarize(values):
    return dict(quantiles=np.quantile(values, [.025, .5, .975]).tolist(),
                log_sd=float(np.std(np.log(values))),
                probability_above_10=float(np.mean(values > 10)))


def null_reuse_experiment(n=1000000):
    f, a, mass = null_draws(n)
    # A separate iid draw supplies independent, same-marginal acceleration.
    _, a_independent, _ = null_draws(n, seed=20260908)
    same_pair_force = mass * a
    detached_force = mass * a_independent
    point_force = a  # m_hat = sigma_F/sigma_a = 1 kg in this fixture.
    log_f, log_a, log_m = np.log(f), np.log(a), np.log(mass)
    out = dict(n=n, seed=20260907, independent_seed=20260908,
               target="positive latent magnitudes under measured zero-zero; unit scales",
               original_force=summarize(f), same_pair_force=summarize(same_pair_force),
               incorrectly_detached_force=summarize(detached_force),
               point_substitution_force=summarize(point_force),
               new_exact_unit_acceleration_force=summarize(mass),
               same_pair_max_relative_error=float(np.max(np.abs(same_pair_force / f - 1))),
               magnitude_correlation_original=float(np.corrcoef(f, a)[0, 1]),
               magnitude_correlation_point_substitution=float(np.corrcoef(point_force, a)[0, 1]),
               log_mass_acceleration_correlation=float(np.corrcoef(log_m, log_a)[0, 1]),
               log_force_variance=float(log_f.var()),
               log_force_variance_from_joint=float(log_m.var() + log_a.var() + 2 * np.cov(log_m, log_a, ddof=0)[0, 1]),
               log_force_variance_if_dependence_omitted=float(log_m.var() + log_a.var()),
               predicted_to_original_log_variance_ratio=float(np.var(np.log(detached_force)) / log_f.var()))
    # Deterministic half-normal quantiles from a standard-normal bisection.
    def half_normal_q(p):
        lo, hi = 0., 12.
        for _ in range(80):
            mid = (lo + hi) / 2
            if math.erf(mid / math.sqrt(2)) < p:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    interval_f = [half_normal_q(.025), half_normal_q(.975)]
    interval_m = [half_cauchy_quantile(.025), half_cauchy_quantile(.975)]
    out["interval_products"] = dict(original_force_95=interval_f, mass_95=interval_m,
                                     marginal_95_endpoint_product=positive_interval_product(interval_m, interval_f),
                                     note="Endpoint products do not encode the dependence or define a central 95% interval.")
    out["analytic_oracles"] = dict(log_force_variance=math.pi ** 2 / 8,
                                   log_mass_variance=math.pi ** 2 / 4,
                                   log_mass_acceleration_covariance=-math.pi ** 2 / 8,
                                   log_mass_acceleration_correlation=-1 / math.sqrt(2),
                                   detached_log_force_variance=3 * math.pi ** 2 / 8,
                                   original_force_probability_above_10=math.erfc(10 / math.sqrt(2)),
                                   detached_force_quantiles=[detached_force_quantile(p) for p in (.025, .5, .975)])
    return out


def informative_summary(force, acceleration, sf=1., sa=1.):
    e = estimate(force, acceleration, isotropic_covariance(2, sf, sa), direction_order=256, ratio_order=2048)
    qmin = isotropic_qmin(force, acceleration, sf, sa)
    force_interval = [3 * e.quantile(.025), 3 * e.quantile(.975)]
    return dict(mass=e.mass, interval=[e.quantile(.025), e.quantile(.975)],
                force_interval_at_exact_acceleration_3_if_model_assumed=force_interval,
                force_interval_after_global_gate=None if qmin > 9 else force_interval,
                fit_status="rejected_at_9" if qmin > 9 else "passes_global_gate_at_9",
                qmin=qmin,
                q_at_point=profile(force, acceleration, isotropic_covariance(2, sf, sa), e.mass))


def paths_experiment():
    fixed_noise, scaled_noise, precision = [], [], []
    for factor in (1., .3, .1, .03, .01, .001):
        f, a = [20 * factor, 0], [10 * factor, 0]
        fixed_noise.append(dict(factor=factor, **informative_summary(f, a)))
        scaled_noise.append(dict(factor=factor, **informative_summary(f, a, factor, factor)))
    for sigma in (2., 1., .5, .25):
        precision.append(dict(sigma=sigma, **informative_summary([4, 0], [2, 0], sigma, sigma)))
    return dict(shrink_readings_keep_noise=fixed_noise, shrink_readings_and_noise=scaled_noise,
                keep_readings_improve_precision=precision)


def angle_transitions():
    # Equal strengths p=q=5: point remains 1; interval topology and fit change.
    out = []
    for degrees in (0, 30, 45, 50, 50.2081805004, 55, 90, 120, 180):
        theta = math.radians(degrees)
        f, a = [5 * math.cos(theta), 5 * math.sin(theta)], [5, 0]
        out.append(dict(degrees=degrees, **informative_summary(f, a)))
    moderate = []
    p = math.sqrt(6)
    for degrees in (0, 60, 90, 119, 120, 121, 150, 180):
        theta = math.radians(degrees)
        f, a = [p * math.cos(theta), p * math.sin(theta)], [p, 0]
        intervals = compatibility_intervals(f, a, 1, 1, 9)
        moderate.append(dict(degrees=degrees, **informative_summary(f, a),
                             intervals=[[lo, hi if math.isfinite(hi) else None] for lo, hi in intervals]))
    return dict(strong_equal_channels=out, moderate_equal_channels=moderate,
                strong_gate_angle_degrees=math.degrees(math.acos(1 - 9 / 25)),
                moderate_point_leaves_set_angle_degrees=120.)


def main():
    result = dict(date="2026-09-07", managed_commit="743f82b2aaa3cf09d6a1be45508b6da1c85bbfc6",
                  model="Gaussian independent isotropic errors; flat positive magnitudes; 2D unknown direction for numerical paths",
                  null_reuse=null_reuse_experiment(), paths=paths_experiment(), angle_transitions=angle_transitions())
    (ROOT / "results" / "uncertainty_walkthrough.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
