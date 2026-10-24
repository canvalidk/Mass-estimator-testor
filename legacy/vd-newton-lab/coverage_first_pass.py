"""Bounded equation-comparison fixtures; no admission gate or literature claim.

Run with the existing NumPy-only estimator. Results record numerical refinement,
not certified error bounds. All vector fixtures are 2D; the report's mathematical
argument, rather than these examples, addresses the full finite-input domain.
"""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
from estimator import estimate

ROOT = Path(__file__).resolve().parent


def readout(force, acceleration, covariance, direction_order, ratio_order):
    result = estimate(force, acceleration, covariance,
                      direction_order=direction_order, ratio_order=ratio_order)
    return dict(mass=result.mass, lower=result.quantile(.025),
                median=result.quantile(.5), upper=result.quantile(.975))


def main():
    started = time.monotonic()
    factors = np.array([[1., 0., 0., 0.], [.4, 1., 0., 0.],
                        [.35, .2, 1., 0.], [-.1, .25, .3, 1.]])
    scales = np.diag([.5, 1., .75, 1.25])
    lower = scales @ factors
    covariances = {
        "small_isotropic": .25**2 * np.eye(4),
        "large_isotropic": 2.**2 * np.eye(4),
        "anisotropic": np.diag([.5, 1., .75, 1.25])**2,
        "correlated": lower @ lower.T,
    }
    cases = {
        "both_zero": ([0., 0.], [0., 0.]),
        "force_zero": ([0., 0.], [1., 0.]),
        "acceleration_zero": ([2., 0.], [0., 0.]),
        "aligned": ([2., 0.], [1., 0.]),
        "perpendicular": ([0., 2.], [1., 0.]),
        "anti_aligned": ([-2., 0.], [1., 0.]),
    }
    records = []
    for covariance_name, covariance in covariances.items():
        for name, (force, acceleration) in cases.items():
            coarse = readout(force, acceleration, covariance, 128, 1024)
            fine = readout(force, acceleration, covariance, 256, 2048)
            initial_refinement = {key: abs(coarse[key]/fine[key]-1) for key in fine}
            refined = max(initial_refinement.values()) > .001
            if refined:
                # Refine the numerical integral when the first grid comparison
                # resolves a tail quantile poorly; no change to the equation.
                coarse = readout(force, acceleration, covariance, 384, 8192)
                fine = readout(force, acceleration, covariance, 512, 16384)
            assert all(math.isfinite(x) and x > 0 for x in fine.values())
            assert fine["lower"] < fine["median"] < fine["upper"]
            denominator = float(np.dot(acceleration, acceleration))
            ols = float(np.dot(force, acceleration) / denominator) if denominator else None
            records.append(dict(case=name, covariance=covariance_name,
                                force=force, acceleration=acceleration, vd=fine,
                                targeted_refinement=refined,
                                initial_refinement_relative=initial_refinement,
                                ordinary_least_squares=ols,
                                ols_status="defined" if denominator else "undefined",
                                refinement_relative={key: abs(coarse[key]/fine[key]-1)
                                                     for key in fine}))
    paths = []
    # One continuous angle path and a signed component path crossing each zero.
    for degrees in np.linspace(0., 180., 13):
        theta = math.radians(degrees)
        force, acceleration = [2*math.cos(theta), 2*math.sin(theta)], [1., 0.]
        paths.append(dict(path="angle_degrees", parameter=float(degrees),
                          vd=readout(force, acceleration, covariances["correlated"], 256, 2048)))
    for t in [-1., -.1, -.01, 0., .01, .1, 1.]:
        for name, force, acceleration in [
                ("force_through_zero", [t, 0.], [1., 0.]),
                ("acceleration_through_zero", [2., 0.], [t, 0.]),
                ("both_through_zero", [2*t, 0.], [t, 0.])]:
            paths.append(dict(path=name, parameter=t,
                              vd=readout(force, acceleration, covariances["correlated"], 256, 2048)))
    for item in paths:
        assert all(math.isfinite(x) and x > 0 for x in item["vd"].values())

    # Independent oracle: null isotropic law is half-Cauchy, point 1.
    null = next(r["vd"] for r in records if r["case"] == "both_zero"
                and r["covariance"] == "small_isotropic")
    assert abs(null["mass"]-1) < 1e-10
    for key, p in [("lower", .025), ("median", .5), ("upper", .975)]:
        assert abs(null[key]/math.tan(math.pi*p/2)-1) < 1e-9

    # Published Aach Eq18, its native equal-isotropic-error subset only.
    aach = []
    for sign in [1., -1.]:
        x1, x2 = np.array([1., 0.]), np.array([sign, 0.])
        A, B, C = float(x1@x1), float(x1@x2), float(x2@x2)
        statistic = .5*(A+C-math.hypot(A-C, 2*B))
        assert statistic == 0.
        aach.append(dict(relative_sign=sign, D_squared=statistic))

    summary = dict(cases=len(records), path_points=len(paths),
                   targeted_refinements=sum(r["targeted_refinement"] for r in records),
                   max_refinement_relative={key:max(r["refinement_relative"][key] for r in records)
                                            for key in null},
                   null_oracle="passed", finite_positive_checks="passed",
                   runtime_seconds=time.monotonic()-started)
    result = dict(date="2026-09-08", managed_commit="adfe14bbaafb29856f0a30da088495b7c9d483d2",
                  model="Gaussian, flat positive magnitudes, uniform common 2D direction; no gate",
                  covariance_order=["F_x", "F_y", "a_x", "a_y"],
                  covariance_matrices={k:v.tolist() for k,v in covariances.items()},
                  resolutions=dict(coarse=[128,1024], fine=[256,2048],
                                   targeted_coarse=[384,8192], targeted_fine=[512,16384]),
                  uncertainty="Conditional latent mass 2.5%, 50%, 97.5% quantiles; not sampling coverage",
                  source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                 for name in ["estimator.py", "coverage_first_pass.py"]},
                  cases=records, paths=paths, aach_sign_counterexample=aach, summary=summary)
    destination = ROOT / "results" / "coverage_first_pass_2026-09-08.json"
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(json.dumps([r for r in records if r["covariance"]=="small_isotropic"], indent=2))


if __name__ == "__main__":
    main()
