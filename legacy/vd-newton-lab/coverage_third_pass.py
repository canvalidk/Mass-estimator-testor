"""Bounded checks of an explicitly instantiated Bayesian EIV competitor.

The Gaussian latent-vector and lognormal positive-scale priors and posterior
mean are this review's choices within published general frameworks; these
numerical results are not attributed to the original papers. A separate 2D
check verifies the exact change of coordinates for the target prior measure.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from estimator import estimate

ROOT = Path(__file__).resolve().parent


def alternative(force, acceleration, covariance, *, count=12001, extent=12.):
    """v~N(0,I), log(m)~N(0,1), y|(m,v)~N([mv,v], covariance).

    All quantities are in fixed numerical units. The prior is held fixed
    between cases, including changes to the measurement error covariance.
    Analytically integrate v; trapezoid-integrate t=log(m).
    """
    f, a = np.asarray(force), np.asarray(acceleration)
    d, y = len(f), np.r_[f, a]
    t = np.linspace(-extent, extent, count)
    m = np.exp(t)
    b = np.concatenate((m[:, None, None]*np.eye(d)[None, :, :],
                        np.broadcast_to(np.eye(d), (count, d, d))), axis=1)
    c = covariance[None, :, :] + b @ np.swapaxes(b, 1, 2)
    signs, logdet = np.linalg.slogdet(c)
    assert np.all(signs > 0)
    solve = np.linalg.solve(c, np.broadcast_to(y, (count, 2*d))[..., None])[..., 0]
    logweight = -.5*(t*t + logdet + np.einsum("i,ni->n", y, solve))
    peak = logweight.max()
    density = np.exp(logweight-peak)
    norm = np.trapezoid(density, t)
    density /= norm
    cdf = np.r_[0., np.cumsum((density[1:]+density[:-1])*np.diff(t)/2)]
    quantiles = np.exp(np.interp([.025, .5, .975], cdf, t))
    mean = np.trapezoid(m*density, t)
    return {"mean": float(mean), "quantiles": quantiles.tolist(),
            "log_normalizer_without_constants": float(peak+math.log(norm)),
            "grid_edge_density": [float(density[0]), float(density[-1])]}


def flat_cartesian_angle_density(force, acceleration, covariance, theta, scale):
    """2D target equals flat dm d^2v: integrate v as a full Gaussian.

    With m=scale*tan(theta), rewrite B_m = D_theta/cos(theta).
    For d=2 the cos^2 from |B'WB|^-1/2 cancels dm/dtheta.
    This stable expression includes both endpoints without any mass cutoff.
    """
    y = np.r_[force, acceleration]
    w = np.linalg.inv(covariance)
    dmat = np.concatenate((scale*np.sin(theta)[:, None, None]*np.eye(2),
                           np.cos(theta)[:, None, None]*np.eye(2)), axis=1)
    aw = np.einsum("nki,kl,nlj->nij", dmat, w, dmat)
    bw = np.einsum("nki,k->ni", dmat, w@y)
    logdet = np.linalg.slogdet(aw)[1]
    quad = np.einsum("ni,ni->n", bw, np.linalg.solve(aw, bw[..., None])[..., 0])
    logs = -.5*logdet+.5*quad
    density = np.exp(logs-logs.max())
    h = theta[1]-theta[0]
    norm = h/3*(density[0]+density[-1]+4*density[1:-1:2].sum()+2*density[2:-1:2].sum())
    return density/norm


def main():
    cases = {
        "both_zero": ([0.,0.,0.], [0.,0.,0.]),
        "force_zero": ([0.,0.,0.], [1.,0.,0.]),
        "acceleration_zero": ([2.,0.,0.], [0.,0.,0.]),
        "aligned": ([2.,0.,0.], [1.,0.,0.]),
        "perpendicular": ([0.,2.,0.], [1.,0.,0.]),
        "anti_aligned": ([-2.,0.,0.], [1.,0.,0.]),
    }
    lower = np.array([[.4,0,0,0,0,0], [.1,.7,0,0,0,0],
                      [-.1,.2,.9,0,0,0], [.2,.1,0,.6,0,0],
                      [0,-.2,.1,.15,.8,0], [.1,0,-.1,0,.2,1.1]])
    covariances = {"small_isotropic": .25**2*np.eye(6),
                   "large_isotropic": 2.**2*np.eye(6),
                   "full_correlated": lower@lower.T}
    rows = []
    for name, covariance in covariances.items():
        for case, (f,a) in cases.items():
            coarse = alternative(f,a,covariance,count=6001,extent=10.)
            fine = alternative(f,a,covariance,count=14401,extent=12.)
            change = max(abs(coarse["mean"]/fine["mean"]-1),
                         np.max(abs(np.array(coarse["quantiles"])/fine["quantiles"]-1)))
            assert fine["mean"] > 0 and math.isfinite(fine["mean"])
            assert 0 < fine["quantiles"][0] < fine["quantiles"][1] < fine["quantiles"][2] < math.inf
            assert change < 1e-4, (name, case, change)
            rows.append({"covariance":name,"case":case,"force":f,"acceleration":a,
                         "result":fine,"refinement_relative":float(change)})

    # Continuous paths, with one prescription even at zeros and adverse angles.
    paths = []
    for degree in np.linspace(0.,180.,13):
        angle = math.radians(degree)
        result = alternative([2*math.cos(angle),2*math.sin(angle),0.], [1.,0.,0.],
                             covariances["full_correlated"],count=6001)
        paths.append({"path":"angle_degrees","at":float(degree),"result":result})
    for value in [-.1,-.01,0.,.01,.1]:
        for name, f, a in [("force_through_zero",[value,0.,0.],[1.,0.,0.]),
                            ("acceleration_through_zero",[2.,0.,0.],[value,0.,0.])]:
            result = alternative(f,a,covariances["full_correlated"],count=6001)
            paths.append({"path":name,"at":value,"result":result})
    for row in paths:
        assert row["result"]["mean"] > 0 and all(q > 0 and math.isfinite(q) for q in row["result"]["quantiles"])
    aligned = next(r for r in rows if r["case"]=="aligned" and r["covariance"]=="small_isotropic")
    anti = next(r for r in rows if r["case"]=="anti_aligned" and r["covariance"]=="small_isotropic")
    assert abs(aligned["result"]["mean"]-anti["result"]["mean"]) > .1
    accurate = alternative([2.,0.,0.],[1.,0.,0.], .03**2*np.eye(6))
    assert abs(accurate["mean"]/2.-1) < .005

    # An independent Cartesian marginalization verifies the target's 2D prior
    # identity on all six configurations, including a correlated covariance.
    indices = [0,1,3,4]
    cov2 = covariances["full_correlated"][np.ix_(indices,indices)]
    identity_rows = []
    for case, (f,a) in cases.items():
        result = estimate(f[:2],a[:2],cov2,direction_order=512,ratio_order=4096)
        density = flat_cartesian_angle_density(f[:2],a[:2],cov2,result.ratio_angle,result.mass_scale)
        error = float(np.max(abs(density-result.angle_density)))
        assert error < 1e-8, (case,error)
        identity_rows.append({"case":case,"max_angle_density_absolute_error":error})
    output = {"meaning":"Own Gaussian-vector/lognormal-scale specialization, not a printed competitor estimate.",
              "prior":"v~N3(0,I), log(m)~N(0,1), independent; fixed across cases",
              "covariances":{k:v.tolist() for k,v in covariances.items()},
              "cases":rows,"paths":paths,"well_measured_aligned":accurate,
              "target_2d_prior_identity":identity_rows,
              "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "limits":"Numerical refinement is not a certified error bound or an everywhere proof."}
    destination = ROOT/"results"/"coverage_third_pass_2026-09-08.json"
    destination.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"cases":len(rows),"paths":len(paths),
                      "max_refinement_relative":max(r["refinement_relative"] for r in rows),
                      "identity_max_absolute_error":max(r["max_angle_density_absolute_error"] for r in identity_rows),
                      "aligned_mean":aligned["result"]["mean"],"anti_aligned_mean":anti["result"]["mean"],
                      "well_measured_mean":accurate["mean"],"output":str(destination)},indent=2))


if __name__ == "__main__":
    main()
