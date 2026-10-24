"""Independent scalar oracles and explicit LS/Fieller boundary witnesses.

This is a bounded equation-comparison experiment, not a new estimator.
The correlated normal-means extension is our derivation, not a theorem
attributed to Liseo. Output is separate from earlier research records.
Only NumPy and the Python standard library are required.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from estimator import estimate
from literature_equivalence_check import katz, liseo_positive_angle

ROOT = Path(__file__).resolve().parent


def normal_cdf(values):
    arr = np.asarray(values)
    return np.fromiter((.5*math.erfc(-float(v)/math.sqrt(2)) for v in arr.ravel()),
                       dtype=float, count=arr.size).reshape(arr.shape)


def normal_pdf(values):
    return np.exp(-.5*np.asarray(values)**2)/math.sqrt(2*math.pi)


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def cartesian_oracle(force, acceleration, sf, sa, rho, masses, *, order, signs=(1., -1.)):
    """Integrate over signed acceleration, analytically over conditional force.

    Uses original Cartesian Gaussian posterior, unlike the reference's
    ratio-angle/radial integration. The acceleration interval truncates
    beyond |observed acceleration|+12 sigma; no mass cutoff is imposed.
    """
    n, w = nodes(order)
    extent = abs(acceleration) + 12*sa
    alpha, weights = (n+1)*extent/2, w*extent/2
    conditional_sd = sf*math.sqrt(1-rho*rho)
    ztotal = ftotal = atotal = 0.
    cumulatives = np.zeros(len(masses))
    for sign in signs:
        marginal = normal_pdf((sign*alpha-acceleration)/sa)/sa
        conditional_mean = sign*force + rho*sf/sa*(alpha-sign*acceleration)
        z = conditional_mean/conditional_sd
        positive_probability = normal_cdf(z)
        weighted_marginal = weights*marginal
        ztotal += float(weighted_marginal @ positive_probability)
        ftotal += float(weighted_marginal @ (conditional_mean*positive_probability
                                           + conditional_sd*normal_pdf(z)))
        atotal += float(weighted_marginal @ (alpha*positive_probability))
        lower_cdf = normal_cdf(-z)
        for index, mass in enumerate(masses):
            cumulatives[index] += float(weighted_marginal @
                (normal_cdf((mass*alpha-conditional_mean)/conditional_sd)-lower_cdf))
    return dict(same_sign_probability=ztotal, mass=ftotal/atotal,
                cdf_at_supplied_masses=(cumulatives/ztotal).tolist())


def main():
    rows = []
    probs = [.025, .5, .975]
    independent_errors = []
    null_cauchy_checks = []
    known_direction_rows = []
    sf, sa = .7, 1.3
    for zf, za in [(0.,0.), (0.,2.), (2.,0.), (3.,2.), (-2.,-3.), (1.,-1.), (-2.,.5)]:
        force, acceleration = sf*zf, sa*za
        for rho in [-.6, 0., .6]:
            covariance = np.array([[sf*sf, rho*sf*sa], [rho*sf*sa, sa*sa]])
            reference = estimate([force], [acceleration], covariance, ratio_order=8192)
            quantiles = [reference.quantile(p) for p in probs]
            coarse = cartesian_oracle(force, acceleration, sf, sa, rho, quantiles, order=300)
            oracle = cartesian_oracle(force, acceleration, sf, sa, rho, quantiles, order=600)
            point_error = abs(reference.mass/oracle["mass"]-1)
            cdf_error = max(abs(np.array(oracle["cdf_at_supplied_masses"])-probs))
            quadrature_change = max(abs(np.array(coarse["cdf_at_supplied_masses"])
                                         -oracle["cdf_at_supplied_masses"]))
            assert point_error < 1e-7, (zf, za, rho, point_error)
            assert cdf_error < 1e-4, (zf, za, rho, cdf_error)
            assert quadrature_change < 1e-6, (zf, za, rho, quadrature_change)
            rows.append(dict(z_force=zf, z_acceleration=za, covariance=covariance.tolist(),
                             rho=rho, point=reference.mass, quantiles=quantiles,
                             oracle=oracle, point_relative_error=point_error,
                             quantile_cdf_absolute_error=float(cdf_error),
                             oracle_quadrature_cdf_change=float(quadrature_change)))
            if zf == 0. and za == 0.:
                # Ratio of centered correlated normals is shifted Cauchy;
                # conditioning its sign gives a positive-truncated Cauchy.
                location = rho*sf/sa
                width = sf/sa*math.sqrt(1-rho*rho)
                cdf_zero = .5 + math.atan(-location/width)/math.pi
                expected = [(.5+math.atan((m-location)/width)/math.pi-cdf_zero)
                            /(1-cdf_zero) for m in quantiles]
                error = max(abs(np.asarray(expected)-probs))
                assert error < 1e-7
                assert abs(oracle["same_sign_probability"]-(.5+math.asin(rho)/math.pi))<1e-12
                null_cauchy_checks.append(dict(rho=rho, location=location, scale=width,
                                               quantile_cdf_absolute_error=float(error)))
            if rho == 0.:
                density = liseo_positive_angle(reference.ratio_angle, zf, za)
                density_error = float(np.max(abs(reference.angle_density/density-1)))
                independent_errors.append(density_error)
                assert density_error < 1e-9
                known = estimate([force], [acceleration], covariance,
                                 known_direction=[1.], ratio_order=8192)
                quotient = katz(force,sf)/katz(acceleration,sa)
                katz_error = abs(known.mass/quotient-1)
                assert katz_error < 1e-9
                known_direction_rows.append(dict(z_force=zf, z_acceleration=za,
                                                 katz_quotient=quotient,
                                                 vd_point=known.mass,
                                                 relative_error=katz_error))

    # Keeping the known sign does not license multiplying marginal posteriors
    # when measurement errors are correlated.
    force, acceleration, rho = 1.4, 1.3, .6
    cov = np.array([[sf*sf,rho*sf*sa],[rho*sf*sa,sa*sa]])
    known = estimate([force],[acceleration],cov,known_direction=[1.],ratio_order=8192)
    corr_oracle = cartesian_oracle(force,acceleration,sf,sa,rho,[],order=600,signs=(1.,))
    assert abs(known.mass/corr_oracle["mass"]-1)<1e-7
    naive = katz(force,sf)/katz(acceleration,sa)
    correlation_witness = dict(force=force,acceleration=acceleration,covariance=cov.tolist(),
                               joint_point=known.mass,cartesian_point=corr_oracle["mass"],
                               independent_katz_quotient=naive,
                               note="The independent quotient discards the supplied correlation.")

    # This is the known-covariance Gaussian Fieller inversion, a specialization
    # of the ratio confidence-set construction, not its sample-covariance t case.
    critical = 1.959963984540054
    fixtures = []
    for name, force, acceleration in [("zero_zero",0.,0.),("force_zero",0.,4.),
                                     ("acceleration_zero",4.,0.),("aligned",4.,4.),
                                     ("opposite_sign",-4.,4.)]:
        # Q(m)<=critical^2 is equivalent to A m^2+B m+C<=0.
        A=acceleration**2-critical**2
        B=-2*force*acceleration
        C=force**2-critical**2
        if name=="zero_zero":
            positive_set={"type":"all_positive"}
        elif name=="force_zero":
            positive_set={"type":"upper_bounded","upper":critical/math.sqrt(A)}
        elif name=="acceleration_zero":
            positive_set={"type":"lower_bounded","lower":math.sqrt(C)/critical}
        elif name=="aligned":
            roots=sorted([(-B-math.sqrt(B*B-4*A*C))/(2*A),
                          (-B+math.sqrt(B*B-4*A*C))/(2*A)])
            positive_set={"type":"bounded","bounds":roots}
        else:
            assert A>0 and B>0 and C>0
            positive_set={"type":"empty"}
        fixtures.append(dict(case=name,force=force,acceleration=acceleration,
                             covariance=np.eye(2).tolist(),quadratic=[A,B,C],
                             positive_fieller_set=positive_set,
                             point_ratio=force/acceleration if acceleration else None))

    # Independent algebra: the positive-constrained TLS profile for equal
    # anti-aligned unit vectors has no finite minimizer.
    masses=[.001,.01,.1,1.,10.,100.,1000.]
    tls=[dict(mass=m,profile=(1+m)**2/(1+m*m)) for m in masses]
    assert all(r["profile"]>1 for r in tls)
    summary=dict(correlated_scalar_cases=len(rows),known_direction_independent_cases=len(known_direction_rows),
                 max_cartesian_point_relative_error=max(r["point_relative_error"] for r in rows),
                 max_quantile_cdf_absolute_error=max(r["quantile_cdf_absolute_error"] for r in rows),
                 max_oracle_quadrature_cdf_change=max(r["oracle_quadrature_cdf_change"] for r in rows),
                 max_liseo_independent_density_relative_error=max(independent_errors),
                 max_katz_quotient_relative_error=max(r["relative_error"] for r in known_direction_rows),
                 max_correlated_null_cauchy_cdf_error=max(r["quantile_cdf_absolute_error"] for r in null_cauchy_checks),
                 all_assertions="passed")
    result=dict(date="2026-09-08",managed_commit="adfe14bbaafb29856f0a30da088495b7c9d483d2",
                scope="Scalar Gaussian only. Correlation extension, positive conditioning and quotient are explicitly our adaptations.",
                oracle="Conditional-normal Cartesian integration; acceleration truncated beyond 12 sigma; 300 vs 600 nodes.",
                quantile_probabilities=probs,scalar_cases=rows,known_direction_cases=known_direction_rows,
                correlated_null_cauchy_checks=null_cauchy_checks,
                correlation_witness=correlation_witness,known_covariance_fieller=fixtures,
                anti_aligned_tls=tls,summary=summary,
                source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                               for name in ["coverage_second_pass.py","estimator.py","literature_equivalence_check.py"]})
    target=ROOT/"results"/"coverage_second_pass_2026-09-08.json"
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))
    print(json.dumps(correlation_witness,indent=2))
    print(json.dumps(fixtures,indent=2))


if __name__=="__main__":
    main()
