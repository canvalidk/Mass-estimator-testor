"""Scratch proof: L11 is impossible for every estimator on its stated domain.

Source and scope
----------------
VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8,
Newton-analysis/_dscn_mass_estimator/mass_estimator_lemmas_2026-09-27.md,
L11/T17. L11 requires mass consistency whenever R_N^2=sum_i r_i^2
diverges, with arbitrary per-reading candidates. This proof concerns that
quantifier; it does not identify mass with a metaphysical hidden parameter.
Here a prescribed m is a coefficient in the candidate experiment used to
check the definition's operational reproducibility.

Proof
-----
1. One noise pair, s=1, is a subcase of ours27. For positive m=tan(theta),
   let u_theta=(sin(theta),cos(theta)). Reading i is
       Z_i = sign_i * r_i * u_theta + E_i,
   where E_i are independent N(0,I_2), and all covariances are known.
   Each deterministic sign sequence is an allowed candidate sequence.

2. Choose epsilon>0 and r_i^2=epsilon*i^(-3/4). For EVERY sign sequence,
       R_N^2=epsilon*sum_{i<=N} i^(-3/4) -> infinity,
   while sum_i r_i^4 = epsilon^2*sum_i i^(-3/2) <= 3 epsilon^2.
   The latter bound follows from 1+integral_1^infinity x^(-3/2)dx=3.

3. Temporarily average the signs independently with probabilities 1/2.
   This is an impossibility-proof device, not an imposed estimator prior.
   Write phi_0 for the standard two-dimensional Gaussian density. The
   resulting one-reading density f_{i,theta} has likelihood ratio
       L_{i,theta}(z)=f_{i,theta}(z)/phi_0(z)
                    =exp(-r_i^2/2)*cosh(r_i*u_theta^T z).
   Since u_theta^T Z ~ N(0,1) under phi_0,
       E_0[L]=1,
       E_0[L^2]=exp(-r_i^2)*[1+exp(2r_i^2)]/2
               =cosh(r_i^2).
   This identity is independent of theta.

4. Let P_{theta,N} be the product mixture experiment for the first N
   readings; Q_N is independent pure standard Gaussian noise. Independence
   gives the EXACT chi-squared divergence
       chi2(P_{theta,N} || Q_N) = product_i cosh(r_i^2)-1.
   The elementary inequality log(cosh x)<=x^2/2 follows by differentiating:
   x-tanh x >=0 for x>=0. Therefore, uniformly in N and theta,
       chi2(P_{theta,N} || Q_N)
          <= exp((1/2)*sum_i r_i^4)-1
          <= exp(3 epsilon^2/2)-1.

5. Total variation is TV(P,Q)=(1/2)integral|p-q|. Cauchy-Schwarz implies
       TV(P,Q)<=(1/2)*sqrt(chi2(P||Q)).
   The triangle inequality through the SAME Q_N then gives, for any two
   distinct positive masses,
       TV(P_{theta1,N},P_{theta2,N})
          <= sqrt(exp(3 epsilon^2/2)-1).
   At epsilon=.1 this is about .123, uniformly in N. Thus every test of
   these two masses has type-I plus type-II error at least about .877:
   a consistent discriminator cannot exist.

6. Suppose an estimator satisfied L11 for every deterministic candidate
   sequence at both masses. Choose disjoint fixed neighborhoods of the
   two masses. Conditional on EACH infinite deterministic sign sequence,
   the estimator's error probability would tend to zero because R_N^2
   diverges. Error probabilities lie in [0,1], so dominated convergence
   over the product measure on infinite signs says the mixture error also
   tends to zero at each mass. Classifying by the estimator's neighborhood
   would give a consistent discriminator, contradicting step 5.

Conclusion
----------
No estimator can satisfy L11 on all its stated sequences. Total signal
energy divergence must be supplemented by a detectability condition or
by additional structure on the candidate sequence. This proof does not
claim a unique replacement condition for every conceivable estimator.

Numerical checks below verify the Gaussian integrals independently by
Gauss-Hermite quadrature, exhibit discrete amplitudes, and check finite
products against the analytic bound. The universal result is the proof,
not an extrapolation from those numerical checks.
"""
from pathlib import Path
import json
import math
import numpy as np
from numpy.polynomial.hermite import hermgauss


def normal_integrals(r_squared, nodes):
    x, w = hermgauss(nodes)
    standard_normal = np.sqrt(2) * x
    probability_weights = w / np.sqrt(np.pi)
    likelihood_ratio = np.exp(-.5 * r_squared) * np.cosh(np.sqrt(r_squared) * standard_normal)
    return float(probability_weights @ likelihood_ratio), float(probability_weights @ likelihood_ratio**2)


def main():
    epsilon = .1
    integral_checks = []
    for r_squared in (.0001, .001, .01, .1, .3, 1., 3.):
        mean, second = normal_integrals(r_squared, 96)
        fine_mean, fine_second = normal_integrals(r_squared, 160)
        exact = math.cosh(r_squared)
        integral_checks.append({
            "r_squared": r_squared,
            "integral_L": mean,
            "integral_L_squared": second,
            "cosh_r_squared": exact,
            "first_moment_error": abs(mean - 1),
            "second_moment_relative_error": abs(second / exact - 1),
            "refinement_relative_error": abs(fine_second / second - 1),
        })
        assert abs(mean - 1) < 2e-13
        assert abs(second / exact - 1) < 2e-13
        assert abs(fine_second / second - 1) < 2e-13
        assert abs(fine_mean - 1) < 2e-13

    universal_chi_squared_bound = math.expm1(1.5 * epsilon**2)
    universal_two_mass_tv_bound = math.sqrt(universal_chi_squared_bound)
    finite_checks = []
    for n in (1, 10, 1000, 1000000):
        indices = np.arange(1, n + 1, dtype=float)
        r_squared = epsilon * indices**(-.75)
        sum_r_fourth = float(r_squared @ r_squared)
        # logcosh(x)=log1p(2*sinh(x/2)^2), accurate even for small x.
        log_product = float(np.log1p(2 * np.sinh(r_squared / 2)**2).sum())
        chi_squared = math.expm1(log_product)
        finite_bound = math.expm1(.5 * sum_r_fourth)
        finite_checks.append({
            "N": n,
            "first_r_squared": float(r_squared[0]),
            "last_r_squared": float(r_squared[-1]),
            "total_signal_energy_R_squared": float(r_squared.sum()),
            "sum_r_fourth": sum_r_fourth,
            "exact_chi_squared_from_product": chi_squared,
            "finite_chi_squared_upper_bound": finite_bound,
            "two_mass_TV_upper_bound_via_pure_noise": math.sqrt(chi_squared),
        })
        assert chi_squared <= finite_bound + 1e-14
        assert finite_bound <= universal_chi_squared_bound
        assert sum_r_fourth <= 3 * epsilon**2

    out = {
        "source_commit": "b4bfc2040031ef010416334d0c5f4aa2ad05cbd8",
        "epsilon": epsilon,
        "amplitude_rule": "r_i^2 = epsilon * i^(-3/4)",
        "total_signal_energy_diverges": True,
        "universal_sum_r_fourth_upper_bound": 3 * epsilon**2,
        "universal_chi_squared_to_pure_noise_upper_bound": universal_chi_squared_bound,
        "universal_TV_between_two_masses_upper_bound": universal_two_mass_tv_bound,
        "universal_lower_bound_on_sum_of_testing_errors": 1 - universal_two_mass_tv_bound,
        "gaussian_integral_checks": integral_checks,
        "finite_discrete_amplitude_checks": finite_checks,
        "all_checks_passed": True,
    }
    Path(__file__).with_name("impossibility_check_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
