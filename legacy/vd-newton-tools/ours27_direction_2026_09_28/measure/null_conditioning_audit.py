"""Scratch derivations and numerical audit, 2026-09-28. No proposed adopted law.

Managed context remains VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.
Continuing ../.. /ours27_bridge_2026_09_28/WORKING_STATE.md.  All observations
below concern candidate measures, not an ontic Gaussian distribution of pushes.

COMMON NOISE-UNIT GAUSSIAN MIXTURE
--------------------------------
Use standardized channel coordinates, mass angle theta, M candidate components,
A_theta=x sin(theta)+y cos(theta), H=||A||^2, T=||x||^2+||y||^2.
For candidate v|lambda ~ N(0,lambda I_M), beta=lambda/(1+lambda),

 L(theta,lambda) = (1+lambda)^(-M/2) exp(-T/2+beta H/2),
 v|theta,lambda,data ~ N(beta A,beta I_M).

The common mixing MEASURE q(dlambda)=dlambda/(1+lambda) is improper but
sigma-finite. It is log observable total variance: t=1+lambda>=1, dt/t.
Its candidate measures mu_M(dv)=integral N_M(dv;0,lambda I)q(dlambda)
are sigma-finite, locally finite for every M>=1, and exactly projective:
integral mu_(M+k)(dv,dw)=mu_M(dv), by Tonelli. This is not a proper
probability process with a finite zero-dimensional marginal: mu_0 has
infinite mass. Once data with M>0 arrive, the joint posterior is proper.

 beta|theta,data has density proportional
   (1-beta)^(M/2-1) exp(beta H/2),  0<beta<1.
 For H>=0 its normalization is finite, bounded between 2/M and exp(H/2)*2/M.
 p(theta|data) is proportional
   pi(theta) exp(-Q/2) gamma(M/2,H/2)/(H/2)^(M/2), Q=T-H,
 with the continuous ratio at H=0 equal to 2/M.
 Conditional candidate moments of every fixed polynomial order are bounded
 uniformly in beta because 0<=beta<=1; finite instrument/angle geometry then
 supplies finite ratio-of-means readout moments. Existence is not uniqueness.

APPENDING MEASURED ZERO COMPONENTS
--------------------------------
k extra zero components multiply the joint law by (1-beta)^(k/2).
They change M to M+k, leaving H,T unchanged. The new beta posterior is
stochastically smaller, strictly so for a nondegenerate old beta posterior.
For evidence F_M(H), d log F_M/dH = E_M[beta|H]/2, and consequently
 d log(F_(M+k)/F_M)/dH < 0.
Thus an added null observation reduces angular odds favoring larger H.
At all-zero data H=0 the mass marginal remains pi; that special test cannot
detect the change. The all-zero beta mean changes from 2/(M+2) to
2/(M+k+2). The candidate conditional law changes as well as the mass law.

NULL-NEUTRALITY THEOREM (one M and one k suffice)
-----------------------------------------------
Let q be ANY positive mixing measure on finite lambda>=0 for which F_M(0)
is finite and positive. Suppose appending k>0 zero components leaves the
normalized mass marginal unchanged for data whose H(theta) spans an open
interval and whose pi is positive there. Then
 F_(M+k)(H)=c F_M(H) on that interval. Uniqueness of Laplace transforms of
finite measures on beta in [0,1] yields
 [(1-beta)^(k/2)-c] (1-beta)^(M/2)q_beta(dbeta)=0.
Since the bracket has a unique zero, q is concentrated at ONE finite beta.
The converse holds. It does NOT force beta=1: every fixed finite variance
is an exact counterexample to that stronger statement. Flat Cartesian
measure is the separate infinite-scale angular limit, not a finite-variance
Gaussian probability mixture component at beta=1.

STRONGER NOISE-LIMIT CONSEQUENCE
-------------------------------
Give each added zero component noise SD epsilon in BOTH channels, with the
same instrument ratio. If candidate covariance is lambda times the LOCAL
conditional-noise covariance K_m, standardization removes epsilon exactly.
The same (1-beta)^(k/2) factor occurs for EVERY epsilon>0. Hence asking that
the normalized mass law converge to the old law as epsilon->0 is equivalent
here to exact null neutrality. The theorem forces a fixed scale, even if
finite-noise measured-null neutrality was explicitly relinquished.

This conclusion assumes the lower-dimensional comparison OMITS the added
coordinates with the same mixing q and retained prior. It is not a universal
claim about every meaning of known-empty, every metric, or every conditioning
operation. For all-zero retained data both mass laws happen to equal pi.

PHYSICAL-METRIC PRIOR AND THE CORRECT ZERO-CONDITIONING LIMIT
----------------------------------------------------------
Instead take candidate variance lambda in fixed physical (channel-rescaled)
coordinates, retained noise=1, and added zero-coordinate noise epsilon^2.
The null likelihood factor is (1+lambda/epsilon^2)^(-k/2).
For fixed lambda>0 this is asymptotic to epsilon^k lambda^(-k/2).
If integral lambda^(-k/2)q(dlambda)L_old(theta,lambda) is finite, dominated
convergence proves that the limit uses q_cond proportional lambda^(-k/2)q.
That is exactly the Gaussian zero-section conditioning factor. It need not
equal the marginal/omission prior q, despite candidate marginal projectivity.
Changing the lower-dimensional prior back to q would be a substantive reset.

For q(dlambda)=dlambda/(1+lambda), k=1 is integrable and the limiting beta
law is tilted Beta(1/2,(M+1)/2). More generally the formal conditional law
has beta powers beta^(-k/2)(1-beta)^((M+k)/2-1), so k>=2 is improper at0.
For k>=2 the epsilon>0 posteriors remain proper but concentrate lambda at0
as epsilon->0, erasing the finite retained-data angular evidence: p->pi.

General boundary statement: if q has density c lambda^(a-1)(1+o(1)) at0,
a>0, and old likelihood is positive continuous at0 with adequate tails,
k<2a gives the finite conditional slice above; k>=2a yields lambda->0.
For k>2a, lambda/epsilon^2 has limiting beta-prime density proportional
t^(a-1)(1+t)^(-k/2). For k=2a the normalizer has the logarithmic boundary
factor epsilon^(2a) log(1/epsilon^2). This is size-direction coupling of a
shared physical-scale model, not proof that the physics requires such a model.

CONSTRUCTIVE ESCAPE AND ITS SCOPE
-------------------------------
In noise coordinates let the candidate covariance be lambda D, with D fixed
by a DECLARED instrument design. Retained block D0>0, added zero block
D_epsilon=epsilon^p I, p>0. Null factor
 |I+lambda D_epsilon|^(-1/2) tends to1 pointwise and is <=1. A proper retained
posterior therefore dominates it even when q is improper. The posterior of
retained variables converges in total variation to the original posterior.
Added physical candidates have posterior covariance bounded by their
measurement-noise covariance, so their norm contribution also vanishes.
Thus projective radial mixing and omission-style known-empty continuity CAN
coexist if the standardized candidate covariance shrinks on those axes.
This statement keeps the retained mass measure fixed (as happens for the
same-ratio uniform-angle example). If a mass measure is recomputed from a
larger instrument geometry, its own continuity must be established separately.

This is not forced or unique. It requires an instrument covariance functional
which makes candidate scale shrink FASTER than local noise on the added axis.
For example a global full-design D based on normalized physical covariance
can do this (prior scratch instrument_covariance_scope.py), but D and D^2
both do it and differ elsewhere. A fixed design's principal submatrices give
projective candidate measures; recomputing the global normalization after
removing an arbitrary coordinate is a different rule and need not be projective.

Intrinsic subexperiment consistency is a stronger demand. For independent
axes with covariance c_i^2 times the SAME 2x2 channel covariance, require
each selected-axis prior covariance to be obtainable using that axis's
instrument information alone, with the same global mixing variable. Physical
unit covariance then makes its physical variance c_i^2 times one fixed
coefficient. Sign-flip symmetry removes cross terms. Its standardized
variance is independent of c_i. Under these EXTRA locality conditions the
shrinking-D escape is excluded and the no-go applies. Regrouping one already
declared full design does not by itself entail this stronger locality.

FULL-COVARIANCE WIDTH LIMIT
--------------------------
The dt/t candidate law is algebraically defined for arbitrary full covariance:
use K_m and H=sum y' Sigma^-1 B K_m B' Sigma^-1 y, M=sum candidate dimensions.
But scalar shared scale learns unweighted candidate energy R^2/M, whereas
the needed regular power is J/(J+G), J=sum v'Tv, G=sum tr(T),
T=K^(1/2)V^(-1)K^(1/2). These coincide under common isotropic geometry or
isotropic empirical candidate allocation, not in general. Existence and
projectivity do not solve that independent spectral-allocation flaw.

The executable section checks the concrete angular changes and the two
distinct epsilon limits. It is not a coverage sweep or a fitted correction.
"""

from __future__ import annotations

import json
import math
import numpy as np


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return lo + (hi - lo) * (x + 1) / 2, w * (hi - lo) / 2


def logsumexp_columns(log_terms):
    peak = np.max(log_terms, axis=0)
    return peak + np.log(np.sum(np.exp(log_terms - peak), axis=0))


def angular_case(m_count=1, added=0, epsilon=1.0, metric="noise", n=500):
    # x=2,y=0 on the retained scalar axis; common channel ratio=1.
    # Integrate the two halves separately so P(m<1) has no binning error.
    theta1, tw1 = gl(96, 0, math.pi / 4)
    theta2, tw2 = gl(96, math.pi / 4, math.pi / 2)
    theta = np.concatenate([theta1, theta2])
    tw = np.concatenate([tw1, tw2])
    h = 4 * np.sin(theta) ** 2
    loglam, lw = gl(n, -90, 50)
    lam = np.exp(loglam)
    beta = lam / (1 + lam)
    # q(dlambda)=dlambda/(1+lambda), including the log-lambda Jacobian.
    logweight = np.log(lw) + loglam - (1 + m_count / 2) * np.log1p(lam)
    if metric == "noise":
        logweight -= added / 2 * np.log1p(lam)
    elif metric == "physical":
        logweight -= added / 2 * np.log1p(lam / epsilon**2)
    elif metric == "shrinking":
        # D_extra=epsilon^2 in noise coordinates; retained D=1.
        logweight -= added / 2 * np.log1p(lam * epsilon**2)
    elif metric == "conditional":
        logweight -= added / 2 * loglam
    else:
        raise ValueError(metric)
    joint = logweight[:, None] + beta[:, None] * h[None, :] / 2
    loge = logsumexp_columns(joint)
    angular = tw * np.exp(loge - np.max(loge))
    total = np.sum(angular)
    local_mean_beta = np.sum(np.exp(joint - loge) * beta[:, None], axis=0)
    return {
        "P_mass_below_one": float(np.sum(angular[: len(theta1)]) / total),
        "posterior_mean_beta": float(np.sum(angular * local_mean_beta) / total),
    }


def main():
    result = {"base": angular_case()}
    result["noise_unit_null"] = {str(k): angular_case(added=k) for k in (1, 2, 3)}
    result["physical_zero_limit"] = {
        str(k): {
            str(eps): angular_case(added=k, epsilon=eps, metric="physical")
            for eps in (1, 0.1, 0.01, 1e-4, 1e-8)
        }
        for k in (1, 2, 3)
    }
    result["physical_k1_exact_conditional_limit"] = angular_case(added=1, metric="conditional")
    result["shrinking_D_k3"] = {
        str(eps): angular_case(added=3, epsilon=eps, metric="shrinking")
        for eps in (1, 0.1, 0.01, 1e-4)
    }
    result["quadrature_check_800"] = {
        "base": angular_case(n=800),
        "physical_k1_eps1e8": angular_case(added=1, epsilon=1e-8, metric="physical", n=800),
        "physical_k2_eps1e8": angular_case(added=2, epsilon=1e-8, metric="physical", n=800),
        "shrinking_k3_eps1e4": angular_case(added=3, epsilon=1e-4, metric="shrinking", n=800),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
