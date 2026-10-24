"""Constructive learning law for repeated diagonal instrument groups.

For group g (same declared force/acceleration noise pair), M_g coordinates
share beta_g; beta_g prior Beta(1/2,1/2). Conditional evidence:
  (1-beta_g)^(M_g/2) exp(beta_g H_g(m)/2).
Posterior is proper at every finite input; conditional candidate coordinates
are N(beta_g A_g,beta_g I). Splitting retains the group labels.

At a fixed candidate mass u0=log m, let w_g=dtheta_g/du and
C_g=lim mean squared fixed candidate coefficients. With finite groups,
regular replication and positive C_g, beta_g->C_g/(1+C_g).
Then:
 expected curvature = sum M_g beta_g C_g w_g^2,
 score variance     = sum M_g beta_g^2 (1+C_g) w_g^2,
which agree exactly in the limit. No Gaussian generating-push assumption.

This is a positive result for repeated diagonal instrument groups, not a
unique full-covariance estimator. A new zero observation within a group
updates its learned scale. Exact grouping of almost-equal noise pairs needs
a continuous extension; that unresolved issue is not hidden here.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import math
import json
import numpy as np


@lru_cache(maxsize=None)
def phi_rule(order=1200):
    z,w=np.polynomial.legendre.leggauss(order)
    phi=(z+1)*np.pi/4
    return np.sin(phi)**2, np.log(w*np.pi/4), np.log(np.cos(phi))


def group_moments(H,M,order=1200):
    beta,lw,lc=phi_rule(order)
    logw=lw+M*lc+.5*H*beta
    shift=logw.max()
    weights=np.exp(logw-shift)
    return float(weights@beta/weights.sum()),float(shift+np.log(weights.sum()))


def exact_replication_result():
    powers=[Fraction(1,10),Fraction(99,10)]
    slopes=[Fraction(1,2),Fraction(10,101)]
    beta=[c/(1+c) for c in powers]
    j=sum(b*c*w*w for b,c,w in zip(beta,powers,slopes))
    v=sum(b*b*(1+c)*w*w for b,c,w in zip(beta,powers,slopes))
    assert j==v
    rows=[]
    for N in (10,100,1000,10000):
        learned=[group_moments(N*(1+float(c)),N)[0] for c in powers]
        curvature=N*sum(b*float(c)*float(w)**2 for b,c,w in zip(learned,powers,slopes))
        rows.append({"N":N,"beta_at_expected_statistics":learned,
                     "law_variance_times_N":N/curvature})
    coarse=group_moments(10000*(1+float(powers[1])),10000)
    fine=group_moments(10000*(1+float(powers[1])),10000,1600)
    assert abs(coarse[0]-fine[0])<1e-11
    return {"noise_SD_F":[1,10],"noise_SD_a":[1,1],
            "candidate_squared_sizes":[float(x) for x in powers],
            "beta_limits":[float(x) for x in beta],
            "information_per_repetition":float(j),
            "information_exact":str(j),
            "asymptotic_variance_times_N":float(1/j),
            "sampling_over_law_variance_limit":float(v/j),
            "finite_expected_statistics":rows,
            "quadrature_beta_refinement":abs(coarse[0]-fine[0])}


def formula_checks():
    # Same group: batch and sequential factors have identical joint exponent.
    beta,lw,lc=phi_rule(400)
    H1,M1,H2,M2=5.2,3,2.1,2
    together=lw+(M1+M2)*lc+.5*(H1+H2)*beta
    sequential=lw+M1*lc+.5*H1*beta+M2*lc+.5*H2*beta
    assert np.max(abs(together-sequential))<1e-13
    # A whole independent zero-only group contributes a mass-independent factor.
    zero_log=group_moments(0.,2)[1]
    assert math.isfinite(zero_log)
    # An independent measured-zero group of any positive absolute noise yields
    # unchanged active mass evidence (same declared mass prior held fixed).
    mgrid=np.exp(np.linspace(-3,3,251))
    x,y=3.,1.
    H=(x*mgrid/np.sqrt(1+mgrid*mgrid)+y/np.sqrt(1+mgrid*mgrid))**2
    logs=np.array([group_moments(float(h),1,400)[1] for h in H])
    padded=logs+zero_log
    assert np.max(abs((padded-padded[0])-(logs-logs[0])))<1e-14
    # The forced exact-group rule is discontinuous if distinct groups suddenly
    # merge at equal noise: quantify a law-shape difference at aligned data.
    merged=np.array([group_moments(float(h),3,400)[1] for h in H])
    log_odds_change=(merged-merged[0])-(logs-logs[0])
    return {"sequential_error":float(np.max(abs(together-sequential))),
            "independent_zero_group_log_factor":zero_log,
            "distinct_vs_merged_group_max_log_odds_change":
                float(np.max(abs(log_odds_change))),
            "scope":"Null-group identity holds with starting mass law held fixed."}


def main():
    out={"source_commit":"592d144531fa4a30d1c2229a4c3b528e84043d77",
         "constructive_regular_result":exact_replication_result(),
         "formula_checks":formula_checks(),
         "status":"Constructive repeated-diagonal-group law; continuous full-covariance grouping not selected."}
    Path(__file__).with_name("grouped_learning_results.json").write_text(
        json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
