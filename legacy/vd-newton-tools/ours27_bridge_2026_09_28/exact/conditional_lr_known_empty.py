"""Exact conditional LR is calibrated, but fails the known-empty-axis limit.

This is a limitation of this confidence-set construction, not a defect in its
finite-sample coverage statement. No candidate probability law is asserted.

Whiten the measurements z and write P_theta for the orthogonal projection onto
the whitened codirectional candidate subspace. Under a tested theta, the mean
is in range(P_theta). P_theta z and (I-P_theta)z are independent, and the latter
is standard normal in the orthocomplement. Conditional simulation

  z_star = P_theta z + W_theta eta, eta ~ N(0,I)

therefore eliminates all nuisance candidate components exactly. Computing the
same global LR statistic for observed and simulated data gives a valid
conditional test. Supremizing over the closed positive-angle domain handles
endpoint limits, ties and higher-order peaks without selecting a local mode.

However, take a one-axis reading F=2,a=0, sigmaF=sigmaa=1. Test m=1, so theta
=pi/4. The observed tested projection has A=sqrt(2), B=sqrt(2); the positive
domain global LR (twice log likelihood ratio) equals 2. Conditional resampling
keeps A=sqrt(2) and sets B=eta. The simulated LR exceeds 2 exactly when
|eta|>=sqrt(2), giving p1=erfc(1).

Now append an axis measured as F=a=0, with sigmaF=sigmaa=epsilon>0. In whitened
coordinates this row is zero, and its projection geometry is independent of
epsilon. Conditional resampling appends A2=0, B2=zeta, zeta~N(0,1). The observed
LR is still 2. The conditional simulated LR exceeds 2 exactly when

  eta^2+zeta^2+2*sqrt(2)*|eta| >= 6.

For eta^2+zeta^2 <= 2 the interior principal direction cannot reach the threshold;
otherwise the maximum on the allowed +/-pi/4 angular interval is an endpoint,
yielding the displayed condition. In polar coordinates (eta,zeta), the required
radius is r_min(phi)=sqrt(6+2*cos(phi)^2)-sqrt(2)*|cos(phi)|. Its standard-normal
tail is exp(-r_min^2/2), so the exact p-value is

  p2=(2/pi)*integral_0^(pi/2) exp(-r_min(phi)^2/2) dphi.

p2 is strictly greater than p1 and independent of epsilon. Thus arbitrarily
precise observed-empty added axes do not recover the lower-dimensional p-curve.
The known-exact-empty submodel (epsilon=0 with the axis removed) is a different
conditional experiment; taking a small-noise limit does not commute with the
profile nuisance fit and conditional resampling.
"""
from pathlib import Path
from math import erfc
import json
import numpy as np
from numpy.polynomial.legendre import leggauss


def polar_integral(order):
    z,w=leggauss(order)
    phi=(z+1)*np.pi/4
    radius=np.sqrt(6+2*np.cos(phi)**2)-np.sqrt(2)*np.cos(phi)
    return float(np.dot(w,np.exp(-radius*radius/2))/2)


def numerical_lr(a,b):
    """Global projection energy on theta0 +/- pi/4, minus tested energy."""
    h=float(a@a);q=float(b@b);c=float(a@b)
    boundary=(h+q)/2+abs(c)
    if h>=q:
        best=max(boundary,(h+q+np.hypot(h-q,2*c))/2)
    else:
        best=max(boundary,h)
    return best-h


def main():
    p1=erfc(1.)
    vals=[polar_integral(n) for n in (32,64,128,256)]
    assert max(abs(np.array(vals)-vals[-1]))<1e-13
    rng=np.random.default_rng(280932)
    eta=rng.normal(size=(20000,2))
    algebra=eta[:,0]**2+eta[:,1]**2+2*np.sqrt(2)*abs(eta[:,0])>=6
    direct=np.array([numerical_lr(np.array([np.sqrt(2),0.]),b)>=2 for b in eta])
    assert np.array_equal(algebra,direct)
    out={'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
         'base_data':{'F':[2.],'a':[0.],'sigmaF':[1.],'sigmaa':[1.]},
         'tested_m':1.,'base_p_value':p1,
         'appended_zero_axis_p_value_all_positive_epsilon':vals[-1],
         'absolute_jump':vals[-1]-p1,'quadrature_orders':[32,64,128,256],
         'quadrature_values':vals,'global_lr_algebra_checks':len(eta),
         'monte_carlo_appended_p_value':float(np.mean(algebra)),
         'status':'Exact conditional coverage benchmark; fails known-empty-axis continuity.'}
    Path(__file__).with_name('conditional_lr_known_empty_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
