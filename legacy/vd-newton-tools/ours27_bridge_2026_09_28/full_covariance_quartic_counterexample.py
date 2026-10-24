"""Exact counterexample to extending one local scalar profile power globally.

Two independent scalar readings, noise ratios s1=sqrt(2)+1, s2=sqrt(2)-1.
Acceleration SDs are one; force SDs are s_i. Standardized readings x_i=y_i=k,
so physical F_i=s_i*k, a_i=k. The unrestricted profile peaks at m=1.

With u=log m and a=log(s1), omitting an additive constant,
 ell(u)=k^2 * [1 + (sech(u-a)+sech(u+a))/2].
At u=0, j=-ell''=0 but ell''''=-3k^2/sqrt(2)<0.
The exact conditional score variance in u is positive:
 Vc=k^2*(1+1/sqrt(2))/4.
Thus beta=j/Vc=0 erases all profile evidence, even as k increases.

This is an allowed finite-input stress test; at large k it is rare under
the declared Gaussian Newton-II experiment. The project explicitly retains
such inputs and asks for their conditional least-impossible account. It is
not evidence of a typical-data consistency failure of the one-pair formula.
"""
from pathlib import Path
import json
import numpy as np


def sech(x): return 1/np.cosh(x)


def halfwidth(k):
    # Exact stable identity: 1/sqrt(2)-mean(sech(u+-a))
    # = (cosh(u)-1)^2/[sqrt(2)*(cosh(u)^2+1)].
    def loss(u): return k*k*(2*np.sinh(u/2)**2)**2/(np.sqrt(2)*(np.cosh(u)**2+1))
    lo=0.;hi=5.
    for _ in range(70):
        mid=(lo+hi)/2
        if loss(mid)>np.log(2): hi=mid
        else: lo=mid
    return (lo+hi)/2


def main():
    a=np.log(np.sqrt(2)+1)
    grid=np.linspace(-5,5,100001)
    shape=.5*(sech(grid-a)+sech(grid+a))
    # Quartic flatness permits equal floating-point maxima near zero.
    assert abs(grid[np.argmax(shape)])<.001
    second=sech(a)*(2*np.tanh(a)**2-1)
    fourth=sech(a)*(24*np.tanh(a)**4-28*np.tanh(a)**2+5)
    assert abs(second)<1e-15
    assert abs(fourth+3/np.sqrt(2))<1e-14
    rows=[]
    for k in [2.,10.,100.,1000.]:
        s=np.array([np.sqrt(2)+1,np.sqrt(2)-1])
        F=s*k;acc=np.full(2,k);m=1.
        V=s*s+m*m
        b=(m*F+s*s*acc)/V
        vc=float(np.sum(b*b/V))
        expected=k*k*(1+1/np.sqrt(2))/4
        assert abs(vc/expected-1)<1e-14
        width=halfwidth(k)
        rows.append({'k':k,'conditional_score_variance_in_log_mass':vc,
            'profile_second_derivative':0.,'profile_fourth_derivative':float(k*k*fourth),
            'scalar_power_beta':0.,'profile_half_max_log_mass_width':width,
            'width_times_sqrt_k':width*np.sqrt(k)})
    out={'status':'Rejects naive full-covariance global use of local curvature ratio; not one-pair regular calibration.',
         'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
         'rows':rows,'quartic_asymptotic_width_constant':float((8*np.sqrt(2)*np.log(2))**.25)}
    Path(__file__).with_name('full_covariance_quartic_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
