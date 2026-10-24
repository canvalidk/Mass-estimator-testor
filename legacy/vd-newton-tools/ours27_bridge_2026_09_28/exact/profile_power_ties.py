"""Probe for discontinuity in global-profile-selected scalar power.

Scalar per-reading independent channel variances are already a permitted
subcase of full-covariance ours27. Q=sum (F-ma)^2/(sF^2+m^2*sA^2).
At a selected stationary point mhat, proposed beta=(-ell''(mhat))/Jhat,
Jhat=sum b_i(mhat)^2/V_i(mhat), ell=-Q/2.
The calculation below uses analytic derivatives in m and verifies chart
invariance by angle-chart finite differences. It does not fit any power.
"""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import Polynomial


def polynomial_stationary_points(f,a,sf,sa):
    """Enumerate all finite stationary points by clearing positive denominators.

    Q' = sum -2(F_i-ma_i)(a_i*sF_i^2+m*sA_i^2*F_i)/V_i^2.
    Every V_i is positive on the real line, so clearing them adds no real roots.
    Check the projective seam separately: Q(infinity)=sum a_i^2/sA_i^2.
    """
    vv=[Polynomial([ss*ss,0.,tt*tt]) for ss,tt in zip(sf,sa)]
    numerator=Polynomial([0.])
    for i in range(len(f)):
        term=-2*Polynomial([f[i],-a[i]])*Polynomial([a[i]*sf[i]**2,sa[i]**2*f[i]])
        for j in range(len(f)):
            if i!=j:
                term*=vv[j]*vv[j]
        numerator+=term
    roots=numerator.roots()
    real=sorted(float(z.real) for z in roots if abs(z.imag)<1e-7)
    rows=[]
    for m in real:
        for _ in range(5):
            pp=power(m,f,a,sf,sa)
            m-=pp['Q_prime']/(2*pp['negative_ell_second'])
        rows.append(power(m,f,a,sf,sa))
    return {'stationary_polynomial_degree':len(numerator.coef)-1,
            'all_finite_real_stationary_points':rows,
            'Q_at_projective_infinity':float(np.sum(a*a/(sa*sa))),
            'stationary_polynomial_coefficients_increasing_order':numerator.coef.tolist()}


def q_value(m, f, a, sf, sa):
    m=np.asarray(m)[...,None]
    return np.sum((f-m*a)**2/(sf*sf+m*m*sa*sa),axis=-1)


def power(m,f,a,sf,sa):
    r=f-m*a;v=sf*sf+m*m*sa*sa
    vp=2*m*sa*sa;vpp=2*sa*sa
    qp=np.sum(-2*a*r/v-r*r*vp/v**2)
    qpp=np.sum(2*a*a/v+4*a*r*vp/v**2-r*r*vpp/v**2+2*r*r*vp*vp/v**3)
    b=(m*sa*sa*f+sf*sf*a)/v
    j=np.sum(b*b/v)
    return {"m":float(m),"Q":float(q_value(m,f,a,sf,sa)),"Q_prime":float(qp),
            "negative_ell_second":float(qpp/2),"J_hat":float(j),"beta":float(qpp/(2*j))}


def golden(fun,left,right):
    r=(np.sqrt(5)-1)/2
    c=right-r*(right-left);d=left+r*(right-left)
    fc=fun(c);fd=fun(d)
    for _ in range(110):
        if fc<fd:
            right,d,fd=d,c,fc;c=right-r*(right-left);fc=fun(c)
        else:
            left,c,fc=c,d,fd;d=left+r*(right-left);fd=fun(d)
    return (left+right)/2


def extrema(f,a,sf,sa,n=4096):
    theta=np.linspace(-np.pi/2,np.pi/2,n,endpoint=False)+np.pi/(2*n)
    values=q_value(np.tan(theta),f,a,sf,sa)
    loc=np.where((values<np.roll(values,1))&(values<np.roll(values,-1)))[0]
    out=[]
    for index in loc:
        center=theta[index];half=np.pi/n
        t=golden(lambda t:float(q_value(np.tan(t),f,a,sf,sa)),center-half,center+half)
        # Positive-domain illustrative examples stay far from the projective seam.
        m=np.tan(t)
        for _ in range(6):
            pp=power(m,f,a,sf,sa)
            m-=pp["Q_prime"]/(2*pp["negative_ell_second"])
        out.append(power(m,f,a,sf,sa))
    return sorted(out,key=lambda x:x["Q"])


def tuned_example():
    sf=np.array([.1,3.,.4]);sa=np.ones(3)
    def at(t):
        a=np.array([1.,np.sqrt(t),0.]);f=np.array([.1,10*np.sqrt(t),.4*np.sqrt(.3)])
        modes=sorted(extrema(f,a,sf,sa),key=lambda x:x["m"])
        if len(modes)!=2:
            raise AssertionError((t,modes))
        return modes,f,a
    low,high=.02,.16
    # Q(low-m mode)-Q(high-m mode) increases with the second reading's weight.
    for _ in range(65):
        middle=(low+high)/2
        modes,_,_=at(middle)
        if modes[0]["Q"]<modes[1]["Q"]:
            low=middle
        else:
            high=middle
    tie=(low+high)/2
    modes,f,a=at(tie)
    z,w=leggauss(1500);theta=(z+1)*np.pi/4;w=w*np.pi/4
    qq=q_value(np.tan(theta),f,a,sf,sa);qq-=min(qq)
    law_jumps=[]
    for scale in (1.,5.,10.):
        probability=[w*np.exp(-.5*row["beta"]*scale*scale*qq) for row in modes]
        probability=[p/sum(p) for p in probability]
        law_jumps.append({"common_observation_amplitude_scale":scale,
                          "law_TV_difference_between_tied_mode_beta_choices":float(.5*sum(abs(probability[0]-probability[1])))})
    # The second derivative transforms covariantly at each stationary point.
    chart_checks=[]
    for row in modes:
        theta0=np.arctan(row["m"]);step=1e-4
        qm=float(q_value(np.tan(theta0-step),f,a,sf,sa))
        qp=float(q_value(np.tan(theta0+step),f,a,sf,sa))
        hessian_theta=(qp-2*row["Q"]+qm)/(2*step*step)
        j_theta=row["J_hat"]*(1+row["m"]**2)**2
        chart_checks.append({"m":row["m"],"beta_analytic_m_chart":row["beta"],
                             "beta_finite_difference_theta_chart":hessian_theta/j_theta})
    assert abs(modes[0]["Q"]-modes[1]["Q"])<1e-12
    assert abs(modes[0]["beta"]-modes[1]["beta"])>.05
    perturbations=[]
    for delta in (-1e-4,-1e-6,-1e-8,0.,1e-8,1e-6,1e-4):
        mm,ff,aa=at(tie*(1+delta))
        chosen=min(mm,key=lambda x:x["Q"])
        perturbations.append({"relative_data_weight_perturbation":delta,"chosen":chosen,"modes":mm})
    return {"force_sd":sf.tolist(),"acceleration_sd":sa.tolist(),"F":f.tolist(),"a":a.tolist(),
            "second_reading_squared_amplitude_tie":tie,"tied_modes":modes,"perturbations":perturbations,
            "positive_domain_law_jump_uniform_angle_reference":law_jumps,"chart_checks":chart_checks,
            "global_check":polynomial_stationary_points(f,a,sf,sa)}


def search_beta_above_one():
    rng=np.random.default_rng(280930)
    best=None
    for attempt in range(400):
        sf=np.exp(rng.uniform(-2.5,2.5,3));sa=np.ones(3)
        f=rng.normal(size=3)*sf*3;a=rng.normal(size=3)*3
        modes=extrema(f,a,sf,sa,1024)
        selected=modes[0]
        if best is None or selected["beta"]>best["selected"]["beta"]:
            best={"attempt":attempt,"force_sd":sf.tolist(),"acceleration_sd":sa.tolist(),
                  "F":f.tolist(),"a":a.tolist(),"selected":selected,"all_local_minima":modes}
        if selected["beta"]>1.2:
            break
    return best


def exact_quartic_example():
    """An exact full-covariance failure: curvature zero at a strict global peak.

    Take two scalar readings, sF=(sqrt(2)+1, sqrt(2)-1), sA=(1,1),
    F_i=k*sF_i and a_i=k. Put u=log(m), alpha=log(sqrt(2)+1).
    Up to a u-independent constant,

      ell(u) = k^2/2 [sech(u-alpha)+sech(u+alpha)]
             = sqrt(2)*k^2*cosh(u)/(cosh(u)^2+1).

    For t=cosh(u)>=1, ell(0)-ell(u)
       = k^2/sqrt(2) * (t-1)^2/(t^2+1).
    Thus u=0 is the unique global positive-mass peak. Negative m gives
    negative cross terms, so it cannot defeat this unrestricted peak.
    The peak is strictly quartic: ell''(0)=0 and ell''''(0)=-3k^2/sqrt(2).
    The denominator J in the u chart at m=1 is positive:
       J = k^2*(1+1/sqrt(2))/4.
    Therefore beta=-ell''/J=0 at every k>0. The proposed law equals its
    reference prior although the profile concentrates at rate k^(-1/2)
    in u as k increases. This is not a removable tie or flat-profile case.
    """
    sf=np.array([np.sqrt(2)+1,np.sqrt(2)-1]);sa=np.ones(2)
    rows=[]
    for k in (1.,10.,100.):
        f=k*sf;a=k*sa
        check=power(1.,f,a,sf,sa)
        analytic_j=k*k*(1+1/np.sqrt(2))/4
        assert abs(check['negative_ell_second'])<1e-8
        assert abs(check['J_hat']/analytic_j-1)<1e-14
        check.update({'k':k,'analytic_J_u':float(analytic_j),
                      'analytic_ell_fourth_u':float(-3*k*k/np.sqrt(2)),
                      'analytic_beta':0.})
        rows.append(check)
    u=np.linspace(-4,4,401);t=np.cosh(u)
    gap=(t-1)**2/(np.sqrt(2)*(t*t+1))
    observed_gap=(q_value(np.exp(u),sf,sa,sf,sa)-q_value(1.,sf,sa,sf,sa))/2
    error=float(np.max(abs(gap-observed_gap)))
    assert error<1e-14
    return {'force_sd':sf.tolist(),'acceleration_sd':sa.tolist(),
            'data_definition':'F=k*force_sd, a=k*ones(2)',
            'rows':rows,'max_global_gap_identity_error':error,
            'conclusion':'Strict unique quartic global peak gives beta=0 for arbitrarily strong observations.'}


def main():
    sf=np.array([1.7,1.2,.13]);sa=np.ones(3)
    f=np.array([2.5,-2.3,.38]);a=np.array([.1,4.1,1.9])
    rounded={"force_sd":sf.tolist(),"acceleration_sd":sa.tolist(),"F":f.tolist(),"a":a.tolist(),
             "all_local_minima":extrema(f,a,sf,sa)}
    assert rounded["all_local_minima"][0]["beta"]>1
    out={"source_commit":"b4bfc2040031ef010416334d0c5f4aa2ad05cbd8",
         "tie_example":tuned_example(),"largest_beta_found":search_beta_above_one(),
         "rounded_beta_above_one_example":rounded,"exact_quartic_example":exact_quartic_example()}
    found=out['largest_beta_found']
    found['global_check']=polynomial_stationary_points(np.array(found['F']),np.array(found['a']),
                                                      np.array(found['force_sd']),np.array(found['acceleration_sd']))
    Path(__file__).with_name("profile_power_ties_results.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
