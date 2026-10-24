"""Numerical cross-checks for analytical readout identities, not sampling calibration."""
import json
import math
import random
from pathlib import Path

rng = random.Random(9052026)

def close(a, b, tol=1e-10):
    assert abs(a-b) <= tol*max(1,abs(a),abs(b)), (a,b)

def simpson(fn, lo, hi, n=40000):
    h=(hi-lo)/n
    return h/3*(fn(lo)+fn(hi)+sum((4 if i%2 else 2)*fn(lo+i*h) for i in range(1,n)))

for trial in range(1000):
    n=12
    raw=[rng.random()+.01 for _ in range(n)]
    p=[x/sum(raw) for x in raw]
    f=[math.exp(rng.uniform(-6,6)) for _ in p]
    a=[math.exp(rng.uniform(-6,6)) for _ in p]
    B=sum(w*x for w,x in zip(p,f))
    A=sum(w*x for w,x in zip(p,a))
    C=sum(w*math.sqrt(x*y) for w,x,y in zip(p,f,a))
    m=B/A
    loss=lambda r: sum(w*(math.sqrt(x)/r**.25-math.sqrt(y)*r**.25)**2 for w,x,y in zip(p,f,a))
    expected=B/math.sqrt(m)+A*math.sqrt(m)-2*C
    close(loss(m), expected)
    for shift in (-3,-.7,.7,3):
        z=m*math.exp(shift)
        assert loss(z)>loss(m)
        close(loss(z), 2*math.sqrt(A*B)*math.cosh(shift/2)-2*C)
    swapped=sum(w*(math.sqrt(y)/(1/m)**.25-math.sqrt(x)*(1/m)**.25)**2 for w,x,y in zip(p,f,a))
    close(loss(m),swapped)
    kappa=C/math.sqrt(A*B)
    assert 0<kappa<=1+1e-12
    qf=[w*x/B for w,x in zip(p,f)]
    qa=[w*x/A for w,x in zip(p,a)]
    close(kappa,sum(math.sqrt(x*y) for x,y in zip(qf,qa)))
    parts=[]
    for indices in (range(0,4),range(4,8),range(8,12)):
        pk=sum(p[i] for i in indices)
        ak=sum(p[i]*a[i] for i in indices)/pk
        fk=sum(p[i]*f[i] for i in indices)/pk
        parts.append((pk,ak,fk/ak))
    regrouped=sum(pk*ak*mk for pk,ak,mk in parts)/sum(pk*ak for pk,ak,mk in parts)
    close(regrouped,m)
    assert min(x[2] for x in parts)-1e-10<=m<=max(x[2] for x in parts)+1e-10
    # Independent unit rescalings preserve the selected physical mass.
    cf,ca=3.7,.23
    r2=(cf*B)/(ca*A)
    l2=sum(w*(math.sqrt(cf*x)/r2**.25-math.sqrt(ca*y)*r2**.25)**2 for w,x,y in zip(p,f,a))
    close(l2,math.sqrt(cf*ca)*loss(m))

c=math.sqrt(2/math.pi)
root_mean=simpson(lambda t:2*c*t*t*math.exp(-t**4/2),0,math.sqrt(12))
kappa0=root_mean**2/c
close(kappa0,math.gamma(.75)**2/math.sqrt(math.pi))
log_density=lambda t:1/(math.pi*math.cosh(t))
close(simpson(log_density,-40,40),1)
log_var=simpson(lambda t:t*t*log_density(t),-40,40)
close(log_var,math.pi**2/4)

out={"status":"passed","random_positive_joint_laws":1000,
"identities":["unique loss minimum","channel exchange","unit covariance","mixture aggregation","Hellinger affinity"],
"null_affinity":kappa0,"null_Hellinger_squared":1-kappa0,
"null_log_ratio_variance":log_var,"null_log_ratio_sd":math.sqrt(log_var),
"limits":"Finite checks supplement the analytical proofs; no admission threshold or sampling coverage calibrated."}
Path(__file__).with_name('check-results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
