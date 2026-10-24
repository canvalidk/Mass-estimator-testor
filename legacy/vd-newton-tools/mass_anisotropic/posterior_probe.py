import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator/investigation'))
sys.path.insert(0,str(ROOT/'03_trace_and_evaluator/mass_estimator'))
from comparison_diagonal_posterior import infer_batch, _infer_chunk
from comparison_posterior import infer_batch as iso
from compare_estimators import refinement_check
from estimator import estimate

cases=[([2.,-1.,.5],[.8,.1,.2],[.7,1,1.4],[1.2,.9,.8]),
       ([0.,0.,0.],[0.,0.,0.],[.5,1,2],[2,1,.5])]
for f,a,sf,sa in cases:
    result=infer_batch(f,a,force_sd=sf,acceleration_sd=sa)
    refs=[estimate(f,a,np.diag(np.r_[np.array(sf)**2,np.array(sa)**2]),direction_order=n,ratio_order=256).mass for n in (24,48)]
    print('fullcov',result['points']['flat_joint'],refs)
f=np.array([[0.,0.,0.],[2.,-1.,3.],[8.,0.,0.],[32.,0.,0.]])
a=np.array([[0.,0.,0.],[.1,.7,-.4],[-8.,0.,0.],[-32.,0.,0.]])
direct=_infer_chunk(f,a,np.ones(3),np.ones(3),2.,96)
reference=iso(f,a,2.,96)
for name,point in direct['points'].items():
    print('iso',name,point,reference['points'][name])
profiles=[([.5,1,2],[2,1,.5]),([.25,.25,4],[4,.25,.25]),([.25,1,4],[1,4,.25])]
direction=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1]])
direction=direction/np.linalg.norm(direction,axis=1)[:,None]
rng=np.random.default_rng(20260919)
for sf,sa in profiles:
    sf,sa=np.array(sf),np.array(sa)
    f=np.vstack([direction*s for s in (.25,2,8,16,32)])
    a=np.vstack([direction*s for s in (8,8,8,8,8)])
    f+=rng.normal(size=f.shape)*sf
    a+=rng.normal(size=a.shape)*sa
    start=time.perf_counter()
    coarse=infer_batch(f,a,force_sd=sf,acceleration_sd=sa,order=96)
    fine=infer_batch(f,a,force_sd=sf,acceleration_sd=sa,order=192)
    print('refine',sf,sa,json.dumps(refinement_check(coarse,fine)),time.perf_counter()-start)
