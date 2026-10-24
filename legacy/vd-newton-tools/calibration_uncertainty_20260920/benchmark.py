from pathlib import Path
import sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'03_trace_and_evaluator/mass_estimator/investigation'))
from calibration_posterior import infer_batch
rng=np.random.default_rng(20260920)
for nu in (9,87):
    sf=.5*np.sqrt(rng.chisquare(nu,64)/nu)
    sa=.5*np.sqrt(rng.chisquare(nu,64)/nu)
    for f,a in ((0,0),(.5,2),(8,8),(32,8)):
        F=rng.normal(0,.5,(64,3))+[f,0,0]
        A=rng.normal(0,.5,(64,3))+[a,0,0]
        start=time.monotonic()
        low=infer_batch(F,A,force_sd=sf,acceleration_sd=sa,degrees=nu)
        elapsed=time.monotonic()-start
        high=infer_batch(F,A,force_sd=sf,acceleration_sd=sa,degrees=nu,order=192,calibration_order=48)
        print(dict(nu=nu,f=f,a=a,seconds64=elapsed,max_relative_point=max(abs(low['points']['calibrated_joint']/high['points']['calibrated_joint']-1)),max_interval_diff=max(float(np.max(abs(low['distributions']['calibrated_joint'][key]-high['distributions']['calibrated_joint'][key]))) for key in ('log_lower_95','log_upper_95'))),flush=True)
