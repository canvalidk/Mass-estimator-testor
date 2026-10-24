"""Compare independent local-information derivations and existing code."""
from pathlib import Path
import json
import math
import sys
import numpy as np

WORKSPACE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(WORKSPACE/'03_trace_and_evaluator'/'mass_estimator'))
from practical_formulas import aligned_local_uncertainty


def main():
    rng=np.random.default_rng(20260915)
    differences=[]
    for dim in (1,2,3):
        for _ in range(20):
            m=math.exp(rng.uniform(-2,2))
            a=rng.normal(size=dim)
            raw=rng.normal(size=(2*dim,2*dim))
            sigma=raw@raw.T+np.eye(2*dim)
            ja=np.vstack((m*np.eye(dim),np.eye(dim)))
            jm=np.r_[a,np.zeros(dim)]
            information=ja.T@np.linalg.solve(sigma,ja)
            cross=ja.T@np.linalg.solve(sigma,jm)
            schur=float(jm@np.linalg.solve(sigma,jm)-cross@np.linalg.solve(information,cross))
            v=sigma[:dim,:dim]+m*m*sigma[dim:,dim:]-m*(sigma[:dim,dim:]+sigma[dim:,:dim])
            contrast=float(a@np.linalg.solve(v,a))
            implemented=aligned_local_uncertainty(m*a,a,sigma)
            differences.extend((abs(schur/contrast-1),abs(implemented*math.sqrt(contrast)-1)))
            assert abs(schur/contrast-1)<1e-11
            assert abs(implemented*math.sqrt(contrast)-1)<1e-11
            sf,sa=np.exp(rng.uniform(-2,2,size=2))
            iso=np.diag([sf*sf]*dim+[sa*sa]*dim)
            info=float(a@a)/(sf*sf+m*m*sa*sa)
            actual=aligned_local_uncertainty(m*a,a,iso)
            assert abs(actual*math.sqrt(info)-1)<1e-11
            rho_f=m*np.linalg.norm(a)/sf
            rho_a=np.linalg.norm(a)/sa
            assert abs(m*m*info*(rho_f**-2+rho_a**-2)-1)<1e-12
    result=dict(status='All assertions passed',random_covariances=60,
                maximum_relative_difference=max(differences),
                scope='Direct Fisher Schur complement, residual-contrast formula, existing aligned local uncertainty, isotropic SNR formula')
    Path(__file__).with_name('local_information_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
