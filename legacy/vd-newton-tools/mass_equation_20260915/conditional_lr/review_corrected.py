"""Recheck corrected CLR numerics; preserve all pre-correction evidence."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import review_conditional as previous

cc=previous.cc

def main():
    cases=[]
    for c in np.geomspace(1e-8,1e12,41):
        for s in np.r_[0.,np.geomspace(1e-8,1e12,41)]:
            ref=previous.closed3(float(c),float(s))
            got=float(cc.conditional_cdf(c,s,3))
            cases.append((abs(got-ref),float(c),float(s),got,ref))
    ordinary=[]
    for k in (1,2,3,24,192,1536,12288):
        q1,qk=cc.chi_quantile(1,.95),cc.chi_quantile(k,.95)
        for L in (qk*(1+1e-8),qk*1.1,2*k+10,10*k+10,1e5,1e8,1e12):
            if L<=qk: continue
            b=float(cc.fitted_threshold(L,k))
            c=L-b
            got=float(cc.conditional_cdf(c,b,k,256))
            ref=previous.closed3(c,b) if k==3 else None
            ordinary.append(dict(k=k,L=L,c=c,cdf_refined=got,closed3=ref,
                                 float_spacing=float(np.spacing(L))))
    old=json.loads(Path(__file__).with_name('review_results.json').read_text())
    critical=[]
    for item in old['conditional_mc']:
        k,s,c=item['k'],item['s'],item['c']
        got=float(cc.conditional_cdf(c,s,k))
        critical.append(dict(k=k,s=s,c=c,cdf=got,error=abs(got-.95)))
    # Existing full-simulation summary is only read, never rewritten.
    result={'script_sha256':hashlib.sha256(Path(cc.__file__).read_bytes()).hexdigest(),
            'k3_grid_count':len(cases),'k3_worst_cases':sorted(cases,reverse=True)[:12],
            'previous_mc_critical_checks':critical,'threshold_grid':ordinary}
    out=Path(__file__).with_name('review_corrected_results.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
