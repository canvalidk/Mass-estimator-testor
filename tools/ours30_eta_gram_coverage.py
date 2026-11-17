"""B22: paired Gaussian coverage of Bartlett, Gram and raw tempering.

Run: python tools/ours30_eta_gram_coverage.py --series 300 --workers 4
     --output results/b22-coverage.json

Seed 2211; independent SeedSequence [2211,k,N,signal_index,replicate].
The required 24 unique cells are augmented by null N=80,200 in 1D and 3D.
All six variants use each same simulated series. Truth c=s=m=1, unit
Gaussian channel noise, fixed total reading signal in a uniform direction.
Width means ln(upper/lower), not a width in mass units. Output contains
binomial Wilson intervals; flags outside [92%,98%] are empirical, not tests.

The study-only fast path reduces one isotropic instrument to exact 2x2
Grams. It uses full.py's adaptive-grid constants, half-Cauchy prior and CDF
inversion. The default run validates it against all six full APIs first;
--validate-only checks this separately. No change to production numerics.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from met import full  # noqa: E402

SIGNALS = (0., .5, 1., 2., 5.)
MODELS = tuple(f'{base}_{eta}' for base in ('ours30', 'ours29')
               for eta in ('bartlett', 'gram', 'raw'))


def cells():
    required = {(k,n,r) for k in (1,3) for n in (5,20) for r in SIGNALS}
    required |= {(3,n,r) for n in (80,200) for r in (1.,2.)}
    extra = {(k,n,0.) for k in (1,3) for n in (80,200)}
    return sorted(required | extra)


def world(k, n, signal, replicate):
    rng = np.random.default_rng([2211,k,n,SIGNALS.index(signal),replicate])
    direction = rng.normal(size=(n,k))
    direction /= np.linalg.norm(direction,axis=1)[:,None]
    latent = signal/math.sqrt(2)*direction
    return latent+rng.normal(size=(n,k)), latent+rng.normal(size=(n,k))


def summaries(force, acceleration, positivity):
    a = np.einsum('ij,ij->i',force,force)
    b = np.einsum('ij,ij->i',acceleration,acceleration)
    cross = np.einsum('ij,ij->i',force,acceleration)
    reading_gap = np.hypot(a-b,2*cross)
    q = np.maximum((a+b-reading_gap)/2,0.)
    q_used = np.where(cross >= 0,q,np.minimum(a,b)) if positivity else q
    numerator = force.shape[1]-(.5 if positivity else 1.)
    g = np.ones(len(a))
    if numerator > 0:
        np.divide(numerator,q_used,out=g,where=q_used>0)
        g = np.minimum(g,1.)
    aa, bb, cc = float(g@a),float(g@b),float(g@cross)
    gap = math.hypot(aa-bb,2*cc)
    phi = .5*math.atan2(2*cc,bb-aa)
    score = g*(.5*(a-b)*math.sin(2*phi)+cross*math.cos(2*phi))
    j = float(score@score)
    bartlett = 1. if j <= 1e-300 or gap <= 0 else min(1.,gap/j)
    hmax = (aa+bb+gap)/2
    misfit = max(aa+bb-hmax,0.)
    gram = 1. if hmax <= 1e-300 else 1.-misfit/hmax
    return (aa,bb,cc), g, {'bartlett':bartlett,'gram':gram,'raw':1.}


def interval(gram, eta, refinement=1):
    aa,bb,cc = gram
    def log_law(lm):
        phi = np.arctan(np.exp(lm))
        sn,cs = np.sin(phi),np.cos(phi)
        return full._log_prior(lm,1.)+.5*eta*(aa*sn*sn+bb*cs*cs+2*cc*sn*cs)
    coarse = np.linspace(-full.COARSE_HALF,full.COARSE_HALF,full.COARSE_NODES)
    lp = log_law(coarse)
    keep = np.flatnonzero(lp > lp.max()-full.SUPPORT_NATS)
    step = coarse[1]-coarse[0]
    grid = np.linspace(coarse[keep[0]]-2*step,coarse[keep[-1]]+2*step,
                       (full.FINE_NODES-1)*refinement+1)
    quantiles,_ = full._quantiles(grid,log_law(grid),[.025,.5,.975])
    return quantiles


def fast_results(force, acceleration, refinement=1):
    output = {}
    for base,positive in [('ours30',True),('ours29',False)]:
        gram,g,etas = summaries(force,acceleration,positive)
        for name,eta in etas.items():
            q = interval(gram,eta,refinement)
            output[f'{base}_{name}'] = (q,eta,g)
    return output


def validate():
    maxima = dict(eta=0.,g=0.,log_quantile=0.,refinement_log_quantile=0.)
    fixtures = [(1,5,0.,0),(1,20,.5,1),(3,5,2.,2),(3,20,5.,3),
                (3,80,1.,4),(3,200,2.,5),(1,200,0.,6)]
    for k,n,r,replicate in fixtures:
        force,acceleration = world(k,n,r,replicate)
        readings = [full.isotropic(f,a,1.,1.) for f,a in zip(force,acceleration)]
        fast = fast_results(force,acceleration)
        refined = fast_results(force,acceleration,refinement=2)
        for model,(q,eta,g) in fast.items():
            base,method = model.split('_')
            options = dict(readout=False,calibrate=method!='raw',
                           eta_form='gram' if method=='gram' else 'bartlett')
            actual = (full.ours30 if base=='ours30' else full.ours29)(readings,1.,**options)
            expected = np.array([actual['interval'][0],actual['median'],actual['interval'][1]])
            maxima['eta'] = max(maxima['eta'],abs(eta-actual['eta']))
            maxima['g'] = max(maxima['g'],float(np.max(np.abs(g-actual['g']))))
            maxima['log_quantile'] = max(maxima['log_quantile'],float(np.max(np.abs(np.log(q/expected)))))
            maxima['refinement_log_quantile'] = max(maxima['refinement_log_quantile'],
                float(np.max(np.abs(np.log(q/refined[model][0])))))
            np.testing.assert_allclose(g,actual['g'],rtol=1e-9,atol=2e-10)
            np.testing.assert_allclose(eta,actual['eta'],rtol=1e-9,atol=2e-10)
            np.testing.assert_allclose(np.log(q),np.log(expected),rtol=0,atol=2e-8)
            np.testing.assert_allclose(np.log(q),np.log(refined[model][0]),rtol=0,atol=2e-6)
    return dict(fixtures=len(fixtures),api_comparisons=6*len(fixtures),max_errors=maxima)


def wilson(successes,n):
    z=1.959963984540054
    p=successes/n
    centre=(p+z*z/(2*n))/(1+z*z/n)
    radius=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [centre-radius,centre+radius]


def run_cell(job):
    k,n,r,series = job
    values={name:[] for name in MODELS}
    for replicate in range(series):
        force,acceleration=world(k,n,r,replicate)
        results=fast_results(force,acceleration)
        for name,(q,eta,_) in results.items():
            covered=q[0] <= 1. <= q[2]
            # Resolve any near-boundary classification with a doubled grid.
            if min(abs(math.log(q[0])),abs(math.log(q[2]))) < 2e-5:
                base,method=name.split('_')
                gram,_,etas=summaries(force,acceleration,base=='ours30')
                refined=interval(gram,etas[method],refinement=2)
                assert covered == (refined[0] <= 1. <= refined[2]), (k,n,r,replicate,name)
            values[name].append((covered,math.log(q[2]/q[0]),eta))
    summary={}
    for name,rows in values.items():
        arr=np.asarray(rows)
        hits=int(arr[:,0].sum())
        cover=hits/series
        summary[name]=dict(covered=hits,coverage=cover,coverage_wilson95=wilson(hits,series),
                          median_log_width=float(np.median(arr[:,1])),mean_eta=float(arr[:,2].mean()),
                          gram_outside_92_98=bool(name.endswith('_gram') and (cover<.92 or cover>.98)))
    comparison={}
    for base in ('ours30','ours29'):
        bartlett=np.asarray(values[base+'_bartlett'])
        gram=np.asarray(values[base+'_gram'])
        comparison[base]=dict(coverage_disagreements=int(np.count_nonzero(bartlett[:,0]!=gram[:,0])),
                             mean_absolute_log_width_difference=float(np.mean(np.abs(bartlett[:,1]-gram[:,1]))))
    return dict(k=k,n=n,signal=r,series=series,models=summary,paired_gram_vs_bartlett=comparison)


def print_cell(result):
    entries=[]
    for model in MODELS:
        m=result['models'][model]
        entries.append(f"{100*m['coverage']:.1f}/{m['median_log_width']:.3f}")
    print(f"{result['k']} {result['n']:3} {result['signal']:3.1f} | "+' | '.join(entries),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--series',type=int,default=300)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--output',type=Path,default=Path('results/b22-coverage.json'))
    parser.add_argument('--validate-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.series<150 and not args.validate_only:
        parser.error('B22 requires at least 150 series per cell')
    validation=validate()
    print('Validated:',json.dumps(validation),flush=True)
    if args.validate_only:
        return
    args.output.parent.mkdir(parents=True,exist_ok=True)
    results=[]
    if args.resume and args.output.exists():
        saved=json.loads(args.output.read_text())
        if saved['series']!=args.series or saved['base_seed']!=2211:
            parser.error('resume seed/series mismatch')
        results=saved['cells']
    completed={(x['k'],x['n'],x['signal']) for x in results}
    jobs=[(*cell,args.series) for cell in cells() if cell not in completed]
    def save():
        payload=dict(base_seed=2211,series=args.series,truth=1.,prior_centre=1.,width='ln(hi/lo)',
                     required_cells=24,total_cells=28,validation=validation,
                     cells=sorted(results,key=lambda x:(x['k'],x['n'],x['signal'])))
        args.output.write_text(json.dumps(payload,indent=2)+'\n')
    print('k   N sig | '+' | '.join(MODELS)+' (coverage % / median log-width)',flush=True)
    if args.workers==1:
        for job in jobs:
            result=run_cell(job); results.append(result); print_cell(result); save()
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures=[executor.submit(run_cell,job) for job in jobs]
            for future in as_completed(futures):
                result=future.result(); results.append(result); print_cell(result); save()
    save()


if __name__=='__main__':
    main()
