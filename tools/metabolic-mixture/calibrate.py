"""Noise sensitivity near an independently fitted solution; no global-search claim."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from prototype import Model,generate,POPS,HERE

def run(seed):
    parent=HERE/'draft/calibration'/str(seed)
    generate(parent/'data',seed)
    m=Model(parent/'data')
    observed=np.array([float(r['fraction']) for r in m.meas]);sd=np.array([float(r['sd']) for r in m.meas])
    initial=np.load(HERE/'draft/reference/parameters.npy')
    def residual(p):
        decoded=m.decode(p)
        r=(m.predict(decoded)-observed)/sd
        rates=[(sum(m.weights[pop]*decoded[pop][0][row['reaction_id']] for pop in POPS)-float(row['rate']))/float(row['sd']) for row in m.rates if float(row['sd'])>0]
        penalty=[1e3*min(v,0) for pop in POPS for rid,v in decoded[pop][0].items() if not m.rev[rid]]
        return np.r_[r,rates,penalty]
    nf=len(m.flux.free);ne=len(m.flux.exch);nu=len(m.unknown['A'])
    lo=np.tile(np.r_[np.full(nf,1e-5),np.zeros(ne),np.full(nu,.5)],2)
    hi=np.tile(np.r_[np.full(nf,300),np.full(ne,500),np.full(nu,200)],2)
    begun=time.monotonic()
    fit=least_squares(residual,initial,bounds=(lo,hi),x_scale='jac',diff_step=1e-4,max_nfev=80,ftol=2e-7,xtol=2e-7,gtol=1e-5)
    decoded=m.decode(fit.x);truth=json.loads((parent/'injected_truth.json').read_text())
    errors={}
    for pop in POPS:
        net,exchange,pools=decoded[pop]
        for rid,tol in [('G6PDH',.45),('PDH',.12),('PC',.12),('ME',.08),('CS',.10),('IDH',.20)]:
            got=exchange[rid]+max(-net[rid],0) if rid=='IDH' else net[rid]
            want=truth[pop]['exchange'][rid] if rid=='IDH' else truth[pop]['net'][rid]
            errors[pop+'/'+rid]={'relative_error':abs(got-want)/abs(want),'tolerance':tol}
        for met in m.unknown[pop]:
            want=truth[pop]['pools'][met]
            errors[pop+'/'+met]={'relative_error':abs(pools[met]-want)/want,'tolerance':.25}
    report={'seed':seed,'ssr':float(2*fit.cost),'elapsed_sec':time.monotonic()-begun,'errors':errors,'all_targets_within_tolerance':all(v['relative_error']<=v['tolerance'] for v in errors.values()),'initialization':'reference estimate from a different noise realization; no injected parameters used for initialization'}
    (parent/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('calibration',seed,'SSR',report['ssr'],'targets',report['all_targets_within_tolerance'],flush=True)
    return report

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,range(1000,1004)))
    (HERE/'draft/calibration_results.json').write_text(json.dumps(results,indent=2)+'\n')
