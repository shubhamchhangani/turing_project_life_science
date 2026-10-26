import sys,copy,time,json,csv
from pathlib import Path
import numpy as np
sys.path.insert(0,'tools/metabolic-flux-analysis-v2')
import dynamic_model as D
import make_dataset as G
M=D.M
# A second, vector-based implementation of the transient RHS; no flattened
# index arrays or compiled summation kernel are used here.
def rhs_independent(y,coeff,loss):
    d={m:y[D.OFF[m]:D.OFF[m]+2**M.NCARB[m]] for m in D.METS}
    d={m:v/v.sum() for m,v in d.items()}
    out=-loss*y
    for k,(_,_,m,tr,_) in enumerate(M.TRANSITIONS):
        out[D.OFF[m]:D.OFF[m]+2**M.NCARB[m]]+=coeff[k]*tr.apply(d)
    return out
r=Path('tasks/life-sciences/biology/metabolic-flux-analysis-v2/environment/data')
known,unknown,schedule,meas,rates=D.read_data(r)
start=time.monotonic()
base=D.predictions(G.T.NET,G.T.EXCH,G.POOLS,schedule,meas,rtol=2e-10,atol=1e-13)
fast=D.rhs_fast
D.rhs_fast=rhs_independent
ind=D.predictions(G.T.NET,G.T.EXCH,G.POOLS,schedule,meas,rtol=1e-9,atol=1e-12,method='RK45')
D.rhs_fast=fast
error=float(np.max(np.abs(ind-base)))
print('Independent RHS / RK45 maximum MID error',error,'seconds',time.monotonic()-start,flush=True)
assert error<1e-7
observed=np.array([float(r['fraction']) for r in meas]);sd=np.array([float(r['sd']) for r in meas])
results={'independent_max_mid_difference':error}
# Incorrect schedule ignores the chase and leaves the tracer on indefinitely.
wrong_schedule=[x for x in schedule if float(x['start_min'])==0]
pred=D.predictions(G.T.NET,G.T.EXCH,G.POOLS,wrong_schedule,meas)
results['ignore_chase_ssr']=float(np.sum(((pred-observed)/sd)**2))
wrong_pools={**G.POOLS,**{m:1 for m in G.UNKNOWN}}
pred=D.predictions(G.T.NET,G.T.EXCH,wrong_pools,schedule,meas)
results['unit_unknown_pools_ssr']=float(np.sum(((pred-observed)/sd)**2))
# Replace all time points by the isotopic steady-state prediction for the pulse.
ss={exp:{f:M.raw_mid(d,f) for f in M.FRAGMENTS} for exp in ['glc12','gln5'] for d in [M.simulate(G.T.NET,G.T.EXCH,M.sources_for(exp))]}
pred=np.array([ss[x['experiment']][x['fragment']][int(x['isotopologue'][1:])] for x in meas])
results['steady_state_ssr']=float(np.sum(((pred-observed)/sd)**2))
print(json.dumps(results,indent=2),flush=True)
Path('tools/metabolic-flux-analysis-v2/physics_checks.json').write_text(json.dumps(results,indent=2)+'\n')
