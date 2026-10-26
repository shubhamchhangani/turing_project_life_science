"""Mass-conserving transient isotopomer balances; contains no fitted answers."""
import csv
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
import forward_model as M

METS = list(M.NCARB)
OFF = {}; size = 0
for met in METS:
    OFF[met] = size
    size += 2 ** M.NCARB[met]
STARTS=np.array([OFF[m] for m in METS])
LENGTHS=np.array([2 ** M.NCARB[m] for m in METS])
REACTIONS = list(M.REACTIONS)
A=[]; B=[]; DEST=[]; EVENT=[]; EVENT_MET=[]; EVENT_RID=[]; EVENT_DIR=[]; EVENT_WEIGHT=[]
for event, (rid, direction, met, transition, weight) in enumerate(M.TRANSITIONS):
    ns=[2 ** M.NCARB[s] for s in transition.subs]
    aa=np.repeat(np.arange(ns[0]),ns[1]) if len(ns)==2 else np.arange(ns[0])
    bb=np.tile(np.arange(ns[1]),ns[0])+OFF[transition.subs[1]] if len(ns)==2 else np.full(ns[0],-1)
    A.extend(aa+OFF[transition.subs[0]]); B.extend(bb)
    DEST.extend(transition.map+OFF[met]); EVENT.extend([event]*len(aa))
    EVENT_MET.append(met); EVENT_RID.append(rid); EVENT_DIR.append(direction); EVENT_WEIGHT.append(weight)
A=np.array(A);B=np.array(B);DEST=np.array(DEST);EVENT=np.array(EVENT)

@njit(cache=True)
def rhs_fast(y, coeff, loss):
    dy=-loss*y
    # Normalize substrate probability vectors when forming reaction products.
    # Exact balances conserve unit sums; this suppresses floating-point drift
    # in repeated bimolecular convolution without changing those balances.
    prob=y.copy()
    for k in range(len(STARTS)):
        start=STARTS[k];n=LENGTHS[k];total=0.0
        for q in range(start,start+n): total+=y[q]
        if total > 0:
            for q in range(start,start+n):prob[q]/=total
    for j in range(len(A)):
        value=prob[A[j]]
        if B[j]>=0: value*=prob[B[j]]
        dy[DEST[j]]+=coeff[EVENT[j]]*value
    return dy


def initial():
    return np.concatenate([M.natural_dist(M.NCARB[m]) for m in METS])


def sources(row):
    return {'Glc':M.positional_dist([float(row[f'glucose_C{i}']) for i in range(1,7)]),
            'Gln':M.positional_dist([float(row[f'glutamine_C{i}']) for i in range(1,6)]),
            'CO2':M.positional_dist([float(row['co2_13c'])])}


def simulate(net, exchange, pools, schedule, times, rtol=2e-7, atol=1e-10, method='DOP853'):
    vd=M.directional_fluxes(net,exchange)
    coeff=np.array([vd[(r,d)]*w/pools[m] for r,d,w,m in zip(EVENT_RID,EVENT_DIR,EVENT_WEIGHT,EVENT_MET)])
    loss=np.zeros(size)
    for c,m in zip(coeff,EVENT_MET): loss[OFF[m]:OFF[m]+2**M.NCARB[m]]+=c
    assert np.all(np.isfinite(coeff)) and np.all(coeff>=0)
    y=initial(); found={}
    schedule=sorted(schedule,key=lambda r:float(r['start_min']))
    end=max(times)
    for k,row in enumerate(schedule):
        start=float(row['start_min'])
        if start>=end: break
        stop=min(float(schedule[k+1]['start_min']) if k+1<len(schedule) else end,end)
        for met,d in sources(row).items(): y[OFF[met]:OFF[met]+len(d)]=d
        requested=[float(t) for t in times if start<t<=stop]
        eval_times=sorted(set(requested+[stop]))
        sol=solve_ivp(lambda t,z:rhs_fast(z,coeff,loss),(start,stop),y,t_eval=eval_times,rtol=rtol,atol=atol,method=method)
        if not sol.success or not np.isfinite(sol.y).all(): raise RuntimeError('transient integration failed')
        for t in requested:
            z=sol.y[:,eval_times.index(t)]
            found[t]={m:z[OFF[m]:OFF[m]+2**M.NCARB[m]] for m in METS}
        y=sol.y[:,-1]
    return found


def read_data(folder):
    folder=Path(folder)
    def rows(name):
        with (folder/name).open() as f:return list(csv.DictReader(f))
    pool_rows=rows('pool_sizes.csv')
    known={r['metabolite']:float(r['pool_size']) for r in pool_rows if r['pool_size']!='unknown'}
    unknown=[r['metabolite'] for r in pool_rows if r['pool_size']=='unknown']
    return known,unknown,rows('labeling_schedule.csv'),rows('mid_measurements.csv'),rows('extracellular_rates.csv')


def predictions(net,exchange,pools,schedules,measurements,**kwargs):
    exps=sorted({r['experiment'] for r in measurements}); pred={}
    for exp in exps:
        rows=[r for r in measurements if r['experiment']==exp]
        ts=sorted({float(r['time_min']) for r in rows})
        ds=simulate(net,exchange,pools,[r for r in schedules if r['experiment']==exp],ts,**kwargs)
        for t in ts:
            for f in M.FRAGMENTS:pred[exp,t,f]=M.raw_mid(ds[t],f)
    return np.array([pred[r['experiment'],float(r['time_min']),r['fragment']][int(r['isotopologue'][1:])] for r in measurements])
