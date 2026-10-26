"""Joint transient flux/pool fit from measured inputs, without true parameters."""
import sys,csv,json,time,os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import network_model as S
import dynamic_model as D

DATA=Path(sys.argv[1] if len(sys.argv)>1 else '/root/data')
OUT=Path(sys.argv[2] if len(sys.argv)>2 else '/root/results')
STARTS=int(sys.argv[3]) if len(sys.argv)>3 else 4
SEED=int(sys.argv[4]) if len(sys.argv)>4 else 17

def main():
    rows=list(csv.DictReader((DATA/'network.tsv').open(),delimiter='\t'))
    rxns=[r['reaction_id'] for r in rows];rev={r['reaction_id']:r['reversible']=='yes' for r in rows}
    fm=S.FluxModel(rxns,rev)
    known,unknown,schedule,meas,rates=D.read_data(DATA)
    observed=np.array([float(r['fraction']) for r in meas]);sd=np.array([float(r['sd']) for r in meas])
    rates={r['reaction_id']:(float(r['rate']),float(r['sd'])) for r in rates}
    nfree=len(fm.free);nx=len(fm.exch)
    def decode(p):
        net=fm.net(p[:nfree]);exchange=dict(zip(fm.exch,p[nfree:nfree+nx]));pools={**known,**dict(zip(unknown,p[nfree+nx:]))}
        return net,exchange,pools
    count=0
    def residual(p):
        nonlocal count
        count+=1
        net,exchange,pools=decode(p)
        physical={r:max(v,0) if not rev[r] else v for r,v in net.items()}
        pred=D.predictions(physical,exchange,pools,schedule,meas,rtol=2e-8,atol=1e-11)
        extras=[(net[r]-v)/s for r,(v,s) in rates.items() if s>0]
        penalty=[1e3*min(net[r],0) for r in rxns if not rev[r]]
        return np.r_[(pred-observed)/sd,extras,penalty]
    rng=np.random.default_rng(SEED);best=None
    lo=np.r_[np.full(nfree,1e-5),np.zeros(nx),np.full(len(unknown),.5)]
    hi=np.r_[np.full(nfree,300),np.full(nx,500),np.full(len(unknown),200)]
    begun=time.monotonic()
    for start in range(STARTS):
        z=[rates[r][0] if r in rates else rng.uniform(8,45) for r in fm.free]
        p=np.r_[z,rng.uniform(5,100,nx),rng.uniform(8,60,len(unknown))]
        p=np.clip(p,lo+1e-4,hi-1e-4)
        fit=least_squares(residual,p,bounds=(lo,hi),x_scale='jac',diff_step=1e-4,max_nfev=180,ftol=2e-7,xtol=2e-7,gtol=1e-5)
        print('start',start,'SSR',2*fit.cost,'nfev',fit.nfev,'sec',time.monotonic()-begun,flush=True)
        if best is None or fit.cost<best.cost:
            best=fit
            OUT.mkdir(parents=True,exist_ok=True)
            net,exchange,pools=decode(fit.x);vf,vb=fm.directional(net,[exchange[r] for r in fm.exch])
            with (OUT/'fluxes.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['reaction_id','net_flux','forward_flux','reverse_flux'])
                for r in rxns:w.writerow([r,f'{net[r]:.6f}',f'{vf[r]:.6f}',f'{vb.get(r,0):.6f}'])
            with (OUT/'pools.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['metabolite','pool_size'])
                for m in unknown:w.writerow([m,f'{pools[m]:.6f}'])
            print('unknown pools', {m:pools[m] for m in unknown},flush=True)
    print('total residual evaluations',count,flush=True)

if __name__=='__main__':main()
