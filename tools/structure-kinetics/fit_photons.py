"""Input-only multistart inference baseline. This does not yet use the cryo-EM images."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from likelihood import load,loglike

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--starts',type=int,default=3);p.add_argument('--iterations',type=int,default=150);args=p.parse_args()
    c,traces=load(args.data);rng=np.random.default_rng(7729);bounds=[]
    for _ in range(4):bounds.extend([c['parameter_bounds']['azimuth'],c['parameter_bounds']['elevation']])
    bounds+=([c['parameter_bounds']['energies']]*3+[c['parameter_bounds']['barriers']]*6)
    start=time.monotonic();runs=[];args.out.parent.mkdir(parents=True,exist_ok=True)
    for run in range(args.starts):
        initial=np.r_[rng.uniform([-2.,-.8],[2.,.8],(4,2)).ravel(),rng.uniform(-.2,.8,3),rng.uniform(3.,5.,6)]
        initial_ll=loglike(c,traces,initial)
        result=minimize(lambda x:-loglike(c,traces,x)/len(traces),initial,method='L-BFGS-B',bounds=bounds,
                        options={'maxiter':args.iterations,'ftol':1e-11,'gtol':1e-5,'maxls':30})
        row=dict(start=run,initial_log_likelihood=initial_ll,log_likelihood=-float(result.fun)*len(traces),
                 parameters=result.x.tolist(),converged=bool(result.success),message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev))
        runs.append(row)
        report=dict(status='photon_only_baseline',complete=(run==args.starts-1),seconds=time.monotonic()-start,runs=runs,
                    data_source=str(args.data),truth_used=False,images_used=False)
        temp=args.out.with_suffix('.tmp');temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(args.out)
        print(json.dumps({k:v for k,v in row.items() if k!='parameters'}),flush=True)
