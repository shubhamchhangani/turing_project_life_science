"""Author-only, AFTER fitting: compare estimates with separately held simulation truth.

Never import this module into an inverse solver. These diagnostics do not define
acceptance tolerances and do not prove statistical identifiability.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from model import coordinates,kinetics
from likelihood import load,decode

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);args=p.parse_args()
    r=args.root;c,_=load(r/'public');fit=json.loads((r/'crossmodal-fit.json').read_text())
    assert fit['complete'],'Do not audit an unfinished fit'
    t=np.load(r/'private/truth.npz',allow_pickle=False);truthxyz=coordinates(c,t['angles']);qt,_=kinetics(c,t['energies'],t['barriers'],300.)
    rows=[]
    for run in fit['runs']:
        a,e,b=decode(run['parameters']);xyz=coordinates(c,a)
        distance=np.sqrt(np.mean(np.sum((truthxyz[:,None,6:]-xyz[None,:,6:])**2,axis=-1),axis=-1))
        ii,jj=linear_sum_assignment(distance);q,_=kinetics(c,e,b,300.);q=q[np.ix_(jj,jj)];mask=~np.eye(4,dtype=bool)
        rows.append(dict(start=run['start'],mobile_domain_rmsd_angstrom=distance[ii,jj].tolist(),
            log_rate_rmse=float(np.sqrt(np.mean(np.log(q[mask]/qt[mask])**2))),
            max_rate_fold_error=float(np.exp(np.max(np.abs(np.log(q[mask]/qt[mask])))))))
    out=dict(scope='Post-fit diagnostic only. No truth used by inverse solver.',results=rows)
    (r/'crossmodal-diagnostics.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
