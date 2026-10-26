"""Initialize geometry from images, then infer kinetics from photons, without truth.

This diagnostic is not a full joint maximum-likelihood solver. It does not use
poorly fitted images or cooling populations in the kinetic fit.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.optimize import least_squares,minimize
from model import coordinates
from likelihood import load,loglike

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--image-fit',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    c,traces=load(args.data);imagefit=json.loads(args.image_fit.read_text());rows=[r for r in imagefit['rows'] if r['reduced_chi_square']<1.5]
    assert len(rows)>=12,'Insufficient noise-level image fits to initialize four conformations'
    observed_angles=np.array([r['parameters'][5:7] for r in rows]);features=coordinates(c,observed_angles)[:,6:].reshape(len(rows),-1)
    centers,labels=kmeans2(features,4,minit='++',seed=219)
    angles=[];populations=[];scatter=[]
    for k in range(4):
        subset=features[labels==k];assert len(subset)>=3,'Unsupported geometry cluster'
        center=np.median(subset,axis=0);seed=observed_angles[np.argmin(np.linalg.norm(features-center,axis=1))]
        fit=least_squares(lambda a:coordinates(c,np.array([a]))[0,6:].ravel()-center,seed)
        angles.append(fit.x);populations.append(len(subset));scatter.append(float(np.sqrt(np.mean((subset-center)**2))))
    # State order is a reporting convention. Choose the most frequently recovered
    # geometry as state 0, while allowing negative relative well energies below.
    order=np.argsort(populations)[::-1];angles=np.array(angles)[order];rng=np.random.default_rng(819);runs=[];started=time.monotonic()
    bounds=[[-1.5,1.5]]*3+[[1.,7.5]]*6
    for attempt in range(3):
        initial=np.r_[rng.uniform(-.2,.8,3),rng.uniform(3.,5.,6)]
        result=minimize(lambda z:-loglike(c,traces,np.r_[angles.ravel(),z])/len(traces),initial,
                        method='L-BFGS-B',bounds=bounds,options={'maxiter':220,'ftol':1e-12,'gtol':1e-5})
        row=dict(start=attempt,parameters=np.r_[angles.ravel(),result.x].tolist(),log_likelihood=-float(result.fun)*len(traces),
                 converged=bool(result.success),iterations=int(result.nit),message=str(result.message))
        runs.append(row);out=dict(status='crossmodal_initialization_only',complete=attempt==2,truth_used=False,
           good_image_count=len(rows),cluster_sizes=np.array(populations)[order].tolist(),
           cluster_coordinate_scatter_angstrom=np.array(scatter)[order].tolist(),seconds=time.monotonic()-started,runs=runs,
           limitations=['Not full joint refinement','Poor image fits excluded','Cooling likelihood not yet included','No uncertainty calibration'])
        args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in row.items() if k!='parameters'}),flush=True)
