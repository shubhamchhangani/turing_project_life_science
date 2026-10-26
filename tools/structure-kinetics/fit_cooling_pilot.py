"""Coupled kinetic objective with photons and all images at input-estimated poses.

Geometry and poses stay fixed at inferred values: this is not full joint
refinement and does not marginalize orientation uncertainty.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.spatial.transform import Rotation
from model import coordinates,projection,frozen
from likelihood import load,loglike

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);a=parser.parse_args();r=a.root
    c,traces=load(r/'public');data=np.load(r/'public/images.npz',allow_pickle=False)
    assigned=json.loads((r/'image-assignments-complete.json').read_text());assert assigned['complete']
    assert len(assigned['rows'])==len(data['images']) and all(row['reduced_chi_square']<1.5 for row in assigned['rows'])
    fits=json.loads((r/'crossmodal-fit.json').read_text());best=max(fits['runs'],key=lambda row:row['log_likelihood']);x=np.array(best['parameters']);fixed=x[:8];xyz=coordinates(c,fixed.reshape(4,2))
    chisq=[];protocol=[]
    for row in assigned['rows']:
        i=row['image'];pose=min(row['fits'],key=lambda v:v['chi_square'])['pose'];rotation=Rotation.from_rotvec(pose[:3]).as_matrix()
        chisq.append([float(np.sum(((projection(c,xx,rotation,pose[3:],float(data['defocus'][i]))-data['images'][i])/c['image_noise_sd'])**2)) for xx in xyz]);protocol.append(int(data['protocol'][i]))
    chisq=np.asarray(chisq);protocol=np.asarray(protocol);started=time.monotonic()
    def components(z):
        weights=np.array([frozen(c,z[:3],z[3:],d) for d in c['cooling_seconds']])
        image_ll=float(np.sum(logsumexp(-chisq/2+np.log(weights[protocol]),axis=1)))
        photon_ll=loglike(c,traces,np.r_[fixed,z])
        return photon_ll,image_ll,weights
    initial=components(x[8:]);bounds=[[-1.5,1.5]]*3+[[1.,7.5]]*6
    result=minimize(lambda z:-sum(components(z)[:2])/len(traces),x[8:],bounds=bounds,method='L-BFGS-B',
        options={'maxiter':100,'ftol':1e-12,'gtol':1e-5})
    final=components(result.x)
    out=dict(status='coupled_kinetic_pilot_completed',complete=True,truth_used=False,parameters=np.r_[fixed,result.x].tolist(),
        photon_log_likelihood=final[0],image_log_likelihood_without_constant=final[1],frozen_populations=final[2].tolist(),
        joint_log_likelihood_gain=float(sum(final[:2])-sum(initial[:2])),converged=bool(result.success),
        iterations=int(result.nit),seconds=time.monotonic()-started,
        limitations=['Geometry and poses fixed at input-derived estimates','No orientation marginalization','No calibrated parameter confidence intervals','No resource-limited or model-failure validation'])
    (r/'cooling-fit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='parameters'},indent=2))
