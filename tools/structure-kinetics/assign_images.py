"""Refine every image against inferred conformations; true poses/classes unavailable."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from model import coordinates,projection

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--geometry-fit',type=Path,required=True)
    p.add_argument('--pose-fit',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--orientations',type=int,default=160);p.add_argument('--refine',type=int,default=3)
    p.add_argument('--indices',default='');a=p.parse_args()
    c=json.loads((a.data/'config.json').read_text());data=np.load(a.data/'images.npz',allow_pickle=False)
    runs=json.loads(a.geometry_fit.read_text())['runs'];best=max(runs,key=lambda r:r['log_likelihood']);angles=np.array(best['parameters'][:8]).reshape(4,2);xyz=coordinates(c,angles)
    previous={r['image']:r['parameters'][:5] for r in json.loads(a.pose_fit.read_text())['rows']}
    rng=np.random.default_rng(8541);frames=Rotation.random(a.orientations,random_state=rng).as_rotvec();start=time.monotonic();rows=[]
    indices=list(map(int,a.indices.split(','))) if a.indices else list(range(len(data['images'])))
    for i in indices:
        img=data['images'][i]
        seeds=[np.r_[rv,0.,0.] for rv in frames];seeds.append(np.array(previous[i]));fits=[]
        for k in range(4):
            def residual(x):return ((projection(c,xyz[k],Rotation.from_rotvec(x[:3]).as_matrix(),x[3:],float(data['defocus'][i]))-img)/c['image_noise_sd']).ravel()
            ranked=np.argsort([np.dot(v:=residual(x),v) for x in seeds])[:a.refine]
            trials=[least_squares(residual,seeds[j],bounds=([-4*np.pi]*3+[-4,-4],[4*np.pi]*3+[4,4]),max_nfev=90,ftol=1e-8,xtol=1e-8) for j in ranked]
            r=min(trials,key=lambda f:f.cost);fits.append(dict(class_id=k,chi_square=float(2*r.cost),pose=r.x.tolist()))
        winner=min(fits,key=lambda f:f['chi_square']);row=dict(image=i,protocol=int(data['protocol'][i]),fits=fits,
            class_id=winner['class_id'],reduced_chi_square=winner['chi_square']/(img.size-5))
        rows.append(row);report=dict(status='fixed_geometry_image_profile',complete=i==indices[-1],truth_used=False,
            geometry_parameters=angles.tolist(),seconds=time.monotonic()-start,rows=rows,
            note='Pose-profile likelihood, not a marginal orientation likelihood or full joint refinement')
        tmp=a.out.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(a.out)
        if (i+1)%16==0:print(json.dumps(dict(images=i+1,noise_level_fits=sum(r['reduced_chi_square']<1.5 for r in rows),seconds=report['seconds'])),flush=True)
