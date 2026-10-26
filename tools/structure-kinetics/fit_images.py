"""Input-only per-image pose/geometry baseline; test this shortcut before joint fitting."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from model import projection,coordinates

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--images',type=int,default=8);p.add_argument('--orientations',type=int,default=128);p.add_argument('--refine',type=int,default=6);args=p.parse_args()
    c=json.loads((args.data/'config.json').read_text());data=np.load(args.data/'images.npz',allow_pickle=False)
    rng=np.random.default_rng(9105);frames=Rotation.random(args.orientations,random_state=rng).as_rotvec();started=time.monotonic();rows=[]
    for i in range(min(args.images,len(data['images']))):
        observed=data['images'][i];defocus=float(data['defocus'][i]);evaluations=[0]
        def residual(x):
            evaluations[0]+=1;xyz=coordinates(c,np.array([x[5:7]]))[0]
            return ((projection(c,xyz,Rotation.from_rotvec(x[:3]).as_matrix(),x[3:5],defocus)-observed)/c['image_noise_sd']).ravel()
        seeds=[np.r_[rv,0.,0.,rng.uniform(-2,2),rng.uniform(-.8,.8)] for rv in frames]
        order=np.argsort([np.dot(r:=residual(s),r) for s in seeds])[:args.refine]
        best=None
        lower=[-4*np.pi]*3+[-4.,-4.,c['parameter_bounds']['azimuth'][0],c['parameter_bounds']['elevation'][0]]
        upper=[4*np.pi]*3+[4.,4.,c['parameter_bounds']['azimuth'][1],c['parameter_bounds']['elevation'][1]]
        for j in order:
            fit=least_squares(residual,seeds[j],bounds=(lower,upper),max_nfev=110,ftol=1e-8,xtol=1e-8,gtol=1e-6)
            if best is None or fit.cost<best.cost:best=fit
        row=dict(image=i,reduced_chi_square=float(2*best.cost/(observed.size-7)),parameters=best.x.tolist(),evaluations=evaluations[0])
        rows.append(row)
        report=dict(status='image_only_shortcut_screen',truth_used=False,photon_data_used=False,
            complete=(i==min(args.images,len(data['images']))-1),seconds=time.monotonic()-started,rows=rows)
        args.out.parent.mkdir(exist_ok=True,parents=True);tmp=args.out.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(args.out)
        print(json.dumps({k:v for k,v in row.items() if k!='parameters'}),flush=True)
