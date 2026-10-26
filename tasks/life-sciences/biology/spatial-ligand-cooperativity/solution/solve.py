"""Input-only maximum-likelihood reconstruction for a sparse-field membrane."""
import argparse
import csv
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from sparse_physics import SparseFields


def fit(data,output):
    model=json.loads((data/'model.json').read_text())
    sites=[(s['x'],s['y']) for s in model['field_sites']]
    physics=SparseFields(model['width'],model['height'],sites)
    with (data/'edge_moments.csv').open() as stream:
        edge_values={(int(r['x']),int(r['y']),r['direction']):float(r['mean_product']) for r in csv.DictReader(stream)}
    with (data/'site_moments.csv').open() as stream:
        field_values={r['receptor_id']:float(r['mean_state']) for r in csv.DictReader(stream)}
    target=np.array([edge_values[p[0],p[1],'E' if d==0 else 'N'] for p,q,d in physics.edges]+[field_values[s['id']] for s in model['field_sites']])
    calls=0;start=time.monotonic();last=None
    def objective(parameters):
        nonlocal calls,last
        value,gradient=physics.evaluate(parameters);residual=gradient-target;calls+=1
        print('evaluation',calls,'objective',value-parameters@target,'max_moment_error',max(abs(residual)),'seconds',time.monotonic()-start,flush=True)
        last=residual
        return value-parameters@target,residual
    initial=np.r_[np.full(len(physics.edges),.35),np.full(len(sites),.05)]
    bounds=[tuple(model['coupling_bounds'])]*len(physics.edges)+[tuple(model['field_bounds'])]*len(sites)
    result=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=bounds,
                    options={'maxiter':250,'ftol':1e-15,'gtol':1e-8,'maxls':30,'maxcor':40})
    assert max(abs(last))<2e-7, 'Moment fit did not converge sufficiently'
    output.mkdir(parents=True,exist_ok=True)
    with (output/'couplings.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['x','y','direction','K'])
        for (p,q,d),value in zip(physics.edges,result.x):writer.writerow([p[0],p[1],'E' if d==0 else 'N',f'{value:.10f}'])
    with (output/'biases.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['receptor_id','h'])
        for site,value in zip(model['field_sites'],result.x[len(physics.edges):]):writer.writerow([site['id'],f'{value:.10f}'])
    print('Fit completed',result.message,'calls',calls,'seconds',time.monotonic()-start,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path('/root/data'))
    parser.add_argument('--output',type=Path,default=Path('/root/results'))
    args=parser.parse_args();fit(args.data,args.output)
