"""Author-only alternative: inverse cooperative receptor interactions.

Analytic finite-torus likelihood derivatives. This does not change the frozen
forward task currently being evaluated.
"""
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import inv
from scipy.optimize import minimize
from prototype import Torus,pfaffian_log,enumerate_spins


class InverseTorus(Torus):
    def logz_gradient(self,couplings):
        t=np.tanh(couplings)
        if np.any(abs(t)<1e-12):
            raise ValueError('Zero weights require a limiting evaluation')
        index=np.arange(len(self.edges))*2
        roots=[];derivatives=[]
        for tx,ty in [(0,0),(0,1),(1,0),(1,1)]:
            w=t*np.where((tx*self.wrap[:,0]+ty*self.wrap[:,1])%2,-1,1)
            b=self.base.copy()
            entries=self.reversal/w
            b[index,index+1]=entries
            b[index+1,index]=-entries
            inverse=inv(b,overwrite_a=False,check_finite=False)
            sign,logabs=pfaffian_log(b)
            sign*=np.prod(np.sign(entries))
            logabs-=np.log(abs(entries)).sum()
            roots.append(sign*np.exp(logabs))
            derivatives.append((1/t-t)*(1-inverse[index+1,index]*entries))
        roots=np.array(roots)*np.array([-1,1,1,1])
        total=roots.sum()
        assert total>0
        z=self.lx*self.ly*np.log(2)+np.log(np.cosh(couplings)).sum()+np.log(total/2)
        gradient=t+(roots[:,None]*np.array(derivatives)).sum(axis=0)/total
        return z,gradient


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fit',action='store_true')
    args=parser.parse_args()
    rng=np.random.default_rng(20261027)
    small=InverseTorus(3,4)
    k=rng.uniform(-.25,.65,len(small.edges))
    z,g=small.logz_gradient(k)
    truth=enumerate_spins(small,k,[[p,q] for p,q,d in small.edges])
    small_error=float(max(abs(g-truth[1])))
    print('small gradient error',max(abs(g-truth[1])),'logZ error',abs(z-truth[0]),flush=True)
    assert max(abs(g-truth[1]))<1e-10
    root=Path(__file__).resolve().parents[2]
    data=root/'tasks/life-sciences/biology/receptor-cluster-thermodynamics/environment/data'
    model=json.loads((data/'model.json').read_text())
    large=InverseTorus(model['width'],model['height'])
    with (data/'couplings.csv').open() as stream:k=np.array([float(r['K']) for r in csv.DictReader(stream)])
    start=time.monotonic()
    z,g=large.logz_gradient(k)
    print('full evaluation seconds',time.monotonic()-start,'logZ',z,'range',g.min(),g.max(),flush=True)
    out=Path(__file__).resolve().parent/'draft/inverse';out.mkdir(exist_ok=True)
    np.save(out/'correlations.npy',g)
    (out/'gradient-status.json').write_text(json.dumps({'small_gradient_error':small_error,'full_logz':z,'full_gradient_min':float(g.min()),'full_gradient_max':float(g.max())},indent=2)+'\n')
    if args.fit:
        calls=[0]
        observed=np.round(g,12)
        start=time.monotonic()
        def objective(candidate):
            value,gradient=large.logz_gradient(candidate)
            calls[0]+=1
            residual=gradient-observed
            print('evaluation',calls[0],'objective',value-np.dot(candidate,observed),'gradient',max(abs(residual)),'seconds',time.monotonic()-start,flush=True)
            return value-np.dot(candidate,observed),residual
        result=minimize(objective,np.full(len(k),.35),jac=True,method='L-BFGS-B',bounds=[(-.4,.95)]*len(k),options={'maxiter':160,'ftol':1e-15,'gtol':2e-10,'maxls':30,'maxcor':25})
        np.save(out/'fitted.npy',result.x)
        report={'success':bool(result.success),'message':str(result.message),'max_coupling_error':float(max(abs(result.x-k))),'calls':calls[0],'seconds':time.monotonic()-start}
        (out/'fit-report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(report,flush=True)
