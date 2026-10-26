"""Inverse receptor cooperativity with three spatially separated ligand biases."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import inv
from scipy.optimize import minimize
from prototype import Torus,pfaffian_log


class FieldTorus(Torus):
    def __init__(self,width,height,field_sites):
        super().__init__(width,height)
        assert len(field_sites)==3
        self.field_sites=field_sites
        self.masks=[0,3,5,6]
        self.edge_permutation=[]
        for p,q,d in self.edges:
            x,y=p
            self.edge_permutation.append(self.edge_id[(x,(y-1)%height,1)] if d==0 else self.edge_id[((x-1)%width,y,0)])
        self.edge_permutation=np.array(self.edge_permutation)

    def evaluate_gradient(self,parameters,hessian=False):
        couplings=parameters[:-3]
        fields=parameters[-3:]
        weights=self.dual_weights(couplings)
        index=np.arange(len(self.edges))*2
        moments=[];derivatives=[];second_derivatives=[]
        for mask in self.masks:
            sites=[p for i,p in enumerate(self.field_sites) if mask>>i&1]
            current=self.disorder(weights,sites)
            total=0.
            derivative=np.zeros(len(self.edges))
            second=np.zeros((len(self.edges),len(self.edges))) if hessian else None
            for tx,ty in [(0,0),(0,1),(1,0),(1,1)]:
                w=current*np.where((tx*self.wrap[:,0]+ty*self.wrap[:,1])%2,-1,1)
                b=self.base.copy()
                entries=self.reversal/w
                b[index,index+1]=entries
                b[index+1,index]=-entries
                inverse=inv(b,overwrite_a=False,check_finite=False)
                sign,logabs=pfaffian_log(b)
                sign*=np.prod(np.sign(entries))
                logabs-=np.log(abs(entries)).sum()
                value=sign*np.exp(logabs)
                total+=value/4
                u=inverse[index+1,index]*entries
                gradient=2*(u-1)
                derivative+=value/4*gradient
                if hessian:
                    cross=inverse[np.ix_(index+1,index)]
                    bb=inverse[np.ix_(index+1,index+1)]
                    aa=inverse[np.ix_(index,index)]
                    log_hessian=-4*np.outer(entries,entries)*(cross*cross.T-bb*aa.T)
                    log_hessian[np.diag_indices_from(log_hessian)]+=4*u
                    second+=value/4*(log_hessian+np.outer(gradient,gradient))
            moments.append(total)
            derivatives.append(derivative[self.edge_permutation])
            if hessian:second_derivatives.append(second[np.ix_(self.edge_permutation,self.edge_permutation)])
        t=np.tanh(fields)
        coefficients=np.array([1.,t[0]*t[1],t[0]*t[2],t[1]*t[2]])
        sum_loops=np.dot(coefficients,moments)
        assert sum_loops>0
        logz=np.log(2)+sum(couplings)+np.log(np.cosh(fields)).sum()+np.log(sum_loops)
        edge_gradient=1+coefficients@np.array(derivatives)/sum_loops
        jacobian=np.array([[0,t[1],t[2],0],[0,t[0],0,t[2]],[0,0,t[0],t[1]]])*(1-t*t)[:,None]
        magnetization=t+jacobian@np.array(moments)/sum_loops
        gradient=np.r_[edge_gradient,magnetization]
        if not hessian:return logz,gradient
        n=len(couplings)
        covariance=np.empty((n+3,n+3))
        covariance[:n,:n]=np.tensordot(coefficients,np.array(second_derivatives),axes=1)/sum_loops-np.outer(edge_gradient-1,edge_gradient-1)
        cross=jacobian@np.array(derivatives)/sum_loops-np.outer(magnetization-t,edge_gradient-1)
        covariance[n:,:n]=cross;covariance[:n,n:]=cross.T
        second_coefficients=np.zeros((3,3,4))
        for term,mask in enumerate(self.masks):
            for i in range(3):
                if not mask>>i&1:continue
                others=[j for j in range(3) if mask>>j&1 and j!=i]
                second_coefficients[i,i,term]=-2*t[i]*(1-t[i]**2)*np.prod(t[others])
                for j in others:second_coefficients[i,j,term]=(1-t[i]**2)*(1-t[j]**2)
        covariance[n:,n:]=np.diag(1-t*t)+second_coefficients@np.array(moments)/sum_loops-np.outer(magnetization-t,magnetization-t)
        covariance=(covariance+covariance.T)/2
        return logz,gradient,covariance


def exhaustive(model,parameters):
    n=model.lx*model.ly
    spins=2*((np.arange(1<<n)[:,None]>>np.arange(n))&1)-1
    index=lambda p:p[1]*model.lx+p[0]
    sufficient=np.column_stack([spins[:,index(p)]*spins[:,index(q)] for p,q,d in model.edges]+[spins[:,index(p)] for p in model.field_sites])
    energy=sufficient@parameters
    shift=max(energy);weights=np.exp(energy-shift)
    return shift+np.log(weights.sum()),weights@sufficient/weights.sum()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fit',action='store_true')
    parser.add_argument('--width',type=int,default=17)
    parser.add_argument('--height',type=int,default=19)
    args=parser.parse_args()
    rng=np.random.default_rng(20261028)
    errors=[]
    for width,height in [(3,3),(3,4),(4,4)]:
        small=FieldTorus(width,height,[(0,0),(width-1,height-1),(1,height-1)])
        for low,high in [(.15,.35),(.35,.65),(-.2,.65)]:
            p=np.r_[rng.uniform(low,high,len(small.edges)),rng.uniform(-.6,.6,3)]
            actual=small.evaluate_gradient(p);direct=exhaustive(small,p)
            error=max(abs(actual[0]-direct[0]),max(abs(actual[1]-direct[1])))
            errors.append({'width':width,'height':height,'range':[low,high],'max_error':float(error)})
            assert error<1e-10
    print('small-case errors',max(r['max_error'] for r in errors),flush=True)
    sites=[(0,0),(args.width//2,args.height//2),(args.width-3,2)]
    model=FieldTorus(args.width,args.height,sites)
    rng=np.random.default_rng(20261028)
    couplings=[]
    for (x,y),q,d in model.edges:
        k=.46+.12*np.sin(.43*x+.61*y+.3*d)+rng.uniform(-.09,.09)
        if (17*x+11*y+29*d)%113==17:k=-rng.uniform(.10,.21)
        couplings.append(round(k,9))
    truth=np.r_[couplings,[.43,-.31,.24]]
    start=time.monotonic();logz,observables=model.evaluate_gradient(truth)
    print('generation seconds',time.monotonic()-start,'logZ',logz,'fields',observables[-3:],flush=True)
    out=Path(__file__).resolve().parent/'draft/field-inverse';out.mkdir(exist_ok=True)
    (out/'generation.json').write_text(json.dumps({'width':args.width,'height':args.height,'field_sites':sites,'seed':20261028,'log_partition':logz,'truth':truth.tolist(),'observables':observables.tolist(),'small_cases':errors},indent=2)+'\n')
    if args.fit:
        observed=np.round(observables,12)
        calls=[0];start=time.monotonic()
        def objective(p):
            value,g=model.evaluate_gradient(p)
            residual=g-observed;calls[0]+=1
            print('evaluation',calls[0],'objective',value-p@observed,'gradient',max(abs(residual)),'seconds',time.monotonic()-start,flush=True)
            return value-p@observed,residual
        result=minimize(objective,np.r_[np.full(len(couplings),.35),[.05,.05,.05]],jac=True,method='L-BFGS-B',bounds=[(-.35,.9)]*len(couplings)+[(-.8,.8)]*3,options={'maxiter':180,'ftol':1e-14,'gtol':1e-7,'maxls':30,'maxcor':30})
        np.save(out/'fitted.npy',result.x)
        report={'success':bool(result.success),'message':str(result.message),'max_coupling_error':float(max(abs(result.x[:-3]-truth[:-3]))),'max_field_error':float(max(abs(result.x[-3:]-truth[-3:]))),'calls':calls[0],'seconds':time.monotonic()-start}
        (out/'fit-report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(report,flush=True)
