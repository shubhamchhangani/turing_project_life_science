"""Validate the finite-field Fisher information and solve the inverse problem."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import solve
from field_inverse_prototype import FieldTorus,exhaustive


def fit(model,observed,lower,upper,max_iterations=30):
    parameters=np.r_[np.full(len(model.edges),.35),[.05,.05,.05]]
    start=time.monotonic();calls=0
    value,gradient,hessian=model.evaluate_gradient(parameters,True);calls+=1
    objective=value-parameters@observed
    for iteration in range(max_iterations):
        residual=gradient-observed
        print('Newton iteration',iteration,'calls',calls,'objective',objective,'gradient',max(abs(residual)),'elapsed',time.monotonic()-start,flush=True)
        if max(abs(residual))<2e-8:
            return parameters,{'iterations':iteration,'calls':calls,'seconds':time.monotonic()-start,'max_moment_error':float(max(abs(residual))),'converged':True,'minimum_fisher_eigenvalue':float(np.linalg.eigvalsh(hessian)[0])}
        step=solve(hessian,-residual,assume_a='pos',check_finite=False)
        alpha=1.
        positive=step>0;negative=step<0
        if positive.any():alpha=min(alpha,.99*min((upper[positive]-parameters[positive])/step[positive]))
        if negative.any():alpha=min(alpha,.99*min((lower[negative]-parameters[negative])/step[negative]))
        for backtrack in range(20):
            candidate=parameters+alpha*step
            nv,ng,nh=model.evaluate_gradient(candidate,True);calls+=1
            no=nv-candidate@observed
            if no<=objective+1e-4*alpha*residual@step+1e-11:
                parameters,value,gradient,hessian,objective=candidate,nv,ng,nh,no
                break
            alpha*=.5
        else:raise RuntimeError('Newton line search did not converge')
    raise RuntimeError('Newton iteration allowance exhausted')


if __name__=='__main__':
    rng=np.random.default_rng(991)
    small=FieldTorus(3,3,[(0,0),(2,2),(1,2)])
    parameters=np.r_[rng.uniform(-.15,.65,len(small.edges)),[.4,-.25,.15]]
    z,g,h=small.evaluate_gradient(parameters,True)
    n=9
    spins=2*((np.arange(1<<n)[:,None]>>np.arange(n))&1)-1
    idx=lambda p:p[1]*3+p[0]
    sufficient=np.column_stack([spins[:,idx(p)]*spins[:,idx(q)] for p,q,d in small.edges]+[spins[:,idx(p)] for p in small.field_sites])
    weights=np.exp(sufficient@parameters-z)
    means=weights@sufficient
    covariance=(sufficient.T*weights)@sufficient-np.outer(means,means)
    error=float(np.max(abs(h-covariance)));print('Fisher exhaustive error',error,flush=True)
    assert error<1e-10
    root=Path(__file__).resolve().parent/'draft/field-inverse'
    data=json.loads((root/'generation.json').read_text())
    model=FieldTorus(data['width'],data['height'],[tuple(p) for p in data['field_sites']])
    observed=np.round(data['observables'],12)
    lower=np.r_[np.full(len(model.edges),-.35),[-.8]*3]
    upper=np.r_[np.full(len(model.edges),.9),[.8]*3]
    fitted,report=fit(model,observed,lower,upper)
    truth=np.array(data['truth'])
    report.update({'max_coupling_error':float(max(abs(fitted[:-3]-truth[:-3]))),'max_field_error':float(max(abs(fitted[-3:]-truth[-3:]))),'small_fisher_error':error})
    np.save(root/'newton-fitted.npy',fitted)
    (root/'newton-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)
