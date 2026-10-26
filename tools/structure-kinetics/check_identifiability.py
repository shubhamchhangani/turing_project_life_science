"""Local likelihood-curvature diagnostic, using public data and completed fits only.

Inverse Hessians are local approximations, not calibrated confidence intervals
or proofs of global identifiability. Compare finite-difference resolutions.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from likelihood import load,loglike

def hessian(fun,x,h):
    n=len(x);H=np.zeros((n,n));f=fun(x)
    for i in range(n):
        ei=np.zeros(n);ei[i]=h;H[i,i]=(fun(x+ei)-2*f+fun(x-ei))/h**2
        for j in range(i):
            ej=np.zeros(n);ej[j]=h
            H[i,j]=H[j,i]=(fun(x+ei+ej)-fun(x+ei-ej)-fun(x-ei+ej)+fun(x-ei-ej))/(4*h*h)
    return H

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).parent/'draft');args=parser.parse_args()
    root=args.root;c,traces=load(root/'public');fits=json.loads((root/'crossmodal-fit.json').read_text())
    assert fits['complete'],'Wait for the multistart fit to complete'
    best=max(fits['runs'],key=lambda r:r['log_likelihood']);theta=np.asarray(best['parameters']);fixed=theta[:8];x=theta[8:]
    started=time.monotonic();fun=lambda z:-loglike(c,traces,np.r_[fixed,z])
    H1=hessian(fun,x,1e-3);H2=hessian(fun,x,5e-4);eigenvalues=np.linalg.eigvalsh(H2)
    result=dict(scope='local photon-only curvature with image-derived geometry fixed; no global identifiability claim',
                hessian_eigenvalues=eigenvalues.tolist(),
                finite_difference_relative_disagreement=float(np.linalg.norm(H1-H2)/np.linalg.norm(H2)),
                positive_definite=bool(eigenvalues[0]>0),truth_used=False)
    if eigenvalues[0]>0:
        cov=np.linalg.inv(H2);grads=[]
        index=3
        for i in range(4):
            for j in range(i+1,4):
                for source in [i,j]:
                    g=np.zeros(9)
                    if source:g[source-1]=1.
                    g[index]=-1.;grads.append(g)
                index+=1
        std=np.sqrt(np.einsum('ij,jk,ik->i',grads,cov,grads))
        result.update(log_rate_local_standard_errors=std.tolist(),condition_number=float(eigenvalues[-1]/eigenvalues[0]),
                      largest_local_rate_uncertainty_factor=float(np.exp(max(std))))
    result['seconds']=time.monotonic()-started
    (root/'identifiability-screen.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
