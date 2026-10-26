"""Independent signed expectations through sequential rank-two skew updates.

Constructs its own skew matrix and paths. Does not import reference code or
consume reference outputs. No square-root sign is guessed.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def calculate(model,progress=None):
    xy=model['coordinates'];edges=model['edges'];m=len(edges)
    directed=[(a,b,k) for k,(a,b,j) in enumerate(edges)]+[(b,a,k) for k,(a,b,j) in enumerate(edges)]
    outgoing=[[] for _ in xy];adj=[[] for _ in xy]
    for e,(a,b,k) in enumerate(directed):outgoing[a].append(e);adj[a].append((b,k))
    vectors=np.array([complex(*xy[b])-complex(*xy[a]) for a,b,k in directed]);r=vectors/abs(vectors)
    t=np.tanh([j for a,b,j in edges]);rr=[];cc=[];vv=[]
    for e,(a,b,k) in enumerate(directed):
        rr.append(e);cc.append((e+m)%(2*m));vv.append(r[e]/t[k])
        for f in outgoing[a]:
            if f==e:continue
            turn=np.angle(vectors[f]/(-vectors[e]))
            rr.append(e);cc.append(f);vv.append(-r[e]*np.exp(.5j*turn))
    matrix=coo_matrix((vv,(rr,cc)),shape=(2*m,2*m)).tocsc()
    assert abs(matrix+matrix.T).max()<1e-12
    lu=splu(matrix)
    # det(skew matrix) differs from the loop polynomial square by unit phase
    # and the inverse-square product of all tanh couplings.
    loop=.5*sum(np.log(abs(lu.U.diagonal())))+sum(np.log(abs(t)))
    logz=len(xy)*np.log(2)+sum(np.logaddexp(j,-j)-np.log(2) for a,b,j in edges)+loop
    sites=model['ligand_sites']+model['reporter_sites'];root=sites[-1]
    parents={root:None};queue=[root]
    for a in queue:
        for b,k in sorted(adj[a],reverse=True):
            if b not in parents:parents[b]=(a,k);queue.append(b)
    paths=[]
    for site in sites:
        p=set()
        while site!=root:site,k=parents[site];p.add(k)
        paths.append(p)
    union=sorted(set().union(*paths));index={e:i for i,e in enumerate(union)}
    selected=np.array([i for e in union for i in [e,e+m]])
    rhs=np.zeros((2*m,len(selected)),complex);rhs[selected,np.arange(len(selected))]=1
    inverse=lu.solve(rhs)[selected,:];inverse=(inverse-inverse.T)/2
    count=1<<len(sites);moments=np.zeros(count)
    for mask in range(count):
        if mask.bit_count()%2:continue
        changed=set()
        for j,path in enumerate(paths):
            if mask>>j&1:changed.symmetric_difference_update(path)
        if not changed:moments[mask]=1;continue
        chosen=sorted(changed)
        idx=np.array([i for e in chosen for i in [2*index[e],2*index[e]+1]])
        b=inverse[np.ix_(idx,idx)].astype(np.clongdouble);ratio=np.clongdouble(1)
        for edge in chosen:
            delta=r[edge]*(t[edge]-1/t[edge]);q=1-delta*b[0,1]
            ratio*=q/t[edge]
            g=b[2:,0];h=b[2:,1]
            b=b[2:,2:]+delta/q*(np.outer(g,h)-np.outer(h,g))
        assert abs(ratio.imag)<1e-7
        moments[mask]=ratio.real
        if progress and mask%1024==0:progress(mask,count)
    nlig=len(model['ligand_sites']);clamp=sum((s==-1)<<i for i,s in enumerate(model['ligand_spins']))
    joint=[]
    for pattern in range(1<<len(model['reporter_sites'])):
        assignment=clamp|(pattern<<nlig)
        joint.append(sum((-1 if (mask&assignment).bit_count()%2 else 1)*v for mask,v in enumerate(moments))/count)
    pc=sum(joint)
    return {'log_partition':float(logz),'log_clamp_probability':float(np.log(pc)),
            'probabilities':[float(x/pc) for x in joint],
            'source':'Independent skew-matrix rank-two update implementation; no oracle outputs.'}


if __name__=='__main__':
    p=Path(__file__).parent/'draft';t0=time.monotonic()
    result=calculate(json.loads((p/'model.json').read_text()),lambda i,n:print(i,n,round(time.monotonic()-t0,2),flush=True))
    result['seconds']=time.monotonic()-t0
    (p/'expectations.json').write_text(json.dumps(result,indent=2)+'\n');print('Completed',result['seconds'])
