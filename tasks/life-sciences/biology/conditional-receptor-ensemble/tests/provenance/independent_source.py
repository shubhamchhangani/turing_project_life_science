"""Independent expectations from star-path high-temperature graph polynomials.

Does not import the solution engine or read its outputs. Uses outgoing-edge
weights (rather than incoming), a common-root tree rather than paired paths,
and explicit character sums for the conditional law.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,time
from pathlib import Path
import numpy as np
from scipy.sparse import eye,coo_matrix
from scipy.sparse.linalg import splu


def calculate(data):
    vertices=data['coordinates'];bonds=data['edges'];n=len(vertices);m=len(bonds)
    directed=[(a,b,k) for k,(a,b,j) in enumerate(bonds)]+[(b,a,k) for k,(a,b,j) in enumerate(bonds)]
    outgoing={a:[] for a in range(n)}
    adjacency={a:[] for a in range(n)}
    for i,(a,b,k) in enumerate(directed):outgoing[a].append(i);adjacency[a].append((b,k))
    rr=[];cc=[];vv=[]
    for e,(a,b,k) in enumerate(directed):
        za=complex(*vertices[a]);zb=complex(*vertices[b])
        for f in outgoing[b]:
            _,c,_=directed[f]
            if c==a:continue
            turn=np.angle((complex(*vertices[c])-zb)/(zb-za))
            rr.append(e);cc.append(f);vv.append(np.exp(.5j*turn))
    rr=np.array(rr);cc=np.array(cc);vv=np.array(vv)
    bond_ids=np.array([k for a,b,k in directed]);identity=eye(2*m,format='csc',dtype=complex)
    weights=np.tanh([j for a,b,j in bonds])
    def loop(w):
        kw=identity-coo_matrix((vv*w[bond_ids[rr]],(rr,cc)),shape=(2*m,2*m)).tocsc()
        return float(np.log(np.abs(splu(kw).U.diagonal())).sum()/2)
    base=loop(weights)
    sites=data['ligand_sites']+data['reporter_sites'];root=sites[0]
    queue=[root];parents={root:None}
    for a in queue:
        for b,k in sorted(adjacency[a],reverse=True):
            if b not in parents:parents[b]=(a,k);queue.append(b)
    paths=[]
    for v in sites:
        path=set()
        while v!=root:v,k=parents[v];path.add(k)
        paths.append(path)
    count=1<<len(sites);moments=np.zeros(count);started=time.monotonic()
    for subset in range(count):
        if subset.bit_count()%2:continue
        changed=set()
        for bit,path in enumerate(paths):
            if subset>>bit&1:changed.symmetric_difference_update(path)
        if changed:
            ids=np.array(sorted(changed));w=weights.copy();w[ids]=1/w[ids]
            moments[subset]=np.exp(np.log(weights[ids]).sum()+loop(w)-base)
        else:moments[subset]=1
        if subset%1024==0:print(subset,count,round(time.monotonic()-started,2),flush=True)
    nlig=len(data['ligand_sites']);clamp=sum((s==-1)<<i for i,s in enumerate(data['ligand_spins']))
    outputs=[]
    for pattern in range(1<<len(data['reporter_sites'])):
        assignment=clamp|(pattern<<nlig)
        value=sum((-1 if (subset&assignment).bit_count()%2 else 1)*v for subset,v in enumerate(moments))/count
        outputs.append(value)
    normalization=sum(outputs)
    return {'log_partition':float(n*np.log(2)+sum(np.logaddexp(j,-j)-np.log(2) for a,b,j in bonds)+base),
            'log_clamp_probability':float(np.log(normalization)),
            'probabilities':[float(v/normalization) for v in outputs],
            'source':'Independent star-path high-temperature graph sums; no oracle outputs consumed.'}


if __name__=='__main__':
    p=Path(__file__).parent/'draft';t=time.monotonic()
    result=calculate(json.loads((p/'model.json').read_text()))
    result['seconds']=time.monotonic()-t
    (p/'expectations.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Completed',result['seconds'],flush=True)
