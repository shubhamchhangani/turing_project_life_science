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
from source_kernel import source_moments,source_invert


def calculate(model,progress=None):
    original=model
    planar=[];nonlocal_edges=[]
    for a,b,j in model['edges']:
        ax,ay=model['coordinates'][a];bx,by=model['coordinates'][b]
        (planar if (bx-ax,by-ay) in [(1,0),(0,1),(1,1)] else nonlocal_edges).append([a,b,j])
    known=model['ligand_sites']+model['reporter_sites']
    extra_sites=sorted(({s for a,b,j in nonlocal_edges for s in [a,b]}|{s for s,h in model['fields']})-set(known))
    model=dict(model,edges=planar,reporter_sites=model['reporter_sites']+extra_sites)
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
    count=1<<len(sites);assert len(union)<=128
    path_words=np.zeros((len(sites),2),np.uint64)
    for bit,path in enumerate(paths):
        for edge in path:
            k=index[edge];path_words[bit,k//64]^=np.uint64(1)<<np.uint64(k%64)
    union=np.array(union);delta=r[union]*(t[union]-1/t[union])
    if progress:progress(0,count)
    joint=source_moments(inverse,path_words,delta,t[union],count)
    if progress:progress(count,count)
    source_invert(joint)
    site_index={s:i for i,s in enumerate(sites)};ratio=0.
    probabilities=np.zeros(1<<len(original['reporter_sites']))
    # Chunked aggregation bounds memory even for a large special-site marginal.
    for start in range(0,count,262144):
        stop=min(count,start+262144);ids=np.arange(start,stop,dtype=np.int64)
        spin={s:1-2*((ids>>site_index[s])&1) for s in sites}
        energy=sum(h*spin[s] for s,h in original['fields'])
        for a,b,j in nonlocal_edges:energy=energy+j*spin[a]*spin[b]
        weighted=joint[start:stop]*np.exp(energy);ratio+=weighted.sum()
        select=np.ones(stop-start,dtype=bool)
        for s,value in zip(original['ligand_sites'],original['ligand_spins']):select&=spin[s]==value
        labels=np.zeros(stop-start,dtype=np.int64)
        for bit,s in enumerate(original['reporter_sites']):labels+=((1-spin[s])//2)<<bit
        probabilities+=np.bincount(labels[select],weights=weighted[select],minlength=len(probabilities))
    pc=probabilities.sum()/ratio
    return {'log_partition':float(logz+np.log(ratio)),'log_clamp_probability':float(np.log(pc)),
            'probabilities':(probabilities/probabilities.sum()).tolist(),
            'source':'Independent skew rank-two updates, tensor inversion and explicit scaffold reweighting; no oracle outputs.'}


if __name__=='__main__':
    p=Path(__file__).parent/'draft';t0=time.monotonic()
    result=calculate(json.loads((p/'model.json').read_text()),lambda i,n:print(i,n,round(time.monotonic()-t0,2),flush=True))
    result['seconds']=time.monotonic()-t0
    (p/'expectations.json').write_text(json.dumps(result,indent=2)+'\n');print('Completed',result['seconds'])
