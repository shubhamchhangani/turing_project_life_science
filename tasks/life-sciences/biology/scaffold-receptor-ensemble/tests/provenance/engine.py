"""Input-only planar zero-field Ising spin moments via even-subgraph sums."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from collections import deque


class Planar:
    def __init__(self, model):
        self.model=model
        xy=np.array(model['coordinates'])
        edges=model['edges']; self.n=len(xy)
        self.t=np.tanh([e[2] for e in edges])
        self.adj=[[] for _ in xy]
        src=[]; dst=[]
        for i,(a,b,_) in enumerate(edges):
            src.extend([a,b]); dst.extend([b,a])
            self.adj[a].append((b,i)); self.adj[b].append((a,i))
        self.src=np.array(src);self.dst=np.array(dst)
        angles=np.angle((xy[self.dst]-xy[self.src])@np.array([1,1j]))
        outgoing=[[] for _ in xy]
        for e,a in enumerate(src):outgoing[a].append(e)
        rows=[];cols=[];phases=[]
        for e,b in enumerate(dst):
            for f in outgoing[b]:
                if f==(e^1):continue
                turn=(angles[f]-angles[e]+np.pi)%(2*np.pi)-np.pi
                rows.append(e);cols.append(f);phases.append(-np.exp(.5j*turn))
        self.rows=np.array(rows);self.cols=np.array(cols)
        self.phases=np.array(phases);self.dim=len(src)
        self.base=self.loop_log(self.t)
        self.logz=self.n*np.log(2)+sum(np.logaddexp(e[2],-e[2])-np.log(2) for e in edges)+self.base

    def loop_log(self, weights):
        diagonal=np.arange(self.dim)
        matrix=csc_matrix((np.r_[self.phases*weights[self.cols//2],np.ones(self.dim)],
                           (np.r_[self.rows,diagonal],np.r_[self.cols,diagonal])),shape=(self.dim,self.dim))
        factor=splu(matrix)
        return .5*float(np.log(np.abs(factor.U.diagonal())).sum())

    def path(self, a,b):
        previous={a:None};todo=deque([a])
        while b not in previous:
            v=todo.popleft()
            for w,e in self.adj[v]:
                if w not in previous:previous[w]=(v,e);todo.append(w)
        edges=[]
        while b!=a:
            b,e=previous[b];edges.append(e)
        return edges

    def correlation(self, sites):
        if len(sites)%2:return 0.0
        changed=set()
        for a,b in zip(sites[::2],sites[1::2]):
            changed.symmetric_difference_update(self.path(a,b))
        if not changed:return 1.0
        index=np.array(sorted(changed));weights=self.t.copy()
        weights[index]=1/weights[index]
        return float(np.exp(np.log(self.t[index]).sum()+self.loop_log(weights)-self.base))


def distribution(model, progress=None):
    engine=Planar(model);sites=model['ligand_sites']+model['reporter_sites']
    count=1<<len(sites);moments=np.zeros(count)
    for mask in range(count):
        if mask.bit_count()%2==0:
            moments[mask]=engine.correlation([s for j,s in enumerate(sites) if mask>>j&1])
        if progress and mask%256==0:progress(mask,count)
    # Walsh inversion: bit 0 is spin +1, bit 1 is spin -1.
    probabilities=moments.copy()
    step=1
    while step<count:
        blocks=probabilities.reshape(-1,2*step)
        left=blocks[:,:step].copy();right=blocks[:,step:].copy()
        blocks[:,:step]=left+right;blocks[:,step:]=left-right
        step*=2
    probabilities/=count
    nlig=len(model['ligand_sites'])
    clamp=sum((s==-1)<<i for i,s in enumerate(model['ligand_spins']))
    joint=probabilities[clamp::1<<nlig]
    clamp_probability=joint.sum()
    return {'log_partition':engine.logz,'log_clamp_probability':float(np.log(clamp_probability)),
            'probabilities':(joint/clamp_probability).tolist(),
            'min_unconditional_probability':float(probabilities.min())}
