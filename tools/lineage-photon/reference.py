"""Sparse augmented-generator solution for path counts and lineage smoothing."""
import json, math, time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.stats import binom

def states(model):
    return [(s,m,p) for s in range(2) for m in range(model['rna_capacity']+1)
            for p in range(model['protein_capacity']+1)]

def generator(model, stimulus):
    st=states(model); lookup={v:i for i,v in enumerate(st)}; n=len(st)
    r=model['rates']; rows=[];cols=[];values=[];br=[];bc=[];bv=[]
    M=model['rna_capacity'];P=model['protein_capacity']
    for j,(s,m,p) in enumerate(st):
        transitions=[((1-s,m,p), r['off'] if s else r['on']*stimulus*(1+r['feedback']*p*p/(r['half']**2+p*p)),False)]
        for b,prob in [(1,1-r['burst_two']),(2,r['burst_two'])]:
            if m+b<=M: transitions.append(((s,m+b,p),r['transcription']*(1 if s else r['leak'])*prob,True))
        if m: transitions.append(((s,m-1,p),r['rna_decay']*m,False))
        if p<P: transitions.append(((s,m,p+1),r['translation']*m,False))
        if p: transitions.append(((s,m,p-1),r['protein_decay']*p,False))
        total=0.
        for dst,rate,marked in transitions:
            i=lookup[dst];rows.append(i);cols.append(j);values.append(rate);total+=rate
            if marked:br.append(i);bc.append(j);bv.append(rate)
        rows.append(j);cols.append(j);values.append(-total)
    Q=sparse.csc_matrix((values,(rows,cols)),shape=(n,n))
    A=sparse.csc_matrix((bv,(br,bc)),shape=(n,n))
    H=sparse.diags([r['photon_background']+r['photon_brightness']*p for s,m,p in st],format='csc')
    return Q,A,H

def exposure(model, interval, vector, backward=False, marked=False):
    Q,A,H=generator(model,interval['stimulus']);K=vector.shape[0]-1;y=interval['photons'];n=Q.shape[0]
    B=Q-(A if marked else 0)-H
    if backward:B=B.T;A=A.T
    T=sparse.kron(sparse.eye((K+1)*(y+1)),B,format='csc')
    if marked and K:
        T+=sparse.kron(sparse.diags(np.ones(K),-1,shape=(K+1,K+1)),sparse.kron(sparse.eye(y+1),A),format='csc')
    if y:
        T+=sparse.kron(sparse.eye(K+1),sparse.kron(sparse.diags(np.ones(y),-1,shape=(y+1,y+1)),H),format='csc')
    seed=np.zeros((K+1,y+1,n));seed[:,0,:]=vector
    out=expm_multiply(T*interval['duration'],seed.ravel()).reshape(seed.shape)
    return out[:,y,:]

def messages(model, marked=False):
    st=states(model);n=len(st);K=model['count_limit'] if marked else 0
    mother=np.zeros((K+1,n));mother[0,st.index(tuple(model['initial']))]=1.
    for interval in model['mother']:mother=exposure(model,interval,mother,marked=marked)
    branches=[]
    for name in ['daughter_a','daughter_b']:
        v=np.zeros_like(mother);v[0,:]=1.
        for interval in reversed(model[name]):v=exposure(model,interval,v,True,marked)
        branches.append(v)
    return mother,*branches

def contract(model, mother, a, b):
    st=states(model);ix={v:i for i,v in enumerate(st)};K=len(mother)-1
    parent=np.zeros(len(st));joint=np.zeros((model['rna_capacity']+1,)*2);counts=np.zeros(K+1)
    for j,(s,m,p) in enumerate(st):
        for ma in range(m+1):
            for pa in range(p+1):
                w=binom.pmf(ma,m,model['partition_rna'])*binom.pmf(pa,p,model['partition_protein'])
                ia=ix[s,ma,pa];ib=ix[s,m-ma,p-pa]
                c=np.convolve(mother[:,j],np.convolve(a[:,ia],b[:,ib])[:K+1])[:K+1]*w
                counts+=c
                if K==0:parent[j]+=c[0];joint[ma,m-ma]+=c[0]
    return counts,parent,joint

def solve(model):
    start=time.monotonic()
    u,parent,joint=contract(model,*messages(model));evidence=float(u[0])
    count,*_=contract(model,*messages(model,True));posterior=count/evidence
    return dict(log_evidence=math.log(evidence),parent=(parent/evidence).tolist(),
                daughter_rna=(joint/evidence).tolist(),counts=posterior.tolist(),
                count_overflow=float(1-posterior.sum()),seconds=time.monotonic()-start)

if __name__=='__main__':
    d=Path(__file__).parent/'draft';result=solve(json.loads((d/'model.json').read_text()))
    (d/'reference.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['parent','daughter_rna','counts']}),flush=True)
