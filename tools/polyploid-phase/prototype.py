"""Original synthetic polyploid MEC experiment; no author truth used by fit()."""
import argparse,json,time
from pathlib import Path
import numpy as np
from numba import njit


def generate(seed=20261031,k=8,n=240,m=2400):
    rng=np.random.default_rng(seed)
    # Mosaic homologs, with both private and shared ancestral tracts.
    ancestors=rng.integers(0,2,(k,n),dtype=np.int8)
    H=ancestors.copy()
    for h in range(k):
        for start in range(0,n,30):
            a=rng.integers(k)
            H[h,start:start+30]=ancestors[a,start:start+30]
    keep=(H.sum(0)>0)&(H.sum(0)<k);H=H[:,keep];n=H.shape[1]
    positions=[];alleles=[];weights=[];starts=[0]
    for r in range(m):
        h=rng.integers(k); length=int(rng.integers(25,66));start=int(rng.integers(max(1,n-length+1)))
        js=np.arange(start,min(n,start+length));js=js[rng.random(len(js))>.18]
        w=rng.choice([1,2,3],len(js),p=[.25,.5,.25])
        # Error odds exp(-w) give integer mismatch weights a likelihood interpretation.
        calls=H[h,js]^(rng.random(len(js))<1/(1+np.exp(w)))
        positions.extend(js);alleles.extend(calls);weights.extend(w);starts.append(len(positions))
    return H,tuple(np.asarray(x,dtype=np.int32) for x in (starts,positions,alleles,weights))

@njit
def costs(H,starts,pos,alleles,weights):
    m=len(starts)-1;k=len(H);D=np.zeros((m,k),np.int32)
    for r in range(m):
        for q in range(starts[r],starts[r+1]):
            for h in range(k):
                if H[h,pos[q]]!=alleles[q]:D[r,h]+=weights[q]
    return D

@njit
def coordinate(H,starts,pos,alleles,weights,dosage,iterations=50):
    k,n=H.shape;m=len(starts)-1;previous=2147483647
    for it in range(iterations):
        D=costs(H,starts,pos,alleles,weights); assignment=np.empty(m,np.int32);obj=0
        for r in range(m):assignment[r]=np.argmin(D[r]);obj+=D[r,assignment[r]]
        if obj>=previous:break
        previous=obj
        diff=np.zeros((k,n),np.int32)
        for r in range(m):
            h=assignment[r]
            for q in range(starts[r],starts[r+1]):diff[h,pos[q]]+=weights[q]*(1-2*alleles[q])
        for j in range(n):
            order=np.argsort(diff[:,j]);H[:,j]=0
            for x in range(dosage[j]):H[order[x],j]=1
    return H

@njit
def anneal(H,starts,pos,alleles,weights,dosage,seed,sweeps=300):
    np.random.seed(seed);k,n=H.shape;m=len(starts)-1
    # Per-variant inverted incidence lists, each read calls a variant at most once.
    counts=np.zeros(n,np.int32)
    for p in pos:counts[p]+=1
    offsets=np.zeros(n+1,np.int32)
    for j in range(n):offsets[j+1]=offsets[j]+counts[j]
    rr=np.empty(len(pos),np.int32);aa=rr.copy();ww=rr.copy();cur=offsets.copy()
    for r in range(m):
        for q in range(starts[r],starts[r+1]):
            p=pos[q];a=cur[p];rr[a]=r;aa[a]=alleles[q];ww[a]=weights[q];cur[p]+=1
    D=costs(H,starts,pos,alleles,weights);bestH=H.copy();score=0
    for r in range(m):score+=np.min(D[r])
    best=score
    for sweep in range(sweeps):
        T= max(.12,3.5*(1-sweep/sweeps)**2)
        for step in range(n*k):
            j=np.random.randint(n);a=np.random.randint(k);b=np.random.randint(k)
            if H[a,j]==H[b,j]:continue
            delta=0
            for q in range(offsets[j],offsets[j+1]):
                r=rr[q];old=np.min(D[r]);v=ww[q]*(1 if H[a,j]==aa[q] else -1)
                D[r,a]+=v;D[r,b]-=v;delta+=np.min(D[r])-old
            if delta<=0 or np.random.random()<np.exp(-delta/T):
                H[a,j],H[b,j]=H[b,j],H[a,j];score+=delta
                if score<best:best=score;bestH=H.copy()
            else:
                for q in range(offsets[j],offsets[j+1]):
                    r=rr[q];v=ww[q]*(1 if H[a,j]==aa[q] else -1);D[r,a]-=v;D[r,b]+=v
        # Swap suffixes to cross phasing barriers while preserving every dosage.
        for move in range(12):
            a=np.random.randint(k);b=np.random.randint(k);j=np.random.randint(n)
            if a==b:continue
            oldA=D[:,a].copy();oldB=D[:,b].copy();oldscore=score
            for r in range(m):
                old=np.min(D[r])
                for q in range(starts[r],starts[r+1]):
                    p=pos[q]
                    if p>=j and H[a,p]!=H[b,p]:
                        v=weights[q]*(1 if H[a,p]==alleles[q] else -1);D[r,a]+=v;D[r,b]-=v
                score+=np.min(D[r])-old
            delta=score-oldscore
            if delta<=0 or np.random.random()<np.exp(-delta/T):
                for p in range(j,n):H[a,p],H[b,p]=H[b,p],H[a,p]
                if score<best:best=score;bestH=H.copy()
            else:D[:,a]=oldA;D[:,b]=oldB;score=oldscore
    return bestH,best


def fit(data,dosage,k,seed=1,restarts=10,sweeps=300):
    rng=np.random.default_rng(seed);best=10**12;out=None;start=time.monotonic()
    for run in range(restarts):
        H=np.zeros((k,len(dosage)),np.int8)
        for j,d in enumerate(dosage):H[rng.permutation(k)[:d],j]=1
        H=coordinate(H,*data,dosage)
        H,score=anneal(H,*data,dosage,seed+run,sweeps)
        if score<best:best=score;out=H.copy()
        print(json.dumps(dict(run=run,score=int(score),best=int(best),seconds=time.monotonic()-start)),flush=True)
    return out,best

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--restarts',type=int,default=6);p.add_argument('--sweeps',type=int,default=150);p.add_argument('--reads',type=int,default=1800);a=p.parse_args()
    H,data=generate(m=a.reads);dosage=H.sum(0).astype(np.int32)
    planted=int(costs(H,*data).min(1).sum());print(json.dumps(dict(planted_cost=planted,variants=H.shape[1],reads=a.reads)),flush=True)
    out,score=fit(data,dosage,len(H),restarts=a.restarts,sweeps=a.sweeps)
    Path('tools/polyploid-phase/draft/prototype-result.json').write_text(json.dumps(dict(planted_cost=planted,reference_cost=int(score),eligible=bool(score<=planted)),indent=2)+'\n')
