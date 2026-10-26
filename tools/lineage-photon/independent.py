"""Independent row-generator assembly and positive uniformization recurrence."""
import json, math, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.stats import poisson

def setup(m, stimulus):
    nr=m['rna_capacity']+1;np_=m['protein_capacity']+1;n=2*nr*np_;r=m['rates']
    coords=np.array([(s,a,b) for s in (0,1) for a in range(nr) for b in range(np_)])
    ix=lambda s,a,b:(s*nr+a)*np_+b
    ii=[];jj=[];vv=[];ti=[];tj=[];tv=[];h=[]
    for i,(s,a,b) in enumerate(coords):
        dest=[]
        if s:dest.append((ix(0,a,b),r['off'],False))
        else:dest.append((ix(1,a,b),stimulus*r['on']*(1+r['feedback']/(1+(r['half']/b)**2)) if b else stimulus*r['on'],False))
        if a<nr-1:dest.append((ix(s,a+1,b),r['transcription']*(s+(1-s)*r['leak'])*(1-r['burst_two']),True))
        if a<nr-2:dest.append((ix(s,a+2,b),r['transcription']*(s+(1-s)*r['leak'])*r['burst_two'],True))
        if a>0:dest.append((ix(s,a-1,b),a*r['rna_decay'],False))
        if b<np_-1 and a>0:dest.append((ix(s,a,b+1),a*r['translation'],False))
        if b>0:dest.append((ix(s,a,b-1),b*r['protein_decay'],False))
        for j,v,tag in dest:
            ii.append(i);jj.append(j);vv.append(v)
            if tag:ti.append(i);tj.append(j);tv.append(v)
        ii.append(i);jj.append(i);vv.append(-sum(x[1] for x in dest))
        h.append(r['photon_background']+b*r['photon_brightness'])
    return csr_matrix((vv,(ii,jj)),shape=(n,n)),csr_matrix((tv,(ti,tj)),shape=(n,n)),np.array(h),coords

def step(m, interval, seed, reverse=False, tagged=False):
    R,D,h,_=setup(m,interval['stimulus']);n=len(h);K=seed.shape[1]-1;y=interval['photons']
    core=R-(D if tagged else 0)-csr_matrix((h,(np.arange(n),np.arange(n))),shape=(n,n))
    rate=float(max(-core.diagonal()))
    P=eye(n,format='csr')+core/rate;D=D/rate
    if not reverse:P=P.T.tocsr();D=D.T.tocsr()
    v=np.zeros((n,K+1,y+1));v[:,:,0]=seed
    pieces=max(1,math.ceil(rate*interval['duration']/20));mu=rate*interval['duration']/pieces
    limit=max(0,int(mu))
    while poisson.sf(limit,mu)>2e-19:limit+=1
    weights=[math.exp(-mu)]
    for j in range(1,limit+1):weights.append(weights[-1]*mu/j)
    for _ in range(pieces):
        term=v;result=weights[0]*term
        for w in weights[1:]:
            nxt=(P@term.reshape(n,-1)).reshape(v.shape)
            if tagged and K:nxt[:,1:,:]+=(D@term[:,:-1,:].reshape(n,-1)).reshape(n,K,y+1)
            if y:nxt[:,:,1:]+=h[:,None,None]/rate*term[:,:,:-1]
            term=nxt;result+=w*term
        v=result
    return v[:,:,y]

def histories(m,tagged):
    _,_,_,c=setup(m,1);K=m['count_limit'] if tagged else 0
    a=np.zeros((len(c),K+1));a[np.nonzero(np.all(c==m['initial'],axis=1))[0][0],0]=1
    for t in m['mother']:a=step(m,t,a,tagged=tagged)
    daughters=[]
    for name in ['daughter_a','daughter_b']:
        v=np.zeros_like(a);v[:,0]=1
        for t in m[name][::-1]:v=step(m,t,v,True,tagged)
        daughters.append(v)
    return a,*daughters,c

def solve(m):
    started=time.monotonic();plain=histories(m,False);tagged=histories(m,True)
    a,u,v,c=plain;ta,tu,tv,_=tagged;K=m['count_limit'];nr=m['rna_capacity']+1;np_=m['protein_capacity']+1
    parents=np.zeros(len(c));joint=np.zeros((nr,nr));counts=np.zeros(K+1)
    for i,(s,r,p) in enumerate(c):
        for x in range(r+1):
            px=math.comb(r,x)*m['partition_rna']**x*(1-m['partition_rna'])**(r-x)
            for y in range(p+1):
                w=px*math.comb(p,y)*m['partition_protein']**y*(1-m['partition_protein'])**(p-y)
                ia=(s*nr+x)*np_+y;ib=(s*nr+r-x)*np_+p-y
                term=w*a[i,0]*u[ia,0]*v[ib,0]
                parents[i]+=term;joint[x,r-x]+=term
                branch=np.polynomial.polynomial.polymul(tu[ia],tv[ib])
                total=np.polynomial.polynomial.polymul(ta[i],branch)[:K+1]
                counts[:len(total)]+=w*total
    z=parents.sum()
    return dict(log_evidence=float(np.log(z)),parent=(parents/z).tolist(),daughter_rna=(joint/z).tolist(),
                counts=(counts/z).tolist(),count_overflow=float(1-counts.sum()/z),seconds=time.monotonic()-started)

if __name__=='__main__':
    d=Path(__file__).parent/'draft';out=solve(json.loads((d/'model.json').read_text()))
    (d/'independent.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['parent','daughter_rna','counts']}),flush=True)
