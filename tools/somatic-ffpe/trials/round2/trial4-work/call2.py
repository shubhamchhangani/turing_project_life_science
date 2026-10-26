import numpy as np, sys
from scipy.stats import binom, poisson
r=np.load('/root/work/ref.npy'); n=len(r); idx=np.arange(n)
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4].astype(np.int64); NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].astype(np.int64)
Tc=C.sum(0); Nc=NC.sum(0); T=Tc.sum(1); N=Nc.sum(1)
SA=Tc.copy(); SD=np.repeat(T[:,None],4,1)
m=(r==1); SA[m,3]=C[1][m,3]; SD[m,3]=C[1][m].sum(1)
m=(r==2); SA[m,0]=C[0][m,0]; SD[m,0]=C[0][m].sum(1)
score=SA+Nc; score[idx[r<4],r[r<4]]=-1; a=score.argmax(1)
sa=SA[idx,a]; sd=SD[idx,a]; an=Nc[idx,a]
isct=((r==1)&(a==3))|((r==2)&(a==0))
p,f=float(sys.argv[1]),float(sys.argv[2])
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(n)
for s0,e0,c in segs: cn[s0:e0]=c
DT=p*cn+2*(1-p); DN=f*cn+2*(1-f)
# empirical error rates (per fragment, specific alt)
valid=r<4
def rate(cnts,deps,mask):
    return cnts[mask].sum()/deps[mask].sum()
low=(sa<=2)&(an<=1)&valid
et_ct=rate(sa,sd,low&isct); et_o=rate(sa,sd,low&~isct)
# for specific-alt rate, non-argmax alts also count; approximate using argmax alt (upward biased slightly -> conservative)
en_ct=rate(an,N,low&isct); en_o=rate(an,N,low&~isct)
print('rates',et_ct,et_o,en_ct,en_o)
et=np.where(isct,et_ct,et_o); en=np.where(isct,en_ct,en_o)
cand=np.where(((sa+an)>=2)&valid)[0]
S,SDc,A,NN,ET,EN,dT,dN,CN=sa[cand],sd[cand],an[cand],N[cand],et[cand],en[cand],DT[cand],DN[cand],cn[cand]
def lik(qt,qn): return binom.logpmf(S,SDc,np.clip(qt,1e-9,1-1e-9))+binom.logpmf(A,NN,np.clip(qn,1e-9,1-1e-9))
l0=lik(ET,EN)
ls=np.full(len(cand),-np.inf); xb=np.zeros(len(cand))
for x in np.concatenate([np.linspace(0.02,1,99),np.linspace(1.05,5,80)]):
    ok=x<=CN
    l=np.where(ok,lik(x*p/dT+ET,x*f/dN+EN),-np.inf); u=l>ls; ls[u]=l[u]; xb[u]=x
lg=np.full(len(cand),-np.inf)
for nk in range(0,6):
    ok=nk<=CN
    l=np.where(ok,lik((1-p+p*nk)/dT,(1-f+f*nk)/dN),-np.inf); lg=np.maximum(lg,l)
llr=ls-l0; gl=ls-lg
np.savez('/root/work/calls2.npz',cand=cand,llr=llr,gl=gl,xb=xb,a=a[cand],sa=S,sd=SDc,an=A,N=NN,isct=isct[cand],cn=CN)
# null expected FP above thresholds: approximate by simulation of noise counts at all valid sites (one replicate)
rng=np.random.default_rng(1)
vs=np.where(valid)[0]
reps=3; fps={}
for rep in range(reps):
    s0=rng.binomial(sd[vs],et[vs]); a0=rng.binomial(N[vs],en[vs])
    # 3 alternative alleles: approximate by max of 3 draws
    for k in range(2):
        s1=rng.binomial(sd[vs],et[vs]); a1=rng.binomial(N[vs],en[vs]); 
        better=(s1+a1)>(s0+a0); s0=np.where(better,s1,s0); a0=np.where(better,a1,a0)
    cc=np.where((s0+a0)>=2)[0]; vv=vs[cc]
    def lk(qt,qn): return binom.logpmf(s0[cc],sd[vv],np.clip(qt,1e-9,1-1e-9))+binom.logpmf(a0[cc],N[vv],np.clip(qn,1e-9,1-1e-9))
    L0=lk(et[vv],en[vv]); LS=np.full(len(cc),-np.inf)
    for x in np.linspace(0.02,1,99):
        LS=np.maximum(LS,lk(x*p/DT[vv]+et[vv],x*f/DN[vv]+en[vv]))
    nl=LS-L0
    for t in [6,8,10,12,14,16,20]: fps.setdefault(t,[]).append((nl>t).sum())
for t in fps: print('thr',t,'nullFP',np.mean(fps[t]),'calls',((llr>t)&(gl>5)).sum())
