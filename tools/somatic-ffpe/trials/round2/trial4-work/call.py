import numpy as np
from scipy.stats import binom
r=np.load('/root/work/ref.npy'); n=len(r); idx=np.arange(n)
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4].astype(np.int64); NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].astype(np.int64)
Tc=C.sum(0); Nc=NC.sum(0); T=Tc.sum(1); N=Nc.sum(1)
# safe counts per alt allele
SA=Tc.copy(); SD=np.repeat(T[:,None],4,1)
m=(r==1); SA[m,3]=C[1][m,3]; SD[m,3]=C[1][m].sum(1)   # C>T use orient1
m=(r==2); SA[m,0]=C[0][m,0]; SD[m,0]=C[0][m].sum(1)   # G>A use orient0
score=SA+Nc; score[idx[r<4],r[r<4]]=-1; a=score.argmax(1)
sa=SA[idx,a]; sd=SD[idx,a]; an=Nc[idx,a]
p,f=0.615,0.2
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(n)
for s0,e0,c in segs: cn[s0:e0]=c
K=(f/p)*(p*cn+2*(1-p))/(f*cn+2*(1-f))
isct=((r==1)&(a==3))|((r==2)&(a==0))
et=np.where(isct,1.5e-4,2.5e-4)/3*3  # per-fragment error to specific alt approx
et=np.where(isct,1.5e-4,1e-4); en=2e-4
vg=np.linspace(0.005,0.8,160)
cand=np.where(((sa>=2)|(an>=2))&(r<4))[0]
print('cand',len(cand))
ll0=binom.logpmf(sa[cand],sd[cand],et[cand])+binom.logpmf(an[cand],N[cand],en)
best=np.full(len(cand),-np.inf); bv=np.zeros(len(cand))
for v in vg:
    l=binom.logpmf(sa[cand],sd[cand],v+et[cand])+binom.logpmf(an[cand],N[cand],np.minimum(K[cand]*v+en,1))
    upd=l>best; best[upd]=l[upd]; bv[upd]=v
llr=best-ll0
vn=an[cand]/np.maximum(N[cand],1)
germ=vn>0.3
np.savez('/root/work/calls.npz',cand=cand,llr=llr,bv=bv,germ=germ,a=a[cand],sa=sa[cand],sd=sd[cand],an=an[cand],N=N[cand],isct=isct[cand],cn=cn[cand])
h=np.histogram(llr[~germ],bins=[-100,0,2,4,6,8,10,12,15,20,30,50,1e9])
print(list(zip(h[1][:-1],h[0])))
