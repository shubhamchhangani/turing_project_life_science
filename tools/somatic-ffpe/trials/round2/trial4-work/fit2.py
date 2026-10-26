import numpy as np
from scipy.stats import binom
r=np.load('/root/work/ref.npy'); idx=np.arange(len(r))
Tc=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4]; Nc=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4]
T=Tc.sum(1); N=Nc.sum(1)
S=Tc+Nc; c2=S.copy(); c2[idx[r<4],r[r<4]]=-1; a=c2.argmax(1)
vt=Tc[idx,a]/np.maximum(T,1); vn=Nc[idx,a]/np.maximum(N,1)
# select hets by tumor in balanced segs and by tumor VAF elsewhere? use combined: both samples alt>=3 and tumor vaf in 0.1-0.9 and normal vaf 0.1-0.9
het=(Tc[idx,a]>=5)&(Nc[idx,a]>=3)&(vt>0.1)&(vt<0.9)&(vn>0.1)&(vn<0.9)&(r<4)
print(het.sum())
segs=[(0,300000,4,2),(300000,550000,3,1),(550000,670143,2,0),(670143,800000,4,2),(800000,1000000,5,2)]
for a0,b0,c,nb in segs:
    h=np.where(het[a0:b0])[0]+a0
    print(c,nb,len(h),'T mean vaf',vt[h].mean().round(3),'sd',vt[h].std().round(3),'binom sd',np.sqrt(0.25/T[h]).mean().round(3),'| N mean',vn[h].mean().round(3),'sd',vn[h].std().round(3),'binsd',np.sqrt(0.25/N[h]).mean().round(3), 'N mean|dev|',np.abs(vn[h]-.5).mean().round(3))
    # Mirror: average of max allele dev using tumour phase (tumour tells which allele is major)
    sgn=np.sign(vt[h]-0.5)
    print('   normal dev phased by tumour:',((vn[h]-0.5)*sgn).mean().round(4))
