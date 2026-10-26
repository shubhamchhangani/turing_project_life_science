import numpy as np
r=np.load('/root/work/ref.npy'); idx=np.arange(len(r))
C=np.load('/root/work/tumor.cnt.npy'); Tc=C.sum(0)[:,:4]; Nc=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4]
T=Tc.sum(1); N=Nc.sum(1)
c2=Tc.copy(); c2[idx[r<4],r[r<4]]=-1; a=c2.argmax(1)
at=Tc[idx,a]; an=Nc[idx,a]; vt=at/np.maximum(T,1)
B='ACGT'
from collections import Counter
for lo,hi in [(1,1),(2,2),(3,4),(5,8),(9,15),(16,1000)]:
    s=np.where((at>=lo)&(at<=hi)&(an==0)&(r<4))[0]
    cnt=Counter(B[r[i]]+'>'+B[a[i]] for i in s)
    o0=C[0][s,a[s]]; o1=C[1][s,a[s]]
    print(lo,hi,len(s),dict(sorted(cnt.items())),'orient0 frac',(o0.sum()/(o0+o1).sum()).round(3))
print('per-site orientation')
for typ in [(1,3),(2,0)]:
    s=np.where((at>=3)&(an==0)&(r==typ[0])&(a==typ[1]))[0]
    f=C[0][s,a[s]]/at[s]; print(B[typ[0]],'>',B[typ[1]],np.histogram(f,bins=np.linspace(0,1,11))[0])
    rf0=C[0][s,r[s]]/(C[0][s,r[s]]+C[1][s,r[s]]); print('  ref orient0 frac mean',rf0.mean())
# likely somatic: an>=3
s=np.where((at>=5)&(an>=3)&(vt<0.45))[0]
f=C[0][s,a[s]]/at[s]; print('somatic-like',len(s),np.histogram(f,bins=np.linspace(0,1,11))[0])
