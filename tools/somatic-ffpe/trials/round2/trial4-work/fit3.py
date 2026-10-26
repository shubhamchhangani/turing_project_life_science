import numpy as np
r=np.load('/root/work/ref.npy'); idx=np.arange(len(r))
Tc=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4]; Nc=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4]
T=Tc.sum(1); N=Nc.sum(1)
S=Tc+Nc; c2=S.copy(); c2[idx[r<4],r[r<4]]=-1; a=c2.argmax(1)
vt=Tc[idx,a]/np.maximum(T,1); vn=Nc[idx,a]/np.maximum(N,1)
m=(np.r_[0:300000]) 
sel=(Tc[idx,a]>=3)|(Nc[idx,a]>=3)
for name,(a0,b0) in [('bal',(0,300000)),('loh',(550000,670143))]:
    s=np.where(sel[a0:b0])[0]+a0
    H,_,_=np.histogram2d(vt[s],vn[s],bins=[np.linspace(0,1,21),np.linspace(0,1,21)])
    print(name,len(s)); print('rows=tumorVAF bins 0.05, cols=normal VAF')
    for i in range(20): print(f'{i*0.05:.2f}',' '.join(f'{int(x):4d}' for x in H[i]))
