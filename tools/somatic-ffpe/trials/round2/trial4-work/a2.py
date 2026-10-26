import numpy as np
r=np.load('/root/work/ref.npy')
T=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4]; Nn=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4]
dT=T.sum(1); dN=Nn.sum(1)
idx=np.arange(len(r))
def alt(c):
    c2=c.copy(); c2[idx[r<4],r[r<4]]=-1; a=c2.argmax(1); return a, c[idx,a]
aN,cN=alt(Nn); aT,cT=alt(T)
vN=cN/np.maximum(dN,1); vT=T[idx,aN]/np.maximum(dT,1)
het=(vN>0.2)&(vN<0.8)&(dN>30)
hom=(vN>=0.8)&(dN>30)
print('het',het.sum(),'hom',hom.sum())
np.save('/root/work/het.npy',np.where(het)[0])
B=5000
for st in range(0,1000000,B):
    m=slice(st,st+B); h=np.where(het[m])[0]+st
    bt=np.abs(vT[h]-0.5); bn=np.abs(vN[h]-0.5)
    print(st, round(dT[m].mean(),1), round(dN[m].mean(),1), round(dT[m].mean()/dN[m].mean(),3), len(h), round(np.median(bt),3) if len(h) else '', round(np.median(bn),3) if len(h) else '')
