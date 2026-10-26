import numpy as np
r=np.load('/root/work/ref.npy'); idx=np.arange(len(r))
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4]; NC=np.load('/root/work/normal.cnt.npy')[:,:,:4]
Tc=C.sum(0); Nc=NC.sum(0)
T=Tc.sum(1); N=Nc.sum(1)
S=Tc+Nc
c2=S.copy(); c2[idx[r<4],r[r<4]]=-1; a=c2.argmax(1)
ct=(r==1)&(a==3); ga=(r==2)&(a==0)
# safe counts
sa=np.where(ct,C[1][idx,a],np.where(ga,C[0][idx,a],Tc[idx,a]))
sd=np.where(ct,C[1].sum(1),np.where(ga,C[0].sum(1),T))
an=Nc[idx,a]; at=Tc[idx,a]
np.savez('/root/work/cand.npz',a=a,sa=sa,sd=sd,an=an,at=at,T=T,N=N,ct=ct|ga)
for lab,m in [('nonCT',~(ct|ga)),('CT',ct|ga)]:
    print(lab)
    H=np.zeros((8,8),int)
    for i in range(8):
        for j in range(8):
            H[i,j]=((np.minimum(sa,7)==i)&(np.minimum(an,7)==j)&m&(r<4)).sum()
    print(' rows tumour safe alt 0..7+, cols normal alt 0..7+'); print(H)
