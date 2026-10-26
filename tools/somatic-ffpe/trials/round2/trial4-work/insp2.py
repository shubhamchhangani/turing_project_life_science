import numpy as np
r=np.load('/root/work/ref.npy'); n=len(r); idx=np.arange(n)
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4]; NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].sum(0)
for ref,alt,safe in [(1,3,1),(2,0,0)]:
    m=np.where((r==ref)&(C[safe][:,alt]>=3)&(NC[:,alt]<=1))[0]
    for i in m:
        print(i,'CT' if ref==1 else 'GA','safe',C[safe][i,alt],'dmg',C[1-safe][i,alt],'safedepth',C[safe][i].sum(),'dmgdepth',C[1-safe][i].sum(),'normal',NC[i,alt])
