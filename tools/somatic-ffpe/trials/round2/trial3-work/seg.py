import numpy as np
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
B=2500
for c in ['chr1','chr2','chr3','chr4']:
    dt=T[c][:,:,:4].sum((1,2)); dn=N[c][:,:,:4].sum((1,2))
    a=dt[:len(dt)//B*B].reshape(-1,B).mean(1); b=dn[:len(dn)//B*B].reshape(-1,B).mean(1)
    print(c,' '.join(f'{x:.0f}' for x in a))
    print(c,' '.join(f'{x:.0f}' for x in b))
