import numpy as np
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
B='ACGT'
for name,X in (('T',T),('N',N)):
    num=np.zeros((4,4,2)); den=np.zeros((4,2))
    for c in ['chr1','chr2','chr3','chr4']:
        x=X[c][:,:,:4].astype(float); ref=X[c+'_ref']; L=len(ref)
        nn=N[c][:,:,:4].sum(1); dn=nn.sum(1); nref=nn[np.arange(L),ref.clip(0,3)]
        ok=(ref<4)&((dn-nref)<=1)
        for r in range(4):
            m=ok&(ref==r)
            den[r]+=x[m].sum(2).sum(0)
            num[r]+=x[m].sum(0).T
    for r in range(4):
        print(name,B[r],' '.join(f'{B[a]}:{num[r,a,0]/den[r,0]:.2e}/{num[r,a,1]/den[r,1]:.2e}' for a in range(4) if a!=r))
