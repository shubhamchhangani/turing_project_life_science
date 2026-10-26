import numpy as np, pandas as pd
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
B='ACGT'; out=[]
for c in ['chr1','chr2','chr3','chr4']:
    t=T[c][:,:,:4].astype(np.int32); n=N[c][:,:,:4].sum(1).astype(np.int32); ref=T[c+'_ref']; L=len(ref)
    Cn=np.full(L,{'chr1':4,'chr2':3,'chr3':4,'chr4':5}[c]); 
    if c=='chr3': Cn[:120150]=2
    dn=n.sum(1); dt0=t[:,0].sum(1); dt1=t[:,1].sum(1)
    for r in range(4):
        for a in range(4):
            if a==r: continue
            m=(ref==r)
            ac0=t[:,0,a]; ac1=t[:,1,a]
            if (r,a)==(1,3): ac=ac1; Dc=dt1; ad=ac0; Dd=dt0
            elif (r,a)==(2,0): ac=ac0; Dc=dt0; ad=ac1; Dd=dt1
            else: ac=ac0+ac1; Dc=dt0+dt1; ad=np.zeros(L,int); Dd=np.zeros(L,int)
            sel=np.nonzero(m&((ac+n[:,a]+ad)>=1))[0]
            out.append(pd.DataFrame(dict(chrom=c,pos=sel+1,ref=B[r],alt=B[a],ac=ac[sel],Dc=Dc[sel],ad=ad[sel],Dd=Dd[sel],an=n[sel,a],Dn=dn[sel],Dt=(dt0+dt1)[sel],C=Cn[sel])))
df=pd.concat(out); df.to_parquet('/root/work/all.parq') if False else df.to_pickle('/root/work/all.pkl'); print(len(df))
