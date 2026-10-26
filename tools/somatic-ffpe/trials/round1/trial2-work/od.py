import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
df['seg']=np.where(df.chrom!='chr3',df.chrom,np.where(df.pos<=120000,'chr3a','chr3b'))
df['nv']=df.nalt/df.ndp.clip(lower=1); df['tv']=df.talt/df.tdp.clip(lower=1)
for s in ['chr1','chr3b','chr2','chr3a','chr4']:
    h=df[(df.seg==s)&(df.nv>0.2)&(df.nv<0.8)&(df.ndp>40)]
    print(s,len(h))
    print(' tumor', np.histogram(h.tv,bins=np.arange(0,1.05,0.05))[0])
    print(' normal',np.histogram(h.nv,bins=np.arange(0,1.05,0.05))[0])
    if s in('chr1','chr3b'):
        print(' tv var',h.tv.var(),'binom',(0.25/h.tdp).mean(),' nv var',h.nv.var(),'binom',(0.25/h.ndp).mean(), 'means',h.tv.mean(),h.nv.mean())
