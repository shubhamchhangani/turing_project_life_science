import numpy as np, pandas as pd
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.post>0.9)].copy()
s['vt']=s.ac/s.Dc; s['vn']=s.an/s.Dn
for C,g in s.groupby('C'):
    hi=g[g.vt>np.median(g.vt)]
    print(C,len(g),'pooled vt',g.ac.sum()/g.Dc.sum(),'vn',g.an.sum()/g.Dn.sum(),'ratio',(g.an.sum()/g.Dn.sum())/(g.ac.sum()/g.Dc.sum()))
    print(np.histogram(g.vt,bins=np.arange(0,0.85,0.025))[0])
print('restricted to high vt')
for C,g in s.groupby('C'):
    thr={2:0.25,3:0.17,4:0.14,5:0.12}[C]
    g=g[g.vt>thr]
    rt=g.ac.sum()/g.Dc.sum(); rn=g.an.sum()/g.Dn.sum()
    print(C,len(g),rt,rn,rn/rt)
    for p in (0.607,):
        # solve f from ratio
        r=rn/rt; x=r*p/(2-2*p+p*C); f=2*x/(1+2*x-x*C); print('  f',f)
