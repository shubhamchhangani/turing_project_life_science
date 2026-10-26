import numpy as np, pandas as pd
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.ac>=4)&(df.an/df.Dn<0.3)&(df.cl!="C>T")].copy()
s['vt']=s.ac/s.Dc
for C,g in s.groupby('C'):
    g=g.assign(b=pd.cut(g.vt,[0,0.04,0.06,0.08,0.1,0.13,0.16,0.2,0.25,0.3,0.4,0.6,1]))
    print('C',C)
    print(g.groupby('b',observed=True).apply(lambda x: pd.Series(dict(n=len(x),vt=x.ac.sum()/x.Dc.sum(),vn=x.an.sum()/x.Dn.sum(),r=(x.an.sum()/x.Dn.sum())/(x.ac.sum()/x.Dc.sum())))).round(3).to_string())
