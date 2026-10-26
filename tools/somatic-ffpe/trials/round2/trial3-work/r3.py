import numpy as np, pandas as pd
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.an/df.Dn<0.3)&(df.cl!="C>T")].copy()
s['vt']=s.ac/s.Dc; s['vn']=s.an/s.Dn
s['mt']=s.vt*(2*(1-0.607)+0.607*s.C)/0.607
s['mn']=s.vn*(2*(1-0.205)+0.205*s.C)/0.205
x=s[(s.ac>=3)|(s.an>=3)]
print(pd.crosstab(x.ac.clip(upper=25),x.an.clip(upper=20)).to_string())
y=s[(s.ac>=5)]
print(pd.crosstab(pd.cut(y.mt,[0,0.15,0.2,0.25,0.3,0.35,0.45,0.6,0.8,1.2,1.6,2.5,4]),pd.cut(y.mn,[-0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2,2.5,3,5])).to_string())
