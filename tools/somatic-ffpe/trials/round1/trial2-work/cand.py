import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
df['nv']=df.nalt/df.ndp.clip(lower=1); df['tv']=df.talt/df.tdp.clip(lower=1)
c=df[(df.talt>=3)]
c=c.assign(sub=c.ref+'>'+c.alt)
print(len(c))
print(pd.crosstab(c['sub'], pd.cut(c.tv,[0,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.6,1.01])))
cc=c[c.nalt==0]
print(pd.crosstab(cc['sub'], pd.cut(cc.tv,[0,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.6,1.01])))
cc=c[(c.nalt>0)&(c.nv<0.3)]
print(pd.crosstab(cc['sub'], pd.cut(cc.tv,[0,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.6,1.01])))
c.to_pickle('/root/work/cand.pkl')
