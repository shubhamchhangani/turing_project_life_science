import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
c=pd.read_pickle('/root/work/cand.pkl')
c['nv']=c.nalt/c.ndp.clip(lower=1); c['tv']=c.talt/c.tdp.clip(lower=1)
bins=[0,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.6,1.01]
print(pd.crosstab(c['sub'], pd.cut(c.tv,bins)))
print(pd.crosstab(c['sub'], pd.cut(c.nv,[-0.01,0,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,1.01])))
