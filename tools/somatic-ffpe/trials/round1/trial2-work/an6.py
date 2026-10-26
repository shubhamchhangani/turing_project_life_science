import pandas as pd, numpy as np
pd.set_option('display.width',250)
c=pd.read_pickle('/root/work/cf2.pkl')
ct=c['sub'].isin(['C>T','G>A'])
d=c[ct&(c.nv<0.3)].copy()
d['anti']=np.where(d['sub']=='C>T',d.tumor_afr-d.a_f1r2,d.a_f1r2)
d['pro']=d.tumor_afr-d.anti
print(pd.crosstab(d.anti.clip(upper=12),d.normal_afr.clip(upper=8)))
x=c[~ct&(c.nv<0.3)]
print('nonCT: talt vs normal'); print(pd.crosstab(x.tumor_afr.clip(upper=20)//2*2,x.normal_afr.clip(upper=8)))
