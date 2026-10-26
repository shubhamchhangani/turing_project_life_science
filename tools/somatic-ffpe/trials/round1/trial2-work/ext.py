import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
df=df[df.alt!='.']
df['sub']=df.ref+'>'+df.alt
ct=df['sub'].isin(['C>T','G>A'])
x=df[(df.talt<3)]
print(pd.crosstab(x.talt.clip(upper=3),x.nalt.clip(upper=6)).to_string())
print('nonct'); x2=x[~ct.loc[x.index]]
print(pd.crosstab(x2.talt,x2.nalt.clip(upper=6)).to_string())
