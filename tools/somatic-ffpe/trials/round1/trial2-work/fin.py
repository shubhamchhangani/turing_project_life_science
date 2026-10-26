import pandas as pd, numpy as np
c=pd.read_pickle('/root/work/cf4.pkl')
s=c[c.call&(c.dll>0)]
print(len(s)); print(s['sub'].value_counts())
print(pd.crosstab(s['sub'].isin(['C>T','G>A']), pd.cut(s.tv,[0,0.03,0.05,0.08,0.12,0.2,0.3,1])))
