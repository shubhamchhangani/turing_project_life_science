import pandas as pd,numpy as np
df=pd.read_csv('/root/work/cand.tsv',sep='\t')
df['mt']=df.ref+'>'+df.alt
for mt in ['C>T','G>A']:
    x=df[(df.mt==mt)&(df.nalt==0)&(df.talt>=4)]
    print(mt,len(x),'F1 sum',x.tF1.sum(),'F2 sum',x.tF2.sum())
    print(pd.crosstab(x.tF1.clip(0,8),x.tF2.clip(0,8)))
