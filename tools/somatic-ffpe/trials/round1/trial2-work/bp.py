import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
x=df[(df.chrom=='chr3')&(df.pos>110000)&(df.pos<130000)].tdp.values
best=None
for i in range(200,len(x)-200):
    s=((x[:i]-x[:i].mean())**2).sum()+((x[i:]-x[i:].mean())**2).sum()
    if best is None or s<best[0]: best=(s,i)
print(110001+best[1])
