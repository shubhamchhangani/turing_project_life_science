import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
B=2000
df['bin']=(df.pos-1)//B
g=df.groupby(['chrom','bin']).agg(t=('tdp','mean'),n=('ndp','mean')).reset_index()
for c in ['chr1','chr2','chr3','chr4']:
    x=g[g.chrom==c]
    print(c,' '.join(f'{v:.0f}' for v in x.t))
    print(c,' '.join(f'{v:.0f}' for v in x.n))
