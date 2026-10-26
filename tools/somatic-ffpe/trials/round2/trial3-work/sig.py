import pandas as pd, numpy as np, pysam
d=pd.read_pickle('/root/work/post3_0.607_0.205.pkl')
fa=pysam.FastaFile('/root/data/reference.fa')
d['ctx']=[fa.fetch(c,p-2,p+1).upper() for c,p in zip(d.chrom,d.pos)]
m=d[d.pmal>0.9]; n=d[d.pnon>0.9]
print(len(m),len(n))
print(pd.concat([m.cl.value_counts().rename('mal'),n.cl.value_counts().rename('non')],axis=1))
def cpg(x): return ((x.ref=='C')&(x.ctx.str[2]=='G'))|((x.ref=='G')&(x.ctx.str[0]=='C'))
print('CpG frac mal',cpg(m).mean(),'non',cpg(n).mean())
print(m.groupby('chrom').size(), n.groupby('chrom').size())
