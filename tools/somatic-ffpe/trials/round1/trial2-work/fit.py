import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
def seg(r):
    if r.chrom=='chr3': return 'chr3a' if r.pos<=120000 else 'chr3b'
    return r.chrom
df['seg']=np.where(df.chrom!='chr3',df.chrom,np.where(df.pos<=120000,'chr3a','chr3b'))
# trim edges
d=df[(df.pos>2000)]
d=d[~((d.chrom=='chr3')&(d.pos>118000)&(d.pos<122000))]
CN={'chr1':4,'chr2':3,'chr3a':2,'chr3b':4,'chr4':5}
m=d.groupby('seg')[['tdp','ndp']].median()
mm=d.groupby('seg')[['tdp','ndp']].mean()
print(m,mm)
for col in ['tdp','ndp']:
  for M in (m,mm):
    x=np.array([CN[s] for s in M.index]); y=M[col].values
    A=np.polyfit(x,y,1); slope,inter=A
    print(col,'slope',slope,'int',inter,'frac',slope/(slope+inter/2))
