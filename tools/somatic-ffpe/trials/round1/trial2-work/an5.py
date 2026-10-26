import pandas as pd, numpy as np
c=pd.read_pickle('/root/work/cf2.pkl')
p,f=0.618,0.205
c['seg']=np.where(c.chrom!='chr3',c.chrom,np.where(c.pos<=120000,'chr3a','chr3b'))
CN={'chr1':4,'chr2':3,'chr3a':2,'chr3b':4,'chr4':5}
c['CN']=c.seg.map(CN)
ct=c['sub'].isin(['C>T','G>A'])
d=c[~ct&(c.nv<0.27)&(c.tv<0.37)|(~ct&(c.seg=='chr3a')&(c.nv<0.27))].copy()
d['tccf']=d.tv*(p*d.CN+2-2*p)/p
d['expn']=d.tccf*f/(f*d.CN+2-2*f)*d.normal_nfr
d['bin']=pd.cut(d.tccf,[0,0.3,0.5,0.7,0.9,1.2,1.6,2.2,4])
g=d.groupby('bin').agg(n=('expn','size'),obs=('normal_afr','sum'),exp=('expn','sum'))
g['ratio']=g.obs/g.exp; print(g)
print(np.histogram(d.tccf,bins=np.arange(0,3.6,0.1)))
