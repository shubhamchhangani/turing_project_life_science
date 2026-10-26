import pandas as pd,numpy as np
from scipy.stats import binom
df=pd.read_csv('/root/work/cand2.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['ct']=((df.ref=='C')&(df.alt=='T'))|((df.ref=='G')&(df.alt=='A'))
df['conc']=df.tconc+df.nconc; df['disc']=df.tdisc+df.ndisc
df['artp']=binom.sf(df.disc-1,df.conc+df.disc,0.03)
df['artok']=df.artp>1e-3
s=df[(~df.ct)&(df.vn<0.3)&df.artok]
print(pd.crosstab(s.talt.clip(0,12),s.nalt.clip(0,10)))
def C(r):
    if r.chrom=='chr3' and r.pos<=120000: return 2
    return {'chr1':4,'chr2':3,'chr3':4,'chr4':5}[r.chrom]
p=0.612
x=df[(df.vn<0.3)&df.artok&(df.talt+df.nalt>=5)&(~df.ct)].copy()
x['C']=x.apply(C,axis=1)
x['mcf']=x.vt*(p*x.C+2*(1-p))/p
x['mcfn']=x.vn*(0.203*x.C+2*(1-0.203))/0.203
print(np.histogram(x.mcf,bins=np.arange(0,3.6,0.1)))
print(x.sort_values('mcf')[['chrom','pos','ref','alt','talt','dt','nalt','dn','C','mcf','mcfn','conc','disc']].to_string())
