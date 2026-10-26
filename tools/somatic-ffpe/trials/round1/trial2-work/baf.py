import pandas as pd, numpy as np
from scipy.stats import binom
df=pd.read_pickle('/root/work/all.pkl')
df['seg']=np.where(df.chrom!='chr3',df.chrom,np.where(df.pos<=120000,'chr3a','chr3b'))
df['nv']=df.nalt/df.ndp.clip(lower=1); df['tv']=df.talt/df.tdp.clip(lower=1)
het=df[(df.ndp>=30)&(df.nv>0.25)&(df.nv<0.75)&(df.tdp>=30)]
AB={'chr1':(2,2),'chr2':(2,1),'chr3a':(2,0),'chr3b':(2,2),'chr4':(3,2)}
def ll(p,col_alt,col_dp,segs):
    t=0
    for s in segs:
        h=het[het.seg==s]; A,B=AB[s]
        v1=(p*A+1-p)/(p*(A+B)+2-2*p); v2=1-v1
        L=0.5*binom.pmf(h[col_alt],h[col_dp],v1)+0.5*binom.pmf(h[col_alt],h[col_dp],v2)
        t+=np.log(L+1e-300).sum()
    return t
ps=np.linspace(0.01,0.99,197)
for segs in [['chr2','chr3a','chr4'],['chr2'],['chr3a'],['chr4']]:
    lt=[ll(p,'talt','tdp',segs) for p in ps]; ln=[ll(p,'nalt','ndp',segs) for p in ps]
    print(segs,'tumor p',ps[np.argmax(lt)],'normal f',ps[np.argmax(ln)])
