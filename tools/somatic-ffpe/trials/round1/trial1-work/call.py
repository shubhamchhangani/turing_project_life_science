import pandas as pd,numpy as np
from scipy.stats import binom
from scipy.optimize import minimize_scalar
p,q=0.612,0.203
df=pd.read_csv('/root/work/cand2.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['C2T']=(df.ref=='C')&(df.alt=='T'); df['G2A']=(df.ref=='G')&(df.alt=='A')
df['conc']=df.tconc+df.nconc; df['disc']=df.tdisc+df.ndisc
df['artp']=binom.sf(df.disc-1,df.conc+df.disc,0.03)
def segC(c,pos):
    if c=='chr3' and pos<=120000: return 2
    return {'chr1':4,'chr2':3,'chr3':4,'chr4':5}[c]
df['C']=[segC(c,x) for c,x in zip(df.chrom,df.pos)]
df['at']=np.where(df.C2T,df.tF2,np.where(df.G2A,df.tF1,df.talt))
df['dtc']=np.where(df.C2T|df.G2A,df.dt/2,df.dt).round().astype(int)
et,en=1.0e-4,1.7e-4
def nll(x,r):
    ft=p*x/(p*r.C+2*(1-p))+et; fn=q*x/(q*r.C+2*(1-q))+en
    return -(binom.logpmf(r['at'],r.dtc,min(ft,1))+binom.logpmf(r.nalt,r.dn,min(fn,1)))
def sep_nll(r):
    # separate fits
    ft=max(r['at']/max(r.dtc,1),et); fn=max(r.nalt/max(r.dn,1),en)
    return -(binom.logpmf(r['at'],r.dtc,ft)+binom.logpmf(r.nalt,r.dn,fn))
pre=df[(df['at']+df.nalt>=2)].copy()
res=[]
for i,r in pre.iterrows():
    o=minimize_scalar(lambda x: nll(x,r),bounds=(0,r.C),method='bounded')
    l0=nll(0,r); l1=o.fun; ls=sep_nll(r)
    res.append((o.x,2*(l0-l1),2*(l1-ls)))
pre[['mcf','LR','incons']]=res
pre.to_csv('/root/work/scored.tsv',sep='\t',index=False)
