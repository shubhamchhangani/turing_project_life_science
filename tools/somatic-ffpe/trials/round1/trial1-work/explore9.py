import pandas as pd,numpy as np
d=pd.read_csv('/root/work/scored.tsv',sep='\t')
d['ct']=d.C2T|d.G2A
d['tot']=d['at']+d.nalt
d['art1']=np.where(d.C2T,d.tF1,np.where(d.G2A,d.tF2,0))
base=d[(d.vn<0.35)&(d.artp>1e-3)&(d.incons<15)]
x=base[base.ct&(base.tot==2)]
print(pd.crosstab([x['at'],x.nalt],x.art1.clip(0,10)))
calls=base[(base.ct&(base.tot>=3))|(~base.ct&(base.tot>=4))]
print(len(calls), calls.ct.sum())
print(np.histogram(calls.mcf,bins=np.arange(0,5.1,0.1))[0])
print(calls[(calls.incons>5)][['chrom','pos','ref','alt','tF1','tF2','at','dtc','nalt','dn','C','mcf','incons','conc','disc']].to_string())
dup=calls.groupby(['chrom','pos']).size(); print('multi',(dup>1).sum())
from scipy.stats import poisson
for ct in [False,True]:
    y=calls[(calls.ct==ct)&(calls.mcf<0.6)]
    print('ct',ct,len(y),np.bincount(y.tot.clip(0,20)))
    # MLE of lambda under truncated poisson
    k=4 if not ct else 3
    ls=np.linspace(1,15,141)
    L=[ (poisson.logpmf(y.tot.clip(0,20),l)-np.log(poisson.sf(k-1,l))).sum() for l in ls]
    lb=ls[int(np.argmax(L))]; print(' lambda',lb,' detect prob',poisson.sf(k-1,lb), 'est missing',len(y)*(1/poisson.sf(k-1,lb)-1))
p,q=0.612,0.203
calls=calls.assign(kt=calls.dtc*p/(p*calls.C+2*(1-p)), kn=calls.dn*q/(q*calls.C+2*(1-q)))
for ct in [False,True]:
  for lo,hi in [(0,0.6),(0.6,1.5),(1.5,5)]:
    y=calls[(calls.ct==ct)&(calls.mcf>=lo)&(calls.mcf<hi)]
    print(ct,lo,len(y),'x_t=%.3f x_n=%.3f'%(y['at'].sum()/y.kt.sum(), y.nalt.sum()/y.kn.sum()), 'nalt0 frac %.2f'%(y.nalt==0).mean(), 'tF1/tF2 %d/%d'%(y.tF1.sum(),y.tF2.sum()))
