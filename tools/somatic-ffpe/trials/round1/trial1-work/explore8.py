import pandas as pd,numpy as np
d=pd.read_csv('/root/work/scored.tsv',sep='\t')
d['ct']=d.C2T|d.G2A
d['tot']=d['at']+d.nalt
base=d[(d.vn<0.35)&(d.artp>1e-3)]
print('incons distribution among tot>=5:',np.histogram(base[base.tot>=5].incons.fillna(99),bins=[0,2,5,10,15,20,30,50,1e9])[0])
b=base[base.incons<15]
for ct in [False,True]:
    x=b[b.ct==ct]
    print('CT' if ct else 'nonCT')
    print(pd.crosstab(x['at'].clip(0,10),x.nalt.clip(0,8)))
print(base[(base.incons>=15)&(base.tot>=4)][['chrom','pos','ref','alt','tF1','tF2','at','dtc','nalt','dn','C','mcf','incons','conc','disc']].to_string())
