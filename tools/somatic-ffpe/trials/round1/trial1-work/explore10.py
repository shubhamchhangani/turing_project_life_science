import pandas as pd,numpy as np
d=pd.read_csv('/root/work/scored.tsv',sep='\t')
d['ct']=d.C2T|d.G2A
d['tot']=d['at']+d.nalt
x=d[(d.vn<0.35)&(d.artp<=1e-3)&(d.incons<15)&(d.tot>=4)]
print(len(x))
print(x.sort_values('conc')[['chrom','pos','ref','alt','tF1','tF2','at','dtc','nalt','dn','C','mcf','incons','tconc','tdisc','nconc','ndisc','tlow']].tail(40).to_string())
