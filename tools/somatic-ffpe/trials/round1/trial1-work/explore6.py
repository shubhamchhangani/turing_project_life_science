import pandas as pd,numpy as np
df=pd.read_csv('/root/work/cand2.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['ct']=((df.ref=='C')&(df.alt=='T'))|((df.ref=='G')&(df.alt=='A'))
s=df[(~df.ct)&(df.vn<0.3)&(df.talt+df.nalt>=4)]
s=s.assign(conc=s.tconc+s.nconc,disc=s.tdisc+s.ndisc)
print(pd.crosstab(s.conc.clip(0,12),s.disc.clip(0,12)))
g=df[(df.vn>0.35)]
g=g.assign(conc=g.tconc+g.nconc,disc=g.tdisc+g.ndisc)
print('germline conc/disc',g.conc.sum(),g.disc.sum())
