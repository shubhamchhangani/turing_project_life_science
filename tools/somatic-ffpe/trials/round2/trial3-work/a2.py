import pandas as pd, numpy as np
d=pd.read_csv('/root/work/cand.tsv',sep='\t')
d['nv']=d.n_alt/d.n_dp
d=d[d.nv<0.3]
ct=d[(d.ref=='C')&(d.alt=='T')]
print(pd.crosstab(ct.t_o1.clip(upper=12),ct.t_o0.clip(upper=15)))
ga=d[(d.ref=='G')&(d.alt=='A')]
print(pd.crosstab(ga.t_o0.clip(upper=12),ga.t_o1.clip(upper=15)))
o=d[~(((d.ref=='C')&(d.alt=='T'))|((d.ref=='G')&(d.alt=='A')))]
print(pd.crosstab(o.t_o1.clip(upper=12),o.t_o0.clip(upper=12)))
print(pd.crosstab((o.t_o0+o.t_o1).clip(upper=30),o.n_alt.clip(upper=10)))
