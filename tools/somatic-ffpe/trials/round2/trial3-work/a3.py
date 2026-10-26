import pandas as pd, numpy as np
d=pd.read_csv('/root/work/cand.tsv',sep='\t')
d['nv']=d.n_alt/d.n_dp
p=0.607;f=0.205
def C(r):
    if r.chrom=='chr1': return 4
    if r.chrom=='chr2': return 3
    if r.chrom=='chr4': return 5
    return 2 if r.pos<=120150 else 4
d=d[(d.nv<0.3)]
o=d[~(((d.ref=='C')&(d.alt=='T'))|((d.ref=='G')&(d.alt=='A')))].copy()
o['ta']=o.t_o0+o.t_o1
o=o[(o.ta>=3)&(o.n_alt>=1) | (o.ta>=6)]
o['C']=o.apply(C,axis=1)
o['tv']=o.ta/o.t_dp
o['mccf']=o.tv*(2*(1-p)+p*o.C)/p   # mutation copies * ccf
o['nmccf']=o.nv*(2*(1-f)+f*o.C)/f
print(len(o))
print(np.histogram(o.mccf,bins=np.arange(0,3.6,0.1)))
print(np.histogram(o.nmccf,bins=np.arange(0,3.6,0.2)))
print(o.groupby('C').size())
