import pandas as pd, numpy as np
d=pd.read_csv('/root/work/cand.tsv',sep='\t')
d['ta']=d.t_o0+d.t_o1; d['tv']=d.ta/d.t_dp; d['nv']=d.n_alt/d.n_dp
d['sub']=d.ref+'>'+d.alt
print(d.groupby('sub').agg(n=('ta','size'),o0=('t_o0','sum'),o1=('t_o1','sum')))
g=d[(d.nv<0.3)]
print('nv<0.3, by ta bins')
g=g.assign(tab=pd.cut(g.ta,[1,2,3,4,5,6,8,10,15,20,30,50,1000]))
print(pd.crosstab(g.tab,g['sub']))
# germline
print('germ', ((d.nv>0.3)).sum())
