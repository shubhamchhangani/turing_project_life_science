import numpy as np
r=np.load('/root/work/ref.npy')
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4]
s=np.where(r[1:-1]==1)[0]+1
alt=C[0][s,3]; dep=C[0][s].sum(1)
B='ACGTN'
ctx=r[s-1]*5+r[s+1]
for c in np.unique(ctx):
    m=ctx==c; fr=alt[m]/np.maximum(dep[m],1)
    print(B[c//5]+'C'+B[c%5], m.sum(), 'mean rate',round(alt[m].sum()/dep[m].sum(),4),'sd of frac',round(fr.std(),4), 'q90', round(np.quantile(fr,0.9),3))
# overdispersion: per-site, expected binomial sd
print('overall', alt.sum()/dep.sum())
