import numpy as np, pandas as pd
from scipy.stats import binom
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.ac>=5)&(df.an/df.Dn<0.3)&(df.post>0.5)]
print(len(s))
grid=np.arange(0.1,3.01,0.05)
ac=s.ac.values[:,None]; Dc=s.Dc.values[:,None]; an=s.an.values[:,None]; Dn=s.Dn.values[:,None]; C=s.C.values[:,None]
def ll(p,f):
    vt=p*grid/(2*(1-p)+p*C); vn=f*grid/(2*(1-f)+f*C)
    LT=binom.logpmf(ac,Dc,vt+1e-4); LN=binom.logpmf(an,Dn,vn+1.6e-4)
    LT=np.where(grid[None,:]<=C+1e-9,LT,-np.inf)
    w=np.ones(len(grid))/len(grid)
    for it in range(300):
        a=LT+LN+np.log(w+1e-300); mx=a.max(1,keepdims=True); R=np.exp(a-mx); tot=R.sum(1,keepdims=True); R/=tot
        w=R.mean(0)
    return (np.log(tot[:,0])+mx[:,0]).sum()
for p in (0.58,0.607,0.63):
    r=[(ll(p,f),f) for f in np.arange(0.16,0.30,0.01)]
    print(p,max(r)[1],[round(x[0]-max(r)[0],1) for x in r])
