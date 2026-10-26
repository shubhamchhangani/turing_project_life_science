import numpy as np, pandas as pd
from scipy.stats import binom
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.post>0.9)]
grid=np.arange(0.1,3.01,0.05)
ac=s.ac.values[:,None]; Dc=s.Dc.values[:,None]; an=s.an.values[:,None]; Dn=s.Dn.values[:,None]; C=s.C.values[:,None]
def ll(p,f):
    vt=p*grid/(2*(1-p)+p*C); vn=f*grid/(2*(1-f)+f*C)
    L=binom.logpmf(ac,Dc,vt+1e-4)+binom.logpmf(an,Dn,vn+1.6e-4); L[:, :]=np.where(grid[None,:]<=C+1e-9,L,-np.inf)
    w=np.ones(len(grid))/len(grid)
    for it in range(200):
        a=L+np.log(w+1e-300); mx=a.max(1,keepdims=True); R=np.exp(a-mx); tot=R.sum(1,keepdims=True); R/=tot
        w=R.mean(0)
    return (np.log(tot[:,0])+mx[:,0]).sum(), w
res=[]
for p in np.arange(0.55,0.67,0.01):
    for f in np.arange(0.16,0.26,0.01):
        l,w=ll(p,f); res.append((l,p,f))
res.sort(reverse=True); print(res[:8])
l,w=ll(res[0][1],res[0][2]); print({round(g,2):round(x,3) for g,x in zip(grid,w) if x>0.005})
