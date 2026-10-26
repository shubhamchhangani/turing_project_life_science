import numpy as np, pandas as pd
from scipy.stats import binom
p=0.607;f=0.205
df=pd.read_pickle('/root/work/post.pkl'); s=df[(df.chrom!='pseudo')&(df.an/df.Dn<0.3)&(df.post>0.5)].copy()
print(len(s))
gt=np.arange(0.1,3.01,0.1); gn=np.arange(0.0,4.01,0.1)
G=np.array([(a,b) for a in gt for b in gn])
ac=s.ac.values[:,None]; Dc=s.Dc.values[:,None]; an=s.an.values[:,None]; Dn=s.Dn.values[:,None]; C=s.C.values[:,None]
vt=p*G[:,0][None]/(2*(1-p)+p*C); vn=f*G[:,1][None]/(2*(1-f)+f*C)
L=binom.logpmf(ac,Dc,np.minimum(vt+1e-4,1))+binom.logpmf(an,Dn,np.minimum(vn+1.6e-4,1))
mx=L.max(1,keepdims=True); E=np.exp(L-mx)
w=np.ones(len(G))/len(G)
for it in range(2000):
    R=E*w; R/=R.sum(1,keepdims=True); w=R.mean(0)
idx=np.argsort(-w)
for i in idx[:40]: print(G[i].round(2), round(w[i]*len(s),1))
np.save('/root/work/G.npy',G); np.save('/root/work/w.npy',w)
