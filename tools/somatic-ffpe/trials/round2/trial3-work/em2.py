import numpy as np, pandas as pd
from scipy.stats import binom
p=0.607; f=0.205; et=1.0e-4; en=1.65e-4
df=pd.read_pickle('/root/work/post.pkl')
W=df.w.values
print('rows',len(df))
gm=np.array([0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,3.0])
ga=np.arange(0.05,1.06,0.1); gb=np.arange(0.4,3.01,0.2)
NG=np.array([(a,b) for a in ga for b in gb if b>=max(0.5,1.8*a)])
mt=np.concatenate([gm,NG[:,0]]); mn=np.concatenate([gm,NG[:,1]]); ismal=np.r_[np.ones(len(gm),bool),np.zeros(len(NG),bool)]
ac=df.ac.values[:,None]; Dc=df.Dc.values[:,None]; an=df.an.values[:,None]; Dn=df.Dn.values[:,None]; C=df.C.values[:,None]
L1=np.empty((len(df),len(mt)),np.float32)
for i in range(0,len(df),10000):
    sl=slice(i,i+10000)
    vt=p*mt[None]/(2*(1-p)+p*C[sl]); vn=f*mn[None]/(2*(1-f)+f*C[sl])
    x=binom.logpmf(ac[sl],Dc[sl],np.minimum(vt+et,1))+binom.logpmf(an[sl],Dn[sl],np.minimum(vn+en,1))
    x[:,ismal]=np.where(gm[None,:]<=C[sl]+1e-9,x[:,ismal],-np.inf)
    L1[sl]=x
L0=binom.logpmf(ac,Dc,et)[:,0]+binom.logpmf(an,Dn,en)[:,0]
clv=df.cl.values; classes=sorted(set(clv))
w=np.ones(len(mt))/len(mt); pi={k:0.01 for k in classes}
for it in range(150):
    pv=np.array([pi[k] for k in clv])
    a=np.log(pv)[:,None]+np.log(w+1e-300)[None,:]+L1; b=np.log1p(-pv)+L0
    mx=np.maximum(a.max(1),b); R=np.exp(a-mx[:,None]); r0=np.exp(b-mx); den=R.sum(1)+r0; R/=den[:,None]
    post=R.sum(1); w=(R*W[:,None]).sum(0); w/=w.sum()
    for k in classes: m_=clv==k; pi[k]=(post[m_]*W[m_]).sum()/W[m_].sum()
df['psom']=post; df['pmal']=R[:,ismal].sum(1); df['pnon']=R[:,~ismal].sum(1)
print('malignant weight',w[ismal].sum(), 'expected n mal',df.pmal.sum(),'non',df.pnon.sum())
for i in np.argsort(-w)[:25]: print('mal' if ismal[i] else 'non',mt[i].round(2),mn[i].round(2),round(w[i]*post.sum(),1))
df=df[df.chrom!='pseudo']; df.to_pickle('/root/work/post2.pkl'); print(pi)
for th in [0.3,0.5,0.7,0.9]:
    s=df[df.pmal>th]; print(th,len(s),'expTP',s.pmal.sum(),'prec',s.pmal.mean(),'recall',s.pmal.sum()/df.pmal.sum())
