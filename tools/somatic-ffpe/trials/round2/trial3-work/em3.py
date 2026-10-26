import numpy as np, pandas as pd, sys
from scipy.stats import binom
from scipy.optimize import minimize
p=float(sys.argv[1]) if len(sys.argv)>1 else 0.607; f=float(sys.argv[2]) if len(sys.argv)>2 else 0.205
et=1.0e-4; en=1.65e-4
df=pd.read_pickle('/root/work/post.pkl')
W=df.w.values; ac=df.ac.values.astype(float); Dc=df.Dc.values.astype(float); an=df.an.values.astype(float); Dn=df.Dn.values.astype(float); C=df.C.values.astype(float)
clv=df.cl.values; classes=sorted(set(clv))
comps=[('mal',0.25,0.25),('mal',0.55,0.55),('mal',1,1),('mal',1.5,1.5),('mal',2,2),('mal',3,3),('non',0.3,1.1),('non',0.6,1.5),('non',0.1,0.8)]
def lik(mt,mn,mal):
    vt=p*mt/(2*(1-p)+p*C); vn=f*mn/(2*(1-f)+f*C)
    l=binom.logpmf(ac,Dc,np.minimum(vt+et,1))+binom.logpmf(an,Dn,np.minimum(vn+en,1))
    if mal: l=np.where(mt<=C+1e-9,l,-np.inf)
    return l
L0=binom.logpmf(ac,Dc,et)+binom.logpmf(an,Dn,en)
w=np.ones(len(comps))/len(comps); pi={k:3e-4 for k in classes}
for it in range(60):
    L=np.stack([lik(a,b,t=='mal') for t,a,b in comps],1)
    pv=np.array([pi[k] for k in clv])
    A=np.log(pv)[:,None]+np.log(w)[None,:]+L; B=np.log1p(-pv)+L0
    mx=np.maximum(A.max(1),B); R=np.exp(A-mx[:,None]); r0=np.exp(B-mx); R/=(R.sum(1)+r0)[:,None]
    post=R.sum(1)
    w=(R*W[:,None]).sum(0); w/=w.sum()
    for k in classes: m_=clv==k; pi[k]=(post[m_]*W[m_]).sum()/W[m_].sum()
    # M-step positions
    new=[]
    for j,(t,a,b) in enumerate(comps):
        r=R[:,j]*W; sel=r>1e-4
        if t=='mal' and a in (1,2,3):  new.append((t,a,b)); continue
        acs,Dcs,ans,Dns,Cs,rs=ac[sel],Dc[sel],an[sel],Dn[sel],C[sel],r[sel]
        def nl(x):
            mt=x[0]; mn=x[0] if t=='mal' else x[1]
            vt=p*mt/(2*(1-p)+p*Cs); vn=f*mn/(2*(1-f)+f*Cs)
            return -(rs*(binom.logpmf(acs,Dcs,np.clip(vt+et,1e-9,1))+binom.logpmf(ans,Dns,np.clip(vn+en,1e-9,1)))).sum()
        x0=[a] if t=='mal' else [a,b]
        res=minimize(nl,x0,method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-3))
        x=res.x; new.append((t,x[0],x[0] if t=='mal' else x[1]))
    comps=new
print('p f',p,f)
tot=(post*W).sum()
for (t,a,b),ww in zip(comps,w): print(t,round(a,3),round(b,3),round(ww*tot,1))
print('pi',pi)
pm=R[:,[i for i,c in enumerate(comps) if c[0]=='mal']].sum(1)
pn=R[:,[i for i,c in enumerate(comps) if c[0]=='non']].sum(1)
df['pmal']=pm; df['pnon']=pn; df['psom']=post
ll=(np.log(np.exp(A-mx[:,None]).sum(1)+np.exp(B-mx))+mx)
print('LL',(ll*W).sum())
d=df[df.chrom!='pseudo']; d.to_pickle(f'/root/work/post3_{p}_{f}.pkl')
print('expected mal',(pm*W).sum(),'non',(pn*W).sum())
for th in [0.3,0.4,0.5,0.6,0.7,0.9]:
    s=d[d.pmal>th]; print(th,len(s),'expTP',round(s.pmal.sum(),1),'prec',round(s.pmal.mean(),3),'recall',round(s.pmal.sum()/(pm*W).sum(),3))
