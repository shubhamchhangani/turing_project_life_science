import numpy as np, sys
from scipy.stats import binom, betabinom
p,f=float(sys.argv[1]),float(sys.argv[2]); usedmg=int(sys.argv[3])
r=np.load('/root/work/ref.npy'); n=len(r); idx=np.arange(n)
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4].astype(np.int64); NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].astype(np.int64)
Tc=C.sum(0); Nc=NC.sum(0); T=Tc.sum(1); N=Nc.sum(1)
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(n)
for s0,e0,c in segs: cn[s0:e0]=c
DT=p*cn+2*(1-p); DN=f*cn+2*(1-f)
d=np.load('/root/work/calls2.npz')
germ=np.zeros(n,bool); from scipy.stats import binom as _b
gm=(d['gl']<5)&(_b.cdf(d['an'],d['N'],0.38)>1e-3)
germ[d['cand'][gm]]=True
print('germ',gm.sum())
valid=(r<4)&~germ&(T>20)
xs0=np.arange(0.1,3.01,0.05); gs=[0.5,1,1.5,2,3,4.5]
xs=[(x,g) for g in gs for x in xs0]
# damage BB params
s=np.where(r==1)[0]; mu=C[0][s,3].sum()/C[0][s].sum()
fr=C[0][s,3]/np.maximum(C[0][s].sum(1),1); nn_=C[0][s].sum(1).mean()
vb=fr.var(); rho=max((vb*nn_/(mu*(1-mu))-1)/(nn_-1),1e-3)
ab=1/rho-1; print('damage mu',mu,'rho',rho)
units=[]
for alt in range(4):
    v=np.where(valid&(r!=alt))[0]
    ct=((r[v]==1)&(alt==3))|((r[v]==2)&(alt==0))
    safe=np.where(r[v]==1,1,0)
    sa=np.where(ct,C[1][v,alt]*(r[v]==1)+C[0][v,alt]*(r[v]==2),Tc[v,alt])
    sd=np.where(ct,np.where(r[v]==1,C[1][v].sum(1),C[0][v].sum(1)),T[v])
    da=np.where(ct,np.where(r[v]==1,C[0][v,alt],C[1][v,alt]),0)
    dd=np.where(ct,np.where(r[v]==1,C[0][v].sum(1),C[1][v].sum(1)),0)
    an=Nc[v,alt]
    keep=(sa+an)>=2
    rng=np.random.default_rng(alt); sub=(~keep)&(rng.random(len(v))<0.01)
    wt=np.where(keep,1.0,100.0)
    k=keep|sub
    units.append((v[k],np.full(k.sum(),alt),ct[k],sa[k],sd[k],da[k],dd[k],an[k],wt[k]))
U=[np.concatenate([u[i] for u in units]) for i in range(9)]
v,al,ct,sa,sd,da,dd,an,wt=U
et=np.where(ct,2.2e-4,1.0e-4); en=np.where(ct,3.8e-4,1.75e-4)
l0=binom.logpmf(sa,sd,et)+binom.logpmf(an,N[v],en)
if usedmg: l0+=np.where(ct,betabinom.logpmf(da,dd,mu*ab,(1-mu)*ab),0)
L=np.empty((len(v),len(xs)),np.float32)
for i,(x,g) in enumerate(xs):
    qt=np.minimum(x*p/DT[v],.99); qn=np.minimum(g*x*f/DN[v],.99)
    l=binom.logpmf(sa,sd,qt+et)+binom.logpmf(an,N[v],qn+en)
    if usedmg:
        m2=np.minimum(mu+qt,.99)
        l+=np.where(ct,betabinom.logpmf(da,dd,m2*ab,(1-m2)*ab),0)
    L[:,i]=np.where(x<=cn[v]+0.01,l-l0,-np.inf)
post=np.zeros(len(v)); W={}
for lab,mask in [('nonCT',~ct),('CT',ct)]:
    Lm=L[mask]; w_=wt[mask]; M=w_.sum(); w=np.full(len(xs),1e-5)
    for it in range(200):
        A=np.log(np.maximum(w,1e-300))[None,:]+Lm
        mx=np.maximum(A.max(1),0); num=np.exp(A-mx[:,None]); den=num.sum(1)+np.exp(-mx)*(1-w.sum())
        pp=num/den[:,None]; w=(pp*w_[:,None]).sum(0)/M
    post[mask]=pp.sum(1)
    print(lab,'est total',round(w.sum()*M,1)); print('  ',' '.join(f'{x[0]:.2f}/{x[1]}:{w[i]*M:.1f}' for i,x in enumerate(xs) if w[i]*M>2))
# one alt per site: keep best
res=[]
for th in [0.3,0.5,0.7,0.9]:
    c=post>th
    print('thr',th,'calls',c.sum(),'expFP',round((1-post[c]).sum(),1),'expFN',round((post[~c]*wt[~c]).sum(),1),'exp total',round((post*wt).sum(),1))
np.savez(f'/root/work/postg2_{usedmg}.npz',v=v,al=al,post=post,ct=ct,sa=sa,sd=sd,da=da,an=an,wt=wt)
