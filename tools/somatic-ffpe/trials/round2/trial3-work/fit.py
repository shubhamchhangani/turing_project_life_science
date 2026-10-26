import numpy as np
from scipy.optimize import least_squares
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
dt3=T['chr3'][:,:,:4].sum((1,2)).astype(float)
# breakpoint
cs=np.cumsum(dt3); L=len(dt3); best=None
for b in range(100000,140000,100):
    m1=cs[b]/b; m2=(cs[-1]-cs[b])/(L-b)
    sse=((dt3[:b]-m1)**2).sum()+((dt3[b:]-m2)**2).sum()
    if best is None or sse<best[0]: best=(sse,b)
bp=best[1]
for b in range(bp-300,bp+300,10):
    m1=cs[b]/b; m2=(cs[-1]-cs[b])/(L-b)
    sse=((dt3[:b]-m1)**2).sum()+((dt3[b:]-m2)**2).sum()
    if sse<best[0]: best=(sse,b)
bp=best[1]; print('bp',bp)
segs=[('chr1',0,300000,2,2),('chr2',0,250000,2,1),('chr3',0,bp,2,0),('chr3',bp,250000,2,2),('chr4',0,200000,3,2)]
mt=[];mn=[];C=[];bafs=[]
for c,s,e,a,b in segs:
    t=T[c][s+200:e-200,:,:4].sum(1); n=N[c][s+200:e-200,:,:4].sum(1); ref=T[c+'_ref'][s+200:e-200]
    dt=t.sum(1); dn=n.sum(1)
    mt.append(np.median(dt)); mn.append(np.median(dn)); C.append(a+b)
    idx=np.arange(len(dt)); nn=n.copy(); nn[idx,ref.clip(0,3)]=-1; al=nn.argmax(1)
    nv=n[idx,al]/np.maximum(dn,1)
    het=(dn>=40)&(nv>0.25)&(nv<0.75)
    # tumor allele fraction of the alt; mirror is ambiguous; use mean of max(baf,1-baf) and also mean of baf^... use 2-component
    tb=t[idx,al][het]/dt[het]
    bafs.append(tb)
    print(c,s,e,'T',np.mean(dt),np.median(dt),'N',np.mean(dn),np.median(dn),'nhet',het.sum(),'meanMajor',np.mean(np.maximum(tb,1-tb)))
mt=np.array(mt);mn=np.array(mn);C=np.array(C)
for name,m in (('tumor',mt),('normal',mn)):
    r=least_squares(lambda x: x[0]*(2*(1-x[1])+x[1]*C)-m,[40,0.5])
    print(name,r.x,'resid',r.fun)
# BAF ML for purity: mixture with symmetric
from scipy.stats import binom
from scipy.special import logsumexp
cnt=[]
for (c,s,e,a,b) in segs:
    t=T[c][s+200:e-200,:,:4].sum(1); n=N[c][s+200:e-200,:,:4].sum(1); ref=T[c+'_ref'][s+200:e-200]
    dt=t.sum(1); dn=n.sum(1); idx=np.arange(len(dt)); nn=n.copy(); nn[idx,ref.clip(0,3)]=-1; al=nn.argmax(1)
    nv=n[idx,al]/np.maximum(dn,1); het=(dn>=40)&(nv>0.2)&(nv<0.8)
    cnt.append((t[idx,al][het],dt[het],n[idx,al][het],dn[het],a,b))
def nll(p,which=0):
    tot=0
    for k,d,kn,dn,a,b in cnt:
        x=p
        bb=(1-x+x*a)/(2*(1-x)+x*(a+b))
        if which==0: kk,dd=k,d
        else: kk,dd=kn,dn
        l1=binom.logpmf(kk,dd,bb); l2=binom.logpmf(kk,dd,1-bb)
        tot+=-np.logaddexp(l1,l2).sum()
    return tot
for w in (0,1):
    g=np.arange(0.02,0.95,0.005); v=[nll(x,w) for x in g]; print('BAF ML',w,g[np.argmin(v)])
