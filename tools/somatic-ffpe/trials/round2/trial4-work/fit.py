import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom
r=np.load('/root/work/ref.npy')
Tc=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4]; Nc=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4]
T=Tc.sum(1).astype(float); N=Nc.sum(1).astype(float)
segs=[(0,300000,4,2),(300000,550000,3,1),(550000,670143,2,0),(670143,800000,4,2),(800000,1000000,5,2)]
het=np.load('/root/work/het.npy'); idx=np.arange(len(r))
c2=Nc.copy(); c2[idx[r<4],r[r<4]]=-1; aN=c2.argmax(1)
for name,C,D in [('tumor',Tc,T),('normal',Nc,N)]:
    dm=[];cn=[]
    hets=[]
    for a,b,c,nb in segs:
        dm.append(D[a+1000:b-1000].mean()); cn.append(c)
        h=het[(het>=a+200)&(het<b-200)]
        alt=C[h,aN[h]]; dep=D[h]
        hets.append((alt,dep,c,nb))
    dm=np.array(dm); cn=np.array(cn)
    # depth fit: dm = k*(2+p*(cn-2))
    A=np.vstack([np.ones(5)*2,cn-2]).T; w=np.array([b-a for a,b,_,_ in segs])
    sol,res,_,_=np.linalg.lstsq(A*np.sqrt(w)[:,None],dm*np.sqrt(w),rcond=None)
    k,kp=sol; p=kp/k
    print(name,'seg depths',np.round(dm,1),'k',k,'p_depth',p, 'pred',np.round(A@sol,1))
    def nll(p):
        L=0
        for alt,dep,c,nb in hets:
            tot=2*(1-p)+p*c
            f1=((1-p)+p*nb)/tot; f2=1-f1
            L-=np.log(0.5*binom.pmf(alt,dep.astype(int),f1)+0.5*binom.pmf(alt,dep.astype(int),f2)+1e-300).sum()
        return L
    ps=np.linspace(0.02,0.98,97); ll=[nll(x) for x in ps]
    i=np.argmin(ll); print(' BAF p',ps[i], 'per-seg:')
    for alt,dep,c,nb in hets:
        if nb*2==c: continue
        def n1(p):
            tot=2*(1-p)+p*c; f1=((1-p)+p*nb)/tot
            return -np.log(0.5*binom.pmf(alt,dep.astype(int),f1)+0.5*binom.pmf(alt,dep.astype(int),1-f1)).sum()
        l=[n1(x) for x in ps]; print('   cn',c,nb,'p',ps[np.argmin(l)],len(alt))
ploidy=sum((b-a)*c for a,b,c,_ in segs)/1e6; print('ploidy',ploidy)
