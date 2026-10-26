import numpy as np, sys
from scipy.stats import binom
r=np.load('/root/work/ref.npy'); n=len(r); idx=np.arange(n)
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4].astype(np.int64); NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].astype(np.int64)
Tc=C.sum(0); Nc=NC.sum(0); T=Tc.sum(1); N=Nc.sum(1)
SA=Tc.copy(); SD=np.repeat(T[:,None],4,1)
m=(r==1); SA[m,3]=C[1][m,3]; SD[m,3]=C[1][m].sum(1)
m=(r==2); SA[m,0]=C[0][m,0]; SD[m,0]=C[0][m].sum(1)
p,f=float(sys.argv[1]),float(sys.argv[2])
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(n)
for s0,e0,c in segs: cn[s0:e0]=c
DT=p*cn+2*(1-p); DN=f*cn+2*(1-f)
d=np.load('/root/work/calls2.npz')
germ=np.zeros(n,bool); germ[d['cand'][d['gl']<5]]=True
valid=(r<4)&~germ&(T>20)
# per site, sum over 3 alt alleles: treat each (site,alt) as unit
units=[]
xs=np.arange(0.05,3.01,0.05)
rates={True:(2.2e-4/1.0,3.8e-4),False:(1.0e-4,1.75e-4)}
LLn_tot=0
Lmat=[];L0=[];info=[]
for alt in range(4):
    v=np.where(valid&(r!=alt))[0]
    ct=((r[v]==1)&(alt==3))|((r[v]==2)&(alt==0))
    s=SA[v,alt]; an=Nc[v,alt]
    rng=np.random.default_rng(alt)
    keep=((s+an)>0)|(rng.random(len(v))<0.05)
    wt=np.where((s+an)>0,1.0,20.0)[keep]
    v=v[keep]; ct=ct[keep]; s=s[keep]; an=an[keep]; sd=SD[v,alt]; nn=N[v]
    et=np.where(ct,7e-5,1e-4); en=np.where(ct,1.3e-4,1.75e-4)
    l0=binom.logpmf(s,sd,et)+binom.logpmf(an,nn,en)
    Lx=np.array([np.where(x<=cn[v],binom.logpmf(s,sd,np.minimum(x*p/DT[v]+et,.999))+binom.logpmf(an,nn,np.minimum(x*f/DN[v]+en,.999)),-np.inf) for x in xs]).T
    Lmat.append((Lx-l0[:,None]).astype(np.float32)); info.append((v,alt,ct,wt))
Lr=np.vstack(Lmat)  # log LR vs null
ctall=np.concatenate([i[2] for i in info]); WT=np.concatenate([i[3] for i in info])
print(Lr.shape)
for lab,mask in [('all',np.ones(len(Lr),bool)),('nonCT',~ctall),('CT',ctall)]:
    L=Lr[mask]; wt=WT[mask]; M=wt.sum()
    w=np.full(len(xs),1e-4)
    for it in range(300):
        A=np.log(np.maximum(w,1e-300))[None,:]+L
        mx=np.maximum(A.max(1),0); num=np.exp(A-mx[:,None]); den=num.sum(1)+np.exp(-mx)*(1-w.sum())
        post=num/den[:,None]
        w=(post*wt[:,None]).sum(0)/M
    print(lab,'est total somatic units',round(w.sum()*M,1))
    print('  x dist:',' '.join(f'{x:.2f}:{w[i]*M:.1f}' for i,x in enumerate(xs) if w[i]*M>0.5))
    np.save(f'/root/work/w_{lab}.npy',w*M)
