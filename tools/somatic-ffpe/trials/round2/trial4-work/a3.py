import numpy as np
T=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4].sum(1).astype(float)
Nn=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4].sum(1).astype(float)
offs=[0,300000,550000,800000,1000000]
def best_split(x):
    n=len(x); cs=np.cumsum(x); tot=cs[-1]
    k=np.arange(200,n-200)
    m1=cs[k-1]/k; m2=(tot-cs[k-1])/(n-k)
    stat=k*(n-k)/n*(m1-m2)**2
    i=stat.argmax(); return k[i], stat[i], m1[i], m2[i]
for ci in range(4):
    a,b=offs[ci],offs[ci+1]; x=T[a:b]
    print('chr',ci+1, best_split(x), x.mean())
# fine around 670k
a=660000;b=680000
x=T[a:b]; print(best_split(x))
print('---')
segs=[(0,300000),(300000,550000),(550000,670143),(670143,800000),(800000,1000000)]
def rec(a,b,depth=0):
    x=T[a:b]
    if b-a<2000: return
    k,st,m1,m2=best_split(x)
    # z-like: stat / var
    v=x.var()
    print('  '*depth,a,b,a+k,round(st/v,1),round(m1,1),round(m2,1))
    if st/v>30 and depth<3:
        rec(a,a+k,depth+1); rec(a+k,b,depth+1)
for a,b in segs: rec(a+300,b-300)
