import numpy as np
from scipy.optimize import brentq
d=np.load('/root/work/post_1.npz'); N=np.load('/root/work/normal.cnt.npy')[:,:,:4].sum(0).sum(1)
p=0.615
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(1000000)
for s0,e0,c in segs: cn[s0:e0]=c
m=(d['post']>0.9)&(~d['ct'])
v=d['v'][m]; vt=d['sa'][m]/d['sd'][m]; an=d['an'][m]
x=vt*(p*cn[v]+2*(1-p))/p
for lab,mm in [('all',x>0),('clonal',x>0.7),('sub',x<0.7)]:
    g=lambda f: (x[mm]*f/(f*cn[v][mm]+2*(1-f))*N[v][mm]).sum()-an[mm].sum()
    print(lab,mm.sum(),brentq(g,0.01,0.9))
print("tumour-selected")
d2=d; m=(~d2["ct"])&(d2["sa"]>=7)&(d2["post"]>0.5)
v=d2["v"][m]; vt=d2["sa"][m]/d2["sd"][m]; an=d2["an"][m]
x=vt*(p*cn[v]+2*(1-p))/p
for lo,hi in [(0,0.4),(0.4,0.7),(0.7,1.5),(1.5,4)]:
    mm=(x>=lo)&(x<hi)
    g=lambda f: (x[mm]*f/(f*cn[v][mm]+2*(1-f))*N[v][mm]).sum()-an[mm].sum()
    print(lo,hi,mm.sum(),brentq(g,0.001,0.95), "mean vn/vt", (an[mm]/N[v][mm]).mean()/vt[mm].mean())
