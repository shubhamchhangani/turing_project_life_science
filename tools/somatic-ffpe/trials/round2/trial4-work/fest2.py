import numpy as np
d=np.load("post_1.npz"); N=np.load("normal.cnt.npy")[:,:,:4].sum(0).sum(1)
p,f=0.615,0.2
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
cn=np.zeros(1000000)
for s0,e0,c in segs: cn[s0:e0]=c
m=(d["post"]>0.5)
v=d["v"][m]; c=cn[v]
xt=(d["sa"][m]/d["sd"][m])*(p*c+2*(1-p))/p
xn=(d["an"][m]/N[v])*(f*c+2*(1-f))/f
H,_,_=np.histogram2d(xt,xn,bins=[np.arange(0,3.01,0.2),np.arange(0,3.01,0.2)])
print("rows xt 0..3 step .2, cols xn")
for i,row in enumerate(H): print(f"{i*0.2:.1f}"," ".join(f"{int(z):3d}" for z in row))
