import numpy as np
d=np.load('/root/work/cand.npz'); r=np.load('/root/work/ref.npy')
sa,sd,an,N,ct=d['sa'],d['sd'],d['an'],d['N'],d['ct']
cand=((sa>=3)|((sa>=2)&(an>=2))|(an>=4))&(r<4)
vt=sa/np.maximum(sd,1); vn=an/np.maximum(N,1)
s=np.where(cand)[0]; print(len(s))
H,_,_=np.histogram2d(vt[s],vn[s],bins=[np.linspace(0,1,21),np.linspace(0,1,21)])
for i in range(20): print(f'{i*0.05:.2f}',' '.join(f'{int(x):4d}' for x in H[i]))
