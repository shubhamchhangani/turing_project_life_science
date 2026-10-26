import numpy as np
d=np.load('/root/work/calls.npz')
sa,sd,an,N,ct,cn,germ=d['sa'],d['sd'],d['an'],d['N'],d['isct'],d['cn'],d['germ']
m=(~germ)&(sa+an>=4)
vt=sa/sd
for c in [2,3,4,5]:
  for lab,mm in [('nonCT',~ct),('CT',ct)]:
    s=m&(cn==c)&mm
    print(c,lab,s.sum(),np.histogram(vt[s],bins=np.arange(0,0.8,0.025))[0])
