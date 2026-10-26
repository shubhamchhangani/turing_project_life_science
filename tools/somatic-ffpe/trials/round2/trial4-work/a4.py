import numpy as np
T=np.load('/root/work/tumor.cnt.npy').sum(0)[:,:4].sum(1).astype(float)
Nn=np.load('/root/work/normal.cnt.npy').sum(0)[:,:4].sum(1).astype(float)
for a,b in [(294000,300000),(300000,303000),(793000,800000),(976000,988000),(0,3000),(547000,553000),(800000,803000),(995000,1000000)]:
    print(a,b)
    for s in range(a,b,500): print('  ',s,round(T[s:s+500].mean(),1),round(Nn[s:s+500].mean(),1))
