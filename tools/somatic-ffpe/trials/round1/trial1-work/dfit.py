import pickle,numpy as np
from scipy.optimize import least_squares
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan.pkl','rb'))
segs=[('chr1',0,300000,4),('chr2',0,250000,3),('chr3',0,120000,2),('chr3',120000,250000,4),('chr4',0,200000,5)]
for name,S in [('T',T),('N',N)]:
    m=[];C=[]
    for c,a,b,cn in segs:
        d=S[c][0][a+500:b-500]; m.append(np.median(d)); C.append(cn)
    m=np.array(m);C=np.array(C)
    f=lambda x: x[0]*(x[1]*C+2*(1-x[1]))-m
    r=least_squares(f,[30,0.5]); print(name,m,r.x, f(r.x))
