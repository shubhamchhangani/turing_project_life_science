import pickle,numpy as np
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan.pkl','rb'))
mt=np.mean(np.concatenate([T[c][0] for c in T])); mn=np.mean(np.concatenate([N[c][0] for c in N]))
B=2000
lv={'t':[0.655,0.857,1.06,1.26]}
for c in T:
    dt=T[c][0]/mt; dn=N[c][0]/mn
    s=''
    for i in range(0,len(dt),B):
        x=dt[i:i+B].mean()
        s+=str(int(np.argmin([abs(x-l) for l in lv['t']])))
    print(c,s)
d=T['chr3'][0]/mt
for i in range(110000,130000,500): print(i,round(d[i:i+500].mean(),3))
