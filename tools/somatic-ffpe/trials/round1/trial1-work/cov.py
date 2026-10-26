import pickle,numpy as np
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan.pkl','rb'))
B=10000
mt=np.mean(np.concatenate([T[c][0] for c in T])); mn=np.mean(np.concatenate([N[c][0] for c in N]))
print('mean depth T',mt,'N',mn)
for c in T:
    dt=T[c][0]; dn=N[c][0]
    for i in range(0,len(dt),B):
        print(c,i,round(dt[i:i+B].mean()/mt,3),round(dn[i:i+B].mean()/mn,3))
