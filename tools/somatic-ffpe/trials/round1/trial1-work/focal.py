import pickle,numpy as np
T=pickle.load(open('/root/work/tumor.scan2.pkl','rb'))
N=pickle.load(open('/root/work/normal.scan2.pkl','rb'))
p=0.612
def segC(c,pos):
    if c=='chr3' and pos<=120000: return 2
    return {'chr1':4,'chr2':3,'chr3':4,'chr4':5}[c]
k=39.92
B=1000
tot=0;n=0
for c in T:
    d=T[c][0]
    for i in range(0,len(d),B):
        C=segC(c,i+1); e=k*(p*C+2*(1-p))
        m=d[i:i+B].mean()
        z=(m-e)/e
        if abs(z)>0.12: print(c,i,round(m,1),round(e,1),round(z,3))
