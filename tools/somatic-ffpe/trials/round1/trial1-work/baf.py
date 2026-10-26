import pickle,numpy as np
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
T=pickle.load(open('/root/work/tumor.scan.pkl','rb'))
import sys
B=int(sys.argv[1]) if len(sys.argv)>1 else 20000
rows=[]
for c in T:
    dt,at=T[c]; dn,an=N[c]
    for (p,a),v in an.items():
        n_alt=v[0]+v[1]
        if dn[p]<30: continue
        f=n_alt/dn[p]
        if 0.25<f<0.75:
            t=at.get((p,a)); ta=(t[0]+t[1]) if t is not None else 0
            rows.append((c,p,f,ta/dt[p],dt[p],dn[p]))
import pandas as pd
df=pd.DataFrame(rows,columns=['c','p','fn','ft','dt','dn'])
df.to_csv('/root/work/hets.tsv',sep='\t',index=False)
print(len(df))
df['bin']=df.p//B
for (c,b),g in df.groupby(['c','bin']):
    mt=np.abs(g.ft-0.5); mn=np.abs(g.fn-0.5)
    print(c,b*B,len(g),'T |baf-.5| med %.3f  N %.3f   T hist'%(mt.median(),mn.median()), np.histogram(g.ft,bins=10,range=(0,1))[0])
