import pandas as pd,numpy as np
from scipy.stats import binom
df=pd.read_csv('/root/work/hets.tsv',sep='\t')
df['at']=(df.ft*df.dt).round().astype(int); df['an']=(df.fn*df.dn).round().astype(int)
def seg(r):
    if r.c=='chr3' and r.p<120000: return (2,0)
    return {'chr1':(2,2),'chr2':(2,1),'chr3':(2,2),'chr4':(3,2)}[r.c]
df['cn']=df.apply(seg,axis=1)
print(df.groupby('c').ft.mean())
def ll(p,a,d,cn):
    A,b=cn; C=A+b
    tot=0
    den=p*C+2*(1-p)
    f1=(p*A+1-p)/den; f2=(p*b+1-p)/den
    return np.logaddexp(binom.logpmf(a,d,f1),binom.logpmf(a,d,f2))+np.log(.5)
for col,dcol,name in [('at','dt','tumor'),('an','dn','normal')]:
    best=[]
    for segname,g in df.groupby(df.cn.astype(str)):
        if segname=='(2, 2)': continue
        ps=np.linspace(0.01,0.99,197)
        L=[ll(p,g[col].values,g[dcol].values,g.cn.iloc[0]).sum() for p in ps]
        i=int(np.argmax(L)); print(name,segname,len(g),round(ps[i],3))
    ps=np.linspace(0.01,0.99,197)
    g=df[df.cn!=(2,2)]
    L=np.array([sum(ll(p,r[col],r[dcol],r.cn) for _,r in g.iterrows()) for p in ps])
    i=int(np.argmax(L)); ci=ps[L>L[i]-1.92]; print(name,'all',round(ps[i],3),ci.min(),ci.max())
