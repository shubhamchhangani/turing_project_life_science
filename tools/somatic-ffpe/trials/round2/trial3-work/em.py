import numpy as np, pandas as pd
from scipy.stats import binom
p=0.607; f=0.205; et=1.0e-4; en=1.65e-4
df=pd.read_pickle('/root/work/all.pkl')
df=df[~((df.an/df.Dn.clip(lower=1)>0.33)&(df.an>=5))]
df=df[df.Dn>0]
cls={('C','T'):'C>T',('G','A'):'C>T',('C','A'):'C>A',('G','T'):'C>A',('C','G'):'C>G',('G','C'):'C>G',('T','A'):'T>A',('A','T'):'T>A',('T','C'):'T>C',('A','G'):'T>C',('T','G'):'T>G',('A','C'):'T>G'}
df['cl']=[cls[(r,a)] for r,a in zip(df.ref,df.alt)]
df['w']=1.0
# pseudo rows for zero-count sites
ref=open('/root/data/reference.fa').read()
import pysam
fa=pysam.FastaFile('/root/data/reference.fa')
pseudo=[]
segs=[('chr1',0,300000,4),('chr2',0,250000,3),('chr3',0,120150,2),('chr3',120150,250000,4),('chr4',0,200000,5)]
med=df.groupby(['C','cl']).agg(Dc=('Dc','median'),Dn=('Dn','median'))
for c,s,e,C in segs:
    seq=fa.fetch(c,s,e).upper()
    for r in 'ACGT':
        nr=seq.count(r)
        for a in 'ACGT':
            if a==r: continue
            k=cls[(r,a)]
            sub=df[(df.chrom==c)&(df.pos>s)&(df.pos<=e)&(df.ref==r)&(df.alt==a)]
            Dc=med.loc[(C,k),'Dc']; Dn=med.loc[(C,k),'Dn']
            if k=='C>T': Dc=med.loc[(C,k),'Dc']
            pseudo.append(dict(chrom='pseudo',pos=0,ref=r,alt=a,ac=0,Dc=Dc,ad=0,Dd=0,an=0,Dn=Dn,Dt=0,C=C,cl=k,w=nr-len(sub)))
df=pd.concat([df,pd.DataFrame(pseudo)],ignore_index=True)
grid=np.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,3.0])
ac=df.ac.values[:,None]; Dc=df.Dc.values[:,None]; an=df.an.values[:,None]; Dn=df.Dn.values[:,None]; C=df.C.values[:,None]
vt=p*grid[None,:]/(2*(1-p)+p*C); vn=f*grid[None,:]/(2*(1-f)+f*C)
valid=grid[None,:]<=C+1e-9
L1=binom.logpmf(ac,Dc,vt+et)+binom.logpmf(an,Dn,vn+en); L1[~valid]=-np.inf
L0=(binom.logpmf(ac,Dc,et)+binom.logpmf(an,Dn,en))[:,0]
w=np.ones(len(grid))/len(grid); classes=sorted(df.cl.unique()); pi={k:1e-4 for k in classes}
W=df.w.values; clv=df.cl.values
for it in range(300):
    pv=np.array([pi[k] for k in clv])
    a=np.log(pv)[:,None]+np.log(w)[None,:]+L1
    b=np.log1p(-pv)+L0
    mx=np.maximum(a.max(1),b)
    den=np.exp(b-mx)+np.exp(a-mx[:,None]).sum(1)
    R=np.exp(a-mx[:,None])/den[:,None]
    post=R.sum(1)
    w=(R*W[:,None]).sum(0); w/=w.sum()
    for k in classes:
        m=clv==k; nk={'C>T':0,'C>A':0}
        pi[k]=(post[m]*W[m]).sum()/W[m].sum()
print('grid weights',dict(zip(grid,np.round(w,3))))
print('pi',pi)
for k in classes:
    m=clv==k; print(k,'expected n somatic',(post[m]*W[m]).sum())
print('total expected',(post*W).sum())
df['post']=post
df.to_pickle('/root/work/post.pkl')
real=df[df.chrom!='pseudo']
for th in [0.3,0.5,0.7,0.9]:
    s=real[real.post>th]; print(th,len(s),'exp TP',s.post.sum(),'prec',s.post.mean(),'recall',s.post.sum()/(post*W).sum())
