import pandas as pd, numpy as np
rows=[]
for line in open('/root/work/all.txt'):
    c,p,r,a,t,n,_=line.split('\t')
    alts=a.split(',')
    tv=list(map(int,t.split(','))); nv=list(map(int,n.split(',')))
    # pick best alt by tumor+normal
    if len(alts)>1 or alts[0]!='<*>':
        best=None;bs=-1
        for i,al in enumerate(alts):
            if al=='<*>': continue
            s=tv[i+1]+nv[i+1]
            if s>bs: bs=s;best=i
        rows.append((c,int(p),r,alts[best],tv[0],tv[best+1],sum(tv),nv[0],nv[best+1],sum(nv)))
    else:
        rows.append((c,int(p),r,'.',tv[0],0,sum(tv),nv[0],0,sum(nv)))
df=pd.DataFrame(rows,columns=['chrom','pos','ref','alt','tref','talt','tdp','nref','nalt','ndp'])
df.to_pickle('/root/work/all.pkl')
print(df.describe())
