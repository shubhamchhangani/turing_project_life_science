import pandas as pd,numpy as np,pickle
N=pickle.load(open('/root/work/normal.scan.pkl','rb'))
df=pd.read_csv('/root/work/cand.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['ct']=((df.ref=='C')&(df.alt=='T'))|((df.ref=='G')&(df.alt=='A'))
def no(r):
    v=N[r.chrom][1].get((r.pos-1,r.alt)); return (v[0],v[1]) if v is not None else (0,0)
s=df[(df.nalt>=3)&(df.vn<0.3)].copy()
s[['nF1','nF2']]=s.apply(lambda r: pd.Series(no(r)),axis=1)
s['grp']=np.where(s.vn>s.vt*0.7,'blob','normalratio')
s['mt']=s.ref+'>'+s.alt
for g,x in s.groupby('grp'):
    print(g,len(x),x.mt.value_counts().to_dict())
    print(' tumor ori',x.tF1.sum(),x.tF2.sum(),' normal ori',x.nF1.sum(),x.nF2.sum(), 'tlow',x.tlow.sum(),'nlow',x.nlow.sum(),'talt',x.talt.sum(),'nalt',x.nalt.sum())
    print(x.chrom.value_counts().to_dict())
s.to_csv('/root/work/s.tsv',sep='\t',index=False)
