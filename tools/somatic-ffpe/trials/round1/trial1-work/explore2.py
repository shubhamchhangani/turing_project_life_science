import pandas as pd,numpy as np
df=pd.read_csv('/root/work/cand.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['ct']=df.ref.isin(['C','G'])&(((df.ref=='C')&(df.alt=='T'))|((df.ref=='G')&(df.alt=='A')))
s=df[(df.vn<0.2)&(df.talt>=3)]
nc=s[~s.ct]
print(nc.sort_values('talt')[['chrom','pos','ref','alt','tF1','tF2','tlow','dt','nalt','dn','vt']].to_string())
