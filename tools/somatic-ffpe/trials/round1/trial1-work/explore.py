import pandas as pd,numpy as np
df=pd.read_csv('/root/work/cand.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['mt']=df.ref+'>'+df.alt
germ=df[df.vn>0.3]
print('germline-like',len(germ))
s=df[(df.vn<0.2)]
print(pd.crosstab(s.talt.clip(0,15),s.nalt.clip(0,8)))
for k in [2,3,4,6,10]:
    x=s[(s.talt>=k)&(s.nalt==0)]
    print(k,'nalt0',x.mt.value_counts().to_dict())
    x=s[(s.talt>=k)&(s.nalt>=2)]
    print(k,'nalt>=2',x.mt.value_counts().to_dict())
