import pandas as pd,numpy as np
df=pd.read_csv('/root/work/cand.tsv',sep='\t')
df['vt']=df.talt/df.dt; df['vn']=df.nalt/df.dn
df['ct']=((df.ref=='C')&(df.alt=='T'))|((df.ref=='G')&(df.alt=='A'))
s=df[(~df.ct)&(df.nalt>=2)&(df.vn<0.3)]
print(len(s))
H,xe,ye=np.histogram2d(s.vt,s.vn,bins=[np.arange(0,0.62,0.04),np.arange(0,0.32,0.02)])
print('rows vt bins, cols vn bins (0.02)')
for i in range(len(xe)-1): print('%.2f'%xe[i],' '.join('%3d'%v for v in H[i]))
