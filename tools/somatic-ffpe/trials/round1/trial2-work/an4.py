import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',40); pd.set_option('display.max_rows',400)
c=pd.read_pickle('/root/work/cf2.pkl')
p,f=0.618,0.205
c['seg']=np.where(c.chrom!='chr3',c.chrom,np.where(c.pos<=120000,'chr3a','chr3b'))
CN={'chr1':4,'chr2':3,'chr3a':2,'chr3b':4,'chr4':5}
c['CN']=c.seg.map(CN)
ct=c['sub'].isin(['C>T','G>A'])
d=c[~ct&(c.nv<0.35)].copy()
d['tccf']=d.tv*(p*d.CN+2-2*p)/p   # = m*ccf
d['nccf']=d.nv*(f*d.CN+2-2*f)/f
print(d[['seg','sub','tumor_afr','tumor_nfr','normal_afr','normal_nfr','tv','nv','tccf','nccf']].sort_values(['seg','tccf']).round(2).to_string())
