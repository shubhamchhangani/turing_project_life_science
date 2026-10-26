import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',40); pd.set_option('display.max_rows',400)
c=pd.read_pickle('/root/work/cf2.pkl')
p,f=0.618,0.205
c['seg']=np.where(c.chrom!='chr3',c.chrom,np.where(c.pos<=120000,'chr3a','chr3b'))
CN={'chr1':4,'chr2':3,'chr3a':2,'chr3b':4,'chr4':5}
c['CN']=c.seg.map(CN)
c['R']=(f/p)*(p*c.CN+2-2*p)/(f*c.CN+2-2*f)
c['ratio']=c.nv/c.tv
d=c[(c.nv>0.15)&(c.nv<0.5)]
print(d[['chrom','pos','sub','tumor_afr','tumor_nfr','tv','normal_afr','normal_nfr','nv','R','ratio','ofrac']].sort_values('nv').round(3).to_string())
