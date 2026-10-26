import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
c=pd.read_pickle('/root/work/cf2.pkl')
c['seg']=np.where(c.chrom!='chr3',c.chrom,np.where(c.pos<=120000,'chr3a','chr3b'))
ct=c['sub'].isin(['C>T','G>A'])
c['anti']=np.where(c['sub']=='C>T',c.tumor_afr-c.a_f1r2,np.where(c['sub']=='G>A',c.a_f1r2,-1))
a=c.anti; n=c.normal_afr; t=c.tumor_afr
callct = ct & (((a>=3)&(n>=1)) | ((a>=2)&(n>=2)) | (a>=4) | ((a>=1)&(n>=3)) | (n>=4))
callnon = ~ct & ((t>=3) | ((t>=2)&(n>=2)))
c['call']=callct|callnon
bord=c[c.call&(c.nv>0.2)&(c.nv<0.5)]
print(bord[['chrom','pos','seg','sub','tumor_afr','tumor_nfr','tv','normal_afr','normal_nfr','nv','anti']].sort_values('nv').round(3).to_string())
c.to_pickle('/root/work/cf3.pkl')
