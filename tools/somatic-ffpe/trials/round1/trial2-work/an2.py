import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
c=pd.read_pickle('/root/work/cf.pkl')
c['tv']=c.tumor_afr/c.tumor_nfr; c['nv']=c.normal_afr/c.normal_nfr.clip(lower=1)
c['ofrac']=c.a_f1r2/c.tumor_afr
ct=c['sub'].isin(['C>T','G>A'])
groups=[('true',~ct&(c.nv<0.3)&(c.normal_afr>=1)),('noct_noN',~ct&(c.normal_afr==0)),('art',ct&(c.normal_afr==0)&(c.tv<0.04)),('ct_high_noN',ct&(c.normal_afr==0)&(c.tv>0.06)&(c.tv<0.3)),('ct_withN',ct&(c.normal_afr>=2)&(c.nv<0.3))]
cols=['tumor_afr','tv','ofrac','abq','aqp','dend','nfam_alt','nfam_both','nfam_mixed','altreads_all','disc']
for name,m in groups:
    d=c[m]
    print(name,len(d)); print(d[cols].describe().loc[['mean','25%','50%','75%']].round(3))
    print(d.groupby('sub').ofrac.agg(['mean','count']).T)
c.to_pickle('/root/work/cf2.pkl')
