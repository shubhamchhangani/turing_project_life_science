import pandas as pd, numpy as np
from scipy.stats import binom
pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
c=pd.read_pickle('/root/work/cf3.pkl')
p,f=0.618,0.205
AB={'chr1':(2,2),'chr2':(2,1),'chr3a':(2,0),'chr3b':(2,2),'chr4':(3,2)}
def lls(r):
    A,B=AB[r.seg]; CN=A+B
    Nt=p*CN+2-2*p; Nn=f*CN+2-2*f
    tg=[]; 
    for k in set([A,B,CN]):  # het on A, het on B, hom (k=CN)
        if k==CN: vt=(p*CN+2-2*p)/Nt; vn=1.0
        else: vt=(p*k+1-p)/Nt; vn=(f*k+1-f)/Nn
        vt=min(vt,0.995); vn=min(vn,0.995)
        tg.append(binom.logpmf(r.tumor_afr,r.tumor_nfr,vt)+binom.logpmf(r.normal_afr,r.normal_nfr,vn))
    xs=np.linspace(0.02,CN,200)
    ls=binom.logpmf(r.tumor_afr,r.tumor_nfr,np.minimum(p*xs/Nt,.995))+binom.logpmf(r.normal_afr,r.normal_nfr,np.minimum(f*xs/Nn,.995))
    return max(tg), ls.max()
m=c.call
res=c[m].apply(lls,axis=1,result_type='expand'); res.columns=['llg','lls']
c.loc[m,'llg']=res.llg; c.loc[m,'lls']=res.lls
c['dll']=c.lls-c.llg
x=c[m]
print((x.dll<0).sum(), 'germline-like of', m.sum())
print(x[(x.dll.abs()<15)][['chrom','pos','seg','sub','tumor_afr','tumor_nfr','tv','normal_afr','normal_nfr','nv','anti','dll']].sort_values('dll').round(3).to_string())
c.to_pickle('/root/work/cf4.pkl')
