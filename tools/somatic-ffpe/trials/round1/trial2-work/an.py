import pandas as pd, numpy as np, pickle
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
c,out=pickle.load(open('/root/work/feat.pkl','rb'))
c=c.reset_index(drop=True)
feats=[]
for i,reads in enumerate(out['tumor']):
    R=pd.DataFrame(reads,columns=['q','alt','bq','qp','f1r2','dup','fs','fe','mq'])
    nd=R[~R.dup]
    # per fragment
    fr=nd.groupby('q').agg(alt=('alt','max'),altmin=('alt','min'),f1r2=('f1r2','first'),bq=('bq','max'),qp=('qp','min'),fs=('fs','first'),fe=('fe','first'))
    a=fr[fr.alt]
    pos=c.pos[i]-1
    # family analysis incl dups
    fam=R.groupby(['fs','fe','f1r2','q']).alt.max().reset_index()
    famg=fam.groupby(['fs','fe'])
    nfam_alt=0; nfam_both=0; nfam_mixed=0
    for k,g in famg:
        if g.alt.any():
            nfam_alt+=1
            if g[g.alt].f1r2.nunique()==2: nfam_both+=1
            if (~g.alt).any(): nfam_mixed+=1
    dend=np.minimum(pos-a.fs, a.fe-1-pos) if len(a) else pd.Series([],dtype=float)
    feats.append(dict(nfr=len(fr),afr=len(a),a_f1r2=int(a.f1r2.sum()),disc=int((fr.alt!=fr.altmin).sum()),
        abq=a.bq.mean() if len(a) else np.nan, aqp=a.qp.mean() if len(a) else np.nan,
        dend=dend.mean() if len(a) else np.nan, nfam_alt=nfam_alt,nfam_both=nfam_both,nfam_mixed=nfam_mixed,
        altdup=int(R[R.dup&R.alt].shape[0]), nfs=len(a)))
F=pd.DataFrame(feats)
nf=[]
for reads in out['normal']:
    R=pd.DataFrame(reads,columns=['q','alt','bq','qp','f1r2','dup','fs','fe','mq'])
    nd=R[~R.dup]; fr=nd.groupby('q').alt.max()
    nf.append((len(fr),int(fr.sum())))
F['nfrN'],F['afrN']=zip(*nf)
c=pd.concat([c,F],axis=1)
c.to_pickle('/root/work/cf.pkl')
c['tv']=c.talt/c.tdp; c['nv']=c.nalt/c.ndp.clip(lower=1)
ct=c['sub'].isin(['C>T','G>A'])
true_like=c[~ct&(c.nv<0.3)&(c.afrN>=1)]
art=c[ct&(c.afrN==0)&(c.tv<0.04)]
for name,d in [('true',true_like),('art',art),('ct_high',c[ct&(c.afrN==0)&(c.tv>0.06)&(c.tv<0.3)]),('ct_withN',c[ct&(c.afrN>=2)&(c.nv<0.3)])]:
    d=d.copy(); d['ofrac']=d.a_f1r2/d.afr
    print(name,len(d)); print(d[['afr','ofrac','abq','aqp','dend','nfam_alt','nfam_both','nfam_mixed','altdup','disc']].describe().loc[['mean','25%','50%','75%']])
    print(d.groupby('sub').ofrac.mean())
