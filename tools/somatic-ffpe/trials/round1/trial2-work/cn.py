import pandas as pd, numpy as np
df=pd.read_pickle('/root/work/all.pkl')
df['nv']=df.nalt/df.ndp.clip(lower=1); df['tv']=df.talt/df.tdp.clip(lower=1)
het=df[(df.ndp>=30)&(df.nv>0.2)&(df.nv<0.8)]
hom=df[(df.ndp>=30)&(df.nv>=0.9)]
print('hets',len(het),'homalt',len(hom))
tot_t=df.tdp.sum(); tot_n=df.ndp.sum()
B=10000
for c in ['chr1','chr2','chr3','chr4']:
    d=df[df.chrom==c]; h=het[het.chrom==c]
    for s in range(0,d.pos.max(),B):
        x=d[(d.pos>s)&(d.pos<=s+B)]
        hh=h[(h.pos>s)&(h.pos<=s+B)]
        r=(x.tdp.sum()/tot_t)/(x.ndp.sum()/tot_n)
        mb=np.abs(hh.tv-0.5).mean()+0.5 if len(hh) else np.nan
        mbn=np.abs(hh.nv-0.5).mean()+0.5 if len(hh) else np.nan
        print(c,s,f'{r:.3f}',len(hh),f'{mb:.3f}',f'{mbn:.3f}', f'{x.tdp.mean():.0f} {x.ndp.mean():.0f}')
