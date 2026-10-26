import numpy as np
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
for c in ['chr1','chr2','chr3','chr4']:
    t=T[c][:,:,:4].sum(1); n=N[c][:,:,:4].sum(1); ref=T[c+'_ref']
    dt=t.sum(1); dn=n.sum(1)
    print(c,'meanT',dt.mean(),'meanN',dn.mean())
    B=10000
    for i in range(0,len(dt),B):
        sl=slice(i,i+B)
        # germline hets: normal VAF 0.3-0.7, dn>=30
        nalt=dn[sl]-n[sl][np.arange(B),ref[sl].clip(0,3)]
        # use major non-ref allele
        nn=n[sl].copy(); nn[np.arange(B),ref[sl].clip(0,3)]=-1
        a=nn.argmax(1); na=n[sl][np.arange(B),a]
        het=(dn[sl]>=30)&(na/np.maximum(dn[sl],1)>0.3)&(na/np.maximum(dn[sl],1)<0.7)
        ta=t[sl][np.arange(B),a]; tb=dt[sl]
        baf=ta[het]/np.maximum(tb[het],1)
        mb=np.median(np.abs(baf-0.5))+0.5 if het.sum() else np.nan
        print(f'{c}:{i:>6} T/N={dt[sl].mean()/dn[sl].mean():.3f} nhet={het.sum():3d} |BAF-.5|med={mb:.3f} hetNvaf={np.median(na[het]/dn[sl][het]) if het.sum() else 0:.3f}')
