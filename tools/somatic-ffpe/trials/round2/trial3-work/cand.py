import numpy as np, pandas as pd
T=np.load('/root/work/tumor.counts.npz'); N=np.load('/root/work/normal.counts.npz')
rows=[]
B='ACGT'
for c in ['chr1','chr2','chr3','chr4']:
    t=T[c][:,:,:4]; n=N[c][:,:,:4]; ref=T[c+'_ref']; L=len(ref)
    ts=t.sum(1); idx=np.arange(L)
    tt=ts.copy(); tt[idx,ref.clip(0,3)]=-1
    for alt in range(4):
        ta=ts[:,alt]
        sel=np.nonzero((ta>=2)&(ref!=alt)&(ref<4))[0]
        for i in sel:
            rows.append((c,i+1,B[ref[i]],B[alt],t[i,0,alt],t[i,1,alt],ts[i].sum(),t[i,0].sum(),t[i,1].sum(),n[i,:,alt].sum(),n[i].sum(),n[i,0,alt],n[i,1,alt]))
df=pd.DataFrame(rows,columns='chrom pos ref alt t_o0 t_o1 t_dp t_dp0 t_dp1 n_alt n_dp n_o0 n_o1'.split())
df.to_csv('/root/work/cand.tsv',sep='\t',index=False)
print(len(df))
