import numpy as np, pysam
fa=pysam.FastaFile('/root/data/reference.fa')
chroms=['chr1','chr2','chr3','chr4']
ref=np.concatenate([np.frombuffer(fa.fetch(c).encode(),np.uint8) for c in chroms])
lut=np.full(256,4,np.uint8)
for i,b in enumerate(b'ACGT'): lut[b]=i
r=lut[ref]; np.save('/root/work/ref.npy',r)
for s in ['tumor','normal']:
    c=np.load(f'/root/work/{s}.cnt.npy'); t=c.sum(0)[:,:4]
    d=t.sum(1); print(s,'depth mean',d.mean(),'median',np.median(d))
    # bins of 10kb
    b=d.reshape(-1,10000).mean(1); print(np.round(b).astype(int))
