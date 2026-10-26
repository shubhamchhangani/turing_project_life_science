import pysam, numpy as np, sys
s=sys.argv[1]; minq=int(sys.argv[2]) if len(sys.argv)>2 else 20
fa=pysam.FastaFile('/root/data/reference.fa')
lut=np.full(256,4,np.uint8)
for i,c in enumerate('ACGT'): lut[ord(c)]=i
bam=pysam.AlignmentFile(f'/root/work/{s}.md.bam')
out={}
for chrom in fa.references:
    L=fa.get_reference_length(chrom)
    ref=lut[np.frombuffer(fa.fetch(chrom).upper().encode(),np.uint8)]
    # counts[pos, orient(0: read1 fwd,1: read1 rev), base5]
    cnt=np.zeros((L,2,5),np.int32)
    # positional: alt count from bases within 10bp of fragment end-of-read
    pend={}
    def proc(r1,r2):
        o=0 if not r1.is_reverse else 1
        a=[]
        for r in (r1,r2):
            st=r.reference_start
            b=lut[np.frombuffer(r.query_sequence.encode(),np.uint8)]
            q=np.asarray(r.query_qualities)
            b=b.copy(); b[q<minq]=4
            a.append((st,b))
        (s1,b1),(s2,b2)=a
        lo=min(s1,s2); hi=max(s1+len(b1),s2+len(b2))
        m=np.full(hi-lo,5,np.int16)
        m[s1-lo:s1-lo+len(b1)]=b1
        seg=m[s2-lo:s2-lo+len(b2)]
        ov=seg!=5
        new=b2.astype(np.int16).copy()
        dis=ov&(seg!=new)
        # where one is N(4) and other valid, take valid
        take=ov&(seg==4)&(new!=4)
        keepold=ov&(new==4)&(seg!=4)
        res=new.copy()
        res[keepold]=seg[keepold]
        res[dis & ~take & ~keepold]=4
        m[s2-lo:s2-lo+len(b2)]=res
        idx=np.nonzero((m!=5)&(m!=4))[0]
        np.add.at(cnt[:,o,:],(idx+lo,m[idx]),1)
    for r in bam.fetch(chrom):
        if r.is_duplicate or r.is_secondary or r.is_supplementary or r.mapping_quality<20: continue
        n=r.query_name
        if n in pend:
            m=pend.pop(n)
            r1,r2=(m,r) if m.is_read1 else (r,m)
            proc(r1,r2)
        else: pend[n]=r
    out[chrom]=cnt; out[chrom+'_ref']=ref
    print(chrom,len(pend),flush=True)
np.savez_compressed(f'/root/work/{s}.counts.npz',**out)
