import pysam, numpy as np, sys
s=sys.argv[1]; minq=int(sys.argv[2]) if len(sys.argv)>2 else 20
fa=pysam.FastaFile('/root/data/reference.fa')
chroms=['chr1','chr2','chr3','chr4']; L={c:fa.get_reference_length(c) for c in chroms}
off={};o=0
for c in chroms: off[c]=o;o+=L[c]
N=o
lut=np.full(256,4,np.uint8)
for i,b in enumerate(b'ACGT'): lut[b]=i
cnt=np.zeros((2,N,5),np.int32)  # orient (read1 fwd=0), pos, base(4=N/conflict)
bam=pysam.AlignmentFile(f'/root/work/{s}.md.bam')
pend={}
for r in bam.fetch(until_eof=True):
    if r.flag & (1024|256|2048): continue
    q=np.frombuffer(r.query_sequence.encode(),np.uint8)
    b=lut[q].copy(); qual=np.array(r.query_qualities,np.uint8)
    b[qual<minq]=4
    st=off[r.reference_name]+r.reference_start
    nm=r.query_name
    if nm in pend:
        st2,b2,ori=pend.pop(nm)
        # merge
        s0=min(st,st2); e0=max(st+150,st2+150)
        arr=np.full(e0-s0,5,np.uint8) # 5=empty
        arr[st2-s0:st2-s0+150]=b2
        seg=arr[st-s0:st-s0+150]
        m=seg==5
        seg2=np.where(m|(seg==4),b,np.where((seg==b)|(b==4),seg,4))
        arr[st-s0:st-s0+150]=seg2
        idx=np.arange(s0,e0)
        ok=arr<5
        np.add.at(cnt[ori],(idx[ok],arr[ok]),1)
    else:
        ori=0 if (r.is_read1 != r.is_reverse) else 1
        pend[nm]=(st,b,ori)
np.save(f'/root/work/{s}.cnt.npy',cnt)
print(s,'done',len(pend))
