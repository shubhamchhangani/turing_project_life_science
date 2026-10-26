import pysam, numpy as np, sys, collections, pickle
sample=sys.argv[1]
bam=pysam.AlignmentFile(f'/root/work/{sample}.md.bam')
ref=pysam.FastaFile('/root/data/reference.fa')
out={}
for chrom in ref.references:
    L=ref.get_reference_length(chrom)
    seq=ref.fetch(chrom).upper()
    depth=np.zeros(L+1,np.int32)
    pending={}
    alts=collections.defaultdict(lambda: np.zeros(6,np.int32))  # key (pos,alt): [F1R2, F2R1, lowq/discordant, dupreads]
    for r in bam.fetch(chrom):
        if r.is_duplicate:
            continue
        mm={}
        q=r.query_qualities; s=r.query_sequence
        for qp,rp in r.get_aligned_pairs(matches_only=True):
            b=s[qp]
            if b!=seq[rp]:
                mm[rp]=(b,q[qp])
        rec=(r.reference_start,r.reference_end,mm)
        n=r.query_name
        if n in pending:
            a=pending.pop(n)
            b1=rec
            s1,e1,m1=a; s2,e2,m2=b1
            # depth over union
            lo=min(s1,s2); hi=max(e1,e2)
            # assume overlapping/contiguous; if gap, add separately
            if s2>e1 or s1>e2:
                depth[s1]+=1;depth[e1]-=1;depth[s2]+=1;depth[e2]-=1
            else:
                depth[lo]+=1;depth[hi]-=1
            r1rev = r.is_reverse if r.is_read1 else (not r.is_reverse)
            ori = 1 if r1rev else 0
            keys=set(m1)|set(m2)
            for p in keys:
                in1 = s1<=p<e1; in2 = s2<=p<e2
                x=m1.get(p); y=m2.get(p)
                if in1 and in2:
                    if x and y and x[0]==y[0]:
                        alts[(p,x[0])][ori]+=1
                        alts[(p,x[0])][4]+=1
                    else:
                        z = x or y
                        if x and y:
                            alts[(p,x[0])][2]+=1; alts[(p,y[0])][2]+=1
                        elif z[1]>=20: alts[(p,z[0])][3]+=1
                        else: alts[(p,z[0])][2]+=1
                else:
                    z=x or y
                    if z[1]>=20: alts[(p,z[0])][ori]+=1
                    else: alts[(p,z[0])][2]+=1
        else:
            pending[n]=rec
    depth=np.cumsum(depth)
    out[chrom]=(depth[:L],dict(alts))
    print(chrom,len(alts),file=sys.stderr)
pickle.dump(out,open(f'/root/work/{sample}.scan2.pkl','wb'))
