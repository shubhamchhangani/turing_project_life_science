import pysam,collections
for s in ['tumor','normal']:
    bam=pysam.AlignmentFile(f'/root/work/{s}.md.bam')
    c=collections.Counter()
    for col in bam.pileup('chr2',45195,45196,truncate=True,min_base_quality=0,stepper='nofilter',ignore_overlaps=False):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip: continue
            r=pr.alignment
            b=r.query_sequence[pr.query_position]; q=r.query_qualities[pr.query_position]
            c[(b,r.is_duplicate,q>=20)]+=1
    print(s,sorted(c.items()))
