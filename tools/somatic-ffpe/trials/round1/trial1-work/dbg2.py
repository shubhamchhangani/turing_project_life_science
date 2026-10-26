import pysam,collections,sys
c,p=sys.argv[1],int(sys.argv[2])-1
for s in ['tumor','normal']:
    bam=pysam.AlignmentFile(f'/root/work/{s}.md.bam')
    fr=collections.defaultdict(list)
    for col in bam.pileup(c,p,p+1,truncate=True,min_base_quality=0,stepper='nofilter',ignore_overlaps=False,max_depth=100000):
        for pr in col.pileups:
            if pr.is_del: continue
            r=pr.alignment
            fr[r.query_name].append((r.query_sequence[pr.query_position],r.query_qualities[pr.query_position],'R1' if r.is_read1 else 'R2','-' if r.is_reverse else '+',pr.query_position,r.is_duplicate,r.reference_start,r.template_length))
    cnt=collections.Counter()
    for n,v in fr.items():
        key=tuple(sorted((x[0],x[1]) for x in v))
        cnt[key]+=1
    print(s); 
    for k,v in sorted(cnt.items(),key=lambda x:-x[1]): print('  ',k,v)
    for n,v in fr.items():
        if any(x[0]!=v[0][0] for x in v) or any(x[1]<20 for x in v) or len(v)==1 and v[0][0]!='G': pass
    alts=[(n,v) for n,v in fr.items() if any(x[0] not in 'ACGT'[0:0]+pysam.FastaFile('/root/data/reference.fa').fetch(c,p,p+1) for x in v)]
    for n,v in alts[:40]: print('   ',n,v)
