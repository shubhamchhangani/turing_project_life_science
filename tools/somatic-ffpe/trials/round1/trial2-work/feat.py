import pandas as pd, numpy as np, pysam, sys, pickle
c=pd.read_pickle('/root/work/cand.pkl')
# also include sites with normal alt support where tumour alt <3 ? keep simple
out={}
for samp in ['tumor','normal']:
    bam=pysam.AlignmentFile(f'/root/work/{samp}.md.bam')
    res=[]
    for r in c.itertuples():
        reads=[]
        for col in bam.pileup(r.chrom, r.pos-1, r.pos, truncate=True, stepper='nofilter', ignore_overlaps=False, min_base_quality=0, max_depth=100000):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip: continue
                a=pr.alignment
                b=a.query_sequence[pr.query_position]
                if b!=r.alt and b!=r.ref: continue
                q=a.query_qualities[pr.query_position]
                qp=pr.query_position if not a.is_reverse else a.query_length-1-pr.query_position
                f1r2 = (a.is_read1 != a.is_reverse)
                # fragment coordinates
                if a.template_length>0: fs=a.reference_start; fe=a.reference_start+a.template_length
                else: fe=a.reference_end; fs=a.reference_end+a.template_length
                reads.append((a.query_name, b==r.alt, q, qp, f1r2, a.is_duplicate, fs, fe, a.mapping_quality))
        res.append(reads)
    out[samp]=res
pickle.dump((c,out),open('/root/work/feat.pkl','wb'))
