import pandas as pd, numpy as np, pysam, sys
c=pd.read_pickle(sys.argv[1]).reset_index(drop=True)
outp=sys.argv[2]
bams={s:pysam.AlignmentFile(f'/root/work/{s}.md.bam') for s in ['tumor','normal']}
rows=[]
def reads_at(bam,chrom,pos,ref,alt):
    L=[]
    for col in bam.pileup(chrom,pos-1,pos,truncate=True,stepper='nofilter',ignore_overlaps=False,min_base_quality=0,max_depth=100000):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip: continue
            a=pr.alignment; qpos=pr.query_position
            b=a.query_sequence[qpos]
            if b!=alt and b!=ref: continue
            qp=qpos if not a.is_reverse else a.query_length-1-qpos
            if a.template_length>0: fs=a.reference_start; fe=a.reference_start+a.template_length
            else: fe=a.reference_end; fs=a.reference_end+a.template_length
            L.append((a.query_name,b==alt,a.query_qualities[qpos],qp,a.is_read1!=a.is_reverse,a.is_duplicate,fs,fe))
    return L
for i,r in enumerate(c.itertuples()):
    if i%2000==0: print(i,flush=True)
    d={}
    for s,bam in bams.items():
        L=reads_at(bam,r.chrom,r.pos,r.ref,r.alt)
        frag={}  # nondup fragments
        fam={}
        for q,al,bq,qp,o,dup,fs,fe in L:
            k=(fs,fe)
            fam.setdefault(k,[]).append((q,al,o))
            if dup: continue
            f=frag.get(q)
            if f is None: frag[q]=[al,al,o,bq,qp,fs,fe]
            else:
                f[0]=f[0] or al; f[1]=f[1] and al; f[3]=max(f[3],bq); f[4]=min(f[4],qp)
        alts=[f for f in frag.values() if f[0]]
        d[s+'_nfr']=len(frag); d[s+'_afr']=len(alts)
        if s=='tumor':
            d['a_f1r2']=sum(1 for f in alts if f[2])
            d['disc']=sum(1 for f in frag.values() if f[0]!=f[1])
            d['abq']=np.mean([f[3] for f in alts]) if alts else np.nan
            d['aqp']=np.mean([f[4] for f in alts]) if alts else np.nan
            d['dend']=np.mean([min(r.pos-1-f[5],f[6]-r.pos) for f in alts]) if alts else np.nan
            nfa=nb=nm=0
            for k,g in fam.items():
                ga=[x for x in g if x[1]]
                if ga:
                    nfa+=1
                    if len(set(x[2] for x in ga))==2: nb+=1
                    if len(ga)<len(g): nm+=1
            d['nfam_alt']=nfa; d['nfam_both']=nb; d['nfam_mixed']=nm
            d['altreads_all']=sum(1 for x in L if x[1])
    rows.append(d)
F=pd.DataFrame(rows)
pd.concat([c,F],axis=1).to_pickle(outp)
print('done')
