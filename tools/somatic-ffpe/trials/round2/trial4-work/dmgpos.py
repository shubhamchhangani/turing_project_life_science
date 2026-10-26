import numpy as np, pysam
from collections import Counter
r=np.load('/root/work/ref.npy')
C=np.load('/root/work/tumor.cnt.npy')[:,:,:4]
NC=np.load('/root/work/normal.cnt.npy')[:,:,:4].sum(0)
sites=np.where((r==1)&(C[0][:,3]>=6)&(NC[:,3]==0)&(C[1][:,3]==0))[0][:300]
offs=[0,300000,550000,800000]; names=['chr1','chr2','chr3','chr4']
bam=pysam.AlignmentFile('/root/work/tumor.md.bam')
distA=[];distR=[];qa=Counter();qr=Counter();dupA=[];
for g in sites:
    ci=np.searchsorted(offs,g,side='right')-1; pos=g-offs[ci]
    for col in bam.pileup(names[ci],pos,pos+1,truncate=True,min_base_quality=0,stepper='nofilter'):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip: continue
            al=pr.alignment
            if al.is_duplicate or not al.is_read1: continue
            b=al.query_sequence[pr.query_position]
            fs=min(al.reference_start,al.next_reference_start); fl=abs(al.template_length)
            d=pos-fs; d=min(d,fl-1-d)
            ori=0 if not al.is_reverse else 1
            if ori!=0: continue
            if b=='T': distA.append(d); qa[al.query_qualities[pr.query_position]]+=1
            elif b=='C': distR.append(d); qr[al.query_qualities[pr.query_position]]+=1
print('alt dist hist',np.histogram(distA,bins=[0,5,10,20,40,80,200])[0]/len(distA))
print('ref dist hist',np.histogram(distR,bins=[0,5,10,20,40,80,200])[0]/len(distR))
print('qual alt',qa,'qual ref',qr)
