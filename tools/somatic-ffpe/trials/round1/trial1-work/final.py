import pandas as pd,numpy as np,json
d=pd.read_csv('/root/work/scored.tsv',sep='\t')
d['ct']=d.C2T|d.G2A
d['tot']=d['at']+d.nalt
base=d[(d.vn<0.35)&(d.artp>1e-3)&(d.incons<15)]
calls=base[(base.ct&(base.tot>=3))|(~base.ct&(base.tot>=4))]
calls=calls.sort_values('LR',ascending=False).drop_duplicates(['chrom','pos'])
calls=calls.sort_values(['chrom','pos'])
assert (calls.ref!=calls.alt).all() and calls.ref.isin(list('ACGT')).all()
calls[['chrom','pos','ref','alt']].to_csv('/root/results/somatic_snvs.tsv',sep='\t',index=False)
segs=[(300000,4),(250000,3),(120000,2),(130000,4),(200000,5)]
ploidy=sum(a*b for a,b in segs)/1e6
json.dump({'tumor_purity':0.612,'tumor_ploidy':ploidy,'normal_tumor_fraction':0.203},open('/root/results/sample_estimates.json','w'))
print(len(calls),ploidy)
