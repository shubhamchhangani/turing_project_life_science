import pandas as pd, numpy as np, json, pysam
c=pd.read_pickle('/root/work/cf4.pkl')
s=c[c.call&(c.dll>0)][['chrom','pos','ref','alt']]
df=pd.read_pickle('/root/work/all.pkl')
df=df[df.alt!='.']
e=df[(df.talt<3)&(df.talt>=1)&(df.nalt>=3)&(df.nalt/df.ndp<0.3)][['chrom','pos','ref','alt']]
print('extra',len(e))
s=pd.concat([s,e]).drop_duplicates(['chrom','pos'])
fa=pysam.FastaFile('/root/data/reference.fa')
for r in s.itertuples():
    assert fa.fetch(r.chrom,r.pos-1,r.pos).upper()==r.ref and r.alt!=r.ref and r.alt in 'ACGT'
s=s.sort_values(['chrom','pos'])
s.to_csv('/root/results/somatic_snvs.tsv',sep='\t',index=False)
json.dump({'tumor_purity':0.618,'tumor_ploidy':3.71,'normal_tumor_fraction':0.205},open('/root/results/sample_estimates.json','w'))
print(len(s))
