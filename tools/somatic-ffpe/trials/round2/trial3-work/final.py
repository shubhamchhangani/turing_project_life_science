import pandas as pd, numpy as np, json
d=pd.read_pickle('/root/work/post4_0.61_0.205.pkl')
s=d[d.pmal>0.5].copy()
s=s.sort_values('pmal',ascending=False).drop_duplicates(['chrom','pos'])
print(len(s), 'exp prec',s.pmal.mean())
ct=s[s.cl=='C>T']
print('C>T calls',len(ct),'clean sum',ct.ac.sum(),'dirty sum',ct.ad.sum(),'Dc',ct.Dc.sum(),'Dd',ct.Dd.sum())
print('low ad with ac>=3:',((ct.ad==0)&(ct.ac>=3)).sum())
s=s.sort_values(['chrom','pos'])
s[['chrom','pos','ref','alt']].to_csv('/root/results/somatic_snvs.tsv',sep='\t',index=False)
json.dump({'tumor_purity':0.61,'tumor_ploidy':3.71,'normal_tumor_fraction':0.205},open('/root/results/sample_estimates.json','w'))
