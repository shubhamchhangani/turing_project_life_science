import numpy as np, json
A=np.load('/root/work/postg2_1.npz'); B=np.load('/root/work/post_1.npz')
def cs(d,t=0.5):
    m=d['post']>t; return {(int(v),int(a)):float(pp) for v,a,pp in zip(d['v'][m],d['al'][m],d['post'][m])}
a=cs(A); b=cs(B)
print(len(a),len(b),len(set(a)&set(b)))
for k in set(a)^set(b): print(k, a.get(k), b.get(k))
# one call per position
best={}
for (v,al),pp in a.items():
    if v not in best or pp>best[v][1]: best[v]=(al,pp)
print('positions',len(best))
r=np.load('/root/work/ref.npy'); B4='ACGT'
offs=[0,300000,550000,800000]; names=['chr1','chr2','chr3','chr4']
import os; os.makedirs('/root/results',exist_ok=True)
with open('/root/results/somatic_snvs.tsv','w') as fh:
    fh.write('chrom\tpos\tref\talt\n')
    for v in sorted(best):
        ci=np.searchsorted(offs,v,side='right')-1
        assert r[v]!=best[v][0]
        fh.write(f'{names[ci]}\t{v-offs[ci]+1}\t{B4[r[v]]}\t{B4[best[v][0]]}\n')
segs=[(0,300000,4),(300000,550000,3),(550000,670143,2),(670143,800000,4),(800000,1000000,5)]
ploidy=sum((e-s)*c for s,e,c in segs)/1e6
json.dump({'tumor_purity':0.615,'tumor_ploidy':round(ploidy,4),'normal_tumor_fraction':0.2},open('/root/results/sample_estimates.json','w'))
