"""Independent scoring and bounded exhaustive checks of the search bookkeeping."""
import itertools,json
from pathlib import Path
import numpy as np
from prototype import costs,coordinate,anneal

def direct(H,data):
 starts,pos,alleles,weights=data
 return sum(min(sum(int(weights[q]) for q in range(int(starts[r]),int(starts[r+1])) if H[h,pos[q]]!=alleles[q]) for h in range(len(H))) for r in range(len(starts)-1))

rng=np.random.default_rng(61);cases=[]
for case in range(12):
 k,n,m=3,5,14
 starts=np.arange(0,m*n+1,n,dtype=np.int32);pos=np.tile(np.arange(n,dtype=np.int32),m)
 alleles=rng.integers(0,2,m*n,dtype=np.int32);weights=rng.integers(1,4,m*n,dtype=np.int32);data=(starts,pos,alleles,weights)
 dosage=rng.integers(1,k,n,dtype=np.int32)
 columns=[list(itertools.combinations(range(k),int(d))) for d in dosage]
 exact=10**9
 for choice in itertools.product(*columns):
  H=np.zeros((k,n),np.int8)
  for j,ones in enumerate(choice):H[list(ones),j]=1
  a=direct(H,data);b=int(costs(H,*data).min(1).sum());assert a==b;exact=min(exact,a)
 H=coordinate(H.copy(),*data,dosage)
 out,score=anneal(H,*data,dosage,case+1,sweeps=80)
 assert np.array_equal(out.sum(0),dosage)
 assert direct(out,data)==score
 assert score>=exact
 cases.append(dict(case=case,exhaustive_min=exact,search_score=int(score)))
Path('tools/polyploid-phase/draft/scoring-validation.json').write_text(json.dumps(dict(cases=cases,passed=True),indent=2)+'\n')
print('All 12 exhaustive scoring/bookkeeping checks passed.')
