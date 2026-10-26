"""Iterated search on public observations; author phase is used only for reporting."""
import json,time
from pathlib import Path
import numpy as np
from prototype import generate,costs,coordinate,anneal
Htrue,data=generate(k=16,n=480,m=2400)
d=Htrue.sum(0).astype(np.int32);k,n=Htrue.shape
threshold=int(costs(Htrue,*data).min(1).sum());del Htrue
rng=np.random.default_rng(20261101);best=10**12;population=[];start=time.monotonic()
for run in range(150):
 if run<12 or run%10==0:
  H=np.zeros((k,n),np.int8)
  for j,v in enumerate(d):H[rng.permutation(k)[:v],j]=1
 else:
  parent=population[int(rng.integers(min(6,len(population))))][1];H=parent.copy()
  for move in range(int(rng.integers(2,14))):
   a,b=rng.choice(k,2,replace=False);lo=int(rng.integers(n));hi=min(n,lo+int(rng.integers(10,110)))
   H[[a,b],lo:hi]=H[[b,a],lo:hi]
  for j in rng.choice(n,size=int(rng.integers(0,30)),replace=False):rng.shuffle(H[:,j])
 H=coordinate(H,*data,d)
 H,score=anneal(H,*data,d,run+1,sweeps=300)
 assert int(costs(H,*data).min(1).sum())==score
 if score<best:
  best=score;np.save('tools/polyploid-phase/draft/best.npy',H)
 population.append((score,H));population.sort(key=lambda x:x[0]);population=population[:16]
 print(json.dumps(dict(run=run,score=int(score),best=int(best),threshold=threshold,seconds=time.monotonic()-start)),flush=True)
 if best<=threshold:
  Path('tools/polyploid-phase/draft/search-result.json').write_text(json.dumps(dict(score=int(best),threshold=threshold,run=run,seconds=time.monotonic()-start),indent=2)+'\n');break
