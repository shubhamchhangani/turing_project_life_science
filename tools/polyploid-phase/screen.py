"""Screen input-only recovery, preserving every result rather than cherry-picking trials."""
import json,time
import numpy as np
from prototype import generate,costs,fit
from pathlib import Path
results=[]
for k,n,m in [(12,360,1800),(16,480,2400),(12,360,3600),(16,480,4800)]:
 H,data=generate(k=k,n=n,m=m);d=H.sum(0).astype(np.int32);planted=int(costs(H,*data).min(1).sum());start=time.monotonic()
 print('CASE',k,n,m,'PLANTED',planted,flush=True)
 out,score=fit(data,d,k,restarts=4,sweeps=240)
 assert int(costs(out,*data).min(1).sum())==score
 r=dict(k=k,n=n,m=m,planted=planted,reference_cost=int(score),seconds=time.monotonic()-start)
 results.append(r);print('RESULT',json.dumps(r),flush=True)
 Path('tools/polyploid-phase/draft/screen-results.json').write_text(json.dumps(results,indent=2)+'\n')
