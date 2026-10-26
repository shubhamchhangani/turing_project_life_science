"""Resource-limited author harness; reads only the public model as scientific input."""
import csv,json,time
from pathlib import Path
from scaffold_engine import distribution

start=time.monotonic()
model=json.loads(Path('/root/data/model.json').read_text())
result=distribution(model,lambda i,n:print('moments',i,n,round(time.monotonic()-start,2),flush=True))
result['seconds']=time.monotonic()-start
out=Path('/root/results');out.mkdir(parents=True,exist_ok=True)
(out/'reference.json').write_text(json.dumps(result,indent=2)+'\n')
with (out/'thermodynamics.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['metric','value'])
    for key in ['log_partition','log_clamp_probability']:w.writerow([key,format(result[key],'.15g')])
with (out/'reporter_distribution.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['pattern','probability'])
    for i,p in enumerate(result['probabilities']):w.writerow([''.join(str((i>>j)&1) for j in range(8)),format(p,'.15g')])
print('Reference complete:',result['seconds'],flush=True)
