import csv,json
from pathlib import Path
from engine import distribution

model=json.loads(Path('/root/data/model.json').read_text())
r=distribution(model)
out=Path('/root/results');out.mkdir(parents=True,exist_ok=True)
with (out/'thermodynamics.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['metric','value'])
 for key in ['log_partition','log_clamp_probability']:w.writerow([key,format(r[key],'.15g')])
with (out/'reporter_distribution.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['pattern','probability'])
 for i,p in enumerate(r['probabilities']):w.writerow([''.join(str((i>>j)&1) for j in range(8)),format(p,'.15g')])
