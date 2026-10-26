"""Check compiled arithmetic on sampled subsets of the actual 26-site instance."""
import json
from pathlib import Path
import numpy as np
from scaffold_engine import split
from signed_engine import Signed
from compiled_moments import moments
from source_kernel import source_moments

p=Path(__file__).parent/'draft'
model=json.loads((p/'model.json').read_text());base,_=split(model);engine=Signed(base)
paths=np.zeros((len(engine.sites),2),np.uint64)
for bit,path in enumerate(engine.paths):
    for edge in path:
        j=engine.lookup[edge];paths[bit,j//64]^=np.uint64(1)<<np.uint64(j%64)
edges=np.array(engine.union);weights=engine.t[edges]
delta=engine.r[2*edges]*(weights-1/weights)
rng=np.random.default_rng(20261026);worst_reference=0.;worst_update=0.;seen=set()
for _ in range(16):
    cohort=rng.choice(len(engine.sites),8,replace=False)
    compiled=moments(engine.inverse,paths[cohort],delta,delta/weights,256)
    updates=source_moments(engine.inverse,paths[cohort],delta,weights,256)
    for mask in range(256):
        original=sum(1<<int(cohort[j]) for j in range(8) if mask>>j&1)
        if mask.bit_count()%2:
            assert compiled[mask]==0 and updates[mask]==0
            continue
        seen.add(original);want=engine.moment(original)
        worst_reference=max(worst_reference,abs(want-compiled[mask]))
        worst_update=max(worst_update,abs(want-updates[mask]))
assert worst_reference<1e-9 and worst_update<1e-9
report={'passed':True,'cohorts':16,'unique_even_subsets':len(seen),
        'max_compiled_vs_python_moment_error':worst_reference,
        'max_rank_update_vs_python_moment_error':worst_update,
        'scope':'Arithmetic spot checks on full graph; not full joint-distribution validation.'}
(p/'large-spot-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
