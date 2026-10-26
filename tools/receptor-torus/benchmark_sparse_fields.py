"""Measure a larger original ligand-field prototype; never alter frozen trials."""
import json
from pathlib import Path
import time
import numpy as np
from sparse_fields_prototype import SparseFields

sites=[(1,1),(6,1),(11,1),(15,4),(2,6),(7,6),(11,9),(15,11),(2,12),(6,15),(11,16),(15,17)]
model=SparseFields(17,19,sites)
rng=np.random.default_rng(20261030)
couplings=[]
for (x,y),q,d in model.edges:
    value=.46+.12*np.sin(.43*x+.61*y+.3*d)+rng.uniform(-.09,.09)
    if (17*x+11*y+29*d)%113==17:value=-rng.uniform(.1,.21)
    couplings.append(round(value,9))
fields=np.array([.62,-.57,.68,-.63,.55,-.71,.59,-.52,.66,-.61,.58,-.54])
parameters=np.r_[couplings,fields]
print('even subsets',len(model.masks),'maximum update rank',max(len(p[1]) for p in model.paths),flush=True)
start=time.monotonic();value,gradient=model.evaluate(parameters);elapsed=time.monotonic()-start
folder=Path(__file__).resolve().parent/'draft/twelve-fields';folder.mkdir(exist_ok=True)
(folder/'generation.json').write_text(json.dumps({'seed':20261030,'width':17,'height':19,'sites':sites,'truth':parameters.tolist(),'log_partition':float(value),'observables':gradient.tolist(),'generation_seconds':elapsed},indent=2)+'\n')
print('generation seconds',elapsed,'logZ',value,'site means',gradient[-12:],flush=True)
