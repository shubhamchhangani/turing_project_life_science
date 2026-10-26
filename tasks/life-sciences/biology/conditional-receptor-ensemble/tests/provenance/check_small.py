"""Validate partition and every joint probability against labelled spin enumeration."""
import itertools,json
from pathlib import Path
import numpy as np
from engine import distribution

records=[]
for seed in range(6):
    rng=np.random.default_rng(seed)
    xy=[[x,y] for y in range(3) for x in range(3)]
    edges=[]
    for a,(x,y) in enumerate(xy):
        for b,(u,v) in enumerate(xy):
            if (u-x,v-y) in [(1,0),(0,1),(1,1)]:edges.append([a,b,float(rng.uniform(.05,.5))])
    m={'coordinates':xy,'edges':edges,'ligand_sites':[0,4],'ligand_spins':[1,-1],
       'reporter_sites':[2,6,8]}
    got=distribution(m)
    spins=np.array(list(itertools.product([1,-1],repeat=9)))
    weights=np.exp(sum(j*spins[:,a]*spins[:,b] for a,b,j in edges))
    z=weights.sum();mask=(spins[:,0]==1)&(spins[:,4]==-1)
    probs=[]
    for pattern in range(8):
        keep=mask.copy()
        for bit,site in enumerate(m['reporter_sites']):keep&=spins[:,site]==(1 if not(pattern>>bit&1) else -1)
        probs.append(weights[keep].sum()/weights[mask].sum())
    error=max(abs(np.array(got['probabilities'])-probs))
    assert error<1e-11
    assert abs(got['log_partition']-np.log(z))<1e-11
    assert abs(got['log_clamp_probability']-np.log(weights[mask].sum()/z))<1e-11
    records.append({'seed':seed,'max_probability_error':float(error)})
report={'passed':True,'cases':records}
(Path(__file__).parent/'draft/exhaustive-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
