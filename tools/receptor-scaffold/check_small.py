"""Independent exhaustive spin sums test both signed implementations."""
import itertools,json
from pathlib import Path
import numpy as np
from scaffold_engine import distribution
from independent_source import calculate

cases=[]
for seed in range(12):
    rng=np.random.default_rng(seed);xy=[[x,y] for y in range(3) for x in range(3)];edges=[]
    for a,(x,y) in enumerate(xy):
        for b,(u,v) in enumerate(xy):
            if (u-x,v-y) in [(1,0),(0,1),(1,1)]:edges.append([a,b,float(rng.uniform(.1,.7)*rng.choice([-1,1]))])
    for a,b in [(0,8),(2,6),(1,7)]:edges.append([a,b,float(rng.uniform(-.5,.5))])
    fields=[[s,float(rng.uniform(-.3,.3))] for s in range(9)]
    model={'coordinates':xy,'edges':edges,'fields':fields,'ligand_sites':[0,4],'ligand_spins':[1,-1],'reporter_sites':[2,6,8]}
    spins=np.array(list(itertools.product([1,-1],repeat=9)))
    weights=np.exp(sum(j*spins[:,a]*spins[:,b] for a,b,j in edges)+sum(h*spins[:,s] for s,h in fields))
    clamp=(spins[:,0]==1)&(spins[:,4]==-1);pc=weights[clamp].sum()/weights.sum();expected=[]
    for i in range(8):
        keep=clamp.copy()
        for j,site in enumerate(model['reporter_sites']):keep&=spins[:,site]==(1 if not(i>>j&1) else -1)
        expected.append(weights[keep].sum()/weights[clamp].sum())
    errors=[]
    for method in [distribution,calculate]:
        got=method(model);err=max(abs(np.array(got['probabilities'])-expected));errors.append(float(err))
        assert err<1e-10
        assert abs(got['log_partition']-np.log(weights.sum()))<1e-10
        assert abs(got['log_clamp_probability']-np.log(pc))<1e-10
    cases.append({'seed':seed,'reference_max_error':errors[0],'independent_source_max_error':errors[1]})
r={'passed':True,'cases':cases}
(Path(__file__).parent/'draft/exhaustive-validation.json').write_text(json.dumps(r,indent=2)+'\n')
print('12 exhaustive cases passed for both signed implementations.')
