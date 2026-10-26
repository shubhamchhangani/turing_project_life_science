"""Exact scaffold reweighting of a signed planar receptor marginal."""
import numpy as np
from signed_engine import Signed
from compiled_moments import enumerate_moments,character_inverse,reweight


def split(model):
    xy=np.array(model['coordinates']);local=[];extra=[]
    for a,b,j in model['edges']:
        dx,dy=xy[b]-xy[a]
        if (int(dx),int(dy)) in [(1,0),(0,1),(1,1)]:local.append([a,b,j])
        else:extra.append([a,b,j])
    sites=model['ligand_sites']+model['reporter_sites']
    auxiliary=sorted(({i for edge in extra for i in edge[:2]}|{i for i,h in model['fields']})-set(sites))
    base=dict(model,edges=local,reporter_sites=model['reporter_sites']+auxiliary)
    return base,extra


def distribution(model,progress=None):
    base,extra=split(model);engine=Signed(base);count=1<<len(engine.sites)
    if progress:progress(0,count)
    values=enumerate_moments(engine)
    if progress:progress(count,count)
    character_inverse(values)
    assert values.min()>-1e-12
    lookup={s:j for j,s in enumerate(engine.sites)}
    nlig=len(model['ligand_sites']);clamp=sum((s==-1)<<i for i,s in enumerate(model['ligand_spins']))
    normalization,clamp_total,probabilities=reweight(values,
        np.array([lookup[s] for s,h in model['fields']]),np.array([h for s,h in model['fields']]),
        np.array([lookup[a] for a,b,j in extra]),np.array([lookup[b] for a,b,j in extra]),
        np.array([j for a,b,j in extra]),nlig,clamp,len(model['reporter_sites']))
    pc=clamp_total/normalization
    return {'log_partition':float(engine.logz+np.log(normalization)),
            'log_clamp_probability':float(np.log(pc)),
            'probabilities':probabilities.tolist()}
