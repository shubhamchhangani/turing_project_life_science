"""Reproduce independent scientific checks without using any submitted answer."""
import csv
import json
from pathlib import Path
import time
import numpy as np
from prototype import Torus,enumerate_spins
from independent import DecoratedTorus

author=Path(__file__).resolve().parent
rng=np.random.default_rng(427819)
cases=[]
for width,height in [(3,3),(3,4),(4,4)]:
    reference=Torus(width,height)
    verifier=DecoratedTorus(width,height)
    probes=[(0,0),(width-1,height-1),(1,0),(2,1)]
    masks=[m for m in range(1,16) if m.bit_count()%2==0]
    queries=[[p for i,p in enumerate(probes) if m>>i&1] for m in masks]
    for low,high in [(.12,.3),(.3,.55),(.65,.9),(-.2,.5)]:
        couplings=rng.uniform(low,high,len(reference.edges))
        truth=enumerate_spins(reference,couplings,queries)
        direct=reference.evaluate(couplings,queries)
        moments=verifier.moments(couplings,probes)
        alternative=verifier.partition(couplings)
        errors={'reference_logz':abs(direct[0]-truth[0]),'reference_moments':float(np.max(np.abs(direct[1]-truth[1]))),'verifier_logz':abs(alternative-truth[0]),'verifier_moments':max(abs(moments[m]-v) for m,v in zip(masks,truth[1]))}
        cases.append({'width':width,'height':height,'coupling_interval':[low,high],'enumerated_assignments':2**(width*height),'max_abs_errors':errors})
        assert max(errors.values())<1e-10

independent=json.loads((author/'draft/independent.json').read_text())
reference_joint={r['pattern']:float(r['probability']) for r in csv.DictReader((author/'draft/reference/joint.csv').open())}
reference_logz=float(list(csv.DictReader((author/'draft/reference/thermodynamics.csv').open()))[0]['value'])
full={'joint_max_abs_difference':max(abs(reference_joint[f'{i:06b}']-value) for i,value in enumerate(independent['joint'])),'logz_abs_difference':abs(reference_logz-independent['log_partition']),'independent_seconds':independent['seconds']}
assert full['joint_max_abs_difference']<1e-10 and full['logz_abs_difference']<1e-10
report={'all_passed':True,'seed':427819,'small_cases':cases,'full_instance':full,'acceptance_absolute_tolerance':1e-7,'rounding_decimals':12}
(author/'draft/physics-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_passed':True,'exhaustive_cases':len(cases),'max_small_error':max(v for c in cases for v in c['max_abs_errors'].values()),'full_instance':full},indent=2))
