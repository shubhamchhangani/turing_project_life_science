"""Trusted forward predictions from submitted numerical values, never submitted code."""
import json,sys
from pathlib import Path
import numpy as np
from mixture_model import Model,POPS

submitted=json.load(sys.stdin)
m=Model(Path('/tests/data'))
decoded={}
for pop in POPS:
    values=submitted['fluxes'][pop]
    net={rid:v['net'] for rid,v in values.items()}
    exchange={rid:min(values[rid]['forward'],values[rid]['reverse']) for rid in m.flux.exch}
    decoded[pop]=(net,exchange,{**m.known[pop],**submitted['pools'][pop]})
pred=m.predict(decoded)
observed=np.array([float(r['fraction']) for r in m.meas])
sd=np.array([float(r['sd']) for r in m.meas])
ssr=float(np.sum(((pred-observed)/sd)**2))
for row in m.rates:
    if float(row['sd'])>0:
        value=sum(m.weights[p]*decoded[p][0][row['reaction_id']] for p in POPS)
        ssr+=((value-float(row['rate']))/float(row['sd']))**2
print(ssr)
