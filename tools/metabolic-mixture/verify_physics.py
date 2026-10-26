"""Independently calculate pooled ion counts before peak normalization."""
import json,re
from pathlib import Path
import numpy as np
from prototype import Model,read_rows,POPS
import dynamic_model as D

root=Path(__file__).resolve().parent/'draft'
m=Model(root/'data')
truth=json.loads((root/'injected_truth.json').read_text())
decoded={p:(v['net'],v['exchange'],v['pools']) for p,v in truth.items()}
abundances={}
for row in read_rows(root/'data','isotope_abundances.csv'):
    abundances.setdefault(row['element'],{})[int(row['mass_shift'])]=float(row['abundance'])
fragments={r['fragment']:r for r in read_rows(root/'data','fragments.csv')}
derivatives={}
for frag,row in fragments.items():
    n=len(row['backbone_carbons'].split(';'))
    polynomial=np.ones(1)
    for element,count in re.findall(r'([A-Z][a-z]?)(\d*)',row['ion_formula']):
        count=int(count or '1')-(n if element=='C' else 0)
        atom=np.zeros(max(abundances[element])+1)
        for shift,probability in abundances[element].items():atom[shift]=probability
        for _ in range(count):polynomial=np.convolve(polynomial,atom)
    derivatives[frag]=polynomial
correct={};wrong={};max_normalization_error=0
for exp in sorted({r['experiment'] for r in m.meas}):
    ts=sorted({float(r['time_min']) for r in m.meas if r['experiment']==exp})
    schedule=[r for r in m.schedule if r['experiment']==exp]
    sims={p:D.simulate(*decoded[p],schedule,ts,method='RK45',rtol=1e-9,atol=1e-12) for p in POPS}
    for t in ts:
        for frag,row in fragments.items():
            met=row['metabolite'];carbons=[int(v)-1 for v in row['backbone_carbons'].split(';')]
            n=len(carbons);pooled=np.zeros(n+3);equal_weight=np.zeros(n+3)
            for p in POPS:
                dist=sims[p][t][met]
                max_normalization_error=max(max_normalization_error,abs(dist.sum()-1))
                backbone=np.zeros(n+1)
                for isotopomer,probability in enumerate(dist):
                    backbone[sum((isotopomer>>c)&1 for c in carbons)]+=probability
                ion=np.convolve(backbone,derivatives[frag])[:n+3]
                pooled+=m.weights[p]*decoded[p][2][met]*ion
                equal_weight+=m.weights[p]*ion
            correct[exp,t,frag]=pooled/pooled.sum()
            wrong[exp,t,frag]=equal_weight/equal_weight.sum()
def assemble(table):
    return np.array([table[r['experiment'],float(r['time_min']),r['fragment']][int(r['isotopologue'][1:])] for r in m.meas])
actual=m.predict(decoded);expected=assemble(correct);bad=assemble(wrong)
observed=np.array([float(r['fraction']) for r in m.meas]);sd=np.array([float(r['sd']) for r in m.meas])
report={'max_abs_difference_independent_pooled_ion_counts':float(np.max(np.abs(actual-expected))),
        'max_pool_normalization_error':float(max_normalization_error),
        'truth_mid_ssr':float(np.sum(((expected-observed)/sd)**2)),
        'omit_pool_amount_weighting_mid_ssr_at_truth':float(np.sum(((bad-observed)/sd)**2)),
        'scope':'Independent observation-model implementation and RK45 integration; reaction atom specification is shared. Negative control is not refitted.'}
assert report['max_abs_difference_independent_pooled_ion_counts']<1e-6,report
assert report['max_pool_normalization_error']<1e-6,report
(root/'physics_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
