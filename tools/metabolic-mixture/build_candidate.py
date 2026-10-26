"""Assemble the unzipped experimental task; author tools/truth stay outside."""
import ast,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DEST=ROOT/'tasks/life-sciences/biology/metabolic-population-deconvolution'
OLD=ROOT/'tasks/life-sciences/biology/metabolic-flux-analysis-v2'
for directory in ('environment/data','solution','tests/data'):(DEST/directory).mkdir(parents=True,exist_ok=True)
for p in (HERE/'draft/data').iterdir():
    for directory in ('environment/data','tests/data'):shutil.copy2(p,DEST/directory/p.name)
for directory in ('solution','tests'):
    for name in ('dynamic_model.py','forward_model.py','network_model.py'):shutil.copy2(HERE/name,DEST/directory/name)
tree=ast.parse((HERE/'prototype.py').read_text())
nodes={n.name:n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef))}
model='''"""Transient population dynamics and molecule-weighted bulk observations.

This module contains model equations only, with no fitted or injected answers.
"""
import csv
from collections import OrderedDict
from pathlib import Path
import numpy as np
import dynamic_model as D
import forward_model as M
import network_model as N
POPS = ('A', 'B')
'''
model+='\n\n'.join(ast.unparse(nodes[n]) for n in ('read_rows','write_rows','Model'))+'\n'
for directory in ('solution','tests'):(DEST/directory/'mixture_model.py').write_text(model)
fit='''"""Fit only the supplied measurements from independent starting points."""
import sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from mixture_model import Model, POPS
'''+ast.unparse(nodes['fit'])+'''
if __name__ == '__main__':
    fit(Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
'''
# The fitter writes a small numerical diagnostic, not an expected-answer file.
fit=fit.replace('import sys,time','import sys,time,json')
(DEST/'solution/mixture_fit.py').write_text(fit)
(DEST/'solution/solve.sh').write_text('#!/bin/bash\nset -euo pipefail\npython /solution/mixture_fit.py /root/data /root/results 2 17\n')
for directory in ('environment','tests'):shutil.copy2(OLD/directory/'Dockerfile',DEST/directory/'Dockerfile')
shutil.copy2(OLD/'tests/test.sh',DEST/'tests/test.sh')
shutil.copy2(HERE/'draft/data/network.tsv',DEST/'tests/network.tsv')
truth=json.loads((HERE/'draft/injected_truth.json').read_text())
expected={'graded':{},'pools':{}}
for pop,value in truth.items():
    for rid,tol in [('G6PDH',.45),('PDH',.12),('PC',.12),('ME',.08),('CS',.10),('IDH',.20)]:
        reverse=rid=='IDH'
        expected['graded'][pop+'_'+rid]={'population':pop,'reaction_id':rid,'column':'reverse' if reverse else 'net','true':value['exchange'][rid]+max(-value['net'][rid],0) if reverse else value['net'][rid],'rel_tol':tol}
    expected['pools'][pop]={m:value['pools'][m] for m in ('OAA','Cit','AKG','Mal')}
(DEST/'tests/expected_fluxes.json').write_text(json.dumps(expected,indent=2)+'\n')
print('Assembled candidate directory; no ZIP created:',DEST)
