"""Build original data and a draft task. This does not approve or ZIP the task."""
import ast,csv,json,shutil
from pathlib import Path
import numpy as np
from prototype import generate
ROOT=Path('tasks/life-sciences/biology/polyploid-haplotype-assembly')
for d in ['environment/data','tests/data','solution']: (ROOT/d).mkdir(parents=True,exist_ok=True)
H,data=generate(k=16,n=480,m=2400);starts,pos,alleles,weights=data
# Independent direct computation on the original synthetic witness, not solver output.
certificate=sum(min(sum(int(weights[q]) for q in range(starts[r],starts[r+1]) if H[h,pos[q]]!=alleles[q]) for h in range(len(H))) for r in range(len(starts)-1))
model=dict(ploidy=len(H),variants=H.shape[1],reads=len(starts)-1,maximum_weighted_corrections=certificate)
for dest in [ROOT/'environment/data',ROOT/'tests/data']:
 (dest/'model.json').write_text(json.dumps(model,indent=2)+'\n')
 with (dest/'variants.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['variant_id','alt_dosage']);w.writerows((f'v{j+1:04d}',int(H[:,j].sum())) for j in range(H.shape[1]))
 with (dest/'read_calls.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['read_id','variant_id','allele','weight'])
  for r in range(len(starts)-1):
   w.writerows((f'r{r+1:05d}',f'v{pos[q]+1:04d}',int(alleles[q]),int(weights[q])) for q in range(starts[r],starts[r+1]))
Path('tools/polyploid-phase/draft/generation.json').write_text(json.dumps(dict(seed=20261031,requested_variants=480,monomorphic_variants_removed=480-H.shape[1],model=model,original_haplotypes=[''.join(map(str,h)) for h in H]),indent=2)+'\n')
base=Path('tasks/life-sciences/biology/spatial-ligand-cooperativity/task.toml').read_text()
base=base.replace('artifacts = ["/root/results/couplings.csv", "/root/results/biases.csv"]','artifacts = ["/root/results/haplotypes.csv"]')
base=base.replace('spatial-ligand-cooperativity','polyploid-haplotype-assembly').replace('Infer heterogeneous receptor cooperativity and twelve spatially separated ligand biases from equilibrium moments on a periodic membrane lattice.','Assemble sixteen homologous chromosome phases from noisy overlapping DNA reads and known allele dosages.').replace('keywords = ["inverse-statistical-mechanics", "membrane-receptors", "ligand-bias", "cooperativity"]','keywords = ["polyploidy", "haplotype-assembly", "DNA-sequencing", "weighted-MEC"]').replace('computational-biophysics','computational-genomics').replace('timeout_sec = 900','timeout_sec = 120')
(ROOT/'task.toml').write_text(base)
(ROOT/'environment/Dockerfile').write_text('''FROM python:3.11.9-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends bash gcc g++ && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==1.26.4 scipy==1.13.1 numba==0.60.0
WORKDIR /root
COPY data /root/data
RUN mkdir -p /root/results
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
''')
(ROOT/'tests/Dockerfile').write_text('''FROM python:3.11.9-slim-bookworm
RUN pip install --no-cache-dir pytest==8.4.1 pytest-json-ctrf==0.3.5
RUN mkdir -p /root/results /logs/verifier
COPY . /tests/
ENV PYTHONDONTWRITEBYTECODE=1
''')
(ROOT/'tests/test.sh').write_text('''#!/bin/bash
set -euo pipefail
mkdir -p /logs/verifier
if timeout --kill-after=5s 100s pytest --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py -rA; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
''')
(ROOT/'solution/solve.sh').write_text('''#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 /solution/solve.py
''')
# Solver kernels exclude generation and all planted chromosome phases.
tree=ast.parse(Path('tools/polyploid-phase/prototype.py').read_text())
functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'costs','coordinate','anneal'}]
(ROOT/'solution/search.py').write_text('import numpy as np\nfrom numba import njit\n\n'+'\n\n'.join(ast.unparse(n) for n in functions)+'\n')
print(json.dumps(model,indent=2))
