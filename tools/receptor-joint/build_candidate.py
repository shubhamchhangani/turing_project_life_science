"""Build a reviewable candidate; never create a release ZIP."""
import json,shutil
from pathlib import Path

root=Path('tasks/life-sciences/biology/conditional-receptor-ensemble')
author=Path('tools/receptor-joint')
template=Path('tasks/life-sciences/biology/polyploid-phase-posterior')
for p in ['environment/data','tests/data','tests/provenance','solution']:(root/p).mkdir(parents=True,exist_ok=True)
for name in ['environment/Dockerfile','tests/Dockerfile','tests/test.sh']:
    shutil.copy2(template/name,root/name)
shutil.copy2(author/'draft/model.json',root/'environment/data/model.json')
shutil.copy2(author/'engine.py',root/'solution/engine.py')
for name in ['generate.py','independent_source.py','check_small.py']:
    shutil.copy2(author/name,root/'tests/provenance'/name)
# The enumeration checker imports the reference only for the explicit cross-check.
shutil.copy2(author/'engine.py',root/'tests/provenance/engine.py')
shutil.copy2(author/'draft/model.json',root/'tests/data/model.json')
if (author/'draft/expectations.json').exists():shutil.copy2(author/'draft/expectations.json',root/'tests/data/expectations.json')
config=(template/'task.toml').read_text().replace('polyploid-phase-posterior','conditional-receptor-ensemble')
config=config.replace('artifacts = ["/root/results/evidence.csv", "/root/results/read_posteriors.csv"]','artifacts = ["/root/results/thermodynamics.csv", "/root/results/reporter_distribution.csv"]')
config=config.replace('Compute exact Bayesian phase uncertainty and read-error posteriors for a synthetic hexaploid chromosome.','Compute the full conditional conformational ensemble of eight reporters in an inhomogeneous cooperative receptor sheet.')
config=config.replace('keywords = ["polyploidy", "phase-uncertainty", "Bayesian-inference", "DNA-sequencing"]','keywords = ["statistical-biophysics", "allostery", "cooperativity", "joint-distribution"]')
config=config.replace('subcategory = "computational-genomics"','subcategory = "statistical-biophysics"')
(root/'task.toml').write_text(config)
(root/'solution/solve.sh').write_text('''#!/bin/bash
set -euo pipefail
python /solution/solve.py
''')
(root/'solution/solve.py').write_text('''import csv,json
from pathlib import Path
from engine import distribution

model=json.loads(Path('/root/data/model.json').read_text())
r=distribution(model)
out=Path('/root/results');out.mkdir(parents=True,exist_ok=True)
with (out/'thermodynamics.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['metric','value'])
 for key in ['log_partition','log_clamp_probability']:w.writerow([key,format(r[key],'.15g')])
with (out/'reporter_distribution.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['pattern','probability'])
 for i,p in enumerate(r['probabilities']):w.writerow([''.join(str((i>>j)&1) for j in range(8)),format(p,'.15g')])
''')
(root/'instruction.md').write_text('''A synthetic receptor sheet contains 1024 two-state subunits. Neighbor interactions
represent cooperative conformational preferences. Six experimentally controlled
subunits are clamped, and eight other subunits are read out simultaneously.
Compute the complete joint conformational law of the eight reporters under this
intervention, retaining the correlations mediated by the rest of the sheet.

The input `/root/data/model.json` defines the model completely:

- `coordinates`: an array of two-dimensional locations; the zero-based array
  index is the site ID. Locations describe the sheet geometry.
- `edges`: records `[i, j, J]`, each specifying one undirected interaction.
  Each listed interaction occurs once. There are no other edges, periodic
  connections, single-site fields, or unlisted interactions.
- `ligand_sites` and `ligand_spins`: the six site IDs and their imposed spins.
- `reporter_sites`: the eight site IDs, in the required output order.

Every subunit has spin $ s_i \\in \\{-1,+1\\} $. Couplings already include the
inverse temperature. For a complete spin assignment the dimensionless weight is
$ w(s)=\\exp(\\sum_{(i,j,J)\\in E} J s_i s_j) $.
The unconstrained partition function is $ Z=\\sum_s w(s) $, with the sum over
all assignments at all 1024 sites. This is a counting sum, not a uniform-prior
average. Let $ C $ denote the event that every ligand site has its specified
spin. Then $ P(C)=\\sum_{s\\in C} w(s)/Z $. The requested reporter law is
$ P(s_R=r\\mid C) $. Clamping conditions on the stated spins; it does not remove
the sites or their incident interactions. This is a finite synthetic equilibrium
model, not a molecular prediction for a specific natural receptor.

Write these two UTF-8 CSV files, with LF or CRLF line endings:

1. `/root/results/thermodynamics.csv`, with literal header `metric,value` and
   exactly two rows, keyed by `log_partition` and `log_clamp_probability`.
   Use natural logarithms of $ Z $ and $ P(C) $, respectively.
2. `/root/results/reporter_distribution.csv`, with literal header
   `pattern,probability` and exactly 256 rows, one for each eight-character
   binary string. Character position $ j $ refers to `reporter_sites[j]`;
   `0` means spin +1 and `1` means spin -1. For example, `10000000` specifies
   only the first reporter in state -1. Each probability is conditional on $ C $.

Rows may appear in any order. Include each key exactly once and no extra
columns or rows. Numeric fields must be finite ASCII decimal numbers, optionally
in scientific notation, without whitespace, grouping separators, or underscores.
Probabilities must lie in [0,1]. Supply every probability with absolute error at
most $ 5\\times10^{-8} $, their sum within $ 10^{-6} $ of 1, and each logarithm
with absolute error at most $ 2\\times10^{-6} $. These accuracy requirements concern
the full finite model, including unobserved spins and all listed edges.

The working directory is `/root`. Python, NumPy, SciPy and a C/C++ compiler are
available. You may use any computational method and install additional tools.
Only the two output CSVs need to be retained in `/root/results`.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
''')
print(root)
