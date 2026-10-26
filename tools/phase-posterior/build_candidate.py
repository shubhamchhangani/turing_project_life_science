"""Create a draft exact-posterior benchmark, never an automatic submission ZIP."""
import json,shutil
from pathlib import Path
ROOT=Path('tasks/life-sciences/biology/polyploid-phase-posterior')
for part in ['environment/data','tests/data','solution']:(ROOT/part).mkdir(parents=True,exist_ok=True)
source=Path('tools/phase-posterior/draft/model-8.json')
for target in ['environment/data/model.json','tests/data/model.json']:shutil.copy2(source,ROOT/target)
base=Path('tasks/life-sciences/biology/polyploid-haplotype-assembly/task.toml').read_text().replace('polyploid-haplotype-assembly','polyploid-phase-posterior').replace('artifacts = ["/root/results/haplotypes.csv"]','artifacts = ["/root/results/evidence.csv", "/root/results/read_posteriors.csv"]')
base=base.replace('Assemble sixteen homologous chromosome phases from noisy overlapping DNA reads and known allele dosages.','Compute exact Bayesian phase uncertainty and read-error posteriors for a synthetic hexaploid chromosome.').replace('keywords = ["polyploidy", "haplotype-assembly", "DNA-sequencing", "weighted-MEC"]','keywords = ["polyploidy", "phase-uncertainty", "Bayesian-inference", "DNA-sequencing"]')
(ROOT/'task.toml').write_text(base)
for name in ['environment/Dockerfile','tests/Dockerfile','tests/test.sh','solution/solve.sh']:
 shutil.copy2(Path('tasks/life-sciences/biology/polyploid-haplotype-assembly')/name,ROOT/name)
for name in ['orbit_graph.cpp','infer.cpp']:shutil.copy2(Path('tools/phase-posterior')/name,ROOT/'solution'/name)
# Build a table-construction helper containing no generation seed or chromosomes.
import ast
nodes=ast.parse(Path('tools/phase-posterior/generate.py').read_text()).body
fn=next(x for x in nodes if isinstance(x,ast.FunctionDef) and x.name=='tables')
(ROOT/'solution/emissions.py').write_text('import numpy as np\n\n'+ast.unparse(fn)+'\n')
print(ROOT)
