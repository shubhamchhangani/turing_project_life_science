"""Assemble the twelve-field candidate without changing any active trial."""
import ast
import json
from pathlib import Path
import shutil

AUTHOR=Path(__file__).resolve().parent
ROOT=AUTHOR.parents[1]
OLD=ROOT/'tasks/life-sciences/biology/inverse-receptor-cooperativity'
TASK=ROOT/'tasks/life-sciences/biology/spatial-ligand-cooperativity'
g=json.loads((AUTHOR/'draft/twelve-fields/generation.json').read_text())
for folder in ['environment/data','tests/data','solution']:(TASK/folder).mkdir(parents=True,exist_ok=True)
for source in (AUTHOR/'draft/twelve-fields/public').iterdir():
    for folder in ['environment/data','tests/data']:shutil.copy2(source,TASK/folder/source.name)
(TASK/'tests/data/injected_parameters.json').write_text(json.dumps({'source':'Original parameter draw, seed 20261030; not reference output.','couplings':g['truth'][:646],'biases':g['truth'][646:]},indent=2)+'\n')
shutil.copy2(OLD/'solution/lattice.py',TASK/'solution/lattice.py')

tree=ast.parse((AUTHOR/'sparse_fields_prototype.py').read_text())
tree.body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.ClassDef))]
for node in tree.body:
    if isinstance(node,ast.ImportFrom) and node.module=='prototype':node.module='lattice'
(TASK/'solution/sparse_physics.py').write_text(ast.unparse(tree)+'\n')
fit=(AUTHOR/'fit_sparse_fields.py').read_text().replace('from sparse_fields_prototype import SparseFields','from sparse_physics import SparseFields')
(TASK/'solution/solve.py').write_text(fit)
(TASK/'solution/solve.sh').write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\npython3 /solution/solve.py\n')
tree=ast.parse((AUTHOR/'independent.py').read_text())
tree.body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.ClassDef))]
for node in tree.body:
    if isinstance(node,ast.ClassDef):node.body=[n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='__init__']
independent=ast.parse((AUTHOR/'sparse_fields_independent.py').read_text())
independent.body=[n for n in independent.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.ClassDef)) and not(isinstance(n,ast.ImportFrom) and n.module=='independent')]
(TASK/'tests/physics.py').write_text(ast.unparse(tree)+'\n'+ast.unparse(independent)+'\n')
tests=(OLD/'tests/test_outputs.py').read_text().replace('IndependentFieldTorus','IndependentSparseFields')
tests=tests.replace('len(submitted) == 649','len(submitted) == len(EDGES) + len(SITES)').replace('submitted[:-3]','submitted[:len(EDGES)]').replace('submitted[-3:]','submitted[len(EDGES):]').replace('All three local ligand biases','All twelve local ligand biases')
(TASK/'tests/test_outputs.py').write_text(tests)
for filename in ['Dockerfile','test.sh']:shutil.copy2(OLD/'tests'/filename,TASK/'tests'/filename)
env=(OLD/'tests/Dockerfile').read_text().replace(' pytest==8.4.1 pytest-json-ctrf==0.3.5','').replace('RUN mkdir -p /root/results /logs/verifier\nCOPY . /tests/','WORKDIR /root\nCOPY data /root/data\nRUN mkdir -p /root/results')
(TASK/'environment/Dockerfile').write_text(env)
toml=(OLD/'task.toml').read_text().replace('inverse-receptor-cooperativity','spatial-ligand-cooperativity').replace('three local ligand biases','twelve spatially separated ligand biases')
(TASK/'task.toml').write_text(toml)
instruction=(OLD/'instruction.md').read_text().replace('Three specified receptors','Twelve specified receptors').replace('the three biases','the twelve biases').replace('a=1}^{3}','a=1}^{12}').replace('the three biased','the twelve biased').replace('the three input receptor IDs','the twelve input receptor IDs').replace('Python, NumPy and SciPy are installed;','Python, NumPy, SciPy and PFAPACK are installed;')
(TASK/'instruction.md').write_text(instruction)
print(TASK)
