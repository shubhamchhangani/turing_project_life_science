"""Build the nonplanar scaffold task, retaining the earlier frozen packages."""
import json,shutil
from pathlib import Path

source=Path('tasks/life-sciences/biology/frustrated-receptor-ensemble')
target=Path('tasks/life-sciences/biology/scaffold-receptor-ensemble')
author=Path('tools/receptor-scaffold')
if not target.exists():shutil.copytree(source,target)
for rel in ['environment/data/model.json','tests/data/model.json']:shutil.copy2(author/'draft/model.json',target/rel)
shutil.copy2(author/'draft/expectations.json',target/'tests/data/expectations.json')
for name in ['engine.py','signed_engine.py','scaffold_engine.py','compiled_moments.py']:shutil.copy2(author/name,target/'solution'/name)
solver=(source/'solution/solve.py').read_text().replace('from signed_engine import distribution','from scaffold_engine import distribution')
(target/'solution/solve.py').write_text(solver)
for name in ['engine.py','signed_engine.py','scaffold_engine.py','compiled_moments.py','source_kernel.py','independent_source.py','generate.py','check_small.py']:shutil.copy2(author/name,target/'tests/provenance'/name)
config=(source/'task.toml').read_text().replace('frustrated-receptor-ensemble','scaffold-receptor-ensemble').replace('receptor sheet with competing positive and negative interactions','receptor sheet with competing interactions, scaffold links and local biases')
(target/'task.toml').write_text(config)
instruction=(source/'instruction.md').read_text()
instruction=instruction.replace('Both signs are intentional.','Both signs are intentional. Additional effective long-range interactions\nrepresent scaffold-mediated coupling, while local fields bias individual\nconformations. The resulting interaction graph need not be planar.')
instruction=instruction.replace('  Each listed interaction occurs once. There are no other edges, periodic\n  connections, single-site fields, or unlisted interactions.','  Each listed interaction occurs once, regardless of its geometric length.\n  There are no other edges, periodic connections, or unlisted interactions.\n- `fields`: records `[i, h]` specifying the dimensionless single-site bias at\n  site i. Unlisted sites have zero bias; each listed site appears once.')
instruction=instruction.replace('\\sum_{(i,j,J)\\in E} J s_i s_j','\\sum_{(i,j,J)\\in E} J s_i s_j + \\sum_{(i,h)\\in H} h s_i')
instruction=instruction.replace('The unconstrained partition function','The two sums range over `edges` and `fields`, respectively.\nThe unconstrained partition function')
(target/'instruction.md').write_text(instruction)
for name in ['run_trials.py','check_regressions.py']:
    text=(Path('tools/receptor-frustration')/name).read_text().replace('frustrated-receptor','scaffold-receptor').replace('tools/receptor-frustration','tools/receptor-scaffold')
    (author/name).write_text(text)
r=json.loads((author/'draft/reference.json').read_text());e=json.loads((author/'draft/expectations.json').read_text())
record={'passed':True,'max_probability_difference':max(abs(a-b) for a,b in zip(r['probabilities'],e['probabilities'])),
        'log_partition_difference':abs(r['log_partition']-e['log_partition']),
        'log_clamp_difference':abs(r['log_clamp_probability']-e['log_clamp_probability']),
        'reference_seconds':r['seconds'],'source_seconds':e['seconds']}
assert record['max_probability_difference']<1e-9 and record['log_partition_difference']<1e-8 and record['log_clamp_difference']<1e-8
(author/'draft/independent-validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
