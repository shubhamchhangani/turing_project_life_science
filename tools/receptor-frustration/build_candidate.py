"""Create a distinct frozen candidate without modifying the previous trial package."""
import json,shutil
from pathlib import Path

source=Path('tasks/life-sciences/biology/conditional-receptor-ensemble')
target=Path('tasks/life-sciences/biology/frustrated-receptor-ensemble')
author=Path('tools/receptor-frustration')
if not target.exists():shutil.copytree(source,target)
for rel in ['environment/data/model.json','tests/data/model.json']:shutil.copy2(author/'draft/model.json',target/rel)
shutil.copy2(author/'draft/expectations.json',target/'tests/data/expectations.json')
for name in ['engine.py','signed_engine.py']:shutil.copy2(author/name,target/'solution'/name)
solver=(source/'solution/solve.py').read_text().replace('from engine import distribution','from signed_engine import distribution')
(target/'solution/solve.py').write_text(solver)
for name in ['engine.py','signed_engine.py','independent_source.py','generate.py','check_small.py']:shutil.copy2(author/name,target/'tests/provenance'/name)
config=(source/'task.toml').read_text().replace('conditional-receptor-ensemble','frustrated-receptor-ensemble').replace('inhomogeneous cooperative receptor sheet','receptor sheet with competing positive and negative interactions')
(target/'task.toml').write_text(config)
instruction=(source/'instruction.md').read_text().replace('represent cooperative conformational preferences.','represent competing conformational preferences: positive couplings favor equal\nspins and negative couplings favor opposite spins. Both signs are intentional.')
(target/'instruction.md').write_text(instruction)
for name in ['run_trials.py','check_regressions.py']:
    text=(Path('tools/receptor-joint')/name).read_text().replace('conditional-receptor','frustrated-receptor').replace('tools/receptor-joint','tools/receptor-frustration')
    (author/name).write_text(text)
r=json.loads((author/'draft/reference.json').read_text());e=json.loads((author/'draft/expectations.json').read_text())
record={'passed':True,'max_probability_difference':max(abs(a-b) for a,b in zip(r['probabilities'],e['probabilities'])),
        'log_partition_difference':abs(r['log_partition']-e['log_partition']),
        'log_clamp_difference':abs(r['log_clamp_probability']-e['log_clamp_probability']),
        'reference_seconds':r['seconds'],'source_seconds':e['seconds']}
assert record['max_probability_difference']<1e-10 and record['log_partition_difference']<1e-8 and record['log_clamp_difference']<1e-9
(author/'draft/independent-validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
