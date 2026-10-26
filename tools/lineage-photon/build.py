"""Build candidate only after independent full-instance and analytic checks."""
import json,shutil
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2];AUTHOR=Path(__file__).parent;D=AUTHOR/'draft'
TASK=ROOT/'tasks/life-sciences/biology/lineage-photon-inference'
OLD=ROOT/'tasks/life-sciences/biology/scaffold-receptor-ensemble'

if __name__=='__main__':
    assert json.loads((D/'small-validation.json').read_text())['passed']
    a=json.loads((D/'reference.json').read_text());b=json.loads((D/'independent.json').read_text())
    differences={key:float(np.max(np.abs(np.asarray(a[key])-np.asarray(b[key])))) for key in ['parent','daughter_rna','counts','count_overflow','log_evidence']}
    assert max(differences.values())<1e-11
    (D/'independent-validation.json').write_text(json.dumps(dict(passed=True,differences=differences,reference_seconds=a['seconds'],source_seconds=b['seconds']),indent=2)+'\n')
    assert not TASK.exists(),'Do not overwrite a frozen package'
    for name in ['environment/data','solution','tests/data','tests/provenance']:(TASK/name).mkdir(parents=True)
    for dest in ['environment/data/model.json','tests/data/model.json']:shutil.copy2(D/'model.json',TASK/dest)
    shutil.copy2(D/'independent.json',TASK/'tests/data/expectations.json')
    shutil.copy2(AUTHOR/'instruction.template.md',TASK/'instruction.md')
    shutil.copy2(AUTHOR/'verifier.template.py',TASK/'tests/test_outputs.py')
    for name in ['test.sh','Dockerfile']:shutil.copy2(OLD/'tests'/name,TASK/'tests'/name)
    shutil.copy2(OLD/'environment/Dockerfile',TASK/'environment/Dockerfile')
    for name in ['reference.py','write_outputs.py']:shutil.copy2(AUTHOR/name,TASK/'solution'/name)
    (TASK/'solution/solve.sh').write_text('#!/bin/bash\nset -euo pipefail\npython /solution/write_outputs.py\n')
    for name in ['generate.py','reference.py','independent.py','check_small.py']:shutil.copy2(AUTHOR/name,TASK/'tests/provenance'/name)
    meta=(OLD/'task.toml').read_text().replace('scaffold-receptor-ensemble','lineage-photon-inference')
    meta=meta.replace('artifacts = ["/root/results/thermodynamics.csv", "/root/results/reporter_distribution.csv"]','artifacts = ["/root/results/evidence.csv", "/root/results/division_posterior.csv", "/root/results/daughter_rna.csv", "/root/results/transcription_counts.csv"]')
    meta=meta.replace('Compute the full conditional conformational ensemble of eight reporters in an receptor sheet with competing interactions, scaffold links and local biases.','Infer hidden gene-expression states and transcription histories from integrated fluorescence across a cell division.')
    meta=meta.replace('["statistical-biophysics", "allostery", "cooperativity", "joint-distribution"]','["systems-biology", "stochastic-gene-expression", "cell-lineage", "bayesian-smoothing"]')
    meta=meta.replace('subcategory = "statistical-biophysics"','subcategory = "systems-biology"')
    (TASK/'task.toml').write_text(meta)
    (TASK/'README.md').write_text('''# Integrated fluorescence lineage inference

This is an original finite, synthetic stochastic gene-expression problem.
Capacities, bursts, feedback, observation timing and partition rules are explicit
model assumptions; no claim is made that the instance predicts a natural cell.

## Independent expectation source

`tests/provenance/generate.py` records the fixed model parameters.
`tests/provenance/independent.py` independently assembles a row generator and
uses positive uniformization recurrences for joint photon and transcription
counts. Conditional division messages enforce complementary inheritance.
`tests/data/expectations.json` comes from this computation, not an oracle run.
The genuine solution separately assembles a column generator and uses SciPy
sparse augmented matrix exponential actions. Both implement the same declared
stochastic model, with different numerical algorithms and separate assemblies.

Eight small instances agree between methods. Analytic limits with
state-independent photon emission and disabled transcription reproduce Poisson
observation masses, zero reaction counts and the unconditional mother prior.
These are cross-algorithm and analytic checks, not exhaustive path enumeration.

## Verification and tolerances

All mother states, all joint daughter RNA pairs, all event-count probabilities
including the tail, and the evidence are checked. Absolute tolerances are 2e-8
for every probability and 2e-7 for log evidence, with normalization within 1e-6.
They comfortably exceed independent numerical disagreement and decimal rounding
error. They are numerical accuracy tolerances, not biological noise estimates.
Contract validation precedes every accuracy check. The separate offline verifier
opens bounded regular output files without following symlinks and never executes
submitted code. Reference and provenance files are not copied to the agent image.

The reference does not fit hidden planted parameters: it computes the full
likelihood and posterior for the public finite model.

## Background

Finite chemical master equations: Munsky and Khammash (2006),
https://cnls.lanl.gov/~munsky/JChemPhys_124_044104.pdf
Sparse matrix exponential action implementation:
https://github.com/scipy/scipy/blob/main/scipy/sparse/linalg/_expm_multiply.py
These describe numerical background, not this task's answers.
''')
    for p in [TASK/'tests/test.sh',TASK/'solution/solve.sh']:p.chmod(0o755)
    worker=(ROOT/'tools/receptor-scaffold/run_trials.py').read_text().replace('scaffold-receptor-ensemble','lineage-photon-inference').replace('scaffold-receptor','lineage-photon')
    (AUTHOR/'run_trials.py').write_text(worker)
    print(json.dumps({'task':str(TASK),'differences':differences}),flush=True)
