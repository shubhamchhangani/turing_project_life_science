"""Finish offline validation after the already-running packaged reference exits."""
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
AUTHOR=Path(__file__).resolve().parent
deadline=time.monotonic()+18000
while True:
    status=subprocess.run(['docker','inspect','--format','{{.State.Running}}','spatial-ligand-reference'],capture_output=True,text=True)
    if status.returncode or status.stdout.strip()!='true':break
    if time.monotonic()>deadline:raise SystemExit('Reference exceeded validation wait allowance.')
    time.sleep(20)
reference_log=(AUTHOR/'spatial-offline-reference.log').read_text()
if 'Fit completed' not in reference_log:raise SystemExit('Reference did not finish successfully; verifier not launched.')
with (AUTHOR/'spatial-offline-verifier.log').open('w') as stream:
    result=subprocess.run(['docker','run','--rm','--network','none','--cpus','1','--memory','2g',
        '-v','/private/tmp/spatial-ligand-reference:/root/results:ro',
        '-v','/private/tmp/spatial-ligand-verifier:/logs/verifier',
        'spatial-ligand-verifier:latest','bash','/tests/test.sh'],stdout=stream,stderr=subprocess.STDOUT,timeout=900)
folder=AUTHOR/'draft/twelve-fields'
shutil.copytree('/private/tmp/spatial-ligand-verifier',folder/'offline-verifier',dirs_exist_ok=True)
generation=json.loads((folder/'generation.json').read_text())
with Path('/private/tmp/spatial-ligand-reference/couplings.csv').open() as stream:
    edges={(int(r['x']),int(r['y']),r['direction']):float(r['K']) for r in csv.DictReader(stream)}
with Path('/private/tmp/spatial-ligand-reference/biases.csv').open() as stream:
    fields={r['receptor_id']:float(r['h']) for r in csv.DictReader(stream)}
parameters=[edges[x,y,d] for y in range(19) for x in range(17) for d in ['E','N']]+[fields[f'F{i+1}'] for i in range(12)]
errors=[abs(a-b) for a,b in zip(parameters,generation['truth'])]
reward=(folder/'offline-verifier/reward.txt').read_text().strip()
final=reference_log.splitlines()[-1]
report={'reference_completed':True,'reference_summary':final,'reference_network':'none','cpus':1,'memory_mb':2048,
    'max_coupling_error':max(errors[:646]),'max_bias_error':max(errors[646:]),'verifier_reward':reward,
    'separate_offline_verifier_passed':result.returncode==0 and reward=='1','zip_created':False}
(folder/'offline-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
raise SystemExit(0 if report['separate_offline_verifier_passed'] else 1)
