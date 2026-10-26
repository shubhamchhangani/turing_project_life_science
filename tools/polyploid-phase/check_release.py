"""Verify the release conditions independently of running subscription jobs."""
import copy,json
from pathlib import Path
from run_trials import release_eligible
base=[dict(status='substantive_attempt_failed',error_type=None,quota_interrupted=False,task_checksum='same',model=m,effort=e) for m,e in [('openai/gpt-6-sol','max'),('anthropic/claude-opus-5-5','low')]]
cases=[]
assert release_eligible(base);cases.append('two_valid_failures')
assert not release_eligible(base,True);cases.append('sleep_interruption')
assert not release_eligible(base[:1]);cases.append('missing_model')
for key,value in [('status','solved_do_not_package'),('status','incomplete'),('status','inconclusive_infrastructure_or_no_attempt'),('error_type','AgentTimeoutError'),('quota_interrupted',True),('task_checksum','other'),('task_checksum',None),('model','anthropic/other'),('effort','high')]:
 changed=copy.deepcopy(base);changed[1][key]=value;assert not release_eligible(changed);cases.append(f'{key}={value}')
Path('tools/polyploid-phase/draft/release-checks.json').write_text(json.dumps(dict(passed=True,cases=cases),indent=2)+'\n')
print(f'{len(cases)} release checks passed')
