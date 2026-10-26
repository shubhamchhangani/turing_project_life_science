"""Run one GPT-6-sol max-effort evaluation with an explicit auth mode.

Credentials remain in the subprocess environment and are never printed.
"""
import os
import argparse
from datetime import datetime, timezone
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--job-name', default='gpt6-sol-max-v2-' + datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S'))
parser.add_argument('--auth', choices=['api', 'chatgpt'], default='api')
parser.add_argument('--task', default='tasks/life-sciences/biology/metabolic-flux-analysis-v2')
parser.add_argument('--jobs-dir', default='jobs/mfa-v2-frontier')
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
env = os.environ.copy()
if args.auth == 'api':
    for line in (root / '.env').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#') and '=' in line:
            key, value = line.split('=', 1)
            env[key.strip()] = value.strip().strip('"').strip("'")
    env.pop('CODEX_FORCE_AUTH_JSON', None)
    env.pop('CODEX_AUTH_JSON_PATH', None)
    agent = 'codex'
else:
    for key in ('OPENAI_API_KEY', 'CODEX_API_KEY', 'OPENAI_BASE_URL', 'OPENAI_API_BASE', 'CODEX_AUTH_JSON_PATH'):
        env.pop(key, None)
    env['CODEX_FORCE_AUTH_JSON'] = '1'
    env['PYTHONPATH'] = str(Path(__file__).resolve().parent)
    agent = 'npx_codex_agent:NpxCodex'
raise SystemExit(subprocess.call([
    'harbor', 'run', '-p', args.task,
    '-a', agent, '-m', 'openai/gpt-6-sol', '--effort', 'max', '-k', '1',
    '--ak', 'version=0.164.0',
    '-o', args.jobs_dir, '--job-name', args.job_name,
], cwd=root, env=env))
