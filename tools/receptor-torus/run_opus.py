"""One isolated Opus 5.5 low-effort trial through the user's Claude subscription.

No API-key billing, credential printing, or credential file export. The native
Claude Code CLI receives an existing OAuth access token in its environment.
The host Keychain lookup is limited to Claude Code's own default service.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[2]


def subscription_token():
    existing=os.environ.get('CLAUDE_CODE_OAUTH_TOKEN')
    if existing:
        return existing
    cli=[shutil.which('claude')] if shutil.which('claude') else ['npx','--yes','@anthropic-ai/claude-code@2.1.300']
    status=subprocess.run(cli+['auth','status'],capture_output=True,text=True)
    try:
        info=json.loads(status.stdout)
    except ValueError:
        raise SystemExit('Claude authentication status could not be read.')
    if not info.get('loggedIn'):
        raise SystemExit('Claude is not logged in. Complete claude auth login --claudeai first.')
    # Never fall back to an API key or a different provider.
    if info.get('authMethod') in {'api_key','apiKey','anthropic_profile'}:
        raise SystemExit('A Claude subscription login is required for this run; API billing is not enabled here.')
    credential_path=Path.home()/'.claude/.credentials.json'
    raw=None
    # Native macOS Claude Code uses Keychain. A legacy plaintext credential
    # file may still exist but must not shadow the currently signed-in account.
    if sys.platform=='darwin':
        result=subprocess.run(['security','find-generic-password','-s','Claude Code-credentials','-w'],capture_output=True,text=True)
        if result.returncode==0:
            raw=result.stdout
    if raw is None and credential_path.is_file():
        raw=credential_path.read_text()
    if raw is None:
        raise SystemExit('The signed-in Claude subscription credential is not available to this runner.')
    try:
        oauth=json.loads(raw).get('claudeAiOauth',{})
    except ValueError:
        raise SystemExit('Claude credential storage could not be decoded; no credential data was printed.')
    del raw
    token=oauth.get('accessToken')
    if not token:
        raise SystemExit('No Claude subscription access token was found; no API fallback was attempted.')
    expires=oauth.get('expiresAt')
    if isinstance(expires,(int,float)) and expires/1000<=time.time()+300:
        raise SystemExit('Claude login needs refreshing before this trial.')
    return token


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-name',default='opus-5-5-low-receptor-torus')
    parser.add_argument('--task',default='tasks/life-sciences/biology/inverse-receptor-cooperativity')
    args=parser.parse_args()
    token=subscription_token()
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('ANTHROPIC_') or key in {'CLAUDE_CODE_USE_BEDROCK','CLAUDE_CODE_USE_VERTEX','CLAUDE_CODE_USE_FOUNDRY','CLAUDE_CONFIG_DIR','AWS_BEARER_TOKEN_BEDROCK'}:
            env.pop(key,None)
    env['CLAUDE_CODE_OAUTH_TOKEN']=token
    env['CLAUDE_FORCE_OAUTH']='1'
    del token
    raise SystemExit(subprocess.call([
        'harbor','run','-p',args.task,
        '-a','claude-code','-m','anthropic/claude-opus-5-5','--effort','low','-k','1',
        '--ak','version=2.1.300','-o','jobs/receptor-torus-frontier','--job-name',args.job_name,
    ],cwd=ROOT,env=env))
