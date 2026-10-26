"""Run the authorized subscription trials; package only two valid normal failures.

No paid API fallback, no quota failures counted, no ZIP for a solved candidate.
The frozen validated task is never edited by this worker.
"""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import select
import shutil
import subprocess
import time
import zipfile

ROOT=Path(__file__).resolve().parents[2]
AUTHOR=Path(__file__).resolve().parent
TASK=ROOT/'tasks/life-sciences/biology/spatial-ligand-cooperativity'
EVIDENCE=AUTHOR/'draft/twelve-fields'
REPORT=ROOT/'handover/spatial-ligand-progress.json'
JOBS=ROOT/'jobs/receptor-torus-frontier'


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


audit=module('task_preflight',ROOT/'tools/preflight_task.py').audit
classify=module('trial_decision',ROOT/'tools/metabolic-flux-analysis-v2/release_decision.py').classify


def capacity():
    binary=Path.home()/'.vscode/extensions/openai.chatgpt-26.1007.21434-darwin-arm64/bin/macos-aarch64/codex'
    p=subprocess.Popen([str(binary),'app-server','--stdio'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
    def send(value):p.stdin.write(json.dumps(value)+'\n');p.stdin.flush()
    try:
        send({'id':1,'method':'initialize','params':{'clientInfo':{'name':'benchmark_capacity_check','version':'1.0'},'capabilities':{'experimentalApi':True}}})
        until=time.monotonic()+25
        while time.monotonic()<until:
            if not select.select([p.stdout],[],[],1)[0]:continue
            line=p.stdout.readline()
            if not line:break
            event=json.loads(line)
            if event.get('id')==1:
                send({'method':'initialized'});send({'id':2,'method':'account/rateLimits/read','params':{}})
            elif event.get('id')==2:return event.get('result',{}).get('rateLimits')
        return None
    finally:p.terminate()


def result_for(job):
    files=list((JOBS/job).glob('*/result.json'))
    if len(files)!=1:return None
    raw=json.loads(files[0].read_text())
    return classify(files[0]) if raw.get('finished_at') else None


def stop_trial(job):
    for config in (JOBS/job).glob('*/config.json'):
        name=config.parent.name.lower()+'__env-main-1'
        if name.startswith('spatial-ligand-cooperativity__'):
            subprocess.run(['docker','stop','--timeout','10',name],capture_output=True,timeout=30)


def release_eligible(records,sleep_interruption=False):
    return (len(records)==2 and not sleep_interruption
        and all(r.get('status')=='substantive_attempt_failed' and r.get('error_type') is None
                and not r.get('quota_interrupted') and r.get('task_checksum') for r in records)
        and {(r.get('model'),r.get('effort')) for r in records}=={('openai/gpt-6-sol','max'),('anthropic/claude-opus-5-5','low')}
        and len({r['task_checksum'] for r in records})==1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gpt-not-before',required=True,type=float)
    args=parser.parse_args()
    lock=(EVIDENCE/'trial-worker.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    frozen=json.loads((EVIDENCE/'package-before-trial.json').read_text())
    state={'status':'VALIDATING_PREREQUISITES','worker_pid':os.getpid(),'zip_created':False,'zip_path':None,
           'package_fingerprint':frozen['package_fingerprint'],
           'gpt_not_before_utc':datetime.fromtimestamp(args.gpt_not_before,timezone.utc).isoformat(),
           'trials':{},'limitation':'Local trials do not guarantee a future Turing-harness result.'}
    def record(status):
        state['status']=status;state['updated_at']=datetime.now(timezone.utc).isoformat()
        temporary=REPORT.with_suffix('.tmp');temporary.write_text(json.dumps(state,indent=2)+'\n');temporary.replace(REPORT)
        print(state['updated_at'],status,flush=True)
    def unchanged():
        check=audit(TASK)
        assert check['static_checks_passed'] and check['package_fingerprint']==frozen['package_fingerprint'],'Frozen task changed or failed preflight'
    processes={};handles=[]
    jobs={'opus':'opus-5-5-low-spatial-ligand','gpt':'gpt6-sol-max-spatial-ligand-npx'}
    try:
        unchanged()
        assert json.loads((EVIDENCE/'offline-validation.json').read_text())['separate_offline_verifier_passed']
        for file in ['regressions.json','independent-validation.json','exhaustive-twelve.json']:
            assert json.loads((EVIDENCE/file).read_text())['all_passed'],file
        def launch(which):
            unchanged()
            assert not (JOBS/jobs[which]).exists(),'Trial name already exists; refusing a duplicate run'
            log=(AUTHOR/f'spatial-{which}-frontier.log').open('w');handles.append(log)
            if which=='opus':
                command=['python3',str(AUTHOR/'run_opus.py'),'--task',str(TASK),'--job-name',jobs[which]]
            else:
                command=['python3',str(ROOT/'tools/metabolic-flux-analysis-v2/run_frontier.py'),'--auth','chatgpt','--task',str(TASK),'--jobs-dir',str(JOBS),'--job-name',jobs[which]]
            processes[which]=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            state['trials'][which]={'status':'running','job':str(JOBS/jobs[which])}
            record(f'{which.upper()}_TRIAL_RUNNING')
        launch('opus')
        state['trials']['gpt']={'status':'queued_for_subscription_reset','job':str(JOBS/jobs['gpt'])}
        record('OPUS_RUNNING_GPT_QUEUED')
        gpt_checked=False;sleep_interruption=False
        previous_wall=time.time();previous_clock=time.monotonic();previous_active=True
        while True:
            now_wall=time.time();now_clock=time.monotonic()
            if previous_active and (now_wall-previous_wall)-(now_clock-previous_clock)>60:
                sleep_interruption=True;state['machine_sleep_detected']=True
            previous_wall,previous_clock=now_wall,now_clock
            for which,process in processes.items():
                job=JOBS/jobs[which]
                if job.is_dir() and not (job/'package-snapshot.json').exists():shutil.copy2(EVIDENCE/'package-before-trial.json',job/'package-snapshot.json')
                if process.poll() is not None and state['trials'][which]['status']=='running':
                    state['trials'][which]=result_for(jobs[which]) or {'status':'inconclusive_launch_failure','job':str(job)}
                    record(f'{which.upper()}_TRIAL_FINISHED')
            if any(r['status']=='solved_do_not_package' for r in state['trials'].values()):
                for which,process in processes.items():
                    if process.poll() is None:
                        stop_trial(jobs[which]);state['trials'][which]={'status':'cancelled_obsolete_candidate','job':str(JOBS/jobs[which])}
                record('DO_NOT_PACKAGE_SOLVED');break
            if not gpt_checked and time.time()>=args.gpt_not_before:
                gpt_checked=True;limits=capacity();state['gpt_capacity_at_launch']=limits
                if limits and (limits.get('primary') or {}).get('usedPercent',100)<=25:
                    launch('gpt')
                else:
                    state['trials']['gpt']={'status':'inconclusive_insufficient_subscription_capacity'}
                    record('HOLD_GPT_CAPACITY')
            if gpt_checked and all(process.poll() is not None for process in processes.values()):
                records=list(state['trials'].values())
                eligible=release_eligible(records,sleep_interruption)
                if not eligible:record('HOLD_INCOMPLETE_OR_INCONCLUSIVE_EVIDENCE');break
                unchanged()
                archive=Path.home()/'Downloads/spatial-ligand-cooperativity-v1.0.0.zip'
                if archive.exists():archive=archive.with_name('spatial-ligand-cooperativity-v1.0.0-'+datetime.now().strftime('%Y%m%d-%H%M%S')+'.zip')
                with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
                    for name in sorted(frozen['files']):z.write(TASK/name,name)
                check=audit(TASK,archive)
                if not check['static_checks_passed']:
                    archive.unlink();raise RuntimeError('Archive failed direct-root or byte-equality verification')
                check['zip_created']=True
                (EVIDENCE/'final-archive-audit.json').write_text(json.dumps(check,indent=2)+'\n')
                state.update(zip_created=True,zip_path=str(archive),zip_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
                record('ZIP_READY_TWO_LOCAL_MODEL_FAILURES');break
            previous_active=any(p.poll() is None for p in processes.values())
            time.sleep(20)
    except Exception as exc:
        state['error']=str(exc);record('HOLD_WORKER_ERROR');raise
    finally:
        for handle in handles:handle.close()
