"""Strictly sequential subscription trials; only two normal failures permit ZIP."""
import fcntl,hashlib,importlib.util,json,os,subprocess,time,zipfile
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
AUTHOR=Path(__file__).resolve().parent
TASK=ROOT/'tasks/life-sciences/biology/frustrated-receptor-ensemble'
EVIDENCE=AUTHOR/'draft'
JOBS=ROOT/'jobs/receptor-torus-frontier'
REPORT=ROOT/'handover/frustrated-receptor-progress.json'


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


audit=module('preflight',ROOT/'tools/preflight_task.py').audit
classify=module('decision',ROOT/'tools/metabolic-flux-analysis-v2/release_decision.py').classify
capacity=module('old_worker_helpers',ROOT/'tools/phase-posterior/run_trials.py').capacity


def permitted_failure(r):
    return (r.get('status')=='substantive_attempt_failed' and r.get('error_type') is None
            and not r.get('quota_interrupted') and bool(r.get('task_checksum')))


if __name__=='__main__':
    lock=(EVIDENCE/'worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    frozen=json.loads((EVIDENCE/'package-before-trial.json').read_text())
    state={'status':'VALIDATING','worker_pid':os.getpid(),'zip_created':False,'zip_path':None,
           'order':'Claude first; GPT only after a valid Claude failure',
           'package_fingerprint':frozen['package_fingerprint'],'trials':{}}
    def save(status):
        state['status']=status;state['updated_at']=datetime.now(timezone.utc).isoformat()
        temp=REPORT.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(REPORT)
        print(state['updated_at'],status,flush=True)
    def unchanged():
        check=audit(TASK)
        assert check['static_checks_passed'] and check['package_fingerprint']==frozen['package_fingerprint']
    try:
        unchanged()
        for file in ['offline-validation.json','regressions.json','exhaustive-validation.json','independent-validation.json']:
            assert json.loads((EVIDENCE/file).read_text())['passed'],file
        for which,job in [('opus','opus-5-5-low-frustrated-receptor'),('gpt','gpt6-sol-max-frustrated-receptor-npx')]:
            if which=='gpt':
                assert permitted_failure(state['trials']['opus'])
                limits=capacity();state['gpt_capacity']=limits
                if not limits or (limits.get('primary') or {}).get('usedPercent',100)>40:
                    save('HOLD_GPT_SUBSCRIPTION_CAPACITY');break
            unchanged();assert not (JOBS/job).exists(),'Refusing duplicate trial'
            if which=='opus':
                command=['python3',str(ROOT/'tools/receptor-torus/run_opus.py'),'--task',str(TASK),'--job-name',job]
            else:
                command=['python3',str(ROOT/'tools/metabolic-flux-analysis-v2/run_frontier.py'),'--auth','chatgpt','--task',str(TASK),'--jobs-dir',str(JOBS),'--job-name',job]
            with (EVIDENCE/f'{which}-trial.log').open('w') as log:
                process=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                state['trials'][which]={'status':'running','job':str(JOBS/job)};save(which.upper()+'_RUNNING')
                wall=time.time();clock=time.monotonic();slept=False
                while process.poll() is None:
                    time.sleep(20)
                    w,c=time.time(),time.monotonic()
                    if (w-wall)-(c-clock)>60:slept=True
                    wall,clock=w,c
            files=list((JOBS/job).glob('*/result.json'))
            result=classify(files[0]) if len(files)==1 else {'status':'inconclusive_launch_failure'}
            state['trials'][which]=result
            if (JOBS/job).exists():(JOBS/job/'package-snapshot.json').write_text(json.dumps(frozen,indent=2)+'\n')
            if result.get('status')=='solved_do_not_package':save('DO_NOT_PACKAGE_SOLVED');break
            if slept or not permitted_failure(result):
                state['machine_sleep_detected']=slept;save('HOLD_INCONCLUSIVE_TRIAL');break
            save(which.upper()+'_VALID_FAILURE')
        else:
            unchanged();a,b=state['trials']['opus'],state['trials']['gpt']
            assert permitted_failure(a) and permitted_failure(b)
            assert a['task_checksum']==b['task_checksum']
            assert {(r.get('model'),r.get('effort')) for r in [a,b]}=={('anthropic/claude-opus-5-5','low'),('openai/gpt-6-sol','max')}
            target=Path.home()/'Downloads/frustrated-receptor-ensemble-v1.0.0.zip'
            if target.exists():target=target.with_name(target.stem+'-'+datetime.now().strftime('%Y%m%d-%H%M%S')+'.zip')
            with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
                for name in sorted(frozen['files']):z.write(TASK/name,name)
            check=audit(TASK,target)
            if not check['static_checks_passed']:
                target.unlink();raise RuntimeError('Archive validation failed')
            (EVIDENCE/'archive-audit.json').write_text(json.dumps(check,indent=2)+'\n')
            state.update(zip_created=True,zip_path=str(target),zip_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
            save('ZIP_READY_TWO_LOCAL_FAILURES')
    except Exception as exc:
        state['error']=str(exc);save('HOLD_WORKER_ERROR');raise
    finally:
        state['worker_finished']=True
        save(state['status'])
