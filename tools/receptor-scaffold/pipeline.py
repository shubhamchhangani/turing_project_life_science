"""Complete reference/source validation, then strictly ordered subscription trials."""
import fcntl,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2];AUTHOR=Path(__file__).resolve().parent;DRAFT=AUTHOR/'draft'
TASK=ROOT/'tasks/life-sciences/biology/scaffold-receptor-ensemble'
REPORT=ROOT/'handover/scaffold-build-progress.json'
PYTHON=ROOT/'tools/.venv-transient/bin/python'
OUTPUT=Path('/private/tmp/scaffold-receptor-reference');LOGS=Path('/private/tmp/scaffold-receptor-verifier')


if __name__=='__main__':
    lock=(DRAFT/'pipeline.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state={'status':'STARTING','worker_pid':os.getpid(),'zip_created':False,
           'model_sha256':hashlib.sha256((DRAFT/'model.json').read_bytes()).hexdigest()}
    processes=[];handles=[]
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def record(status):
        state['status']=status;state['updated_at']=datetime.now(timezone.utc).isoformat()
        temp=REPORT.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(REPORT)
        print(state['updated_at'],status,flush=True)
    def run(args,log):
        with (DRAFT/log).open('w') as f:subprocess.run(list(map(str,args)),cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    try:
        assert json.loads((DRAFT/'exhaustive-validation.json').read_text())['passed']
        OUTPUT.mkdir(exist_ok=True);LOGS.mkdir(exist_ok=True)
        assert not (OUTPUT/'reference.json').exists(),'Refusing to reuse a stale reference result'
        reference=['docker','run','--rm','--name','scaffold-reference-author','--network','none','--cpus','1','--memory','2g',
                   '-v',str(AUTHOR)+':/solution:ro','-v',str(DRAFT/'model.json')+':/root/data/model.json:ro',
                   '-v',str(OUTPUT)+':/root/results','frustrated-receptor-agent:latest',
                   'timeout','--kill-after=10s','18000','python','/solution/run_reference.py']
        source=[str(PYTHON),str(AUTHOR/'independent_source.py')]
        for command,name in [(reference,'reference-container.log'),(source,'independent.log')]:
            log=(DRAFT/name).open('w');handles.append(log)
            processes.append(subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT))
        started=time.monotonic();record('REFERENCE_AND_INDEPENDENT_SOURCE_RUNNING')
        while any(p.poll() is None for p in processes):
            if time.monotonic()-started>18100:raise TimeoutError('Scientific validation exceeded the task time budget')
            for p in processes:
                if p.poll() is not None and p.returncode!=0:raise RuntimeError('A scientific validation process failed; inspect its log')
            state['elapsed_seconds']=round(time.monotonic()-started)
            state['reference_finished']=processes[0].poll() is not None
            state['independent_source_finished']=processes[1].poll() is not None
            record('REFERENCE_AND_INDEPENDENT_SOURCE_RUNNING');time.sleep(30)
        assert all(p.returncode==0 for p in processes)
        state.update(reference_finished=True, independent_source_finished=True,
                     elapsed_seconds=round(time.monotonic()-started))
        shutil.copy2(OUTPUT/'reference.json',DRAFT/'reference.json')
        run([PYTHON,AUTHOR/'build_candidate.py'],'build-candidate.log')
        validation=json.loads((DRAFT/'independent-validation.json').read_text());assert validation['passed']
        assert validation['reference_seconds']<18000
        run([PYTHON,AUTHOR/'check_regressions.py'],'regressions.log')
        run(['docker','build','-t','scaffold-receptor-verifier:latest',TASK/'tests'],'verifier-build.log')
        run(['docker','run','--rm','--network','none','--cpus','1','--memory','2g','-v',str(OUTPUT)+':/root/results:ro',
             '-v',str(LOGS)+':/logs/verifier','scaffold-receptor-verifier:latest','bash','/tests/test.sh'],'verifier-container.log')
        assert (LOGS/'reward.txt').read_text().strip()=='1'
        (DRAFT/'offline-validation.json').write_text(json.dumps({'passed':True,'separate_offline_verifier_passed':True,
             'tests_passed':4,'cpus':1,'memory_mb':2048,'network':'none','reference_seconds':validation['reference_seconds']},indent=2)+'\n')
        readme=(AUTHOR/'README.template.md').read_text()
        readme+='\n## Completed full-instance validation\n\n'
        readme+=f"Reference: {validation['reference_seconds']:.2f} seconds in a one-CPU, 2 GiB offline container; all four separate-verifier tests passed.\n\n"
        readme+=f"Maximum probability difference between implementations: {validation['max_probability_difference']:.3g}. Log-partition difference: {validation['log_partition_difference']:.3g}; log-clamp difference: {validation['log_clamp_difference']:.3g}.\n"
        (TASK/'README.md').write_text(readme)
        run(['python3',ROOT/'tools/preflight_task.py',TASK,'--report',DRAFT/'package-before-trial.json'],'preflight.log')
        assert hashlib.sha256((TASK/'environment/data/model.json').read_bytes()).hexdigest()==state['model_sha256']
        record('VALIDATED_STARTING_CLAUDE_FIRST_TRIALS')
        run(['python3',AUTHOR/'run_trials.py'],'trial-worker.log')
        trial=json.loads((ROOT/'handover/scaffold-receptor-progress.json').read_text())
        state['trial_status']=trial['status'];state['zip_created']=trial.get('zip_created',False);state['zip_path']=trial.get('zip_path')
        record('TRIAL_WORKER_FINISHED')
    except Exception as exc:
        state['error']=str(exc);record('HOLD_PIPELINE_ERROR')
        for p in processes:
            if p.poll() is None:p.terminate()
        subprocess.run(['docker','stop','--timeout','5','scaffold-reference-author'],capture_output=True)
        raise
    finally:
        for h in handles:h.close()
        state['worker_finished']=True;record(state['status'])
