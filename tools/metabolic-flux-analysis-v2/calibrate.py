"""Author-side noise sensitivity check; never included in the task archive."""
import json,subprocess,sys,csv,shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]/'tasks/life-sciences/biology/metabolic-flux-analysis-v2'
BASE=Path('/private/tmp/mfa-v2-calibration')

def one(seed):
    folder=BASE/str(seed);data=folder/'data';out=folder/'output';folder.mkdir(parents=True,exist_ok=True)
    shutil.copytree(ROOT/'environment/data',data,dirs_exist_ok=True)
    subprocess.run([sys.executable,str(HERE/'make_dataset.py'),str(seed),str(data)],check=True,capture_output=True)
    fit=subprocess.run([sys.executable,str(ROOT/'solution/dynamic_fit.py'),str(data),str(out),'2',str(seed)],check=True,capture_output=True,text=True)
    (folder/'fit.log').write_text(fit.stdout)
    with (out/'pools.csv').open() as f:pools={r['metabolite']:float(r['pool_size']) for r in csv.DictReader(f)}
    with (out/'fluxes.csv').open() as f:fluxes={r['reaction_id']:{k:float(v) for k,v in r.items() if k!='reaction_id'} for r in csv.DictReader(f)}
    return {'seed':seed,'pools':pools,'fluxes':fluxes,'log':fit.stdout}

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=2) as pool:
        results=[]
        for res in pool.map(one,range(1000,1006)):
            results.append(res);print(res['seed'],res['pools'],flush=True)
    (HERE/'calibration_results.json').write_text(json.dumps(results,indent=2)+'\n')
