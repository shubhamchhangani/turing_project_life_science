"""Assemble original receptor task; author evidence stays outside the package."""
import ast
import csv
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
TASK=ROOT/'tasks/life-sciences/biology/receptor-cluster-thermodynamics'
AUTHOR=Path(__file__).resolve().parent
for relative in ['environment/data','solution','tests/data']:
    (TASK/relative).mkdir(parents=True,exist_ok=True)

width,height=17,19
rng=np.random.default_rng(20261026)
model={'width':width,'height':height,'probes':[{'id':f'R{i+1}','x':p[0],'y':p[1]} for i,p in enumerate([(0,0),(16,18),(3,13),(10,2),(8,9),(15,7)])]}
for folder in ['environment/data','tests/data']:
    (TASK/folder/'model.json').write_text(json.dumps(model,indent=2)+'\n')
    with (TASK/folder/'couplings.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['x','y','direction','K'])
        # Reset RNG so the independently mounted verifier has identical inputs.
        local=np.random.default_rng(20261026)
        for y in range(height):
            for x in range(width):
                for direction in ['E','N']:
                    k=.49+.125*np.sin((x+2*y)*.53)+local.uniform(-.10,.10)
                    if (x*13+y*7+(direction=='N')*19)%97==11:
                        k=-local.uniform(.10,.22)
                    w.writerow([x,y,direction,f'{k:.9f}'])

reference=ast.parse((AUTHOR/'prototype.py').read_text())
reference.body=[n for n in reference.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in {'pfaffian_log','Torus'}]
source=ast.unparse(reference)+'\n'
source+='''
def solve(data_dir, output_dir):
    import csv, json
    from pathlib import Path
    model=json.loads((data_dir/'model.json').read_text())
    simulator=Torus(model['width'],model['height'])
    rows=list(csv.DictReader((data_dir/'couplings.csv').open()))
    coupling_map={(int(r['x']),int(r['y']),r['direction']):float(r['K']) for r in rows}
    k=np.array([coupling_map[p[0],p[1],'E' if d==0 else 'N'] for p,q,d in simulator.edges])
    sites=[(r['x'],r['y']) for r in model['probes']]
    n=len(sites)
    masks=[m for m in range(1,1<<n) if m.bit_count()%2==0]
    queries=[[p for i,p in enumerate(sites) if m>>i&1] for m in masks]
    logz,values=simulator.evaluate(k,queries)
    moments={0:1.,**dict(zip(masks,values))}
    output_dir.mkdir(parents=True,exist_ok=True)
    with (output_dir/'joint.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\\n');writer.writerow(['pattern','probability'])
        for bits in range(1<<n):
            pattern=format(bits,f'0{n}b')
            spins=[1 if c=='1' else -1 for c in pattern]
            probability=sum(value*np.prod([spins[i] for i in range(n) if mask>>i&1]) for mask,value in moments.items())/(1<<n)
            writer.writerow([pattern,f'{probability:.12f}'])
    with (output_dir/'thermodynamics.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\\n');writer.writerow(['quantity','value'])
        writer.writerow(['log_partition',f'{logz:.12f}'])
    print(json.dumps({'log_partition':logz,'moments_computed':len(moments)},indent=2),flush=True)

if __name__=='__main__':
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path('/root/data'))
    parser.add_argument('--output',type=Path,default=Path('/root/results'))
    args=parser.parse_args()
    solve(args.data,args.output)
'''
(TASK/'solution/solve.py').write_text(source)
(TASK/'solution/solve.sh').write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\npython3 /solution/solve.py\n')
independent=ast.parse((AUTHOR/'independent.py').read_text())
independent.body=[n for n in independent.body if not isinstance(n,ast.If)]
(TASK/'tests/physics.py').write_text(ast.unparse(independent)+'\n')
print(TASK)
