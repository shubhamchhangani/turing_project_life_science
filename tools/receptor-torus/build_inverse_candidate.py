"""Build a new inverse biophysics task from independently planted parameters."""
import ast
import csv
import json
from pathlib import Path

AUTHOR=Path(__file__).resolve().parent
ROOT=AUTHOR.parents[1]
TASK=ROOT/'tasks/life-sciences/biology/inverse-receptor-cooperativity'
generation=json.loads((AUTHOR/'draft/field-inverse/generation.json').read_text())
for folder in ['environment/data','tests/data','solution']:(TASK/folder).mkdir(parents=True,exist_ok=True)
model={'width':generation['width'],'height':generation['height'],'field_sites':[{'id':f'F{i+1}','x':p[0],'y':p[1]} for i,p in enumerate(generation['field_sites'])],'coupling_bounds':[-.35,.9],'field_bounds':[-.8,.8]}
for folder in ['environment/data','tests/data']:
    (TASK/folder/'model.json').write_text(json.dumps(model,indent=2)+'\n')
    with (TASK/folder/'edge_moments.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['x','y','direction','mean_product'])
        index=0
        for y in range(model['height']):
            for x in range(model['width']):
                for direction in ['E','N']:
                    writer.writerow([x,y,direction,f'{generation["observables"][index]:.12f}']);index+=1
    with (TASK/folder/'site_moments.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['receptor_id','mean_state'])
        for site,mean in zip(model['field_sites'],generation['observables'][-3:]):writer.writerow([site['id'],f'{mean:.12f}'])
(TASK/'tests/data/injected_parameters.json').write_text(json.dumps({'source':'Original parameter draw, seed 20261028; not fitted reference output.','couplings':generation['truth'][:-3],'biases':generation['truth'][-3:]},indent=2)+'\n')

tree=ast.parse((AUTHOR/'prototype.py').read_text())
tree.body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in {'pfaffian_log','Torus'}]
for node in tree.body:
    if isinstance(node,ast.ClassDef):node.body=[n for n in node.body if isinstance(n,ast.FunctionDef) and n.name in {'__init__','dual_weights','disorder'}]
(TASK/'solution/lattice.py').write_text(ast.unparse(tree)+'\n')
field_tree=ast.parse((AUTHOR/'field_inverse_prototype.py').read_text())
field_tree.body=[n for n in field_tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,ast.ClassDef)]
for node in field_tree.body:
    if isinstance(node,ast.ImportFrom) and node.module=='prototype':node.module='lattice'
source=ast.unparse(field_tree)+'\n'
newton=ast.parse((AUTHOR/'newton_fit.py').read_text())
newton.body=[n for n in newton.body if isinstance(n,ast.FunctionDef) and n.name=='fit']
newton.body[0].name='newton_fit'
source+='from scipy.linalg import solve\n'+ast.unparse(newton)+'\n'
source+='''
def fit(data,output):
    model=json.loads((data/'model.json').read_text())
    sites=[(p['x'],p['y']) for p in model['field_sites']]
    physics=FieldTorus(model['width'],model['height'],sites)
    with (data/'edge_moments.csv').open() as stream:
        edges={(int(r['x']),int(r['y']),r['direction']):float(r['mean_product']) for r in csv.DictReader(stream)}
    with (data/'site_moments.csv').open() as stream:
        states={r['receptor_id']:float(r['mean_state']) for r in csv.DictReader(stream)}
    target=np.array([edges[p[0],p[1],'E' if d==0 else 'N'] for p,q,d in physics.edges]+[states[p['id']] for p in model['field_sites']])
    lower=np.r_[np.full(len(physics.edges),model['coupling_bounds'][0]),[model['field_bounds'][0]]*3]
    upper=np.r_[np.full(len(physics.edges),model['coupling_bounds'][1]),[model['field_bounds'][1]]*3]
    fitted,report=newton_fit(physics,target,lower,upper)
    output.mkdir(parents=True,exist_ok=True)
    with (output/'couplings.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\\n');writer.writerow(['x','y','direction','K'])
        for value,(p,q,d) in zip(fitted,physics.edges):writer.writerow([p[0],p[1],'E' if d==0 else 'N',f'{value:.10f}'])
    with (output/'biases.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\\n');writer.writerow(['receptor_id','h'])
        for site,value in zip(model['field_sites'],fitted[-3:]):writer.writerow([site['id'],f'{value:.10f}'])
    print(report,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path('/root/data'))
    parser.add_argument('--output',type=Path,default=Path('/root/results'))
    args=parser.parse_args();fit(args.data,args.output)
'''
(TASK/'solution/inverse_fit.py').write_text(source)
(TASK/'solution/solve.sh').write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\npython3 /solution/inverse_fit.py\n')

independent=ast.parse((AUTHOR/'independent.py').read_text())
independent.body=[n for n in independent.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.ClassDef))]
for node in independent.body:
    if isinstance(node,ast.ClassDef):node.body=[n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='__init__']
field=ast.parse((AUTHOR/'field_independent.py').read_text())
field.body=[n for n in field.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.ClassDef)) and not(isinstance(n,ast.ImportFrom) and n.module=='independent')]
(TASK/'tests/physics.py').write_text(ast.unparse(independent)+'\n'+ast.unparse(field)+'\n')
print(TASK)
