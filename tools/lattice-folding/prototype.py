"""Original HP sequence generator and input-only core-threading experiment."""
import argparse
import itertools
import json
import random
import time
from pathlib import Path
from ortools.sat.python import cp_model


def snake(size):
    plane = [(x, y) for y in range(size) for x in (range(size) if y % 2 == 0 else reversed(range(size)))]
    return [(x, y, z) for z in range(size) for x, y in (plane if z % 2 == 0 else reversed(plane))]


def generate(size, core, seed, steps=30000):
    rng = random.Random(seed)
    path = snake(size)
    for _ in range(steps):
        if rng.randrange(2):
            path.reverse()
        head = path[0]
        options = [j for j in range(2, len(path)) if sum(abs(a-b) for a,b in zip(head,path[j])) == 1]
        if options:
            j = rng.choice(options)
            path[:j] = reversed(path[:j])
    sequence = ''.join('H' if all(0 <= v < core for v in xyz) else 'P' for xyz in path)
    assert len(set(path)) == len(path)
    assert all(sum(abs(a-b) for a,b in zip(path[i],path[i+1]))==1 for i in range(len(path)-1))
    return sequence, path


def thread(sequence, core, shell, seconds=120, seed=17, log=False):
    sites = list(itertools.product(range(-shell, core+shell), repeat=3))
    lookup = {xyz:j for j,xyz in enumerate(sites)}
    label = ['H' if all(0 <= v < core for v in xyz) else 'P' for xyz in sites]
    model = cp_model.CpModel()
    variables = [model.new_int_var_from_domain(cp_model.Domain.from_values([j for j,x in enumerate(label) if x == aa]), f'p{i}') for i,aa in enumerate(sequence)]
    model.add_all_different(variables)
    tables = {s:[] for s in ['HH','HP','PH','PP']}
    for j,xyz in enumerate(sites):
        for axis in range(3):
            for step in [-1,1]:
                other = list(xyz)
                other[axis] += step
                k = lookup.get(tuple(other))
                if k is not None:
                    tables[label[j]+label[k]].append((j,k))
    for i in range(len(sequence)-1):
        model.add_allowed_assignments(variables[i:i+2], tables[sequence[i:i+2]])
    h0 = sequence.index('H')
    representatives = [j for j,xyz in enumerate(sites) if all(0 <= v < core for v in xyz) and xyz[0] <= xyz[1] <= xyz[2] <= (core-1)//2]
    model.add_allowed_assignments([variables[h0]], [(j,) for j in representatives])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = seed
    solver.parameters.log_search_progress = log
    start = time.monotonic()
    status = solver.solve(model)
    result = {'status':solver.status_name(status), 'seconds':time.monotonic()-start, 'branches':solver.num_branches, 'conflicts':solver.num_conflicts}
    if status in [cp_model.OPTIMAL,cp_model.FEASIBLE]:
        result['coordinates']=[sites[solver.value(v)] for v in variables]
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--size',type=int,default=5)
    p.add_argument('--core',type=int,default=4)
    p.add_argument('--seed',type=int,default=20261025)
    p.add_argument('--seconds',type=int,default=120)
    p.add_argument('--shell',type=int,default=1)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--log',action='store_true')
    a=p.parse_args()
    seq,witness=generate(a.size,a.core,a.seed)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/'generated.json').write_text(json.dumps({'sequence':seq,'witness':witness},indent=2)+'\n')
    print(f'N={len(seq)}, H={seq.count("H")}, HH={sum(seq[i:i+2]=="HH" for i in range(len(seq)-1))}',flush=True)
    result=thread(seq,a.core,a.shell,a.seconds,log=a.log)
    (a.out/'threaded.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='coordinates'},flush=True)
