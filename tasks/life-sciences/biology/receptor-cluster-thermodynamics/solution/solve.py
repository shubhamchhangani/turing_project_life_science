import itertools
import math
import time
import numpy as np
from scipy.linalg.blas import dger

def pfaffian_log(matrix):
    """Pivoted skew Gaussian elimination; retain sign instead of sqrt(det)."""
    a = np.array(matrix, dtype=float, order='F', copy=True)
    sign = 1.0
    logabs = 0.0
    for k in range(0, len(a) - 1, 2):
        j = k + 1 + np.argmax(np.abs(a[k, k + 1:]))
        if j != k + 1:
            a[[k + 1, j], :] = a[[j, k + 1], :]
            a[:, [k + 1, j]] = a[:, [j, k + 1]]
            sign = -sign
        pivot = a[k, k + 1]
        if pivot == 0:
            return (0.0, -np.inf)
        sign *= np.sign(pivot)
        logabs += np.log(abs(pivot))
        if k + 2 < len(a):
            u = a[k, k + 2:].copy()
            v = a[k + 1, k + 2:].copy()
            a[k + 2:, k + 2:] += (np.outer(v, u) - np.outer(u, v)) / pivot
    return (sign, logabs)

class Torus:

    def __init__(self, lx, ly):
        self.lx, self.ly = (lx, ly)
        self.edges = []
        self.edge_id = {}
        for y in range(ly):
            for x in range(lx):
                for d in [0, 1]:
                    q = ((x + 1) % lx, y) if d == 0 else (x, (y + 1) % ly)
                    self.edge_id[x, y, d] = len(self.edges)
                    self.edges.append(((x, y), q, d))
        size = 2 * len(self.edges)
        origin = []
        target = []
        theta = []
        for p, q, d in self.edges:
            origin.extend([p, q])
            target.extend([q, p])
            theta.extend([0, math.pi] if d == 0 else [math.pi / 2, -math.pi / 2])
        theta = np.array(theta)
        incoming = {p: [] for p in origin}
        for i, p in enumerate(target):
            incoming[p].append(i)
        b = np.zeros((size, size))
        for group in incoming.values():
            for i, j in itertools.permutations(group, 2):
                turn = (theta[j ^ 1] - theta[i] + math.pi) % (2 * math.pi) - math.pi
                value = -1j * np.exp(1j * (theta[i] + turn / 2 - (theta[i] + theta[j]) / 2))
                assert abs(value.imag) < 1e-12
                b[i, j] = round(value.real)
        self.base = b
        self.reversal = []
        self.wrap = []
        for e, (p, q, d) in enumerate(self.edges):
            i = 2 * e
            j = i + 1
            z = 1j * np.exp(1j * (theta[i] - (theta[i] + theta[j]) / 2))
            self.reversal.append(round(z.real))
            self.wrap.append((d == 0 and p[0] == lx - 1, d == 1 and p[1] == ly - 1))
        self.reversal = np.array(self.reversal)
        self.wrap = np.array(self.wrap, dtype=int)
        assert np.max(np.abs(b + b.T)) == 0

    def polynomials(self, weights):
        values = []
        weights = np.asarray(weights)
        for a, b in itertools.product([0, 1], repeat=2):
            w = weights * np.where((a * self.wrap[:, 0] + b * self.wrap[:, 1]) % 2, -1, 1)
            matrix = self.base.copy()
            index = np.arange(len(self.edges)) * 2
            diagonal = self.reversal / w
            matrix[index, index + 1] = diagonal
            matrix[index + 1, index] = -diagonal
            sign, logabs = pfaffian_log(matrix)
            sign *= np.prod(np.sign(diagonal))
            logabs -= np.sum(np.log(np.abs(diagonal)))
            values.append(sign * np.exp(logabs))
        return np.array(values)

    def dual_weights(self, couplings):
        weights = np.empty(len(self.edges))
        for e, (p, q, d) in enumerate(self.edges):
            x, y = p
            if d == 0:
                dual = self.edge_id[x, (y - 1) % self.ly, 1]
            else:
                dual = self.edge_id[(x - 1) % self.lx, y, 0]
            weights[dual] = np.exp(-2 * couplings[e])
        return weights

    def disorder(self, weights, sites):
        result = weights.copy()
        for start, stop in zip(sites[::2], sites[1::2]):
            x, y = start
            xx, yy = stop
            while x != xx:
                e = self.edge_id[x, (y - 1) % self.ly, 1]
                result[e] *= -1
                x = (x + 1) % self.lx
            while y != yy:
                e = self.edge_id[(x - 1) % self.lx, y, 0]
                result[e] *= -1
                y = (y + 1) % self.ly
        return result

    def evaluate(self, couplings, queries):
        w = self.dual_weights(couplings)
        denominator = self.polynomials(w).sum() / 4
        logz = math.log(2) + sum(couplings) + math.log(denominator)
        correlations = [self.polynomials(self.disorder(w, q)).sum() / 4 / denominator for q in queries]
        return (logz, np.array(correlations))

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
        writer=csv.writer(f,lineterminator='\n');writer.writerow(['pattern','probability'])
        for bits in range(1<<n):
            pattern=format(bits,f'0{n}b')
            spins=[1 if c=='1' else -1 for c in pattern]
            probability=sum(value*np.prod([spins[i] for i in range(n) if mask>>i&1]) for mask,value in moments.items())/(1<<n)
            writer.writerow([pattern,f'{probability:.12f}'])
    with (output_dir/'thermodynamics.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(['quantity','value'])
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
