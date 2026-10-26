import argparse
import csv
import itertools
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import inv
from scipy.optimize import minimize
from lattice import Torus, pfaffian_log

class FieldTorus(Torus):

    def __init__(self, width, height, field_sites):
        super().__init__(width, height)
        assert len(field_sites) == 3
        self.field_sites = field_sites
        self.masks = [0, 3, 5, 6]
        self.edge_permutation = []
        for p, q, d in self.edges:
            x, y = p
            self.edge_permutation.append(self.edge_id[x, (y - 1) % height, 1] if d == 0 else self.edge_id[(x - 1) % width, y, 0])
        self.edge_permutation = np.array(self.edge_permutation)

    def evaluate_gradient(self, parameters, hessian=False):
        couplings = parameters[:-3]
        fields = parameters[-3:]
        weights = self.dual_weights(couplings)
        index = np.arange(len(self.edges)) * 2
        moments = []
        derivatives = []
        second_derivatives = []
        for mask in self.masks:
            sites = [p for i, p in enumerate(self.field_sites) if mask >> i & 1]
            current = self.disorder(weights, sites)
            total = 0.0
            derivative = np.zeros(len(self.edges))
            second = np.zeros((len(self.edges), len(self.edges))) if hessian else None
            for tx, ty in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                w = current * np.where((tx * self.wrap[:, 0] + ty * self.wrap[:, 1]) % 2, -1, 1)
                b = self.base.copy()
                entries = self.reversal / w
                b[index, index + 1] = entries
                b[index + 1, index] = -entries
                inverse = inv(b, overwrite_a=False, check_finite=False)
                sign, logabs = pfaffian_log(b)
                sign *= np.prod(np.sign(entries))
                logabs -= np.log(abs(entries)).sum()
                value = sign * np.exp(logabs)
                total += value / 4
                u = inverse[index + 1, index] * entries
                gradient = 2 * (u - 1)
                derivative += value / 4 * gradient
                if hessian:
                    cross = inverse[np.ix_(index + 1, index)]
                    bb = inverse[np.ix_(index + 1, index + 1)]
                    aa = inverse[np.ix_(index, index)]
                    log_hessian = -4 * np.outer(entries, entries) * (cross * cross.T - bb * aa.T)
                    log_hessian[np.diag_indices_from(log_hessian)] += 4 * u
                    second += value / 4 * (log_hessian + np.outer(gradient, gradient))
            moments.append(total)
            derivatives.append(derivative[self.edge_permutation])
            if hessian:
                second_derivatives.append(second[np.ix_(self.edge_permutation, self.edge_permutation)])
        t = np.tanh(fields)
        coefficients = np.array([1.0, t[0] * t[1], t[0] * t[2], t[1] * t[2]])
        sum_loops = np.dot(coefficients, moments)
        assert sum_loops > 0
        logz = np.log(2) + sum(couplings) + np.log(np.cosh(fields)).sum() + np.log(sum_loops)
        edge_gradient = 1 + coefficients @ np.array(derivatives) / sum_loops
        jacobian = np.array([[0, t[1], t[2], 0], [0, t[0], 0, t[2]], [0, 0, t[0], t[1]]]) * (1 - t * t)[:, None]
        magnetization = t + jacobian @ np.array(moments) / sum_loops
        gradient = np.r_[edge_gradient, magnetization]
        if not hessian:
            return (logz, gradient)
        n = len(couplings)
        covariance = np.empty((n + 3, n + 3))
        covariance[:n, :n] = np.tensordot(coefficients, np.array(second_derivatives), axes=1) / sum_loops - np.outer(edge_gradient - 1, edge_gradient - 1)
        cross = jacobian @ np.array(derivatives) / sum_loops - np.outer(magnetization - t, edge_gradient - 1)
        covariance[n:, :n] = cross
        covariance[:n, n:] = cross.T
        second_coefficients = np.zeros((3, 3, 4))
        for term, mask in enumerate(self.masks):
            for i in range(3):
                if not mask >> i & 1:
                    continue
                others = [j for j in range(3) if mask >> j & 1 and j != i]
                second_coefficients[i, i, term] = -2 * t[i] * (1 - t[i] ** 2) * np.prod(t[others])
                for j in others:
                    second_coefficients[i, j, term] = (1 - t[i] ** 2) * (1 - t[j] ** 2)
        covariance[n:, n:] = np.diag(1 - t * t) + second_coefficients @ np.array(moments) / sum_loops - np.outer(magnetization - t, magnetization - t)
        covariance = (covariance + covariance.T) / 2
        return (logz, gradient, covariance)
from scipy.linalg import solve
def newton_fit(model, observed, lower, upper, max_iterations=30):
    parameters = np.r_[np.full(len(model.edges), 0.35), [0.05, 0.05, 0.05]]
    start = time.monotonic()
    calls = 0
    value, gradient, hessian = model.evaluate_gradient(parameters, True)
    calls += 1
    objective = value - parameters @ observed
    for iteration in range(max_iterations):
        residual = gradient - observed
        print('Newton iteration', iteration, 'calls', calls, 'objective', objective, 'gradient', max(abs(residual)), 'elapsed', time.monotonic() - start, flush=True)
        if max(abs(residual)) < 2e-08:
            return (parameters, {'iterations': iteration, 'calls': calls, 'seconds': time.monotonic() - start, 'max_moment_error': float(max(abs(residual))), 'converged': True, 'minimum_fisher_eigenvalue': float(np.linalg.eigvalsh(hessian)[0])})
        step = solve(hessian, -residual, assume_a='pos', check_finite=False)
        alpha = 1.0
        positive = step > 0
        negative = step < 0
        if positive.any():
            alpha = min(alpha, 0.99 * min((upper[positive] - parameters[positive]) / step[positive]))
        if negative.any():
            alpha = min(alpha, 0.99 * min((lower[negative] - parameters[negative]) / step[negative]))
        for backtrack in range(20):
            candidate = parameters + alpha * step
            nv, ng, nh = model.evaluate_gradient(candidate, True)
            calls += 1
            no = nv - candidate @ observed
            if no <= objective + 0.0001 * alpha * residual @ step + 1e-11:
                parameters, value, gradient, hessian, objective = (candidate, nv, ng, nh, no)
                break
            alpha *= 0.5
        else:
            raise RuntimeError('Newton line search did not converge')
    raise RuntimeError('Newton iteration allowance exhausted')

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
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['x','y','direction','K'])
        for value,(p,q,d) in zip(fitted,physics.edges):writer.writerow([p[0],p[1],'E' if d==0 else 'N',f'{value:.10f}'])
    with (output/'biases.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['receptor_id','h'])
        for site,value in zip(model['field_sites'],fitted[-3:]):writer.writerow([site['id'],f'{value:.10f}'])
    print(report,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,default=Path('/root/data'))
    parser.add_argument('--output',type=Path,default=Path('/root/results'))
    args=parser.parse_args();fit(args.data,args.output)
