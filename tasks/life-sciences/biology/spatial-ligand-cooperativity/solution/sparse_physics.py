import itertools
import time
import numpy as np
from scipy.linalg import inv
from pfapack.pfaffian import pfaffian
from lattice import Torus, pfaffian_log

class SparseFields(Torus):

    def __init__(self, width, height, sites):
        super().__init__(width, height)
        self.sites = sites
        self.nfields = len(sites)
        self.masks = [m for m in range(1 << len(sites)) if m.bit_count() % 2 == 0]
        self.selected = [[j for j in range(len(sites)) if m >> j & 1] for m in self.masks]
        self.paths = []
        for selected in self.selected:
            modified = self.disorder(np.ones(len(self.edges)), [sites[j] for j in selected])
            flipped = np.flatnonzero(modified < 0)
            self.paths.append((flipped, np.column_stack([2 * flipped, 2 * flipped + 1]).ravel()))
        self.permutation = np.array([self.edge_id[p[0], (p[1] - 1) % height, 1] if d == 0 else self.edge_id[(p[0] - 1) % width, p[1], 0] for p, q, d in self.edges])

    def evaluate(self, parameters):
        n = len(self.edges)
        k = parameters[:n]
        h = parameters[n:]
        weights = self.dual_weights(k)
        t = np.tanh(h)
        coeff = np.array([np.prod(t[s]) for s in self.selected])
        dc = np.zeros((self.nfields, len(self.masks)))
        for term, selected in enumerate(self.selected):
            for j in selected:
                dc[j, term] = (1 - t[j] ** 2) * np.prod(t[[z for z in selected if z != j]])
        totals = np.zeros(len(self.masks))
        derivatives = np.zeros((len(self.masks), n))
        a = np.arange(n) * 2
        b = a + 1
        for tx, ty in itertools.product([0, 1], repeat=2):
            entries = self.reversal / (weights * np.where((tx * self.wrap[:, 0] + ty * self.wrap[:, 1]) % 2, -1, 1))
            matrix = self.base.copy()
            matrix[a, b] = entries
            matrix[b, a] = -entries
            inverse = inv(matrix, check_finite=False)
            sign, logabs = pfaffian_log(matrix)
            base = sign * np.prod(np.sign(entries)) * np.exp(logabs - np.log(abs(entries)).sum()) / 4
            for term, (flipped, indices) in enumerate(self.paths):
                if not len(flipped):
                    value = base
                    diagonal = inverse[b, a]
                    new_entries = entries
                else:
                    update = inverse[np.ix_(indices, indices)].copy()
                    ii = np.arange(len(flipped)) * 2
                    c = 1 / (2 * entries[flipped])
                    update[ii, ii + 1] += c
                    update[ii + 1, ii] -= c
                    p = pfaffian(update, method='P')
                    ratio = p / np.prod(c)
                    value = base * ratio * (-1) ** len(flipped)
                    correction = inverse[np.ix_(b, indices)] @ inv(update, check_finite=False)
                    diagonal = inverse[b, a] - np.einsum('ij,ji->i', correction, inverse[np.ix_(indices, a)])
                    new_entries = entries.copy()
                    new_entries[flipped] *= -1
                totals[term] += value
                derivatives[term] += value * 2 * (diagonal * new_entries - 1)
        total = coeff @ totals
        assert total > 0
        logz = np.log(2) + sum(k) + sum(np.log(np.cosh(h))) + np.log(total)
        gradient = np.r_[1 + (coeff @ derivatives)[self.permutation] / total, t + dc @ totals / total]
        return (logz, gradient)
