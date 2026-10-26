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
