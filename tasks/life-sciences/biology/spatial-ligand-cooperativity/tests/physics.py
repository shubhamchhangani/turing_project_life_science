import itertools
import math
import numpy as np
from pfapack.pfaffian import pfaffian

class DecoratedTorus:

    def __init__(self, width, height):
        self.width, self.height = (width, height)
        local = np.array([[0, -1, -1, 1], [1, 0, -1, 1], [1, 1, 0, 1], [-1, -1, -1, 0]], dtype=float)
        self.matrix = np.kron(np.eye(width * height), local)
        self.links = []
        permutation = []
        self.lookup = {}
        for y in range(height):
            for x in range(width):
                src = 4 * (y * width + x)
                for direction in [0, 1]:
                    if direction == 0:
                        a = 4 * (y * width + (x + 1) % width)
                        b = src + 2
                        sign = 1
                    else:
                        a = 4 * ((y + 1) % height * width + x) + 1
                        b = src + 3
                        sign = -1
                    self.lookup[x, y, direction] = len(self.links)
                    self.links.append((a, b, sign, int(direction == 0 and x == width - 1), int(direction == 1 and y == height - 1)))
                    permutation.extend([a, b])
        inversions = sum((permutation[i] > permutation[j] for i in range(len(permutation)) for j in range(i + 1, len(permutation))))
        self.match_sign = (-1.0) ** inversions
import itertools
import math
import numpy as np
from scipy.linalg import inv, solve
from pfapack.pfaffian import pfaffian

class IndependentSparseFields(DecoratedTorus):

    def __init__(self, width, height, sites):
        super().__init__(width, height)
        self.sites = sites
        self.nfields = len(sites)
        self.permutation = [self.lookup[(x, (y - 1) % height, 1) if d == 0 else ((x - 1) % width, y, 0)] for x, y, d in self.lookup]
        self.subsets = [s for count in range(0, len(sites) + 1, 2) for s in itertools.combinations(range(len(sites)), count)]
        self.paths = []
        for subset in self.subsets:
            flips = np.zeros(len(self.links), dtype=bool)
            selected = [sites[j] for j in subset][::-1]
            for (x, y), (xx, yy) in zip(selected[::2], selected[1::2]):
                while y != yy:
                    y = (y - 1) % height
                    flips[self.lookup[(x - 1) % width, y, 0]] ^= True
                while x != xx:
                    x = (x - 1) % width
                    flips[self.lookup[x, (y - 1) % height, 1]] ^= True
            edge_ids = np.flatnonzero(flips)
            ids = np.array([v for edge in edge_ids for v in self.links[edge][:2]], dtype=int)
            self.paths.append((edge_ids, ids))

    def observables(self, parameters):
        n = len(self.links)
        k = parameters[:n]
        h = parameters[n:]
        weights = np.empty(n)
        for primal, dual in enumerate(self.permutation):
            weights[dual] = math.exp(-2 * k[primal])
        polynomials = np.zeros(len(self.subsets))
        derivatives = np.zeros((len(self.subsets), n))
        a = np.array([edge[0] for edge in self.links])
        b = np.array([edge[1] for edge in self.links])
        for tx, ty in itertools.product([0, 1], repeat=2):
            matrix = self.matrix.copy()
            baseline_sign = self.match_sign
            baseline_log = 0.0
            entries = []
            for weight, (i, j, sign, wx, wy) in zip(weights, self.links):
                entry = sign / weight * (-1) ** (tx * wx + ty * wy)
                matrix[i, j] = entry
                matrix[j, i] = -entry
                entries.append(entry)
                baseline_sign *= np.sign(entry)
                baseline_log += math.log(abs(entry))
            entries = np.array(entries)
            inverse = inv(matrix, check_finite=False)
            p = pfaffian(matrix * 0.5, method='H', overwrite_a=True)
            baseline = np.sign(p) * baseline_sign * math.exp(math.log(abs(p)) + len(matrix) / 2 * math.log(2) - baseline_log) / 4
            for term, (edge_ids, ids) in enumerate(self.paths):
                if not len(ids):
                    value = baseline
                    diagonal = inverse[b, a]
                    current = entries
                else:
                    q = inverse[np.ix_(ids, ids)].copy()
                    r = np.arange(len(edge_ids)) * 2
                    c = 1 / (2 * entries[edge_ids])
                    q[r, r + 1] += c
                    q[r + 1, r] -= c
                    value = baseline * (-1) ** len(edge_ids) * pfaffian(q, method='H') / np.prod(c)
                    correction = solve(q, inverse[np.ix_(ids, a)], check_finite=False)
                    diagonal = inverse[b, a] - np.einsum('ij,ji->i', inverse[np.ix_(b, ids)], correction)
                    current = entries.copy()
                    current[edge_ids] *= -1
                polynomials[term] += value
                derivatives[term] += value * 2 * (diagonal * current - 1)
        t = np.tanh(h)
        coefficients = np.array([np.prod(t[list(s)]) for s in self.subsets])
        total = coefficients @ polynomials
        edges = 1 + (coefficients @ derivatives)[self.permutation] / total
        fields = []
        for i in range(self.nfields):
            change = sum((p * np.prod(t[[j for j in s if j != i]]) for s, p in zip(self.subsets, polynomials) if i in s))
            fields.append(t[i] + (1 - t[i] ** 2) * change / total)
        logz = math.log(2) + sum(k) + sum(np.log(np.cosh(h))) + math.log(total)
        return (logz, np.r_[edges, fields])
