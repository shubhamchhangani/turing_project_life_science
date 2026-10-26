"""Independent vertex-ordered decorated-graph evaluator using Householder Pfaffians.

The reference uses edge-ordered skew elimination. Here both the indexing and
the elimination algorithm differ; spin insertions also use different paths.
"""
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

    def sectors(self, weights):
        values = []
        for tx, ty in itertools.product([0, 1], repeat=2):
            matrix = self.matrix.copy()
            baseline_sign = self.match_sign
            baseline_log = 0.0
            for weight, (i, j, sign, wx, wy) in zip(weights, self.links):
                value = sign / weight * (-1) ** (tx * wx + ty * wy)
                matrix[i, j] = value
                matrix[j, i] = -value
                baseline_sign *= np.sign(value)
                baseline_log += math.log(abs(value))
            answer = pfaffian(matrix * 0.5, method='H', overwrite_a=True)
            if answer == 0:
                values.append(0.0)
            else:
                values.append(np.sign(answer) * baseline_sign * math.exp(math.log(abs(answer)) + len(matrix) / 2 * math.log(2) - baseline_log))
        return np.array(values)

    def partition(self, couplings):
        p = self.sectors(np.tanh(couplings))
        return self.width * self.height * math.log(2) + sum(np.log(np.cosh(couplings))) + math.log((-p[0] + p[1] + p[2] + p[3]) / 2)

    def moments(self, couplings, probe_sites):
        w = np.empty(len(self.links))
        for (x, y, d), index in self.lookup.items():
            other = (x, (y - 1) % self.height, 1) if d == 0 else ((x - 1) % self.width, y, 0)
            w[self.lookup[other]] = math.exp(-2 * couplings[index])
        base = self.sectors(w).sum()
        moments = {0: 1.0}
        count = len(probe_sites)
        for mask in range(1, 1 << count):
            selected = [probe_sites[i] for i in range(count) if mask >> i & 1]
            if len(selected) % 2:
                moments[mask] = 0.0
                continue
            modified = w.copy()
            selected.reverse()
            for (x, y), (xx, yy) in zip(selected[::2], selected[1::2]):
                while y != yy:
                    y = (y - 1) % self.height
                    modified[self.lookup[(x - 1) % self.width, y, 0]] *= -1
                while x != xx:
                    x = (x - 1) % self.width
                    modified[self.lookup[x, (y - 1) % self.height, 1]] *= -1
            moments[mask] = self.sectors(modified).sum() / base
        return moments

    def joint(self, couplings, probe_sites):
        moments = self.moments(couplings, probe_sites)
        n = len(probe_sites)
        result = []
        for configuration in range(1 << n):
            spins = [1 if configuration >> n - 1 - i & 1 else -1 for i in range(n)]
            result.append(sum((m * np.prod([spins[i] for i in range(n) if mask >> i & 1]) for mask, m in moments.items())) / (1 << n))
        return np.array(result)
