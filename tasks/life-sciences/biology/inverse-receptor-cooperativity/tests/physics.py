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
from scipy.linalg import inv
from pfapack.pfaffian import pfaffian

class IndependentFieldTorus(DecoratedTorus):

    def __init__(self, width, height, sites):
        super().__init__(width, height)
        self.sites = sites
        self.primal_to_dual = []
        for x, y, d in self.lookup:
            other = (x, (y - 1) % height, 1) if d == 0 else ((x - 1) % width, y, 0)
            self.primal_to_dual.append(self.lookup[other])

    def sector_derivatives(self, weights):
        polynomial = 0.0
        derivative = np.zeros(len(weights))
        for tx, ty in itertools.product([0, 1], repeat=2):
            matrix = self.matrix.copy()
            baseline_sign = self.match_sign
            baseline_log = 0.0
            entries = []
            for weight, (i, j, sign, wx, wy) in zip(weights, self.links):
                value = sign / weight * (-1) ** (tx * wx + ty * wy)
                matrix[i, j] = value
                matrix[j, i] = -value
                baseline_sign *= np.sign(value)
                baseline_log += math.log(abs(value))
                entries.append(value)
            inverse = inv(matrix, check_finite=False)
            answer = pfaffian(matrix * 0.5, method='H', overwrite_a=True)
            value = np.sign(answer) * baseline_sign * math.exp(math.log(abs(answer)) + len(matrix) / 2 * math.log(2) - baseline_log)
            gradient = np.array([2 * (inverse[j, i] * entry - 1) for entry, (i, j, _, _, _) in zip(entries, self.links)])
            polynomial += value / 4
            derivative += value / 4 * gradient
        return (polynomial, derivative[self.primal_to_dual])

    def observables(self, parameters):
        k = parameters[:-3]
        h = parameters[-3:]
        w = np.empty(len(k))
        for index, dual in enumerate(self.primal_to_dual):
            w[dual] = math.exp(-2 * k[index])
        masks = [0, 3, 5, 6]
        values = []
        derivatives = []
        for mask in masks:
            selected = [p for i, p in enumerate(self.sites) if mask >> i & 1][::-1]
            modified = w.copy()
            for (x, y), (xx, yy) in zip(selected[::2], selected[1::2]):
                while y != yy:
                    y = (y - 1) % self.height
                    modified[self.lookup[(x - 1) % self.width, y, 0]] *= -1
                while x != xx:
                    x = (x - 1) % self.width
                    modified[self.lookup[x, (y - 1) % self.height, 1]] *= -1
            value, gradient = self.sector_derivatives(modified)
            values.append(value)
            derivatives.append(gradient)
        t = np.tanh(h)
        coeff = np.array([1, t[0] * t[1], t[0] * t[2], t[1] * t[2]])
        total = coeff @ values
        logz = math.log(2) + sum(k) + sum(np.log(np.cosh(h))) + math.log(total)
        edges = 1 + coeff @ np.array(derivatives) / total
        magnetization = []
        for i in range(3):
            change = sum((values[j] * np.prod([t[a] for a in range(3) if masks[j] >> a & 1 and a != i]) for j in range(4) if masks[j] >> i & 1))
            magnetization.append(t[i] + (1 - t[i] ** 2) * change / total)
        return (logz, np.r_[edges, magnetization])
