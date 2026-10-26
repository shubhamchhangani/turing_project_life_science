import numpy as np
from numba import njit

@njit
def costs(H, starts, pos, alleles, weights):
    m = len(starts) - 1
    k = len(H)
    D = np.zeros((m, k), np.int32)
    for r in range(m):
        for q in range(starts[r], starts[r + 1]):
            for h in range(k):
                if H[h, pos[q]] != alleles[q]:
                    D[r, h] += weights[q]
    return D

@njit
def coordinate(H, starts, pos, alleles, weights, dosage, iterations=50):
    k, n = H.shape
    m = len(starts) - 1
    previous = 2147483647
    for it in range(iterations):
        D = costs(H, starts, pos, alleles, weights)
        assignment = np.empty(m, np.int32)
        obj = 0
        for r in range(m):
            assignment[r] = np.argmin(D[r])
            obj += D[r, assignment[r]]
        if obj >= previous:
            break
        previous = obj
        diff = np.zeros((k, n), np.int32)
        for r in range(m):
            h = assignment[r]
            for q in range(starts[r], starts[r + 1]):
                diff[h, pos[q]] += weights[q] * (1 - 2 * alleles[q])
        for j in range(n):
            order = np.argsort(diff[:, j])
            H[:, j] = 0
            for x in range(dosage[j]):
                H[order[x], j] = 1
    return H

@njit
def anneal(H, starts, pos, alleles, weights, dosage, seed, sweeps=300):
    np.random.seed(seed)
    k, n = H.shape
    m = len(starts) - 1
    counts = np.zeros(n, np.int32)
    for p in pos:
        counts[p] += 1
    offsets = np.zeros(n + 1, np.int32)
    for j in range(n):
        offsets[j + 1] = offsets[j] + counts[j]
    rr = np.empty(len(pos), np.int32)
    aa = rr.copy()
    ww = rr.copy()
    cur = offsets.copy()
    for r in range(m):
        for q in range(starts[r], starts[r + 1]):
            p = pos[q]
            a = cur[p]
            rr[a] = r
            aa[a] = alleles[q]
            ww[a] = weights[q]
            cur[p] += 1
    D = costs(H, starts, pos, alleles, weights)
    bestH = H.copy()
    score = 0
    for r in range(m):
        score += np.min(D[r])
    best = score
    for sweep in range(sweeps):
        T = max(0.12, 3.5 * (1 - sweep / sweeps) ** 2)
        for step in range(n * k):
            j = np.random.randint(n)
            a = np.random.randint(k)
            b = np.random.randint(k)
            if H[a, j] == H[b, j]:
                continue
            delta = 0
            for q in range(offsets[j], offsets[j + 1]):
                r = rr[q]
                old = np.min(D[r])
                v = ww[q] * (1 if H[a, j] == aa[q] else -1)
                D[r, a] += v
                D[r, b] -= v
                delta += np.min(D[r]) - old
            if delta <= 0 or np.random.random() < np.exp(-delta / T):
                H[a, j], H[b, j] = (H[b, j], H[a, j])
                score += delta
                if score < best:
                    best = score
                    bestH = H.copy()
            else:
                for q in range(offsets[j], offsets[j + 1]):
                    r = rr[q]
                    v = ww[q] * (1 if H[a, j] == aa[q] else -1)
                    D[r, a] -= v
                    D[r, b] += v
        for move in range(12):
            a = np.random.randint(k)
            b = np.random.randint(k)
            j = np.random.randint(n)
            if a == b:
                continue
            oldA = D[:, a].copy()
            oldB = D[:, b].copy()
            oldscore = score
            for r in range(m):
                old = np.min(D[r])
                for q in range(starts[r], starts[r + 1]):
                    p = pos[q]
                    if p >= j and H[a, p] != H[b, p]:
                        v = weights[q] * (1 if H[a, p] == alleles[q] else -1)
                        D[r, a] += v
                        D[r, b] -= v
                score += np.min(D[r]) - old
            delta = score - oldscore
            if delta <= 0 or np.random.random() < np.exp(-delta / T):
                for p in range(j, n):
                    H[a, p], H[b, p] = (H[b, p], H[a, p])
                if score < best:
                    best = score
                    bestH = H.copy()
            else:
                D[:, a] = oldA
                D[:, b] = oldB
                score = oldscore
    return (bestH, best)
