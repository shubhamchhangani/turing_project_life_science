"""Authoring tool (not shipped to the agent): simulate the hierarchical genotype survey.

Model: for each locus, an ancestral allele frequency p0 drifts to regional
frequencies (beta with F_CT), regional frequencies drift to deme frequencies
(beta with F_SC), and genotypes are drawn with within-deme inbreeding F_IS
(p^2 + F pq, 2pq(1 - F), q^2 + F pq), n = 100 diploid individuals per deme.

Usage: python make_counts.py <seed> <out_csv>
"""
import csv
import sys

import numpy as np

REGIONS = {"North": ["N1", "N2"], "Central": ["C1", "C2"], "South": ["S1", "S2"]}
N_PER_DEME = 100

# locus: (ancestral p0, F_CT, F_SC, F_IS) used to simulate
LOCI = {
    "Locus_A": (0.50, 0.15, 0.010, 0.03),
    "Locus_B": (0.30, 0.010, 0.12, 0.00),
    "Locus_C": (0.60, 0.10, 0.08, 0.08),
    "Locus_D": (0.40, 0.005, 0.005, -0.04),
}


def beta_drift(rng, p, F):
    k = (1 - F) / F
    return rng.beta(p * k, (1 - p) * k)


def main():
    seed = int(sys.argv[1])
    out = sys.argv[2]
    rng = np.random.default_rng(seed)
    rows = []
    for locus, (p0, fct, fsc, fis) in LOCI.items():
        for region, demes in REGIONS.items():
            pr = beta_drift(rng, p0, fct)
            for deme in demes:
                p = beta_drift(rng, pr, fsc)
                q = 1 - p
                probs = np.array([p * p + fis * p * q, 2 * p * q * (1 - fis), q * q + fis * p * q])
                probs = np.clip(probs, 0, None)
                probs /= probs.sum()
                aa_counts = rng.multinomial(N_PER_DEME, probs)
                rows.append([locus, region, deme, *map(int, aa_counts)])
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["locus", "region", "deme", "AA", "Aa", "aa"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
