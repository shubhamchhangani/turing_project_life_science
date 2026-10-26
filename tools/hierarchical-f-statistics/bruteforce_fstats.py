"""Authoring tool: independent gene-level implementation used to cross-check the oracle.

Expands every individual into its two gene copies (indicator of allele A) and
computes the nested-ANOVA sums of squares directly from group means, then uses
the textbook balanced expected mean squares
    E[MS_G] = G
    E[MS_I] = G + 2 I
    E[MS_S] = G + 2 I + 2n S
    E[MS_C] = G + 2 I + 2n S + 2 s n C
(n individuals per deme, s demes per region). Valid only for balanced data,
which is asserted.

Usage: python bruteforce_fstats.py <counts_csv>   -> prints JSON {locus: {F_IS,...}}
"""
import json
import sys

import numpy as np
import pandas as pd


def expand(df_locus):
    region, deme, ind, gene = [], [], [], []
    iid = 0
    for _, r in df_locus.iterrows():
        for genotype, alleles in (("AA", (1, 1)), ("Aa", (1, 0)), ("aa", (0, 0))):
            for _ in range(int(r[genotype])):
                for a in alleles:
                    region.append(r["region"]); deme.append(r["region"] + "/" + r["deme"])
                    ind.append(iid); gene.append(a)
                iid += 1
    return pd.DataFrame({"region": region, "deme": deme, "ind": ind, "x": np.array(gene, float)})


def fstats(g):
    sizes = g.groupby("deme").size() // 2
    demes_per_region = g.groupby("region")["deme"].nunique()
    assert sizes.nunique() == 1 and demes_per_region.nunique() == 1, "balanced design required"
    n, s, r = int(sizes.iloc[0]), int(demes_per_region.iloc[0]), g["region"].nunique()

    x = g["x"].to_numpy()
    m_ind = g.groupby("ind")["x"].transform("mean").to_numpy()
    m_deme = g.groupby("deme")["x"].transform("mean").to_numpy()
    m_reg = g.groupby("region")["x"].transform("mean").to_numpy()
    m_tot = x.mean()

    ss_G = np.sum((x - m_ind) ** 2)
    ss_I = np.sum((m_ind - m_deme) ** 2)
    ss_S = np.sum((m_deme - m_reg) ** 2)
    ss_C = np.sum((m_reg - m_tot) ** 2)
    ms_G = ss_G / (r * s * n)
    ms_I = ss_I / (r * s * (n - 1))
    ms_S = ss_S / (r * (s - 1))
    ms_C = ss_C / (r - 1)

    G = ms_G
    I = (ms_I - ms_G) / 2
    S = (ms_S - ms_I) / (2 * n)
    C = (ms_C - ms_S) / (2 * s * n)
    tot = C + S + I + G
    return {"F_IS": I / (I + G), "F_SC": S / (S + I + G), "F_CT": C / tot, "F_IT": (C + S + I) / tot}


def main():
    df = pd.read_csv(sys.argv[1])
    out = {loc: fstats(expand(d)) for loc, d in sorted(df.groupby("locus"))}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
