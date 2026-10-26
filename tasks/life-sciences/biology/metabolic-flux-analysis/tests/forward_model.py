"""Independent data-generator forward model of central carbon metabolism.

Full-isotopomer steady-state simulation (every metabolite carries a probability
vector over its 2^n carbon isotopomers; bit i set = carbon i+1 is 13C). Each
reaction maps substrate carbons to product carbons with an atom-transition
string. Steady state is found by fixed-point iteration of the isotopomer
balances. Used to generate the task data from known ("true") fluxes.
"""
import itertools

import numpy as np

# ----------------------------------------------------------------------------- network
NCARB = {
    "Glc": 6, "G6P": 6, "F6P": 6, "T3P": 3, "P5P": 5, "E4P": 4, "S7P": 7,
    "PG3": 3, "Pyr": 3, "Lac": 3, "Ala": 3, "Ser": 3, "AcCoA": 2, "OAA": 4,
    "Cit": 6, "AKG": 5, "Suc": 4, "Fum": 4, "Mal": 4, "Asp": 4, "Glu": 5,
    "Gln": 5, "CO2": 1,
}
SOURCES = {"Glc", "Gln", "CO2"}          # fixed input labeling
SYMMETRIC = {"Suc", "Fum"}               # rotationally symmetric 4-carbon dicarboxylates

# id: (substrates [(met, atoms)], products [(met, atoms)], reversible)
REACTIONS = {
    "GLCup":  ([("Glc", "abcdef")], [("G6P", "abcdef")], False),
    "PGI":    ([("G6P", "abcdef")], [("F6P", "abcdef")], True),
    "PFK":    ([("F6P", "abcdef")], [("T3P", "cba"), ("T3P", "def")], False),
    "G6PDH":  ([("G6P", "abcdef")], [("P5P", "bcdef"), ("CO2", "a")], False),
    "TK1":    ([("P5P", "abcde"), ("P5P", "fghij")], [("T3P", "cde"), ("S7P", "abfghij")], True),
    "TA":     ([("S7P", "abcdefg"), ("T3P", "hij")], [("E4P", "defg"), ("F6P", "abchij")], True),
    "TK2":    ([("P5P", "abcde"), ("E4P", "fghi")], [("T3P", "cde"), ("F6P", "abfghi")], True),
    "GAPDH":  ([("T3P", "abc")], [("PG3", "abc")], False),
    "PK":     ([("PG3", "abc")], [("Pyr", "abc")], False),
    "LDH":    ([("Pyr", "abc")], [("Lac", "abc")], False),
    "ALT":    ([("Pyr", "abc")], [("Ala", "abc")], False),
    "PHGDH":  ([("PG3", "abc")], [("Ser", "abc")], False),
    "PDH":    ([("Pyr", "abc")], [("AcCoA", "bc"), ("CO2", "a")], False),
    "PC":     ([("Pyr", "abc"), ("CO2", "d")], [("OAA", "abcd")], False),
    "CS":     ([("AcCoA", "ab"), ("OAA", "cdef")], [("Cit", "fedbac")], False),
    "IDH":    ([("Cit", "abcdef")], [("AKG", "abcde"), ("CO2", "f")], True),
    "AKGDH":  ([("AKG", "abcde")], [("Suc", "bcde"), ("CO2", "a")], False),
    "SDH":    ([("Suc", "abcd")], [("Fum", "abcd")], False),
    "FH":     ([("Fum", "abcd")], [("Mal", "abcd")], True),
    "MDH":    ([("Mal", "abcd")], [("OAA", "abcd")], True),
    "ME":     ([("Mal", "abcd")], [("Pyr", "abc"), ("CO2", "d")], False),
    "GOT":    ([("OAA", "abcd")], [("Asp", "abcd")], False),
    "GLNup":  ([("Gln", "abcde")], [("Glu", "abcde")], False),
    "GLUDH":  ([("Glu", "abcde")], [("AKG", "abcde")], True),
    # effluxes / biomass drains (products leave the system)
    "LACout": ([("Lac", "abc")], [], False),
    "ALAout": ([("Ala", "abc")], [], False),
    "SERout": ([("Ser", "abc")], [], False),
    "ASPout": ([("Asp", "abcd")], [], False),
    "GLUout": ([("Glu", "abcde")], [], False),
    "P5Pout": ([("P5P", "abcde")], [], False),
    "ACCOAout": ([("AcCoA", "ab")], [], False),
}

BALANCED = [m for m in NCARB if m not in SOURCES]


def stoichiometry():
    mets = BALANCED
    S = np.zeros((len(mets), len(REACTIONS)))
    for j, (rid, (subs, prods, _)) in enumerate(REACTIONS.items()):
        for m, _a in subs:
            if m in mets:
                S[mets.index(m), j] -= 1
        for m, _a in prods:
            if m in mets:
                S[mets.index(m), j] += 1
    return mets, S


# ----------------------------------------------------------------------------- isotopomer machinery
def natural_dist(n, p13=0.0107):
    idx = np.arange(2 ** n)
    k = np.array([bin(i).count("1") for i in idx])
    return p13 ** k * (1 - p13) ** (n - k)


def positional_dist(probs):
    """Independent per-carbon 13C probabilities -> isotopomer distribution."""
    n = len(probs)
    d = np.ones(1)
    for i in range(n):       # carbon i -> bit i
        d = np.concatenate([d * (1 - probs[i]), d * probs[i]])
    return d


class Transition:
    """Precomputed map from joint substrate isotopomer index -> product isotopomer index."""

    def __init__(self, subs, prod_atoms, prod_met):
        self.subs = [m for m, _ in subs]
        sizes = [2 ** NCARB[m] for m in self.subs]
        letters = {}
        for s_i, (m, atoms) in enumerate(subs):
            assert len(atoms) == NCARB[m], (m, atoms)
            for c_i, ch in enumerate(atoms):
                letters[ch] = (s_i, c_i)
        assert len(prod_atoms) == NCARB[prod_met], (prod_met, prod_atoms)
        grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing="ij")
        flat = [g.ravel() for g in grids]
        out = np.zeros(flat[0].size, dtype=np.int64)
        for p_i, ch in enumerate(prod_atoms):
            s_i, c_i = letters[ch]
            out |= ((flat[s_i] >> c_i) & 1) << p_i
        self.map = out
        self.nout = 2 ** NCARB[prod_met]

    def apply(self, dists):
        joint = dists[self.subs[0]]
        for m in self.subs[1:]:
            joint = np.multiply.outer(joint, dists[m])
        return np.bincount(self.map, weights=joint.ravel(), minlength=self.nout)


def reverse_reaction(subs, prods):
    return prods, subs


def build_transitions():
    """List of (rid, direction, product_met, Transition, weight)."""
    out = []
    for rid, (subs, prods, rev) in REACTIONS.items():
        dirs = [("f", subs, prods)]
        if rev:
            dirs.append(("b", prods, subs))
        for d, s, p in dirs:
            for pm, pa in p:
                if pm in SOURCES:
                    continue
                variants = [pa]
                if pm in SYMMETRIC:
                    variants = [pa, pa[::-1]]
                for v in variants:
                    out.append((rid, d, pm, Transition(s, v, pm), 1.0 / len(variants)))
    return out


TRANSITIONS = build_transitions()


def directional_fluxes(net, exch):
    """net, exch: dicts by reaction id. Returns dict (rid, 'f'|'b') -> flux."""
    vd = {}
    for rid, (_s, _p, rev) in REACTIONS.items():
        v = net[rid]
        if rev:
            x = exch.get(rid, 0.0)
            vd[(rid, "f")] = x + max(v, 0.0)
            vd[(rid, "b")] = x + max(-v, 0.0)
        else:
            vd[(rid, "f")] = max(v, 0.0)
    return vd


def simulate(net, exch, sources, tol=1e-13, max_iter=20000, init=None):
    """Steady-state isotopomer distributions for all metabolites."""
    vd = directional_fluxes(net, exch)
    dists = {m: natural_dist(NCARB[m]) for m in NCARB} if init is None else {k: v.copy() for k, v in init.items()}
    dists.update(sources)
    influx = {m: 0.0 for m in BALANCED}
    for rid, d, pm, tr, w in TRANSITIONS:
        influx[pm] += vd[(rid, d)] * w
    order = BALANCED
    for it in range(max_iter):
        delta = 0.0
        for m in order:
            acc = np.zeros(2 ** NCARB[m])
            for rid, d, pm, tr, w in TRANSITIONS_BY_PRODUCT[m]:
                f = vd[(rid, d)] * w
                if f > 0:
                    acc += f * tr.apply(dists)
            # renormalize: with bimolecular reactions the normalized state is
            # numerically unstable and round-off would otherwise grow each sweep
            new = acc / acc.sum()
            delta = max(delta, np.max(np.abs(new - dists[m])))
            dists[m] = new
        if delta < tol:
            break
    else:
        raise RuntimeError(f"isotopomer iteration did not converge (delta={delta:.2e})")
    return dists


TRANSITIONS_BY_PRODUCT = {m: [t for t in TRANSITIONS if t[2] == m] for m in BALANCED}


# ----------------------------------------------------------------------------- measurement model
ISOTOPES = {  # element: [(mass shift, abundance)]
    "C": [(0, 0.9893), (1, 0.0107)],
    "H": [(0, 0.999885), (1, 0.000115)],
    "N": [(0, 0.99636), (1, 0.00364)],
    "O": [(0, 0.99757), (1, 0.00038), (2, 0.00205)],
    "Si": [(0, 0.92223), (1, 0.04685), (2, 0.03092)],
}

# fragment id: (metabolite, backbone carbons (1-based), formula of the whole ion)
FRAGMENTS = {
    "Lac261": ("Lac", (1, 2, 3), {"C": 11, "H": 25, "O": 3, "Si": 2}),
    "Ala260": ("Ala", (1, 2, 3), {"C": 11, "H": 26, "N": 1, "O": 2, "Si": 2}),
    "Ala232": ("Ala", (2, 3), {"C": 10, "H": 26, "N": 1, "O": 1, "Si": 2}),
    "Ser390": ("Ser", (1, 2, 3), {"C": 17, "H": 40, "N": 1, "O": 3, "Si": 3}),
    "Asp418": ("Asp", (1, 2, 3, 4), {"C": 18, "H": 40, "N": 1, "O": 4, "Si": 3}),
    "Asp390": ("Asp", (2, 3, 4), {"C": 17, "H": 40, "N": 1, "O": 3, "Si": 3}),
    "Glu432": ("Glu", (1, 2, 3, 4, 5), {"C": 19, "H": 42, "N": 1, "O": 4, "Si": 3}),
    "Glu330": ("Glu", (2, 3, 4, 5), {"C": 16, "H": 36, "N": 1, "O": 2, "Si": 2}),
    "Fum287": ("Fum", (1, 2, 3, 4), {"C": 12, "H": 23, "O": 4, "Si": 2}),
    "Mal419": ("Mal", (1, 2, 3, 4), {"C": 18, "H": 39, "O": 5, "Si": 3}),
    "Cit459": ("Cit", (1, 2, 3, 4, 5, 6), {"C": 20, "H": 39, "O": 6, "Si": 3}),
}


def element_dist(el, count):
    d = np.ones(1)
    base = np.zeros(max(s for s, _ in ISOTOPES[el]) + 1)
    for s, a in ISOTOPES[el]:
        base[s] = a
    for _ in range(count):
        d = np.convolve(d, base)
    return d


def nonbackbone_dist(formula, n_backbone):
    d = np.ones(1)
    for el, cnt in formula.items():
        if el == "C":
            cnt -= n_backbone
        d = np.convolve(d, element_dist(el, cnt))
    return d


NB_DIST = {f: nonbackbone_dist(form, len(carbs)) for f, (_m, carbs, form) in FRAGMENTS.items()}


def backbone_mid(dist, carbons):
    n = int(np.log2(dist.size))
    idx = np.arange(dist.size)
    k = np.zeros(dist.size, dtype=int)
    for c in carbons:
        k += (idx >> (c - 1)) & 1
    return np.bincount(k, weights=dist, minlength=len(carbons) + 1)


def raw_mid(dists, frag):
    met, carbons, _form = FRAGMENTS[frag]
    full = np.convolve(backbone_mid(dists[met], carbons), NB_DIST[frag])
    m = len(carbons) + 3          # report M+0 ... M+(n+2)
    out = full[:m]
    return out / out.sum()


# ----------------------------------------------------------------------------- tracer experiments
def sources_for(experiment):
    if experiment == "glc12":
        glc = positional_dist([0.99, 0.99, 0.0107, 0.0107, 0.0107, 0.0107])
        gln = natural_dist(5)
    elif experiment == "gln5":
        glc = natural_dist(6)
        gln = positional_dist([0.99] * 5)
    else:
        raise ValueError(experiment)
    return {"Glc": glc, "Gln": gln, "CO2": natural_dist(1)}


EXPERIMENTS = ("glc12", "gln5")


def all_combos():
    return list(itertools.product(EXPERIMENTS, FRAGMENTS))
