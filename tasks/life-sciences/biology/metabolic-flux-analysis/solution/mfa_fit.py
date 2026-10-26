"""Reference solution: steady-state 13C-MFA from two parallel tracer experiments.

1. Carbon atom transitions for every reaction in /root/data/network.tsv
   (standard biochemistry; succinate and fumarate are rotationally symmetric).
2. Full-isotopomer steady-state simulation of both experiments.
3. Measurement model: backbone MID of each fragment convolved with the natural
   isotope distribution of every non-backbone atom of the ion, truncated to
   M+0..M+(n+2) and renormalized (fitting the raw data directly is equivalent to
   correcting it first, and keeps the noise model honest).
4. Weighted least squares over MIDs and extracellular rates, net fluxes on the
   stoichiometric null space, exchange fluxes bounded, multi-start.

Usage: python mfa_fit.py [data_dir] [out_csv] [n_starts] [seed]
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/data")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "/root/results/fluxes.csv")
N_STARTS = int(sys.argv[3]) if len(sys.argv) > 3 else 24
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 7

NCARB = {"Glc": 6, "G6P": 6, "F6P": 6, "T3P": 3, "P5P": 5, "E4P": 4, "S7P": 7, "PG3": 3,
         "Pyr": 3, "Lac": 3, "Ala": 3, "Ser": 3, "AcCoA": 2, "OAA": 4, "Cit": 6, "AKG": 5,
         "Suc": 4, "Fum": 4, "Mal": 4, "Asp": 4, "Glu": 5, "Gln": 5, "CO2": 1}
SOURCES = {"Glc", "Gln", "CO2"}
SYMMETRIC = {"Suc", "Fum"}

# Atom transitions (substrate atoms -> product atoms). Effluxes have no products.
ATOMS = {
    "GLCup": ([("Glc", "abcdef")], [("G6P", "abcdef")]),
    "PGI": ([("G6P", "abcdef")], [("F6P", "abcdef")]),
    "PFK": ([("F6P", "abcdef")], [("T3P", "cba"), ("T3P", "def")]),
    "G6PDH": ([("G6P", "abcdef")], [("P5P", "bcdef"), ("CO2", "a")]),
    "TK1": ([("P5P", "abcde"), ("P5P", "fghij")], [("T3P", "cde"), ("S7P", "abfghij")]),
    "TA": ([("S7P", "abcdefg"), ("T3P", "hij")], [("E4P", "defg"), ("F6P", "abchij")]),
    "TK2": ([("P5P", "abcde"), ("E4P", "fghi")], [("T3P", "cde"), ("F6P", "abfghi")]),
    "GAPDH": ([("T3P", "abc")], [("PG3", "abc")]),
    "PK": ([("PG3", "abc")], [("Pyr", "abc")]),
    "LDH": ([("Pyr", "abc")], [("Lac", "abc")]),
    "ALT": ([("Pyr", "abc")], [("Ala", "abc")]),
    "PHGDH": ([("PG3", "abc")], [("Ser", "abc")]),
    "PDH": ([("Pyr", "abc")], [("AcCoA", "bc"), ("CO2", "a")]),
    "PC": ([("Pyr", "abc"), ("CO2", "d")], [("OAA", "abcd")]),
    "CS": ([("AcCoA", "ab"), ("OAA", "cdef")], [("Cit", "fedbac")]),
    "IDH": ([("Cit", "abcdef")], [("AKG", "abcde"), ("CO2", "f")]),
    "AKGDH": ([("AKG", "abcde")], [("Suc", "bcde"), ("CO2", "a")]),
    "SDH": ([("Suc", "abcd")], [("Fum", "abcd")]),
    "FH": ([("Fum", "abcd")], [("Mal", "abcd")]),
    "MDH": ([("Mal", "abcd")], [("OAA", "abcd")]),
    "ME": ([("Mal", "abcd")], [("Pyr", "abc"), ("CO2", "d")]),
    "GOT": ([("OAA", "abcd")], [("Asp", "abcd")]),
    "GLNup": ([("Gln", "abcde")], [("Glu", "abcde")]),
    "GLUDH": ([("Glu", "abcde")], [("AKG", "abcde")]),
    "LACout": ([("Lac", "abc")], []), "ALAout": ([("Ala", "abc")], []),
    "SERout": ([("Ser", "abc")], []), "ASPout": ([("Asp", "abcd")], []),
    "GLUout": ([("Glu", "abcde")], []), "P5Pout": ([("P5P", "abcde")], []),
    "ACCOAout": ([("AcCoA", "ab")], []),
}
FREE_NET = ["G6PDH", "P5Pout", "SERout", "ALAout", "ASPout", "GLUout", "ACCOAout", "GLNup", "ME", "CS"]


# ----------------------------------------------------------------------------- inputs
def read_inputs():
    with open(DATA / "network.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rxns = [r["reaction_id"] for r in rows]
    rev = {r["reaction_id"]: r["reversible"] == "yes" for r in rows}
    assert set(rxns) == set(ATOMS), "network.tsv does not match the atom-transition model"

    iso = {}
    with open(DATA / "isotope_abundances.csv") as fh:
        for r in csv.DictReader(fh):
            iso.setdefault(r["element"], []).append((int(r["mass_shift"]), float(r["abundance"])))

    frags = {}
    with open(DATA / "fragments.csv") as fh:
        for r in csv.DictReader(fh):
            carbons = tuple(int(c) for c in r["backbone_carbons"].split(";"))
            frags[r["fragment"]] = (r["metabolite"], carbons, parse_formula(r["ion_formula"]))

    meas = {}
    with open(DATA / "mid_measurements.csv") as fh:
        for r in csv.DictReader(fh):
            meas.setdefault((r["experiment"], r["fragment"]), []).append(
                (int(r["isotopologue"][1:]), float(r["fraction"]), float(r["sd"])))
    for k in meas:
        meas[k].sort()

    rates = {}
    with open(DATA / "extracellular_rates.csv") as fh:
        for r in csv.DictReader(fh):
            rates[r["reaction_id"]] = (float(r["rate"]), float(r["sd"]))
    return rxns, rev, iso, frags, meas, rates


def parse_formula(s):
    import re
    return {el: int(n) if n else 1 for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", s)}


# ----------------------------------------------------------------------------- isotopomer simulation
def natural(n, p=0.0107):
    out = np.ones(1)
    for _ in range(n):
        out = np.concatenate([out * (1 - p), out * p])
    return out


def positional(probs):
    out = np.ones(1)
    for p in probs:
        out = np.concatenate([out * (1 - p), out * p])
    return out


class Transition:
    def __init__(self, subs, prod_met, prod_atoms):
        self.subs = [m for m, _ in subs]
        where = {}
        for i, (m, atoms) in enumerate(subs):
            for j, ch in enumerate(atoms):
                where[ch] = (i, j)
        grids = np.meshgrid(*[np.arange(2 ** NCARB[m]) for m in self.subs], indexing="ij")
        flat = [g.ravel() for g in grids]
        idx = np.zeros(flat[0].size, dtype=np.int64)
        for k, ch in enumerate(prod_atoms):
            i, j = where[ch]
            idx |= ((flat[i] >> j) & 1) << k
        self.idx, self.n = idx, 2 ** NCARB[prod_met]

    def __call__(self, d):
        joint = d[self.subs[0]]
        for m in self.subs[1:]:
            joint = np.multiply.outer(joint, d[m])
        return np.bincount(self.idx, weights=joint.ravel(), minlength=self.n)


class Simulator:
    def __init__(self, rxns, rev):
        self.rxns, self.rev = rxns, rev
        self.balanced = [m for m in NCARB if m not in SOURCES]
        self.by_product = {m: [] for m in self.balanced}
        for rid in rxns:
            subs, prods = ATOMS[rid]
            dirs = [("f", subs, prods)] + ([("b", prods, subs)] if rev[rid] else [])
            for d, s, p in dirs:
                for pm, pa in p:
                    if pm in SOURCES:
                        continue
                    variants = [pa, pa[::-1]] if pm in SYMMETRIC else [pa]
                    for v in variants:
                        self.by_product[pm].append((rid, d, Transition(s, pm, v), 1 / len(variants)))

    def run(self, vf, vb, sources, init=None, tol=1e-12):
        d = {m: natural(NCARB[m]) for m in NCARB} if init is None else dict(init)
        d.update(sources)
        for _ in range(50000):
            delta = 0.0
            for m in self.balanced:
                acc = np.zeros(2 ** NCARB[m])
                for rid, dr, tr, w in self.by_product[m]:
                    f = (vf if dr == "f" else vb)[rid] * w
                    if f > 0:
                        acc += f * tr(d)
                new = acc / acc.sum()
                delta = max(delta, float(np.max(np.abs(new - d[m]))))
                d[m] = new
            if delta < tol:
                return d
        raise RuntimeError("isotopomer balances did not converge")


# ----------------------------------------------------------------------------- measurement model
def nonbackbone(formula, n_backbone, iso):
    out = np.ones(1)
    for el, cnt in formula.items():
        cnt = cnt - n_backbone if el == "C" else cnt
        base = np.zeros(max(s for s, _ in iso[el]) + 1)
        for s, a in iso[el]:
            base[s] = a
        for _ in range(cnt):
            out = np.convolve(out, base)
    return out


def fragment_mid(dist, carbons, nb, n_report):
    idx = np.arange(dist.size)
    k = np.zeros(dist.size, dtype=int)
    for c in carbons:
        k += (idx >> (c - 1)) & 1
    bb = np.bincount(k, weights=dist, minlength=len(carbons) + 1)
    full = np.convolve(bb, nb)[:n_report]
    return full / full.sum()


# ----------------------------------------------------------------------------- flux parameterization
def stoich(rxns):
    mets = [m for m in NCARB if m not in SOURCES]
    S = np.zeros((len(mets), len(rxns)))
    for j, rid in enumerate(rxns):
        subs, prods = ATOMS[rid]
        for m, _ in subs:
            if m in mets:
                S[mets.index(m), j] -= 1
        for m, _ in prods:
            if m in mets:
                S[mets.index(m), j] += 1
    return S


class FluxModel:
    def __init__(self, rxns, rev):
        self.rxns, self.rev = rxns, rev
        S = stoich(rxns)
        self.fixed = {"GLCup": 100.0}
        self.free = FREE_NET
        indep = list(self.fixed) + self.free
        dep = [r for r in rxns if r not in indep]
        Si = S[:, [rxns.index(r) for r in indep]]
        Sd = S[:, [rxns.index(r) for r in dep]]
        assert np.linalg.matrix_rank(Sd) == len(dep), "free flux set does not determine the network"
        self.dep, self.indep = dep, indep
        self.M = -np.linalg.lstsq(Sd, Si, rcond=None)[0]        # v_dep = M @ v_indep
        self.exch = [r for r in rxns if rev[r]]

    def net(self, z):
        vi = np.r_[[self.fixed["GLCup"]], z]
        vd = self.M @ vi
        v = dict(zip(self.indep, vi))
        v.update(zip(self.dep, vd))
        return v

    def directional(self, v, x):
        vf, vb = {}, {}
        ex = dict(zip(self.exch, x))
        for r in self.rxns:
            if self.rev[r]:
                vf[r] = ex[r] + max(v[r], 0.0)
                vb[r] = ex[r] + max(-v[r], 0.0)
            else:
                vf[r] = max(v[r], 0.0)
        return vf, vb


# ----------------------------------------------------------------------------- fit
def tracer_sources():
    """Labeling of the carbon inputs: unlabeled carbon carries natural 13C (1.07 %)."""
    return {
        "glc12": {"Glc": positional([0.99, 0.99] + [0.0107] * 4), "Gln": natural(5), "CO2": natural(1)},
        "gln5": {"Glc": natural(6), "Gln": positional([0.99] * 5), "CO2": natural(1)},
    }


def main():
    rxns, rev, iso, frags, meas, rates = read_inputs()
    sim = Simulator(rxns, rev)
    fm = FluxModel(rxns, rev)
    sources = tracer_sources()
    nb = {f: nonbackbone(form, len(c), iso) for f, (_m, c, form) in frags.items()}
    keys = sorted(meas)
    y = np.concatenate([[v for _, v, _ in meas[k]] for k in keys])
    sd = np.concatenate([[s for _, _, s in meas[k]] for k in keys])
    rate_ids = [r for r in rates if rates[r][1] > 0]
    cache = {}

    def residuals(p):
        z, x = p[: len(fm.free)], p[len(fm.free):]
        v = fm.net(z)
        vf, vb = fm.directional(v, x)
        pen = [1e3 * min(v[r], 0.0) for r in rxns if not rev[r]]
        sims = {}
        for exp, src in sources.items():
            sims[exp] = sim.run(vf, vb, src, init=cache.get(exp))
            cache[exp] = sims[exp]
        pred = np.concatenate([
            fragment_mid(sims[e][frags[f][0]], frags[f][1], nb[f], len(meas[(e, f)])) for e, f in keys])
        r_mid = (pred - y) / sd
        r_rate = [(v[r] - rates[r][0]) / rates[r][1] for r in rate_ids]
        return np.concatenate([r_mid, r_rate, pen])

    rng = np.random.default_rng(SEED)
    lo = np.r_[np.zeros(len(fm.free)), np.zeros(len(fm.exch))]
    hi = np.r_[np.full(len(fm.free), 400.0), np.full(len(fm.exch), 500.0)]
    best = None
    for s in range(N_STARTS):
        z0 = np.array([rates[r][0] if r in rates else rng.uniform(1, 60) for r in fm.free])
        x0 = rng.uniform(1, 150, len(fm.exch))
        p0 = np.clip(np.r_[z0, x0], lo + 1e-6, hi - 1e-6)
        cache.clear()
        try:
            res = least_squares(residuals, p0, bounds=(lo, hi), x_scale="jac", max_nfev=400, diff_step=1e-6)
        except RuntimeError:
            continue
        print(f"start {s:2d}: SSR = {2 * res.cost:9.2f}", flush=True)
        if best is None or res.cost < best.cost:
            best = res
    v = fm.net(best.x[: len(fm.free)])
    n_data = len(y) + len(rate_ids)
    print(f"best SSR = {2 * best.cost:.2f} with {n_data} measurements and {best.x.size} free parameters")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    vf, vb = fm.directional(v, best.x[len(fm.free):])
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction_id", "net_flux", "forward_flux", "reverse_flux"])
        for r in rxns:
            w.writerow([r, f"{v[r]:.4f}", f"{vf[r]:.4f}", f"{vb.get(r, 0.0):.4f}"])
    for r in ("G6PDH", "PDH", "PC", "ME", "CS", "IDH", "GLUDH", "LDH"):
        print(f"  {r:6s} {v[r]:8.3f}")


if __name__ == "__main__":
    main()
