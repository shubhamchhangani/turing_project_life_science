"""Atom transitions and a stoichiometric null-space flux parameterization."""
import numpy as np
NCARB = {'Glc': 6, 'G6P': 6, 'F6P': 6, 'T3P': 3, 'P5P': 5, 'E4P': 4, 'S7P': 7, 'PG3': 3, 'Pyr': 3, 'Lac': 3, 'Ala': 3, 'Ser': 3, 'AcCoA': 2, 'OAA': 4, 'Cit': 6, 'AKG': 5, 'Suc': 4, 'Fum': 4, 'Mal': 4, 'Asp': 4, 'Glu': 5, 'Gln': 5, 'CO2': 1}
SOURCES = {'Glc', 'Gln', 'CO2'}
ATOMS = {'GLCup': ([('Glc', 'abcdef')], [('G6P', 'abcdef')]), 'PGI': ([('G6P', 'abcdef')], [('F6P', 'abcdef')]), 'PFK': ([('F6P', 'abcdef')], [('T3P', 'cba'), ('T3P', 'def')]), 'G6PDH': ([('G6P', 'abcdef')], [('P5P', 'bcdef'), ('CO2', 'a')]), 'TK1': ([('P5P', 'abcde'), ('P5P', 'fghij')], [('T3P', 'cde'), ('S7P', 'abfghij')]), 'TA': ([('S7P', 'abcdefg'), ('T3P', 'hij')], [('E4P', 'defg'), ('F6P', 'abchij')]), 'TK2': ([('P5P', 'abcde'), ('E4P', 'fghi')], [('T3P', 'cde'), ('F6P', 'abfghi')]), 'GAPDH': ([('T3P', 'abc')], [('PG3', 'abc')]), 'PK': ([('PG3', 'abc')], [('Pyr', 'abc')]), 'LDH': ([('Pyr', 'abc')], [('Lac', 'abc')]), 'ALT': ([('Pyr', 'abc')], [('Ala', 'abc')]), 'PHGDH': ([('PG3', 'abc')], [('Ser', 'abc')]), 'PDH': ([('Pyr', 'abc')], [('AcCoA', 'bc'), ('CO2', 'a')]), 'PC': ([('Pyr', 'abc'), ('CO2', 'd')], [('OAA', 'abcd')]), 'CS': ([('AcCoA', 'ab'), ('OAA', 'cdef')], [('Cit', 'fedbac')]), 'IDH': ([('Cit', 'abcdef')], [('AKG', 'abcde'), ('CO2', 'f')]), 'AKGDH': ([('AKG', 'abcde')], [('Suc', 'bcde'), ('CO2', 'a')]), 'SDH': ([('Suc', 'abcd')], [('Fum', 'abcd')]), 'FH': ([('Fum', 'abcd')], [('Mal', 'abcd')]), 'MDH': ([('Mal', 'abcd')], [('OAA', 'abcd')]), 'ME': ([('Mal', 'abcd')], [('Pyr', 'abc'), ('CO2', 'd')]), 'GOT': ([('OAA', 'abcd')], [('Asp', 'abcd')]), 'GLNup': ([('Gln', 'abcde')], [('Glu', 'abcde')]), 'GLUDH': ([('Glu', 'abcde')], [('AKG', 'abcde')]), 'LACout': ([('Lac', 'abc')], []), 'ALAout': ([('Ala', 'abc')], []), 'SERout': ([('Ser', 'abc')], []), 'ASPout': ([('Asp', 'abcd')], []), 'GLUout': ([('Glu', 'abcde')], []), 'P5Pout': ([('P5P', 'abcde')], []), 'ACCOAout': ([('AcCoA', 'ab')], [])}
FREE_NET = ['G6PDH', 'P5Pout', 'SERout', 'ALAout', 'ASPout', 'GLUout', 'ACCOAout', 'GLNup', 'ME', 'CS']

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
        self.rxns, self.rev = (rxns, rev)
        S = stoich(rxns)
        self.fixed = {'GLCup': 100.0}
        self.free = FREE_NET
        indep = list(self.fixed) + self.free
        dep = [r for r in rxns if r not in indep]
        Si = S[:, [rxns.index(r) for r in indep]]
        Sd = S[:, [rxns.index(r) for r in dep]]
        assert np.linalg.matrix_rank(Sd) == len(dep), 'free flux set does not determine the network'
        self.dep, self.indep = (dep, indep)
        self.M = -np.linalg.lstsq(Sd, Si, rcond=None)[0]
        self.exch = [r for r in rxns if rev[r]]

    def net(self, z):
        vi = np.r_[[self.fixed['GLCup']], z]
        vd = self.M @ vi
        v = dict(zip(self.indep, vi))
        v.update(zip(self.dep, vd))
        return v

    def directional(self, v, x):
        vf, vb = ({}, {})
        ex = dict(zip(self.exch, x))
        for r in self.rxns:
            if self.rev[r]:
                vf[r] = ex[r] + max(v[r], 0.0)
                vb[r] = ex[r] + max(-v[r], 0.0)
            else:
                vf[r] = max(v[r], 0.0)
        return (vf, vb)
